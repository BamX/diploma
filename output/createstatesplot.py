import sys
import re
import subprocess

stateRe = re.compile('\\[(\\d+)\\]\\((\\d+)\\):direction (\\d):row (\\d+)-(\\d+):pass (\\d)')
infoRe = re.compile('\\[(\\d+)\\]:w (\\d+):h (\\d+):t (-?\\d+):b (-?\\d+)')
minTimeDelta = 30

valueZ = 0
valueF = 1
valueS = 2

def readData():
    f = open('states.txt', 'r')
    data = f.read()
    f.close()
    return data

def parseInfo(data):
    results = dict()
    fromW = 0
    fromH = 0
    for m in infoRe.finditer(data):
        (pid, width, height, topN, bottomN) = [int(x) for x in m.groups()]
        results[pid] = (fromW, fromH, width, height, topN, bottomN)
        fromH += height
    return results

def parseStates(data):
    return [[int(x) for x in m.groups()] for m in stateRe.finditer(data)]

def processTimeline(states):
    results = dict()
    currentStates = map(lambda x: None, info)
    for (pid, time, direction, fRow, tRow, pas) in sorted(states, key=lambda v:v[1]):
        currentStates[pid] = (direction, fRow, tRow, pas)
        if time in results:
            results[time][pid] = (direction, fRow, tRow, pas)
        else:
            results[time] = currentStates[:]
    return sorted(list(results.iteritems()))

def proceessMatrixTimeline(info, timeline):
    allHeight = 0
    allWidth = 0
    for pid, (fromW, fromH, width, height, topN, bottomN) in info.iteritems():
        allWidth = width
        allHeight += height

    nextTime = 0
    results = dict()
    for time, states in timeline:
        matrix = [[valueZ for x in range(allWidth)] for x in range(allHeight)]
        for pid in range(len(states)):
            if (states[pid]):
                (direction, fRow, tRow, pas) = states[pid]
                (fromW, fromH, width, height, topN, bottomN) = info[pid]
                value = (valueF, valueS)[pas - 1]
                if direction == 1:
                    for row in range(fromH + fRow, fromH + tRow + 1):
                        for col in range(allWidth):
                            matrix[row][col] = value
                else:
                    for col in range(fromW + fRow, fromW + tRow + 1):
                        for row in range(fromH, fromH + height):
                            matrix[row][col] = value
        while nextTime < time:
            results[nextTime] = matrix
            nextTime += minTimeDelta
    return sorted(list(results.iteritems()))

def outputMatrixTimeline(matrixTimeline):
    f = open('states.out', 'w')
    for time, matrix in matrixTimeline:
        for row in range(len(matrix)):
            for col in range(len(matrix[row])):
                s = '%.1f' % float(matrix[row][col])
                if col < len(matrix[row]) - 1:
                    s = s + ' '
                f.write(s)
            f.write('\n')
        f.write('\n')
    f.close()

data = readData()
info = parseInfo(data)
states = parseStates(data)
timeline = processTimeline(states)
matrixTimeline = proceessMatrixTimeline(info, timeline)
outputMatrixTimeline(matrixTimeline)
print 'Done'
