import sys
import subprocess
import re
import time
import datetime
import shutil
import os
import os.path

BUILDSCRIPT = './build.sh'
RUNSCRIPT = 'debug.run'
SLEEPINTERVAL = 10
MAXRUNNINGJOBS = 5

lineRegex = re.compile('([^ ]+) = (.*)')
traceExitCodeRegex = re.compile('Exit_status=([^ ]*)')
traceTimeRegex = re.compile('resources_used.walltime=(\d+:\d+:\d+)')
traceCpuRegex = re.compile('resources_used.cput=(\d+:\d+:\d+)')

key_walltime = 'resources_used.walltime'
key_cputime = 'resources_used.cput'
key_jobstate = 'job_state'
key_nodes = 'Resource_List.nodes'
key_exit_status = 'exit_status'
value_jobstate_done = 'C'
value_jobstate_queued = 'Q'
value_jobstate_running = 'R'
value_jobstate_none = 'N'

def run(args):
    out = str(subprocess.Popen(args, stdout=subprocess.PIPE).stdout.read())
    #print ' '.join(args), ':', out
    return [s.strip() for s in out.split('\n') if len(s) > 0]

def runTask(nodes, ppn):
    taskName = run(['qsub', '-l', 'nodes=%s:ppn=%s' % (nodes, ppn), RUNSCRIPT])[0]
    return taskName.split('.')[0]

def monitorTask(taskId):
    out = run(['qstat', '-f', taskId])
    params = dict()
    if len(out) > 0:
        for line in out:
            match = lineRegex.search(line)
            if match:
                params[match.group(1)] = match.group(2)
    return params


def parseTime(timeStr):
    splits = [int(x) for x in timeStr.split(':')]
    return splits[0] * 3600 + splits[1] * 60 + splits[2]

def parseNodes(nodes):
    splits = [int(s.strip('=')) for s in nodes.split(':ppn')]
    return ( splits[0], splits[1] )

def taskStatus(taskId):
    data = monitorTask(taskId)
    jobState = value_jobstate_none
    wallTime = 0
    cpuTime = 0
    nodes = 0
    procs = 0
    exitCode = ''

    if key_jobstate in data:
        jobState = data[key_jobstate]
        (nodes, procs) = parseNodes(data[key_nodes])
        if jobState == value_jobstate_done:
            wallTime = parseTime(data[key_walltime])
            cpuTime = parseTime(data[key_cputime])
            exitCode = data[key_exit_status]

    return {
        'state': jobState,
        'time': wallTime,
        'cpu': cpuTime,
        'nodes': nodes,
        'ppn': procs,
        'code': exitCode
    }

def taskTrace(taskId, np):
    out = str(subprocess.Popen(['tracejob', '-q', taskId], stdout=subprocess.PIPE).stdout.read())

    jobState = value_jobstate_none
    wallTime = 0
    cpuTime = 0
    (nodes, procs) = np
    exitCode = ''

    if 'Job Queued' in out:
        jobState = value_jobstate_queued
        if 'Job Run' in out:
            jobState = value_jobstate_running
            if 'Exit_status' in out:
                jobState = value_jobstate_done
                match = traceExitCodeRegex.search(out)
                if match:
                    exitCode = match.group(1)
                match = traceTimeRegex.search(out)
                if match:
                    wallTime = parseTime(match.group(1))
                match = traceCpuRegex.search(out)
                if match:
                    cpuTime = parseTime(match.group(1))
    return {
        'state': jobState,
        'time': wallTime,
        'cpu': cpuTime,
        'nodes': nodes,
        'ppn': procs,
        'code': exitCode
    }

def createOutputFolder():
    currTime = datetime.datetime.now()
    folderName = currTime.strftime('%Y-%m-%d_%H%M%S')
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    return folderName

def copyInputs(outDir):
    if os.path.isfile("config.ini"):
        shutil.copyfile("config.ini", "%s/config.ini" % outDir)

def copyTaskOutputs(taskId, procCount, outDir):
    for prefix in ['e', 'o']:
        fname = "%s.%s%s" % (RUNSCRIPT, prefix, taskId)
        if not os.path.isfile(fname):
            time.sleep(1)
        if os.path.isfile(fname):
            shutil.move(fname, "%s/%s%sx%s" % (outDir, prefix, taskId, procCount))

def waitState(state):
    return state == value_jobstate_queued or state == value_jobstate_running

def saveTaskStatus(taskId, procCount, outDir):
    ftaskStatus = open('%s/s%sx%s' % (outDir, taskId, procCount), 'w')
    #status = str(subprocess.Popen(['qstat', '-f', taskId], stdout=subprocess.PIPE).stdout.read())
    status = str(subprocess.Popen(['tracejob', '-q', taskId], stdout=subprocess.PIPE).stdout.read())
    ftaskStatus.write(status)
    ftaskStatus.close()

def saveResults(times, outDir):
    fstatus = open('%s/status.txt' % outDir, 'w')
    fresults = open('%s/results.txt' % outDir, 'w')

    zeroTime = -1
    for task, result in times:
        fstatus.write('%s(x%s)\t: [%s:%s] w: %s\tc: %s\tnodes=%s:ppn=%s\n' % (task, result['nodes'] * result['ppn'], result['state'], result['code'], result['time'], result['cpu'], result['nodes'], result['ppn']))
        time = result['time']
        procs = result['nodes'] * result['ppn']
        if zeroTime == -1:
            zeroTime = time
        speedUp = 0
        eff = 0
        if zeroTime > 0:
            speedUp = float(time) / zeroTime
            eff = float(time) / zeroTime / procs
        fresults.write("%s\t%s\t%s\t%s\n" % (procs, time, speedUp, eff))

    fstatus.close()
    fresults.close()

def processTasks(tasksParams):
    buildResult = run(BUILDSCRIPT)
    if len(buildResult) > 0:
        print buildResult
        return []

    outDir = createOutputFolder()
    copyInputs(outDir)

    jobsRun = 0
    tasks = []
    taskParamsDict = dict()
    for (nodes, ppn) in tasksParams:
        taskId = runTask(nodes, ppn)
        tasks.append(taskId)
        taskParamsDict[taskId] = (nodes, ppn)
        print "Run %s" % (nodes * ppn)
        time.sleep(2)
        jobsRun += 1
        if jobsRun >= MAXRUNNINGJOBS:
            break

    results = dict()
    anyRuning = True

    while anyRuning:
        anyRuning = False
        sys.stderr.write("\x1b[2J\x1b[H")
        print 'Status:'
        runningJobs = 0

        for t in tasks:
            result = None
            if not t in results or waitState(results[t]['state']):
                result = taskTrace(t, taskParamsDict[t])
                results[t] = result
                anyRuning = anyRuning or waitState(result['state'])
                if not waitState(result['state']):
                    copyTaskOutputs(t, result['nodes'] * result['ppn'], outDir)
                    saveTaskStatus(t, result['nodes'] * result['ppn'], outDir)
                else:
                    runningJobs += 1
            else:
                result = results[t]

            print '%s(x%s)\t: [%s:%s] w: %s\tc: %s\tnodes=%s:ppn=%s' % (t, result['nodes'] * result['ppn'], result['state'], result['code'], result['time'], result['cpu'], result['nodes'], result['ppn'])
        if runningJobs < MAXRUNNINGJOBS and len(tasks) < len(tasksParams):
            exJobsRun = jobsRun
            for (nodes, ppn) in tasksParams[exJobsRun:]:
                taskId = runTask(nodes, ppn)
                tasks.append(taskId)
                taskParamsDict[taskId] = (nodes, ppn)
                print "Run %s" % (nodes * ppn)
                time.sleep(2)
                runningJobs += 1
                jobsRun += 1
                if runningJobs >= MAXRUNNINGJOBS:
                    break

        if anyRuning:
            time.sleep(SLEEPINTERVAL)

    times = sorted(results.iteritems())
    saveResults(times, outDir)
    return times

def printResults(times):
    zeroTime = -1
    for task, result in times:
        time = result['time']
        procs = result['nodes'] * result['ppn']
        if zeroTime == -1:
            zeroTime = time
        speedUp = 0
        eff = 0
        if zeroTime > 0:
            speedUp = float(time) / zeroTime
            eff = float(time) / zeroTime / procs
        print "%s\t%s\t%s\t%s" % (procs, time, speedUp, eff)


#print runTask(4, 4)
#print taskStatus('186314')

#times = processTasks([(5, 6), (3, 7)])
#print taskTrace('188032', (3, 4))
#times = processTasks([(1, 1), (1, 2), (1, 4), (2, 4), (3, 4), (3, 5), (3, 6), (3, 7)])
times = processTasks([(1, 1), (1, 2), (2, 2), (2, 4), (3, 4), (3, 5), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (7, 7), (8, 7), (9, 7)])
#times = processTasks([(1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (3, 4)])
# 6 8 10 12 16 20 24 28 30 36 42
#times = processTasks([(2, 3), (2, 4), (2, 5), (3, 4), (4, 4), (4, 5), (4, 6), (4, 7), (5, 6), (6, 6), (6, 7)])

printResults(times)