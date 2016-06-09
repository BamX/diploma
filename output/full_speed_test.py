import sys
import subprocess
import re
import time
import datetime
import shutil
import os
import os.path
import json
from random import shuffle

BUILDSCRIPT = './build.sh'
BINARYNAME = 'debug'
SLEEPINTERVAL = 10
MAXRUNNINGJOBS = 10
PBSQUEUE = "max"
QSUBALLOWED = True

lineRegex = re.compile('([^ ]+) = (.*)')
additionalLineRegex = re.compile('\t(.*)')
nodesListRegex = re.compile('hadoop2-(\d+)/(\d+)')
traceExitCodeRegex = re.compile('Exit_status=([^ ]*)')
traceTimeRegex = re.compile('resources_used.walltime=(\d+:\d+:\d+)')
traceCpuRegex = re.compile('resources_used.cput=(\d+:\d+:\d+)')

key_walltime = 'resources_used.walltime'
key_cputime = 'resources_used.cput'
key_jobstate = 'job_state'
key_nodes = 'Resource_List.nodes'
key_nodes_list = 'exec_host'
key_exit_status = 'exit_status'
key_mem = 'resources_used.mem'
key_vmem = 'resources_used.vmem'
key_start_time = 'start_time'
value_jobstate_done = 'C'
value_jobstate_queued = 'Q'
value_jobstate_running = 'R'
value_jobstate_none = 'N'
configTemplate = '''
X1      0.15
X2      0.125

Speed   0.0125

InitT   1768
EnvT    303

X1SplitCount %d
X2SplitCount %d
TimeSplitCount %d
Epsilon %f
TMax %d
Repeats %d

MinimumBundle %d
BalanceFactor %f
EnableBalancing %d
TransposeBalanceFactor %f
TransposeBalanceIterationsInterval %d
TransposeBalanceTimeFactor %f
StaticBalanceThresholdFactor %f

# 0 for transpose
# 1 for static
Algorithm %d

EnableConsole %d
EnablePlot %d
EnableMatrix %d
EnableBuckets %d
FramesCount %d

PlotFilename %s
BucketsFilename %s

ViewCount 5
DebugView 4

View0X1  0.0
View0X2  0.115

View1X1  0.14
View1X2  0.115

View2X1  0.14
View2X2  0.0

View3X1  0.0
View3X2  0.0

View4X1  0.15
View4X2  0.125
'''

PBSScriptTemplate = ''':
#PBS -l nodes=%s:ppn=%s
#PBS -d %s
#PBS -e %s
#PBS -o %s
#PBS -q %s
#
mpirun %s %s
'''

def formatConfig(params, taskName, workingDir):
    return configTemplate % (
            params.get('x_split_count', 100), params.get('x_split_count', 100),
            params.get('time_split_count', 20000),
            params.get('epsilon', 0.001),
            params.get('time_max', 600),
            params.get('repeats', 1),
            params.get('bundle_min', 5),
            params.get('balance_factor', 0.3),
            (0, 1)[params.get('enable_balancing', True)],
            params.get('transpose_balance_factor', 0.92),
            params.get('transpose_balance_interval', 15),
            params.get('transpose_balance_time_factor', 1),
            params.get('static_balance_threshold_factor', 0.1),
            params.get('algorithm', 0),
            (0, 1)[params.get('enable_console', False)],
            (0, 1)[params.get('enable_plot', False)],
            (0, 1)[params.get('enable_matrix', False)],
            (0, 1)[params.get('enable_buckets', True)],
            params.get('frames_count', 200),
            '%s/view-%s.csv' % (workingDir, taskName),
            '%s/buckets-%s.csv' % (workingDir, taskName)
        )

def saveConfig(params, taskName, workingDir):
    filename = 'config-%s.ini' % taskName
    with open('%s/%s' % (workingDir, filename), 'w') as configFile:
        configFile.write(formatConfig(params, taskName, workingDir))
    return filename

def formatPBSScript(params, taskName, workingDir):
    return PBSScriptTemplate % (
            params.get('nodes', 3), params.get('ppn', 3),
            workingDir,
            'error-%s.txt' % taskName,
            'output-%s.txt' % taskName,
            PBSQUEUE,
            '../' + BINARYNAME, saveConfig(params, taskName, workingDir)
        )

def savePBSScript(params, taskName, workingDir):
    filename = '%s/pbs-%s.run' % (workingDir, taskName)
    with open('%s' % filename, 'w') as pbsFile:
        pbsFile.write(formatPBSScript(params, taskName, workingDir))
    return filename

def run(args):
    out = str(subprocess.Popen(args, stdout=subprocess.PIPE).stdout.read())
    #print ' '.join(args), ':', out
    return [s.strip() for s in out.split('\n') if len(s) > 0]

def runTask(params, taskName, workingDir):
    pbsFilename = savePBSScript(params, taskName, workingDir)
    taskId = run(['qsub', pbsFilename])[0]
    return taskId.split('.')[0]

def monitorTask(taskId):
    out = run(['qstat', '-f', taskId])
    params = dict()
    lastKey = None
    if len(out) > 0:
        for line in out:
            additional = additionalLineRegex.search(line)
            if additional and lastKey:
                params[lastKey] += additional.group(1)
            else:
                match = lineRegex.search(line)
                if match:
                    lastKey = match.group(1)
                    params[lastKey] = match.group(2)
    return params

def parseTime(timeStr):
    splits = [int(x) for x in timeStr.split(':')]
    return splits[0] * 3600 + splits[1] * 60 + splits[2]

def parseNodes(nodes):
    splits = [int(s.strip('=')) for s in nodes.split(':ppn')]
    return ( splits[0], splits[1] )

def parseNodesList(nodesList):
    nodes = nodesListRegex.findall(nodesList)
    return [ {'node': int(node), 'process': int(process)} for (node, process) in nodes ]

def taskStatus(taskId):
    wallTime = 0
    cpuTime = 0
    memUsed = 0
    vmemUsed = 0
    startTime = None
    nodesList = []

    try:
        data = monitorTask(taskId)
        if key_start_time in data:
            startTime = datetime.datetime.strptime(data[key_start_time], '%a %b %d %H:%M:%S %Y').isoformat()
        if key_nodes_list in data:
            nodesList = parseNodesList(data[key_nodes_list])
        if key_walltime in data:
            wallTime = parseTime(data[key_walltime])
        if key_cputime in data:
            cpuTime = parseTime(data[key_cputime])
        if key_mem in data:
            memUsed = int(data[key_mem][:-2]) * 1024
        if key_vmem in data:
            vmemUsed = int(data[key_vmem][:-2]) * 1024
    except Exception, e:
        pass

    return {
        'time': wallTime,
        'cpu': cpuTime,
        'nodes_list': nodesList,
        'mem': memUsed,
        'vmem': vmemUsed,
        'start_time': startTime
    }

def waitState(state):
    return state != value_jobstate_done

def taskTrace(taskId, case):
    out = str(subprocess.Popen(['tracejob', '-q', taskId], stdout=subprocess.PIPE).stdout.read())

    jobState = value_jobstate_none
    wallTime = 0
    cpuTime = 0
    exitCode = ''
    case = case.copy()

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

    if QSUBALLOWED and waitState(case.get('result_state', 'Q')):
        status = taskStatus(taskId)
        for key in status:
            if status[key]:
                case['result_' + key] = status[key]

    case['result_state'] = jobState
    if not waitState(jobState) or not 'result_time' in case or not 'result_cpu' in case:
        case['result_time'] = wallTime
        case['result_cpu'] = cpuTime
    case['result_code'] = exitCode
    return case

def createWorkingDirectory():
    currTime = datetime.datetime.now()
    folderName = currTime.strftime('%Y-%m-%d_%H%M%S')
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    return os.getcwd() + '/' + folderName

def saveTaskStatus(taskId, taskName, workingDir):
    with open('%s/status-%s.txt' % (workingDir, taskName), 'w') as ftaskStatus:
        status = str(subprocess.Popen(['tracejob', '-q', taskId], stdout=subprocess.PIPE).stdout.read())
        ftaskStatus.write(status)

def readTimeFromOutput(case, workingDir):
    try:
        with open('%s/error-%s.txt' % (workingDir, case['name']), 'r') as errorFile:
            time = float(errorFile.read())
            case['result_time'] = time
    except Exception, e:
        pass

def saveResults(cases, workingDir):
    with open('%s/cases.json' % workingDir, 'w') as fstatus:
        fstatus.write(json.dumps(cases))
    with open('%s/results.txt' % workingDir, 'w') as fresults:
        zeroTime = -1
        for result in cases:
            time = result['result_time'] / result['repeats']
            procs = result['nodes'] * result['ppn']
            if zeroTime == -1:
                zeroTime = time
            speedUp = 0
            eff = 0
            if zeroTime > 0:
                speedUp = float(time) / zeroTime
                eff = float(time) / zeroTime / procs
            fresults.write("%s\t%s\t%s\t%s\n" % (procs, time, speedUp, eff))

def readTestCases(filename):
    testCases = []
    with open(filename, 'r') as fin:
        plain = fin.read()
        data = json.loads(plain)
        configuration = data['configuration']
        templates = data['templates']
        defaults = data['default']
        for caseTag in data['cases']:
            caseData = data['cases'][caseTag]
            caseDefaults = dict()
            if 'default' in caseData:
                caseDefaults = caseData['default']
            cases = []
            if 'casesTemplate' in caseData:
                cases = templates[caseData['casesTemplate']]
            elif 'cases' in caseData:
                cases = caseData['cases']
            for case in cases:
                for repeatNum in xrange(configuration['cases_repeats']):
                    testCase = dict()
                    for key in defaults:
                        if key in case:
                            testCase[key] = case[key]
                        elif key in caseDefaults:
                            testCase[key] = caseDefaults[key]
                        else:
                            testCase[key] = defaults[key]
                    testCase['case_tag'] = caseTag
                    testCase['name'] = 'n%d-p%d-x%d-r%d' % \
                        (testCase['nodes'], testCase['ppn'], testCase['x_split_count'], testCase['repeats'])
                    testCases += [testCase]
    return testCases

def caseProcesses(case):
    return case['nodes'] * case['ppn']

def reorderCases(cases):
    newCases = cases[:]
    shuffle(newCases)
    newCases.sort(key=caseProcesses)
    buckets = []
    limitBucket = 13
    while len(newCases) > 0:
        bucket = []
        bSum = 0
        while len(newCases) > 0 and bSum + newCases[-1]['nodes'] < limitBucket:
            bucket += [newCases[-1]]
            bSum += newCases[-1]['nodes']
            newCases = newCases[:-1]
        while len(newCases) > 0 and bSum + newCases[0]['nodes'] < limitBucket:
            bucket += [newCases[0]]
            bSum += newCases[0]['nodes']
            newCases = newCases[1:]
        buckets += [bucket]
    shuffle(buckets)
    return [case for bucket in buckets for case in bucket]

def tagCasesNumbers(cases):
    caseNum = 0
    for case in cases:
        case['name'] = '%s-%s' % (str(caseNum).zfill(5), case['name'])
        caseNum += 1

def dumpResultCases(keys, cases, workingDir):
    results = [cases[t] for t in keys if t in cases]
    saveResults(results, workingDir)
    return results

def processTasks(testCases):
    testCases = reorderCases(testCases)
    tagCasesNumbers(testCases)

    buildResult = run(BUILDSCRIPT)
    if len(buildResult) > 0:
        print buildResult
        return []

    workingDir = createWorkingDirectory()

    jobsRun = 0
    tasks = []
    taskParamsDict = dict()
    for case in testCases[:MAXRUNNINGJOBS]:
        taskName = case['name']
        taskId = runTask(case, taskName, workingDir)
        case['id'] = taskId
        case['result_start_time'] = datetime.datetime.now().isoformat()
        tasks += [taskId]
        taskParamsDict[taskId] = case
        print "Run %s" % taskName
        time.sleep(3)
        jobsRun += 1

    results = dict()

    runningJobs = jobsRun
    while runningJobs > 0:
        runningJobs = 0
        sys.stderr.write("\x1b[2J\x1b[H")
        print 'Status:'

        for t in tasks:
            if not t in results or waitState(results[t]['result_state']):
                results[t] = taskTrace(t, taskParamsDict[t])
                if waitState(results[t]['result_state']):
                    runningJobs += 1
                else:
                    saveTaskStatus(t, taskParamsDict[t]['name'], workingDir)
                    readTimeFromOutput(results[t], workingDir)
            result = results[t]
            print '%s(x%s)\t: [%s:%s] w: %.3f\tc: %s\t%s' % (t, result['nodes'] * result['ppn'], result['result_state'], result['result_code'], result['result_time'], result['result_cpu'], result['name'])

        if runningJobs < MAXRUNNINGJOBS and len(tasks) < len(testCases):
            exJobsRun = jobsRun
            for case in testCases[exJobsRun:]:
                taskName = case['name']
                taskId = runTask(case, taskName, workingDir)
                case['id'] = taskId
                case['result_start_time'] = datetime.datetime.now().isoformat()
                tasks += [taskId]
                taskParamsDict[taskId] = case
                print "Run %s" % taskName
                time.sleep(3)
                runningJobs += 1
                jobsRun += 1
                if runningJobs >= MAXRUNNINGJOBS:
                    break

        if runningJobs > 0:
            print 'Progress: %d/%d' % (jobsRun, len(testCases))
            dumpResultCases(tasks, results, workingDir)
            time.sleep(SLEEPINTERVAL)

    return dumpResultCases(tasks, results, workingDir)

def printResults(results):
    zeroTime = -1
    zeroProcs = 0
    for result in results:
        time = result['result_time'] / result['repeats']
        procs = result['nodes'] * result['ppn']
        if zeroTime == -1:
            zeroTime = time
            zeroProcs = procs
        speedUp = 0
        eff = 0
        if zeroTime > 0:
            speedUp = float(time) / zeroTime
            eff = float(time) / zeroTime / (procs / zeroProcs)
        print "%s\t%s\t%s\t%s" % (procs, time, speedUp, eff)

#print [caseProcesses(case) for case in reorderCases(readTestCases('test_cases.json'))]
results = processTasks(readTestCases('test_cases.json'))
#printResults(results)
