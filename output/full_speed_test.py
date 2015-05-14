import sys
import subprocess
import re
import time

RUNSCRIPT = 'debug.run'
SLEEPINTERVAL = 10
lineRegex = re.compile('([^ ]+) = (.*)')

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

def taskStatus(taskId):
    data = monitorTask(taskId)
    jobState = value_jobstate_none
    wallTime = 0
    cpuTime = 0
    nodes = ''
    exitCode = ''

    if key_jobstate in data:
        jobState = data[key_jobstate]
        nodes = data[key_nodes]
        if jobState == value_jobstate_done:
            wallTime = parseTime(data[key_walltime])
            cpuTime = parseTime(data[key_cputime])
            exitCode = data[key_exit_status]
        
    return ( jobState, wallTime, cpuTime, nodes, exitCode )

def waitState(state):
    return state == value_jobstate_queued or state == value_jobstate_running

def procCount(nodes):
    splits = [int(s.strip('=')) for s in nodes.split(':ppn')]
    return splits[0] * splits[1]

def processTasks(tasksParams):
    tasks = [runTask(nodes, ppn) for (nodes, ppn) in tasksParams]

    results = dict()
    anyRuning = True

    while anyRuning:
        anyRuning = False
        sys.stderr.write("\x1b[2J\x1b[H")
        print 'Status:'
        for t in tasks:
            if not t in results or waitState(results[t][0]):
                results[t] = taskStatus(t)
                anyRuning = anyRuning or waitState(results[t][0])
            print '%s(x%s)\t: [%s:%s] w: %s\tc: %s\tnodes=%s' % (t, procCount(results[t][3]), results[t][0], results[t][4], results[t][1], results[t][2], results[t][3])
        if anyRuning:
            time.sleep(SLEEPINTERVAL)

    return [(results[t][1], procCount(results[t][3])) for t in tasks]

#print runTask(4, 4)
#print taskStatus('186314')

#times = processTasks([(1, 1), (1, 2), (1, 4), (1, 8), (2, 6), (2, 8), (3, 6), (3, 8)])
# 6 8 10 12 16 20 24 28 30 36 42
times = processTasks([(2, 3), (2, 4), (2, 5), (3, 4), (4, 4), (4, 5), (4, 6), (4, 7), (5, 6), (6, 6), (6, 7)])
for (time, proc) in times:
    print "%s\t%s" % (proc, time)
