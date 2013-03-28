"""
Tools for high level profiling and event logging
"""
import gc
import os
import inspect

global memHardlimit,memLast,memLog,logFile,logLevel,verbose,procID,logAllProcesses,flushBuffer, preInitBuffer

memHardlimit=0.0
memMax=0.0
memLast=0.0
memLog=False
memList=[]
logFile=None
logLevel=0
logAllProcesses=False
verbose=False
procID=None
flushBuffer=False
preInitBuffer=[]

def memProfOn():
    global memLog
    memLog = True

def memHardlimitOn(value):
    global memHardlimit
    memHardlimit = value*1024.0

#these are so the optparser can turn these on and off
def memProfOn_callback(option,opt,value,parser):
    memProfOn()

def memHardlimitOn_callback(option,opt,value,parser):
    memHardlimitOn(value)

def verboseOn_callback(option,opt,value,parser):
    global verbose
    verbose=True

def openLog(filename,level):
    global logFile
    global logLevel
    global preInitBuffer
    logFile=open(filename,'w')
    logLevel = level
    for (string,level,data) in preInitBuffer:
        logEvent(string,level,data)

def closeLog():
    global logFile
    try:
        logFile.close()
    except:
        pass

def logEvent(stringIn,level=1,data=None):
    global logLevel,procID,logAllProcesses,flushBuffer,preInitBuffer
    if procID != None:
        if logAllProcesses or procID==0:
            if level < logLevel:
                if logAllProcesses and procID != None and stringIn != None:
                    string = "Proc %d : " % procID
                    string += stringIn
                else:
                    string = stringIn
                if string!=None:
                    if data != None:
                        string += `data`
                    string+='\n'
                    global logFile,verbose
                    logFile.write(string)
                    if flushBuffer:
                        logFile.flush()
                    if verbose:
                        import sys
                        sys.stdout.write(string)
                        if flushBuffer:
                            sys.stdout.flush()
    elif procID==None:
        preInitBuffer.append((stringIn,level,data))

def memory(message=None,className='',memSaved=None):
    global memLog
    global memLast
    global memList
    global memMax
    if memLog:
        gc.collect()
        try:
            import memory_profiler
            mem = memory_usage(-1,interval=0.1,timeout=0.1)[0]
        except:
            try:
                mem = int(os.popen('ps -p %d -o rss|tail -1' % os.getpid()).read())/1024.0
            except:
                try:
                    #CLE doesn't have full ps functionality
                    #try to find first occurence of pid, then get memory usage slot
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                    plist = os.popen('ps ').read()
                    pr = plist.split('%s' % os.getpid())
                    #print 'pr = %s ' % pr
                    #print 'pr[1] = %s ' % pr[1]
                    mem = int(pr[1].split()[1])/1024.0
                except:
                    logEvent("memory function doesn't work on this platform\n")
                    mem = 0
        if mem > memMax:
            memMax = mem
        if memSaved != None:
            memInc = mem - memSaved
        else:
            memInc = mem-memLast
        memList.append((message,memInc))
        memLast = mem
        caller = inspect.stack()[1][3]
        line = inspect.stack()[1][2]
        filename = inspect.stack()[1][1]
        if className:
            className += '.'
        if memHardlimit:
            if mem > memHardlimit:
                import sys
                if className == None:
                    className = "UnknownClass"
                if caller == None:
                    caller="UnknownCaller"
                if line == None:
                    line = "Unknown Line"
                if message == None:
                    message = ''
                print filename.split("/")[-1],className,caller,`line`,message,memInc,mem,memHardlimit
                sys.exit("ERROR: MEMORY HARDLIMIT REACHED, EXIT From "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %d MB in routine, %d MB in program, %d MB is Hardlimit" % (memInc,mem,memHardlimit))
        if message:
            return "In "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %d MB in routine, %d MB in program" % (memInc,mem)

def memorySummary():
    global memLog
    global memList
    global memMax
    if memLog:
        gc.collect()
        try:
            import memory_profiler
            mem = memory_usage(-1,interval=0.1,timeout=0.1)[0]
        except:
            try:
                mem = int(os.popen('ps -p %d -o vsz|tail -1' % os.getpid()).read())/1024.0
            except:
                try:
                    #CLE doesn't have full ps functionality
                    #try to find first occurence of pid, then get memory usage slot
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                    plist = os.popen('ps ').read()
                    pr = plist.split('%s' % os.getpid())
                    #print 'pr = %s ' % pr
                    #print 'pr[1] = %s ' % pr[1]
                    mem = int(pr[1].split()[1])/1024.0
                except:
                    logEvent("memory function doesn't work on this platform\n")
                    mem = 0
        if mem > 0:
            for pair in memList:
                logEvent(`pair[0]`+"  %"+`100.0*pair[1]/memMax`)
