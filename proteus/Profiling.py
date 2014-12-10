"""
Tools for high level profiling and event logging
"""
import gc
import inspect
import pstats
from time import time

try:
    from cProfile import Profile
except:
    from profile import Profile

global memHardLimit,memLast,memLog,logFile,logLevel,verbose,procID,logAllProcesses,flushBuffer, preInitBuffer

memHardLimit=0.0
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
logDir = '.'

startTime = time()

def memProfOn():
    global memLog
    memLog = True

def memHardLimitOn(value):
    global memHardLimit
    if not memLog:
        memProfOn()
    memHardLimit = value*1024.0

#these are so the optparser can turn these on and off
def memProfOn_callback(option,opt,value,parser):
    memProfOn()

def memHardLimitOn_callback(option,opt,value,parser):
    memHardLimitOn(value)

def verboseOn_callback(option,opt,value,parser):
    global verbose
    verbose=True

def openLog(filename,level,logLocation=None):
    global logFile
    global logLevel
    global preInitBuffer
    global logDir
    global procID
    global logAllProcesses
    filename_full = filename
    import os
    if logLocation != None:
        logDir = logLocation
        filename_full = os.path.join(logDir,filename)
    if  procID == None or procID == 0:
        logFile=open(filename_full,'w')
    elif logAllProcesses:
        logFile=open(filename_full+`procID`,'w')
    logLevel = level
    for string,level,data in preInitBuffer:
        logEvent(string,level,data)

def closeLog():
    global logFile
    try:
        logFile.close()
    except:
        pass

def logEvent(stringIn, level=1, data=None):
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
                    if data!=None:
                        string += repr(data)
                    string +='\n'
                    string = ("[%8d] " % (time() - startTime)) + string
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
        from memory_profiler import memory_usage
        memList = memory_usage(-1)
        if len(memList) < 1:
            memList = memory_usage(-1,timeout=10)
            assert(len(memList) > 0)
        mem = memList[-1]
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
        if memHardLimit:
            if mem > memHardLimit:
                import sys
                from mpi4py import MPI
                if className == None:
                    className = "UnknownClass"
                if caller == None:
                    caller="UnknownCaller"
                if line == None:
                    line = "Unknown Line"
                if message == None:
                    message = ''
                logEvent("PROTEUS ERROR: MEMORY HARDLIMIT REACHED, EXIT From "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %d MB in routine, %d MB in program, %d MB is hard limit" % (memInc,mem,memHardLimit))
                MPI.COMM_WORLD.Abort(1)
                sys.exit("MPI.COMM_WORLD.Abort(1); exit(1)")
        if message:
            return "In "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %d MB in routine, %d MB in program" % (memInc,mem)

def memorySummary():
    global memLog
    global memList
    global memMax
    if memLog:
        gc.collect()
        from memory_profile import memory_usage
        mem = memory_usage(-1)[0]
        if mem > 0:
            for pair in memList:
                logEvent(`pair[0]`+"  %"+`100.0*pair[1]/memMax`)


class Dispatcher():
    """
    Profiles function calls.  Must be enabled like so:

    dispatch = Dispatcher(comm, True)

    A non-profiling dispatcher can be created like:

    dispatch = Dispatcher()

    Then, you can execute function calls using the following syntax:

    ret = dispatch(func, args, kwargs)
    """

    def __init__(self, comm=None, on=False):
        self.comm = comm
        self.on = on

    def __call__(self, func, func_args, func_kwargs, profile_name=''):
        if not self.on:
            return func(*func_args, **func_kwargs)
        else:
            return self.profile_function(func, func_args, func_kwargs, profile_name)

    def profile_function(self, func, args, kwargs, profile_name):
        """
        Profile a Proteus function using the call: func(*args, **kwargs)

        Returns the output of the function call.
        """

        comm = self.comm
        if type(comm) is None:
            raise ValueError("The Dispatcher does not have a valid Comm object")

        prof = Profile()
        func_return = prof.runcall(func, *args, **kwargs)
        profile_rank_name = profile_name + str(comm.rank())
        stripped_profile_name = profile_name + '_c' + str(comm.rank())

        prof.dump_stats(profile_rank_name)
        comm.beginSequential()
        stats = pstats.Stats(profile_rank_name)
        stats.strip_dirs()
        stats.dump_stats(stripped_profile_name)
        stats.sort_stats('cumulative')
        if verbose and comm.isMaster():
            stats.print_stats(30)
            stats.sort_stats('time')
            if verbose and comm.isMaster():
                stats.print_stats(30)
        comm.endSequential()

        return func_return
