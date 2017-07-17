"""
Tools for high level profiling and event logging

.. inheritance-diagram:: proteus.Profiling
   :parts: 1
"""
import gc
import inspect
import pstats
from time import time
import atexit

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
    assert procID is not None, "Initialize Comm and set Profiling.procID before opening log"
    filename_full = filename
    import os
    if logLocation is not None:
        logDir = logLocation
        filename_full = os.path.join(logDir,filename)
    if  procID == 0:
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
    if procID is not None:
        if logAllProcesses or procID==0:
            if level < logLevel:
                if logAllProcesses and procID is not None and stringIn is not None:
                    string = "Proc %d : " % procID
                    string += stringIn
                else:
                    string = stringIn
                if string is not None:
                    if data is not None:
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
    elif procID is None:
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
        if memSaved is not None:
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
                if className is None:
                    className = "UnknownClass"
                if caller is None:
                    caller="UnknownCaller"
                if line is None:
                    line = "Unknown Line"
                if message is None:
                    message = ''
                logEvent("PROTEUS ERROR: MEMORY HARDLIMIT REACHED, EXIT From "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %f MB in routine, %f MB in program, %f MB is hard limit" % (memInc,mem,memHardLimit))
                MPI.COMM_WORLD.Abort(1)
                sys.exit("MPI.COMM_WORLD.Abort(1); exit(1)")
        if message:
            return "In "+filename.split("/")[-1]+", "+className+caller+", line "+`line`+": "+message+", %f MB in routine, %f MB in program" % (memInc,mem)

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
        comm.barrier()#ensure files are ready for master
        if comm.isMaster():
            import copy
            import StringIO
            profilingLog = StringIO.StringIO()
            stats = pstats.Stats(profile_rank_name, stream=profilingLog)
            stats.__dict__['files']=['Maximum times across MPI tasks for',
                                     stats.__dict__['files'][0]]
            statsm = stats.stats
            for i in range(1,comm.size()):
                pstatsi = pstats.Stats(profile_name+str(i))
                statsi = pstatsi.stats
                stats.__dict__['files'].append(pstatsi.__dict__['files'][0])
                for f,c in statsi.iteritems():
                    if f in statsm:
                        if c[2] > statsm[f][2]:
                            statsm[f] = c
                    else:
                        statsm[f] = c
            stats.sort_stats('cumulative')
            stats.print_stats(30)
            stats.sort_stats('time')
            stats.print_stats(30)
            logEvent(profilingLog.getvalue())
            msg = r"""
Wall clock percentage of top 20 calls
-------------------------------------
"""
            total=0.0
            for f in stats.__dict__['fcn_list'][0:20]:
                if f[0] == '~':
                    fname=f[-1].strip("<").strip(">")
                else:
                    fname="function '{2:s}' at {0:s}:{1:d}".format(*f)
                msg+=("{0:11.1%} {1:s}\n".format(statsm[f][2]/stats.__dict__['total_tt'],str(fname)))
                total += statsm[f][2]/stats.__dict__['total_tt']
            logEvent(msg)
            logEvent("Representing "+`total*100.`+"%")
        return func_return

@atexit.register
def ProfilingDtor():
    global procID, verbose
    if procID is None:
        verbose=True
        logEvent(
            "Proteus.Profiling never initialized. Doing it at exit.")
        procID = 0
        openLog("proteus_default.log",level=11,logLocation=".")
    closeLog()


def print_petsc_commandline_options(petsc_options):
    """ Returns a formated string with PETSc command-line options 
    
    petsc_options : list
        A list of the Petsc command line options and values.
    """
    str = ""
    for i,j in zip(petsc_options[1::2], petsc_options[2::2]):
        str += "NAHeader PETScOptions " + i + " " + j +  "\n"
    return str
    
