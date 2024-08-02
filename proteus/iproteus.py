#!/usr/bin/env python

##
# \addtogroup scripts
# Simulation scripts and other tools for pre and postprocessing files
#
#
#
# \file parun
# @{
#   \ingroup scripts
#   \brief A driver for both single and mult-model simulations
#
#
from importlib import reload
import os
import proteus

#remove blanket import statements until after Comm initialized for petsc4py
from proteus import Profiling, Comm
from warnings import *
import optparse
import sys
try:
    from IPython.display import display,Markdown
except:
    pass
import inspect

def display_src(python_object):
    display(Markdown("```Python\n"+inspect.getsource(python_object)+"\n```\n"))

def load_macros(filename,verbose=False):
    with open(filename,'r') as macros:
        markdown = ""
        for line in macros.readlines():
            if len(line.strip('\n')) and not line.startswith('%'):
                markdown += "${0}$\n".format(line.strip('\n'))
        if verbose:
            print(markdown)
        display(Markdown(markdown))

usage = "usage: %prog [options] soModule.py [[soFile.sso] [pModule.py nModule.py]]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-I", "--inspect",
                  help="Inspect namespace at 't0','user_step'",
                  action="store",
                  dest="inspect",
                  default='')
parser.add_option("-i", "--interactive",
                  help="Read input from stdin",
                  action="store_true",
                  dest="interactive",
                  default='')
parser.add_option("-d", "--debug",
                  help="start the python debugger",
                  action="store_true",
                  dest="debug",
                  default=False)
parser.add_option("-V", "--viewer",
                  help="Set the method to use for runtime viewing. Can be vtk or gnuplot",
                  action="store",
                  type="string",
                  dest="viewer",
                  default=False)
parser.add_option("-C", "--plot-coefficients",
                  help="Plot the coefficients of the transport models",
                  action="store_true",
                  dest="plotCoefficients",
                  default=False)
parser.add_option("-P", "--petsc-options",
                  help="Options to pass to PETSc",
                  action="store",
                  type="string",
                  dest="petscOptions",
                  default=None)
parser.add_option("-O", "--petsc-options-file",
                  help="Text file of options to pass to PETSc",
                  action="store",
                  type="string",
                  dest="petscOptionsFile",
                  default=None)
parser.add_option("-D", "--dataDir",
                  help="Data directory",
                  action="store",
                  type="string",
                  dest="dataDir",
                  default='')
parser.add_option("-b", "--batchFile",
                  help="Text file of commands to execute",
                  action="store",
                  type="string",
                  dest="batchFileName",
                  default="")
parser.add_option("-p", "--profile",
                  help="Generate a profile of the run",
                  action="store_true",
                  dest="profile",
                  default=False)
parser.add_option("-T", "--useTextArchive",
                  help="Archive data in ASCII text files",
                  action="store_true",
                  dest="useTextArchive",
                  default=False)
parser.add_option("-m", "--memory",
                  help="Track memory usage of the run",
                  action="callback",
                  callback=Profiling.memProfOn_callback)
parser.add_option("-M", "--memoryHardLimit",
                  help="Abort program if you reach the per-MPI-process memory hardlimit (in GB)",
                  action="callback",
                  type="float",
                  callback=Profiling.memHardLimitOn_callback,
                  default = -1.0,
                  dest = "memHardLimit")
parser.add_option("-l", "--log",
                  help="Store information about what the code is doing, 0=none, 10=everything",
                  action="store",
                  type="int",
                  dest="logLevel",
                  default=1)
parser.add_option("-A", "--logAllProcesses",
                  help="Log events from every MPI process",
                  action="store_true",
                  dest="logAllProcesses",
                  default=False)
parser.add_option("-v", "--verbose",
                  help="Print logging information to standard out",
                  action="callback",
                  callback=Profiling.verboseOn_callback)
parser.add_option("-E", "--ensight",
                  help="write data in ensight format",
                  action="store_true",
                  dest="ensight",
                  default=False)
parser.add_option("-L", "--viewLevels",
                  help="view solution on every level",
                  action="store_true",
                  dest="viewLevels",
                  default=False)
parser.add_option("--viewMesh",
                  help="view mesh",
                  action="store_true",
                  dest="viewMesh",
                  default=False)
parser.add_option("-w", "--wait",
                  help="stop after each nonlinear solver call",
                  action="store_true",
                  dest="wait",
                  default=False)
parser.add_option('--probDir',
                  default='.',
                  help="""where to find problem descriptions""")
parser.add_option("-c","--cacheArchive",
                  default=False,
                  dest="cacheArchive",
                  action="store_true",
                  help="""don't flush the data files after each save, (fast but may leave data unreadable)""")
parser.add_option("-G","--gatherArchive",
                  default=False,
                  dest="gatherArchive",
                  action="store_true",
                  help="""collect data files into single file at end of simulation (convenient but slow on big run)""")
parser.add_option("-s","--subdomainArchives",
                  default=False,
                  dest="subdomainArchives",
                  action="store_true",
                  help="""write xmf files for every processor's subdomain""")
parser.add_option("-n","--no_global_sync",
                  default=True,
                  dest="global_sync",
                  action="store_false",
                  help="""don't use a single hdf5 archive""")
parser.add_option("-H","--hotStart",
                  default=False,
                  dest="hotStart",
                  action="store_true",
                  help="""Use the last step in the archive as the initial condition and continue appending to the
                  archive""")
parser.add_option("-B","--writeVelocityPostProcessor",
                  default=False,
                  dest="writeVPP",
                  action="store_true",
                  help="""Use the last step in the archive as the initial condition and continue appending to the
                  archive""")
parser.add_option("-F","--generatePartitionedMeshFromFiles",
                  default=False,
                  dest="generatePartitionedMeshFromFiles",
                  action="store_true",
                  help="""generate a parallel mesh directly from files (saves memory overhead of storing global mesh)""")

parser.add_option("-S","--save_dof",
                  default=False,
                  dest="save_dof",
                  action="store_true",
                  help="""save the solution degrees of freedom after they have been calculated (for interactive plotting)""")
parser.add_option("--setupOnly",
                  default=False,
                  dest="setupOnly",
                  action="store_true",
                  help="Don't actually compute a solution, just initialize the Numerical Solution")


(opts,args) = parser.parse_args(args=[])

if opts.debug:
    import pdb
    pdb.set_trace()

log = Profiling.logEvent
logDir = None
if opts.dataDir != '':
    if not os.path.exists(opts.dataDir):
        try:
             os.mkdir(opts.dataDir)
        except os.error:
             print("Failed to create: " + opts.dataDir)
             opts.dataDir=''
    logDir = opts.dataDir

log("Initializing Proteus")

log("Initializing MPI")
if opts.petscOptions is not None:
    petsc_argv = sys.argv[:1]+opts.petscOptions.split()
    log("PETSc options from commandline")
    log(str(petsc_argv))
else:
    petsc_argv=sys.argv[:1]
if opts.petscOptionsFile is not None:

    petsc_argv=[sys.argv[0]]
    with open(opts.petscOptionsFile) as petsc_file:
        data = petsc_file.readlines()
    def strip_comments(line):
        if '#' in line:
            line = line[:line.index('#')]
        return line

    stripped_data = [strip_comments(line) for line in data]
    petsc_argv += '\n'.join(stripped_data).split()
    log("PETSc options from commandline")
    log(str(petsc_argv))

Comm.argv = petsc_argv
comm = Comm.init()
Profiling.procID=comm.rank()
if opts.logLevel > 0:
    #mwf looks like just gets .py.log if no sso file given
    #Profiling.openLog(args[0][-3:]+".log",opts.logLevel)
    if len(args) == 0:
        Profiling.openLog("proteus.log",opts.logLevel)
    elif len(args[0].split(' ')) == 1:
        Profiling.openLog(args[0].split('.')[0]+".log",opts.logLevel)
    else:
        Profiling.openLog(args[0][-3:]+".log",opts.logLevel)
if opts.logAllProcesses:
    Profiling.logAllProcesses = True

#blanket import statements can go below here now that petsc4py should be initialized
from proteus import *
try:
    from proteusGraphical import *
except:
    pass

log("Adding "+str(opts.probDir)+" to path for loading modules")
probDir = str(opts.probDir)
if probDir not in sys.path:
    sys.path.insert(0,probDir)

log("Loading problem-specific modules")
pList = []
nList = []
sList = []

if __name__ == '__main__':
    log("Loading problem-specific modules")
    if len(args) < 1:
        raise RuntimeError("No input file specified")
    if len(args) == 1:#arg should be an sso file or an so module
        if args[0][-3:] == '.py':
            log("Loading so module = "+args[0])
            so = __import__(args[0][:-3])
            if so.name is None:
                so.name = args[0][:-3]
            for (pModule,nModule) in so.pnList:
                log("Loading p module = "+pModule)
                pList.append(__import__(pModule))
                if pList[-1].name is None:
                    pList[-1].name = pModule
                log("Loading n module = "+nModule)
                nList.append(__import__(nModule))
            #
            if so.sList == []:
                for i in range(len(so.pnList)):
                    s = default_s
                    sList.append(s)
            else:
                sList = so.sList
        else: #so text file is provided; will be deprecated
            so = default_so
            so.name = args[0][0:-4]
            log("Loading old style sso text file = "+args[0])
            soFile = open(args[0],'r')
            lines = soFile.readlines()
            for mLine in lines:
                filenames = mLine.split()
                #read the so file as text for now
                if filenames[0] == "systemStepControllerType":
                    so.systemStepControllerType = eval(filenames[2])
                elif filenames[0] == "useOneMesh":
                    so.useOneMesh = eval(filenames[1])
                elif filenames[0] == "name":
                    so.name = filenames[1]
                elif filenames[0] == "needEBQ_GLOBAL":
                    so.needEBQ_GLOBAL = filenames[1]
                elif filenames[0] == "needEBQ":
                    so.needEBQ = filenames[1]
                elif len(filenames) == 2:
                    log("Loading p module = "+filenames[0][:-3])
                    pList.append(__import__(filenames[0][:-3]))
                    log("Loading n module = "+filenames[1][:-3])
                    nList.append(__import__(filenames[1][:-3]))
                    if pList[-1].name is None:
                        pList[-1].name = filenames[0][:-3]
    elif len(args) == 2: #p and n modules provided
        so=default_so
        s = default_s
        so.pnList=[(args[0][:-3],args[1][:-3])]
        log("Using default so module")
        log("Loading p module = "+args[0][:-3])
        pList.append(__import__(args[0][:-3]))
        if pList[-1].name is None:
            pList[-1].name = args[0][:-5]
        so.name = pList[-1].name
        log("Loading n module = "+args[1][:-3])
        nList.append(__import__(args[1][:-3]))
        sList.append(s)
        try:
            so.systemStepControllerType = nList[0].systemStepControllerType
        except:
            pass
        try:
            so.tnList = nList[0].tnList
            so.archiveFlag = nList[0].archiveFlag
        except:
            pass
    else:
        raise RuntimeError("Too many input files specified")
    if opts.batchFileName != "":
        log("Attaching input stream to batch file = "+opts.batchFileName)
        inputStream = open(opts.batchFileName,'r')
    else:
        log("Attaching input stream to stdin")
        inputStream = sys.stdin
    simFlagsList = []
    for i in range(len(pList)):
        simFlagsList.append({})
    if sList[0].logAllProcesses or opts.logAllProcesses:
        Profiling.logAllProcesses = True
    Profiling.flushBuffer=sList[0].flushBuffer
    log("Running Proteus version "+proteus.__version__,level=0)
    if opts.viewer:
        log("Starting viewer")
        Viewers.viewerOn(so.name+repr(comm.rank()),opts.viewer)
    log("Setting simulation processing defaults")
    for p,simFlags in zip(pList,simFlagsList):
        pNameProc = p.name + repr(comm.rank())
        simFlags['simulationName'] = p.name
        simFlags['simulationNameProc'] = pNameProc
        simFlags['dataFile']       = pNameProc+'.dat'
        simFlags['components'] = [ci for ci in range(p.coefficients.nc)]
        if opts.viewer != False or opts.ensight:
            simFlags['plotTimes'] = ['All'] #plot solution at each macro time level (also 'Last')
            simFlags['plotQuantities'] = ['u']  #plot soln and exact soln if exists
            p.coefficients.plotCoefficients = opts.plotCoefficients
        if opts.ensight:
            if 'plotOptions' in simFlags:
                if 'ensight' in simFlags['plotOptions']:
                    simFlags['plotOptions']['ensight']['on'] = True
                else:
                    simFlags['plotOptions']['ensight'] = {'on':True}
            else:
                simFlags['plotOptions'] = {'ensight':{'on':True}}
            #control convention on case file for multiple models, doesn't effect parun right now
            simFlags['plotOptions']['ensight']['caseFileName']=simFlags['simulationName']+'master'+ repr(comm.rank())
    if opts.batchFileName != "":
        log("Reading batch file = ",opts.batchFileName)
        #try to split batch file into executable blocks delimited by
        #s or start on line by itself
        #terminate batch file with a q or quit on a line by itself
        batchfile = inputStream.readlines()
        batchBlocks = []; curBlock = 0
        batchBlocks.append([])
        curline = 0; done = False
        while curline < len(batchfile) and not done:
            line = batchfile[curline]
            newBlock = False
            if not line.isspace():
                newBlock = (line.split()[0] == 's' or
                            line.split()[0] == 'start')
            if newBlock:
                batchBlocks[-1].append('#end of block\n')
                batchBlocks.append([])
            if not line.isspace():
                done = (line.split()[0] == 'q' or
                        line.split()[0] == 'quit')
            if not done and not newBlock:
                batchBlocks[-1].append(line)
            curline += 1
        #always pop the last block so have to have a start at the end
        batchBlocks.pop()
        log("Batch file contains %i runs" % (len(batchBlocks),))

    running = True
    runNumber=0
    while running:
        if opts.interactive:
            userInput = True
            sys.stdout.write("Enter python commands or (s)tart/(q)quit\n>>>")
        elif opts.batchFileName != "":
            userInput = False
            lines = '\n'.join(batchBlocks.pop(0))
            exec(lines)
            run     = True
            running = len(batchBlocks) > 0
        else:
            userInput = False
            run = True
            running = False
        while (userInput):
            line = inputStream.readline()
            if line:
                if (line.split()[0] == 's' or
                    line.split()[0] == 'start'):
                    userInput = False
                    run = True
                elif (line.split()[0] == 'q' or
                      line.split()[0] == 'quit'):
                    userInput = False
                    run = False
                    running = False
                else:
                    userInput = True
                    exec(line)
                    sys.stdout.write(">>>")
            else:
                userInput = False
                run = False
                running = False
        if run:
            log("Starting %s run number %i" %(so.name,runNumber))
            runName = so.name+repr(runNumber)
            if opts.profile:
                profiler.run('ns = NumericalSolution.NS_base(so,pList,nList,sList,opts,simFlagsList)',so.name+'_init_prof'+repr(comm.rank()))
                stats = pstats.Stats(so.name+'_init_prof'+repr(comm.rank()))
                stats.strip_dirs()
                stats.dump_stats(so.name+'_init_prof_c'+repr(comm.rank()))
                stats.sort_stats('cumulative')
                if Profiling.verbose:
                    stats.print_stats(30)
                stats.sort_stats('time')
                if Profiling.verbose:
                    stats.print_stats(30)
                profiler.run('ns.calculateSolution(runName)',so.name+'_run_prof'+repr(comm.rank()))
                runNumber+=1
                stats = pstats.Stats(so.name+'_run_prof'+repr(comm.rank()))
                stats.strip_dirs()
                stats.dump_stats(so.name+'_run_prof_c'+repr(comm.rank()))
                stats.sort_stats('cumulative')
                if Profiling.verbose:
                    stats.print_stats(30)
                stats.sort_stats('time')
                if Profiling.verbose:
                    stats.print_stats(30)
            else:
                ns = NumericalSolution.NS_base(so,pList,nList,sList,opts,simFlagsList)
                ns.calculateSolution(runName)
            log("Completed %s run number %i" %(so.name,runNumber))
            import gc
            log("Collecting garbage")
            gc.collect()
            log("Done with garbage")
            if comm.isMaster():
                Profiling.closeLog()
    if opts.viewer:
        if comm.isMaster():
            try:
                if vtkViewers.hasQt:
                    vtkViewers.g.app.exec_()
                for w in list(vtkViewers.windowDict.values()):
                    w.compManager.StopServices()
                sys.exit()
            except:
                pass
        else:
            try:
                for w in list(vtkViewers.windowDict.values()):
                    w.compManager.StartServices()
                sys.exit()
            except:
                pass
        if 'viewerType' in dir(Viewers) and Viewers.viewerType == 'matlab':
            Viewers.viewerPipe.write("quit \n")
