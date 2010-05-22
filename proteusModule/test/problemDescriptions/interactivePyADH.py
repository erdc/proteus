#!/usr/bin/env sage -python
"""
collect basic setup phase commands like import statements and
option setting to be able to calculate numerical solution in 
sage interactive session
"""
import numpy
from pyadh import *
import optparse
import sys


def init():
    comm = Comm.init()
    return comm

def getDefaultOptions():
    usage = "usage: %prog [options] soModule.py [[soFile.sso] [pModule.py nModule.py]]"
    parser = optparse.OptionParser(usage=usage)
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
                      help="Text file of ptions to pass to PETSc",
                      action="store",
                      type="string",
                      dest="petscOptionsFile",
                      default=None)
    parser.add_option("-D", "--dataDir",
                      help="Options to pass to PETSc",
                      action="store",
                      type="string",
                      dest="dataDir",
                      default='')
    parser.add_option("-b", "--batchFile",
                      help="Read input from a file",
                      action="store",
                      type="string",
                      dest="batchFileName",
                      default="")
    parser.add_option("-p", "--profile",
                      help="Generate a profile of the  run",
                      action="store_true",
                      dest="profile",
                      default=True)
    parser.add_option("-T", "--useTextArchive",
                      help="Archive data in ASCII text files",
                      action="store_true",
                      dest="useTextArchive",
                      default=False)
    parser.add_option("-m", "--memory",
                      help="Track memory usage of the  run",
                      action="callback",
                      callback=Profiling.memProfOn_callback)
    parser.add_option("-l", "--log",
                      help="Store information about what the code is doing,0=none,10=everything",
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

    parser.add_option("-H","--hotStart",
                      default=False,
                      dest="hotStart",
                      action="store_true",
                      help="""Use the last step in the archive as the intial condition and continue appending to the archive""")
    (opts,args) = parser.parse_args()
    return opts

def getDefaultDescriptions(pFile=None,nFile=None):
    if pFile != None:
        pList = [__import__(pFile)]
    else:
        pList = [default_p] 
    if nFile != None:
        nList = [__import__(nFile)]
    else:
        nList = [default_n]; 
    sList = [default_s]
    so = default_so
    if pList[-1].name == None:
        pList[-1].name = 'sage-interactive'
    so.name = pList[-1].name
    return (so,pList,nList,sList)

########################################################################
comm = init()
opts = getDefaultOptions()
(so,pList,nList,sList) = getDefaultDescriptions()
