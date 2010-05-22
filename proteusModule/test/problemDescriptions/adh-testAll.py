#!/usr/bin/env python
"""
A script for running the cADH test set through the various interfaces with various viewing and archiving options.
"""
from glob import *
import os,subprocess,optparse,shutil
import optparse
from pyadh import Profiling,Viewers,ADH
log = Profiling.logEvent

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-l", "--level",
                  help="which interface to test",
                  type="int",
                  action="store",
                  dest="level",
                  default=0)
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
parser.add_option("-p", "--profile",
                  help="Generate a profile of the  run",
                  action="store_true",
                  dest="profile",
                  default=False)
parser.add_option("-T", "--useTextArchive",
                  help="Archive data in ASCII text files",
                  action="store_true",
                  dest="useTextArchive",
                  default=False)
parser.add_option("--viewMesh",
                  help="view mesh",
                  action="store_true",
                  dest="viewMesh",
                  default=False)
parser.add_option("-w", "--wait",
                  help="stop after each model run",
                  action="store_true",
                  dest="wait",
                  default=False)
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


if opts.debug:
    import pdb
    pdb.set_trace()

if opts.level > 0:
    import sys
    if opts.petscOptions != None:
        petsc_argv = sys.argv[:1]+opts.petscOptions.split()
        log("PETSc options from commandline")
        log(str(petsc_argv))
    else:
        petsc_argv=sys.argv[:1]
    if opts.petscOptionsFile != None:
        petsc_argv=[sys.argv[0]]
        petsc_argv += open(opts.petscOptionsFile).read().split()
        log("PETSc options from commandline")
        log(str(petsc_argv))

    from pyadh import Comm
    comm = Comm.init(argv=petsc_argv)

if opts.viewer:
    log("Starting viewer")
    Viewers.viewerOn("cADH"+`comm.rank()`,opts.viewer)

bcFiles = glob(os.getenv('PYADH_PACKAGES')+"/adh-test/*/*.bc")
for bc in bcFiles:
    testName = bc.split("/")[-1][:-3]
    if testName != 'het_column_8':#in ['angle_sup','first_flume','ns_cube','pool8']:
        continue#skip this test
    testDir = bc[:-(3+len(testName))]
    os.chdir(testDir)
    print "Changed to directory ",os.getcwd()
    if opts.level == 0:
        try:
            ADH.cadhRun(testName,"cadhRunTest_"+testName)
        except subprocess.CalledProcessError:
            print "cadh failed on "+`testName`+" continuing"
            continue
    if opts.level == 1:
        adh_ns = ADH.ADH_NumericalSolution(adhInput=testName,runname="ADH_NumericalSolution_test"+testName,opts=opts,petscMatrix=False)
        adh_ns.calculateSolution()
        adh_ns.__del__()
    if opts.level == 2:
        adh_transport = ADH.ADH_OneLevelTransport(adhInput=testName,runname="ADH_OneLevelTransport_test"+testName,options=opts,petscMatrix=True)
        adh_transport.__del__()
