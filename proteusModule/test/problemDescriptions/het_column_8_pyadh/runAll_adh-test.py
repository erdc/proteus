from glob import *
import os
from shutil import *
from pyadh.cadh import *
from pyadh import Comm
from pyadh import ADH
import sys
print sys.argv
comm = Comm.init(sys.argv)
log = Profiling.logEvent

bcFiles = glob(os.getenv('PYADH_PACKAGES')+"/adh-test/*/*.bc")
# meshFiles = glob(os.getenv('PYADH_PACKAGES')+"/adh-test/*/*.*dm")
# hotFiles = glob(os.getenv('PYADH_PACKAGES')+"/adh-test/*/*.hot")
# for bc,mesh,hot in zip(bcFiles,meshFiles,hotFiles):
#     bcLocal = bc.split("/")[-1]
#     meshLocal = mesh.split("/")[-1]
#     hotLocal = hot.split("/")[-1]
#     copyfile(bc,bcLocal)
#     copyfile(mesh,meshLocal)
#     copyfile(hot,hotLocal)
#     test_name = bcLocal[:-3]
#     print bcLocal,test_name
#     cadhRun(test_name,"level0_test_run_"+test_name)
models = []
for bc in bcFiles+bcFiles[::-1]+bcFiles+bcFiles[::-1]:
    testName = bc.split("/")[-1][:-3]
    if testName in ['angle_sup','first_flume','ns_cube_snd','pool8','riprap2d']:
        continue#skip this broken test problem
    testDir = bc[:-(3+len(testName))]
    os.chdir(testDir)
    model = ADH.ADH_OneLevelTransport(adhInput=testName,runname="test1")
    #del model.adh_ns
    #del model.adh_transport
    model.__del__()
    #cadhRun(testName,"level0_test_run_"+testName)
    #model = ADH.ADH_NumericalSolution(adhInput=testName,runname="test1",comm=comm)
    #model.calculateSolution("stringNotYetUsed")
