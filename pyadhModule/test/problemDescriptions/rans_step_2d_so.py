from pyadh.default_so import *
from rans_step_2d import *
if useSmagorinsky:
    pnList = [("rans_step_2d_p","rans_step_2d_c0p1c0p1_n")]
else:
    #pnList = [("rans_step_2d_p","rans_step_2d_c0p1c0p1_n"),
    #          ("kEpsilon_step_2d_p","kEpsilon_step_2d_dgp0_n")]
    pnList = [("rans_step_2d_p","rans_step_2d_c0p1c0p1_n"),
              ("kEpsilon_k_step_2d_p","kEpsilon_k_step_2d_c0p1_n"),
              ("kEpsilon_epsilon_step_2d_p","kEpsilon_epsilon_step_2d_c0p1_n")]
#pnList = [("rans_step_2d_p","rans_step_2d_c0p1c0p1_n"),
#          ("kEpsilon_step_2d_p","kEpsilon_step_2d_c0p1_n")]
#modelSpinUpList = [0]
name = "rans_step_2d"
useOneArchive = True
if useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL  = True
needEBQ = True

nDTout = 10
DT = T/nDTout #None#T/100.0
tnList = [0,0.01]
tnList.extend([i*DT for i  in range(1,nDTout+1)])
