from proteus.default_so import *
from miscible_displacement import *
pnList = [("md_flow_p","md_flow_ncp1_n"),
          ("md_trans_p","md_trans_c0p1_n")]
          
modelSpinUpList=[0]
name = "md_example_ncp1_visc_a_5"
systemStepControllerType = Sequential_MinAdaptiveModelStep
#for numerical fluxes
needEBQ_GLOBAL  = True
needEBQ = True

tnList = [i*T/float(nDTout) for i  in range(nDTout+1)]
