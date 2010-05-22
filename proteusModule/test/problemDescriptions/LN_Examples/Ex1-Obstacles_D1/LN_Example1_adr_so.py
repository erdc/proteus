from pyadh.default_so import *
from LN_Example1_adr_ellam_c0p1_n import *
pnList = [("LN_Example1_p","LN_Example1_c0p1_n"),
          ("LN_Example1_adr_p","LN_Example1_adr_ellam_c0p1_n")]
modelSpinUpList=[0]

#systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True

tnList = [i*T/nDTout for i  in range(nDTout+1)]
