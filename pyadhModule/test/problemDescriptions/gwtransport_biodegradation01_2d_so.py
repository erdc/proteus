from pyadh.default_so import *
from gw_biodegradation01 import *
pnList = [("gw_biodegradation01_2d_p","gw_biodegradation01_2d_c0p1_n"),
          ("transport_biodegradation01_2d_p","transport_biodegradation01_2d_c0p1_n")]
modelSpinUpList=[0]
name = "gwbio01"

systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
#systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True

tnList = [i*T/nDTout for i  in range(nDTout+1)]
