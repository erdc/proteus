from pyadh.default_so import *
from impes_modelParams import *
pnList = [("impes_pres_p","impes_pres_n"),
          ("impes_sat_p","impes_sat_dg_n")]
name="impes_twophase_dg"
systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
nOutput=100
DT = T/float(nOutput)
tnList = [i*DT for i  in range(101)]
