from pyadh.default_so import *
from twp_darcy_riem_1d import *
pnList = [("twp_darcy_split_riem_1d_pres_p","twp_darcy_split_riem_1d_pres_c0p1_n"),
          ("twp_darcy_split_riem_1d_sat_p","twp_darcy_split_riem_1d_sat_c0p1_n")]
name="twp_darcy_split_riem_1d"
modelSpinUpList = [0]
systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
#systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
#nDTout
DT = T/100.0
tnList = [i*DT for i  in range(101)]
