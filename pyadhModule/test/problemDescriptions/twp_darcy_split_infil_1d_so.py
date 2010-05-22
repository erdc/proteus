from pyadh.default_so import *
from twp_darcy_split_infil_1d_sat_p import T
pnList = [("twp_darcy_split_infil_1d_pres_p","twp_darcy_split_infil_1d_pres_c0p1_n"),
          ("twp_darcy_split_infil_1d_sat_p","twp_darcy_split_infil_1d_sat_c0p1_n")]
name="twp_darcy_split_infil_1d_bcb"
modelSpinUpList = [0]
systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
#systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
#nDTout
DT = T/100.0
tnList = [i*DT for i  in range(101)]
