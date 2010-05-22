from pyadh.default_so import *
from twp_darcy_split_het_column_8_3d_pres_p import T
pnList = [("twp_darcy_split_het_column_8_3d_pres_p","twp_darcy_split_het_column_8_3d_pres_c0p1_n"),
          ("twp_darcy_split_het_column_8_3d_sat_p","twp_darcy_split_het_column_8_3d_sat_c0p1_n")]
name="twp_darcy_split_het_column_8_3d_cgs"
modelSpinUpList = [0]
systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
#systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
nDTout = 11#51
DT = T/(nDTout-1.)
tnList = [i*DT for i  in range(nDTout)]
