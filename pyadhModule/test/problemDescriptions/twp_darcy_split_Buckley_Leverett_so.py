from pyadh.default_so import *
from darcy_simp_modelParams import *
#pnList = [("twp_darcy_split_Buckley_Leverett_pres_p","twp_darcy_split_Buckley_Leverett_pres_c0p1_n"),
#          ("twp_darcy_split_Buckley_Leverett_sat_p","twp_darcy_split_Buckley_Leverett_sat_c0p1_n")]
pnList = [("twp_darcy_split_Buckley_Leverett_pres_p","twp_darcy_split_Buckley_Leverett_pres_c0p1_n"),
          ("twp_darcy_split_Buckley_Leverett_sat_p","twp_darcy_split_Buckley_Leverett_sat_dgp1_n")]
modelSpinUpList = [0]
name="twp_darcy_split_Buckley_Leverett"
#systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
#nDTout
T=0.5
nDTout = 1
DT = T/nDTout #None#T/100.0
tnList = [i*DT for i  in range(nDTout+1)]
