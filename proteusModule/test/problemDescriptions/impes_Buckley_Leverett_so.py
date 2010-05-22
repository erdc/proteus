from pyadh.default_so import *
from darcy_simp_modelParams import *
pnList = [("impes_Buckley_Leverett_pres_p","impes_Buckley_Leverett_pres_c0p1_n"),
          ("impes_Buckley_Leverett_sat_p","impes_Buckley_Leverett_sat_c0p2_n")]
#pnList = [("impes_Buckley_Leverett_pres_p","impes_Buckley_Leverett_pres_c0p1_n"),
#          ("impes_Buckley_Leverett_sat_p","impes_Buckley_Leverett_sat_dgp2_n")]
modelSpinUpList = [0]
name="impes_Buckley_Leverett"
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
