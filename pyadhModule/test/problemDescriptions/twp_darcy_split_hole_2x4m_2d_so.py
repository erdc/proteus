from pyadh.default_so import *
from twp_darcy_split_hole_2x4m_2d_pres_p import T,tnList
pnList = [("twp_darcy_split_hole_2x4m_2d_pres_p","twp_darcy_split_hole_2x4m_2d_pres_c0p1_n"),
          ("twp_darcy_split_hole_2x4m_2d_sat_p","twp_darcy_split_hole_2x4m_2d_sat_c0p1_n")]
name="twp_darcy_split_hole_2x4m_2d_clay_cgs"
modelSpinUpList = [0]
systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
#systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True
nDTout = 2
#DT = 1#T/100.0
#tnList = [0,1.]
#T=1.
#tnList = [i*DT for i  in range(101)]
