from pyadh.default_so import *
from gw_het_fivespot import *
pnList = [("gw_het_fivespot_2d_p","gw_het_fivespot_2d_c0p1_n"),
          ("transport_het_fivespot_2d_p","transport_het_fivespot_2d_c0p1_n")]
modelSpinUpList=[0]

#systemStepControllerType = Sequential_MinFLCBDFModelStep
#systemStepControllerType = Sequential_FixedStep
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = True
needEBQ = True

tnList = [i*T/nDTout for i  in range(nDTout+1)]
