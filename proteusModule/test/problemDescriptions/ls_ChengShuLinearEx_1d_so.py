from pyadh.default_so import *
from ls_ChengShuLinearEx_1d_p import T
pnList = [("ls_ChengShuLinearEx_1d_p","ls_ChengShuLinearEx_1d_dgpk_n")]
name="ChengShuEx1_1d_LeSaintRaviartdgp3nn161"


systemStepControllerType = Sequential_NonUniformFixedStep

needEBQ_GLOBAL  = True
needEBQ = True
nDTout = 1
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
