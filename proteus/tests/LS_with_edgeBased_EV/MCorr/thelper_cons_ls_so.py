from proteus.default_so import *
import thelper_cons_ls
from thelper_cons_ls import *

thelper_cons_ls.LS_model = 0
thelper_cons_ls.V_model = 0
if applyCorrection:    
    if ct.STABILIZATION_TYPE_ncls>0:
        pnList = [("thelper_ncls_p","thelper_ncls_n"),
                  ("thelper_vof_p","thelper_vof_n"), 
                  ("thelper_MCorr_p","thelper_MCorr_n")]
        thelper_cons_ls.RD_model = None
        thelper_cons_ls.VOF_model = 1
        thelper_cons_ls.MCorr_model = 2
    else:
        pnList = [("thelper_ncls_p","thelper_ncls_n"),
                  ("thelper_redist_p","thelper_redist_n"), 
                  ("thelper_vof_p","thelper_vof_n"), 
                  ("thelper_MCorr_p","thelper_MCorr_n")]
        thelper_cons_ls.RD_model = 1
        thelper_cons_ls.VOF_model = 2
        thelper_cons_ls.MCorr_model =3
else:
    if ct.STABILIZATION_TYPE_ncls>0:
        pnList = [("thelper_ncls_p","thelper_ncls_n"),
                  ("thelper_vof_p","thelper_vof_n")]
        thelper_cons_ls.RD_model = None
        thelper_cons_ls.VOF_model = 1
        thelper_cons_ls.MCorr_model = None
    else:
        pnList = [("thelper_ncls_p","thelper_ncls_n"),
                  ("thelper_redist_p","thelper_redist_n"), 
                  ("thelper_vof_p","thelper_vof_n")]
        thelper_cons_ls.RD_model = 1
        thelper_cons_ls.VOF_model = 2
        thelper_cons_ls.MCorr_model = None
    
systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepExact = True

name=soname

needEBQ_GLOBAL  = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
useOneArchive = True
