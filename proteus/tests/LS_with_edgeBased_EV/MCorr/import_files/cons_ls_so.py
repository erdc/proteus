from proteus.default_so import *
import cons_ls
from cons_ls import *

cons_ls.LS_model = 0
cons_ls.V_model = 0
if applyCorrection:    
    if ct.STABILIZATION_TYPE_ncls>0:
        pnList = [("ncls_p","ncls_n"),
                  ("vof_p","vof_n"), 
                  ("MCorr_p","MCorr_n")]
        cons_ls.RD_model = None
        cons_ls.VOF_model = 1
        cons_ls.MCorr_model = 2
    else:
        pnList = [("ncls_p","ncls_n"),
                  ("redist_p","redist_n"), 
                  ("vof_p","vof_n"), 
                  ("MCorr_p","MCorr_n")]
        cons_ls.RD_model = 1
        cons_ls.VOF_model = 2
        cons_ls.MCorr_model =3
else:
    if ct.STABILIZATION_TYPE_ncls>0:
        pnList = [("ncls_p","ncls_n"),
                  ("vof_p","vof_n")]
        cons_ls.RD_model = None
        cons_ls.VOF_model = 1
        cons_ls.MCorr_model = None
    else:
        pnList = [("ncls_p","ncls_n"),
                  ("redist_p","redist_n"), 
                  ("vof_p","vof_n")]
        cons_ls.RD_model = 1
        cons_ls.VOF_model = 2
        cons_ls.MCorr_model = None
    
systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepExact = True

name=soname

needEBQ_GLOBAL  = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
useOneArchive = True
