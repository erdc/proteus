from pyadh.default_so import *
from square import *
if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_square_1d_p","ls_square_1d_c0p1_n"),
                  ("redist_square_1d_p","redist_square_1d_c0p1_n"), 
                  ("vof_square_1d_p","vof_square_1d_c0p1_n"), 
                  ("ls_consrv_square_1d_p","ls_consrv_square_1d_c0p1_n")]
    else:
        pnList = [("ls_square_1d_p","ls_square_1d_c0p1_n"),
                  ("redist_square_1d_p","redist_square_1d_c0p1_n"),
                  ("vof_square_1d_p","vof_square_1d_c0p1_n")]
else:
    pnList = [("ls_square_1d_p","ls_square_1d_c0p1_n")] 
#mwf hack test NCLS              ("vof_square_1d_p","vof_square_1d_c0p1_n")]#also dgp1,dgpk
name="ls_square_1d_c0p1"
if "FLCBDF" in [timeIntegration_vof,timeIntegration_ls]:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    #systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL  = False#True
needEBQ = False#True
nDTout = 1
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
