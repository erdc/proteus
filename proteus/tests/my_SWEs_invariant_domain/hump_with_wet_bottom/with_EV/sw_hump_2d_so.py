from proteus.default_so import *

pnList = [("ls_vortex_2d_p","ls_vortex_2d_n"),
          ("vof_vortex_2d_p","vof_vortex_2d_n")]

systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepExact = True

name=soname

needEBQ_GLOBAL  = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
#archiveFlag = ArchiveFlags.EVERY_MODEL_STEP
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
useOneArchive = True
