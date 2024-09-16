from proteus.default_so import *
from . import vortex2D
from .vortex2D import *

pnList = [("ncls_p","ncls_n"),
          ("rdls_p","rdls_n")]
#pnList = [("ncls_p","ncls_n")]
 

systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepExact = True

name=soname

needEBQ_GLOBAL  = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
useOneArchive = True

#tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
tnList=[0.,1E-6]
