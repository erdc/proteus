"""
The default values for so-files describing split operator formulations
"""
from SplitOperator import *

name = None

systemStepControllerType = Sequential_MinModelStep

systemStepExact = True

useOneMesh = True

tnList = None

needEBQ_GLOBAL = False

needEBQ = False

modelSpinUpList = []

useOneArchive=True#False

sList = []

from Archiver import ArchiveFlags

archiveFlag = ArchiveFlags.EVERY_USER_STEP
#CEK CHANGED DEFAULT FROM EVERY_SEQUENCE_STEP

dt_system_fixed = None
