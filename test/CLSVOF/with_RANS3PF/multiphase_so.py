from proteus.default_so import *
try:
    from . import multiphase
except:
    import multiphase

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

multiphase.VOS_model=None
multiphase.SED_model=None
if multiphase.useCLSVOF:
    pnList = [("clsvof_p",               "clsvof_n"),#0
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#1
              ("pressureincrement_p", "pressureincrement_n"),#2
              ("pressure_p", "pressure_n"),#3
              ("pressureInitial_p", "pressureInitial_n")]#4
    multiphase.VOF_model=None
    multiphase.LS_model=None
    multiphase.RD_model=None
    multiphase.MCORR_model=None
    multiphase.CLSVOF_model=0
    multiphase.V_model=1
    multiphase.PINC_model=2
    multiphase.PRESSURE_model=3
    multiphase.PINIT_model=4
else:
    pnList = [("vof_p",               "vof_n"),#0
              ("ls_p",                "ls_n"),#1
              ("redist_p",            "redist_n"),#2
              ("ls_consrv_p",         "ls_consrv_n"),#3
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
              ("pressureincrement_p", "pressureincrement_n"),#5
              ("pressure_p", "pressure_n"),#6
              ("pressureInitial_p", "pressureInitial_n")]#7
    multiphase.CLSVOF_model=None
    multiphase.VOF_model=0
    multiphase.LS_model=1
    multiphase.RD_model=2
    multiphase.MCORR_model=3
    multiphase.V_model=4
    multiphase.PINC_model=5
    multiphase.PRESSURE_model=6
    multiphase.PINIT_model=7

name = "multiphase"

modelSpinUpList = [multiphase.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,multiphase.dt_init]+[i*multiphase.dt_fixed for i in range(1,multiphase.nDTout+1)]

info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
