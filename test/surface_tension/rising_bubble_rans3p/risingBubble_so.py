from proteus.default_so import *
try:
    from . import risingBubble
except:
    import risingBubble

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("vof_p",               "vof_n"),#0
          ("ls_p",                "ls_n"),#1
          ("redist_p",            "redist_n"),#2
          ("ls_consrv_p",         "ls_consrv_n"),#3
          ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
          ("pressureincrement_p", "pressureincrement_n"),#5
          ("pressure_p", "pressure_n"),#6
          ("pressureInitial_p", "pressureInitial_n")]#7
risingBubble.VOS_model=None
risingBubble.SED_model=None
risingBubble.VOF_model=0
risingBubble.LS_model=1
risingBubble.RD_model=2
risingBubble.MCORR_model=3
risingBubble.V_model=4
risingBubble.PINC_model=5
risingBubble.PRESSURE_model=6
risingBubble.PINIT_model=7

name = "risingBubble"

#modelSpinUpList = [risingBubble.VOF_model, risingBubble.LS_model, risingBubble.V_model, risingBubble.PINIT_model]
modelSpinUpList = [risingBubble.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,risingBubble.dt_init]+[i*risingBubble.dt_fixed for i in range(1,risingBubble.nDTout+1)]

info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
