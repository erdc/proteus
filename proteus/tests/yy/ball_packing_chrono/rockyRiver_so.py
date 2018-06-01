from proteus.default_so import *
import rockyRiver

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

# if rockyRiver.sedimentDynamics:
#     pnList = [("vos_p",               "vos_n"),#0
#               ("vof_p",               "vof_n"),#1
#               ("ls_p",                "ls_n"),#2
#               ("redist_p",            "redist_n"),#3
#               ("ls_consrv_p",         "ls_consrv_n"),#4
#               ("threep_navier_stokes_sed_p", "threep_navier_stokes_sed_n"),#5
#               ("twp_navier_stokes_p", "twp_navier_stokes_n"),#6
#               ("pressureincrement_p", "pressureincrement_n"),#7
#               ("pressure_p", "pressure_n"),#8
#               ("pressureInitial_p", "pressureInitial_n")]#9
#     cylinder.VOS_model=0
#     cylinder.VOF_model=1
#     cylinder.LS_model=2
#     cylinder.RD_model=3
#     cylinder.MCORR_model=4
#     cylinder.SED_model=5
#     cylinder.V_model=6
#     cylinder.PINC_model=7
#     cylinder.PRESSURE_model=8
#     cylinder.PINIT_model=9
# else:
#     pnList = [("vof_p",               "vof_n"),#0
#               ("ls_p",                "ls_n"),#1
#               ("redist_p",            "redist_n"),#2
#               ("ls_consrv_p",         "ls_consrv_n"),#3
#               ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
#               ("pressureincrement_p", "pressureincrement_n"),#5
#               ("pressure_p", "pressure_n"),#6
#               ("pressureInitial_p", "pressureInitial_n")]#7
#     rockyRiver.VOS_model=None
#     rockyRiver.SED_model=None
#     rockyRiver.VOF_model=0
#     rockyRiver.LS_model=1
#     rockyRiver.RD_model=2
#     rockyRiver.MCORR_model=3
#     rockyRiver.V_model=4
#     rockyRiver.PINC_model=5
#     rockyRiver.PRESSURE_model=6
#     rockyRiver.PINIT_model=7

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
      ("pressureincrement_p", "pressureincrement_n"),#1
      ("pressure_p", "pressure_n"),#2
      ("pressureInitial_p", "pressureInitial_n")]#3

rockyRiver.VOF_model=None
rockyRiver.VOS_model=None
rockyRiver.SED_model=None
rockyRiver.V_model=0
rockyRiver.PINC_model=1
rockyRiver.PRESSURE_model=2
rockyRiver.PINIT_model=3

if rockyRiver.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "rockyRiver"

#modelSpinUpList = [cylinder.VOF_model, cylinder.LS_model, cylinder.V_model, cylinder.PINIT_model]
modelSpinUpList = [rockyRiver.PINIT_model]

# class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
#     def __init__(self,modelList,system=defaultSystem,stepExact=True):
#         Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
#         self.modelList = modelList[:len(pnList)-1]


# systemStepControllerType = Sequential_MinAdaptiveModelStepPS

#There is a bug here that needs to be fixed; 
# the Sequential_FixedStep_Simple ignores the pressure gradient terms?
# systemStepControllerType = Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = 0.01; 
# systemStepExact=False;


class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]

systemStepControllerType = Sequential_MinAdaptiveModelStepPS



needEBQ_GLOBAL = False
needEBQ = False

tnList = rockyRiver.tnList
dt_system_fixed = rockyRiver.dt_fixed

systemStepExact = True


# info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP


