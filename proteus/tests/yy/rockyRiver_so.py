from proteus.default_so import *
import rockyRiver

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
      ]#3

# rockyRiver.VOF_model=None
# rockyRiver.VOS_model=None
# rockyRiver.SED_model=None
# rockyRiver.V_model          =0
# rockyRiver.PINC_model       =1
# rockyRiver.PRESSURE_model   =2
# rockyRiver.PINIT_model      =3

# if rockyRiver.useRANS > 0:
#     pnList.append(("kappa_p",
#                    "kappa_n"))
#     pnList.append(("dissipation_p",
#                    "dissipation_n"))
name = "rockyRiver"

#modelSpinUpList = [cylinder.VOF_model, cylinder.LS_model, cylinder.V_model, cylinder.PINIT_model]
#modelSpinUpList = [rockyRiver.PINIT_model]###############

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
        self.modelList = modelList[:len(pnList)]###############-1

systemStepControllerType = Sequential_MinAdaptiveModelStepPS



needEBQ_GLOBAL = False
needEBQ = False

tnList = rockyRiver.tnList
dt_system_fixed = rockyRiver.dt_fixed

systemStepExact = True


# info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP


