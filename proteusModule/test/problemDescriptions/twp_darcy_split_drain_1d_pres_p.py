from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
#medium type and physical parameters are set here
from twp_darcy_modelParams import *
name = "twp_darcy_split_drain_sand_1d_pres"

#psk model type
model= 'VGM'

se2pc = None
pc2se = None
if model == 'VGM':
    se2pc = lambda se : pcVGM(se,mvg_alpha,mvg_n,mvg_m)
    pc2se = lambda pc : seVGM(pc,mvg_alpha,mvg_n,mvg_m)
elif model == 'BCB':
    se2pc = lambda se : pcBCB(se,bc_pd,bc_lambda)
    pc2se = lambda pc : seBCB(pc,bc_pd,bc_lambda)
#

#spatial domain
L = (10.0,1.0,1.0)
#auxiliary parameters and problem setup 
#where get a 'break' in ic's from bottom to top value 
xICbreak = L[0]


g    = [-gmag]                      # gravity  with direction
#boundary conditions
#air set to 0 pressure head at top
#no flow for water at top
#bottom at equilibrium with water table at some depth
psin_top = 0.0
Se_top = 0.1#0.95#1.0                        # effective saturation top
Sw_top= Se_top*(sw_max-sw_min) + sw_min
psic_top= se2pc(Se_top)
psiw_top = psin_top - psic_top     #psi_w at top

#bottom 
waterTable  = -0.1 #elevation of waterTable
psiwTable    = 0.0
#psi_w at z=0
psiw_bottom = psiwTable + g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
psic_bottom  = psin_bottom-psiw_bottom
Se_bottom =  pc2se(psic_bottom)#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom= Se_bottom*(sw_max-sw_min) + sw_min


#for initial conditions
waterTableIC =  L[0]#-elevation of water table 

hour = 3600.0 #[s]
T = 1200*hour#24.0*hour #sand 12*hour, clay 24*hour                          # time [s]


analyticalSolutions = None

coefficients = TwophaseDarcy_split_pressure(g=g, 
                                            rhon=rhon,
                                            rhow=rhow,
                                            mun    = mun,
                                            muw    = muw,
                                            Ksw=Ksw,
                                            psk_model=model,
                                            vg_alpha = mvg_alpha,
                                            vg_m  = mvg_m,
                                            bc_pd  = bc_pd, 
                                            bc_lambda = bc_lambda,
                                            omega  = omega,
                                            Sw_max = sw_max,
                                            Sw_min = sw_min,
                                            diagonalHet = True,
                                            sd = sd)

#now define the Dirichlet boundary conditions

def getDBC_psiw(x,flag):
    #fixed head on inflow
    if x[0] == L[0]:
        return (lambda x,t: psin_top,2)
    if x[0] == 0.0:
     	return (lambda x,t: psiw_bottom,1)

dirichletConditions = {0:getDBC_psiw}

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiwTable + g[0]/gmag*(x[0]-waterTableIC)
		

initialConditions  = {0:psiw_IC()}

fluxBoundaryConditions = {0:'outFlow'}


advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

