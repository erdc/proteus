from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
#medium type and physical parameters are set here
from twp_darcy_modelParams import *
name = "twp_darcy_fc_drain_sand_3d"


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
nd = 3
meshfile = "het_column_8"
L = (1.0,1.0,10.0)
bottom = -10.0; top = 0.0
#auxiliary parameters and problem setup 
#where get a 'break' in ic's from bottom to top value 
xICbreak = top


g    = [0.0,0.0,-gmag]                      # gravity  with direction
#boundary conditions
#air set to 0 pressure head at top
#no flow for water at top
#bottom at equilibrium with water table at some depth
psin_top = 0.0
Se_top = 0.01#0.95#1.0                        # effective saturation top
Sw_top= Se_top*(sw_max-sw_min) + sw_min
psic_top= se2pc(Se_top)
psiw_top = psin_top - psic_top     #psi_w at top

#bottom 
waterTable  = bottom-0.1 #elevation of waterTable
psiwTable    = 0.0
#psi_w at z=0
psiw_bottom = psiwTable + g[2]/gmag*(bottom-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[2]/gmag*(bottom-waterTable)#waterTable - 0.0
psic_bottom  = psin_bottom-psiw_bottom
Se_bottom =  pc2se(psic_bottom)#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom= Se_bottom*(sw_max-sw_min) + sw_min


#for initial conditions
waterTableIC =  top#-elevation of water table 

hour = 3600.0 #[s]
T = 8.0*hour #24.0*hour #sand 12*hour, clay 24*hour                          # time [s]


analyticalSolutions = None

coefficients = TwophaseDarcy_fc(g=g, 
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
                                density_w_parameters = density_w_exponential,
                                density_n_parameters = density_n_ideal)

#now define the Dirichlet boundary conditions

def getDBC_sw(x,flag):
    #if x[2] >= top-1.0e-10:
    #    return lambda x,t: Se_top*(sw_max - sw_min)+sw_min
    if x[2] <= bottom+1.0e-10:
        return lambda x,t: Se_bottom*(sw_max - sw_min)+sw_min
def getDBC_psiw(x,flag):
    #fixed air pressure head on inflow
    if x[2] >= L[2]-1.0e-10:
        return (lambda x,t: psin_top,2)
    if x[2] <= bottom+1.0e-10:
     	return (lambda x,t: psiw_bottom,1)

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        #if x[2] >= xICbreak:
        #    return Sw_top
        return Sw_bottom
class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiwTable + g[2]/gmag*(x[2]-waterTableIC)
		

initialConditions  = {0:sw_IC(),1:psiw_IC()}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}


def get_w_AFBC(x):
    #if x[0] == 0.0:
    #   return lambda x,t: -0.9*Ksw
    pass
#    if x[0] == 0.0:
#	return lambda x,t: -q
#    if x[0] == L[0]:
#	return lambda x,t: 0.0

def get_n_AFBC(x):
    #fixed total inflow
    #if x[0] == 0.0:
    #    return lambda x,t: -0.9*Ksw
    pass
#    if x[0] == 0.0:
#	return lambda x,t: 0.0
#    if x[0] == L[0]:
#        return lambda x,t: 0.0
#advectiveFluxBoundaryConditions =  {1:get_n_AFBC}
#diffusiveFluxBoundaryConditions = {0:{1:get_w_AFBC},1:{}}

#fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{},1:{}}

