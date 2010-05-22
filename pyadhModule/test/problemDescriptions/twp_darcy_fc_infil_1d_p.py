from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
#medium type and physical parameters are set here
from twp_darcy_modelParams import *
name = "twp_darcy_fc_infil_1d_incomp"

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
Se_top = 0.8#0.95#1.0                        # effective saturation top
Sw_top= Se_top*(sw_max-sw_min) + sw_min
psi_top = 0.1     #psi_w at top
pc_top  = se2pc(Se_top)

psin_top= psi_top + pc_top
waterTable =  -1.0#-1.0#elevation of water table 
psiTable= 0.0     #psi_w at water table
#psi_w at z=0
psi_bottom = psiTable + g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
pc_bottom  = psin_bottom-psi_bottom
Se_bottom =  pc2se(pc_bottom)#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom= Se_bottom*(sw_max-sw_min) + sw_min

hour = 3600.0 #[s]
T = 24.0*hour #sand 12*hour, clay 24*hour                          # time [s]


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
                                density_n_parameters = density_n_ideal,
                                diagonalHet = True,
                                sd = sd)

#now define the Dirichlet boundary conditions

def getDBC_sw(x,flag):
    if x[0] == L[0]:
        return lambda x,t: Se_top*(sw_max - sw_min)+sw_min
    if x[0] == 0.0:
        return lambda x,t: Se_bottom*(sw_max - sw_min)+sw_min
def getDBC_psiw(x,flag):
    #fixed head on inflow
    if x[0] == L[0]:
        return lambda x,t: psi_top
    if x[0] == 0.0:
     	return lambda x,t: psi_bottom

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[0] >= xICbreak:
            return Sw_top
        pc = pc_bottom + (rhon/rhow - 1.0)*g[0]/gmag*x[0]
        se = pc2se(pc)
        return se*(sw_max - sw_min) + sw_min
class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiTable + g[0]/gmag*(x[0]-waterTable)
		

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

fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{},1:{}}

