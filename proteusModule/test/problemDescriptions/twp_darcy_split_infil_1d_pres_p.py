from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_modelParams import *
name = "twp_darcy_split_infil_1d_pres_bcb"

model = 'BCB'#'VGM'
#for convenience , to get bcs straight
def seVGM(psic,alVG,nVG,mVG):
    if psic <= 0: return 1.0
    tmp1 = pow(alVG*psic,nVG)
    tmp2 = pow(1.+tmp1,-mVG)
    return min(max(tmp2,0.0),1.0)
def pcVGM(se,alVG,nVG,mVG):
    if se >= 1.: return 0.0
    tmp1 = pow(se,-1./mVG)
    tmp2 = pow(tmp1-1.0,1.0/nVG)/alVG
    return tmp2
def seBCB(psic,pdBC,lamBC):
    if psic <= pdBC: return 1.0
    tmp1 = pow(pdBC/psic,lamBC)
    return min(max(tmp2,0.0),1.0)
def pcBCB(se,pdBC,lamBC):
    if se >= 1.0: return 0.0
    tmp1 = pow(se,-1.0/lamBC)
    return pdBC*tmp1

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
psi_top = 0.1     #psi_w at top

waterTable =  -1.0#elevation of water table 
psiTable= 0.0     #psi_w at water table
#psi_w at z=0
psi_bottom = psiTable + g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[0]/gmag*(0.0-waterTable)#waterTable - 0.0
pc_bottom  = psin_bottom-psi_bottom

hour = 3600.0 #[s]
T = 24.0*hour #clay 24.0*hour# sand 12*hour                          # time [s]


phase = 'potential' 

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
                                            Sw_min = sw_min)

#now define the Dirichlet boundary conditions

def getDBC_psiw(x,flag):
    #fixed head on inflow
    if x[0] == L[0]:
        return lambda x,t: psi_top
    if x[0] == 0.0:
     	return lambda x,t: psi_bottom


dirichletConditions = {0:getDBC_psiw}

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiTable + g[0]/gmag*(x[0]-waterTable)
	
initialConditions  = {0:psiw_IC()}

fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x):
    pass
#    if x[0] == 0.0:
#	return lambda x,t: -q

advectiveFluxBoundaryConditions =  {0:getAFBC}
diffusiveFluxBoundaryConditions = {0:{}}


