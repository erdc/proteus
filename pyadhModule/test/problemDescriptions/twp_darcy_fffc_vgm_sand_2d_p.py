from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_modelParams import *
model = 'VGM'
name = 'twp_darcy_fffc_vgm_sand_2d_se_1_sgs'

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
L=(3.0,3.0,1.0)
nd = 2
g    = [0.0,-gmag] #   # gravity              (Earth 9.8 m/s^2) 
xICbreak = L[1]

analyticalSolutions = None

coefficients = TwophaseDarcy_fc_ff(g=g, 
                                   rhon=rhon,
                                   rhow=rhow,
                                   mun=mun,
                                   muw=muw,
                                   Ksw=Ksw,
                                   psk_model=model,
                                   vg_alpha = mvg_alpha,
                                   vg_m  = mvg_m,
                                   bc_pd  = bc_pd, 
                                   bc_lambda = bc_lambda,
                                   omega  = omega,
                                   Sw_max = sw_max,
                                   Sw_min = sw_min)

#for bc's
#Q_m_per_d = 0.1 
#q    = Q_m_per_d*m_per_d_to_m_per_s

#now define the Dirichlet boundary conditions
Se_top    = 1.0
Se_bottom = 0.2
Sw_top    = Se_top*(sw_max-sw_min) + sw_min
psi_top   = 0.1

waterTable = 0.0#-1.0#elevation of water table 
psiTable   = 0.0     #psi_w at water table
#psi_w at z=0
psi_bottom = psiTable + g[1]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[1]/gmag*(0.0-waterTable)#waterTable - 0.0
pc_bottom  = psin_bottom-psi_bottom
Se_bottom  = pc2se(pc_bottom)#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom  = Se_bottom*(sw_max-sw_min) + sw_min

slit_left = L[0]/3.0
slit_right= L[0]*2.0/3.0

hour = 3600.0 #[s]
T = 6470#1.2*hour                          # time [s]


def getDBC_sw(x,flag):
    if x[1] == L[1] and slit_left <= x[0] and x[0] <= slit_right:
        return lambda x,t: Sw_top
    if x[1] == 0.0:
        return lambda x,t: Sw_bottom
def getDBC_psiw(x,flag):
    #fixed head on inflow
    if x[1] == L[1] and slit_left <= x[0] and x[0] <= slit_right:
        return lambda x,t: psi_top
    if x[1] == 0.0:
     	return lambda x,t: psi_bottom


dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}
class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= xICbreak and slit_left <= x[0] and x[0] <= slit_right:
            return Sw_top
        pc = pc_bottom + (rhon/rhow - 1.0)*g[1]/gmag*x[1]
        se = pc2se(pc)
        return se*(sw_max - sw_min) + sw_min

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiTable + g[1]/gmag*(x[1]-waterTable)
		
initialConditions  = {0:sw_IC(),1:psiw_IC()}




V = -Ksw*0.5e2

def get_w_psiw_DFBC(x):
#     return lambda x,t: 0.0
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    #if x[1] == 0.0:
    #    return lambda x,t: -V
def get_w_psic_DFBC(x):
#     return lambda x,t: 0.0
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0

def get_n_AFBC(x):
#     return lambda x,t: 0.0
    if x[0] in [0.0,L[0]]:
       return lambda x,t: 0.0
    if x[1] == 0.0:
       return lambda x,t: -V

def get_n_psiw_DFBC(x):
#     return lambda x,t: 0.0
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    if x[1] == 0.0:
        return lambda x,t: 0.0

def get_n_psic_DFBC(x):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0:{},1:{}}
fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}
#advectiveFluxBoundaryConditions =  {1:get_n_AFBC}
#diffusiveFluxBoundaryConditions = {0:{0:get_w_psic_DFBC,
#                                      1:get_w_psiw_DFBC},
#                                   1:{0:get_n_psic_DFBC,
#                                      1:get_n_psiw_DFBC}}


