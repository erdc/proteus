from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
#medium type and physical parameters are set here
from twp_darcy_modelParams import *

#for convenience, to get bcs straight, need to put these somewhere
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


name = "twp_darcy_fcff_vgm_sand_1d"
#spatial domain
L = (1.0,1.0,1.0)
#psk model type
model= 'VGM'
#auxiliary parameters and problem setup 
#where get a 'break' in ic's from bottom to top value 
xICbreak = L[0]*1.0

g    = [-gmag]                      # gravity  with direction
Se_top = 0.9                        # effective saturation top
psi_top   = 0.1                     # pressure head at top [m]
waterTableLocation = -5.0           # where is the water table
psi_waterTable     = 0.0#-10.0#0.0            # pressure head at water table
                                    # hydrostatic ic psi(x) = psi_0 + (z-z_0)rho g
psi_bottom         = psi_waterTable + (0.0 - waterTableLocation)*g[0]/gmag
psic_bottom        = 0.0 - psi_bottom # assume constant air phase pressure head = 0
Se_bottom          = seVGM(psic_bottom,mvg_alpha,mvg_n,mvg_m)                  # effective saturation bottom

#test if can maintain hydrostatic equilibrium
#psi_top = psi_waterTable + (L[0] - waterTableLocation)*g[0]/gmag
#psic_top = 0.0 - psi_top # assume constant air phase pressure head = 0
#Se_top  = seVGM(psic_top,mvg_alpha,mvg_n,mvg_m)    
T = 2000.0                          # time [s]


analyticalSolutions = None

coefficients = TwophaseDarcy_fc_ff(g=g, 
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

class psiw_IC:
    def __init__(self):
        pass
    #def uOfXT(self,x,t):
    #    return psi_waterTable +  g[0]/gmag*(x[0] - waterTableLocation)
    def uOfXT(self,x,t):
        if x[0] <= xICbreak:
            return psi_bottom
        else:
            return psi_top
class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[0] <= xICbreak:
            #psic = 0.0 - psiw_IC().uOfXT(x,t)
            #Se = seVGM(psic,mvg_alpha,mvg_n,mvg_m)
            #mwf debub
            #print "sw_IC x=%s Se=%s " % (x[0],Se)
            Se = Se_bottom
            return Se*(sw_max - sw_min)+sw_min
        else:
            return Se_top*(sw_max - sw_min)+sw_min
	
initialConditions  = {0:sw_IC(),1:psiw_IC()}

#fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}


def get_w_AFBC(x):
    pass
def get_n_AFBC(x):
    pass

advectiveFluxBoundaryConditions =  {1:get_n_AFBC}
diffusiveFluxBoundaryConditions = {0:{1:get_w_AFBC},1:{}}

fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

