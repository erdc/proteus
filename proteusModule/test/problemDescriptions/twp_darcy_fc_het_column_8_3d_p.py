from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
#medium type and physical parameters are set here
from twp_darcy_modelParams import *
#adh column with two-phase flow
name = "twp_darcy_fc_het_col"
model = 'VGM'


#spatial domain
nd = 3
L = (1.0,1.0,10.0)
bottom = -10.0; top = 0.0

meshfile = "het_column_8"
#temporal domain
hour = 3600.0
T = 240.0*hour#1000.0

g = [0.0,0.0,-gmag]

#block heterogeneity
useHet = True
top_mat_index= 0
if useHet:
    bot_mat_index= 2
else:
    bot_mat_index= 0
thetaS_types = numpy.array([0.5,0.4,0.33])
thetaR_types = numpy.array([0.105*0.5,0.074*0.4,0.179*0.33])
sw_max_types = numpy.ones((3,),'d')
sw_min_types = numpy.array([0.105,0.074,0.179])
mvg_alpha_types = numpy.array([4.42*0.3048, 0.487*0.3048, 0.244*0.3048])
mvg_n_types     = numpy.array([2.68,1.37,1.09])
mvg_m_types     = 1.0 - 1.0/mvg_n_types
bc_pd_types     = 1.0/mvg_alpha_types
bc_lambda_types = mvg_n_types - 1.0
#m/d
Ksw_types_md = numpy.array([1.1808e-1,1.1808e-3,1.1808e-5])
#m/s
Ksw_types    = Ksw_types_md*m_per_d_to_m_per_s


se2pc = None
pc2se = None
if model == 'VGM':
    se2pc = lambda se : pcVGM(se,mvg_alpha[top_mat_index],mvg_n[top_mat_index],mvg_m[top_mat_index])
    pc2se = lambda pc : seVGM(pc,mvg_alpha[top_mat_index],mvg_n[top_mat_index],mvg_m[top_mat_index])
elif model == 'BCB':
    se2pc = lambda se : pcBCB(se,bc_pd[top_mat_index],bc_lambda[top_mat_index])
    pc2se = lambda pc : seBCB(pc,bc_pd[top_mat_index],bc_lambda[top_mat_index])
#


coefficients = TwophaseDarcy_fc(g=g, 
                                rhon=rhon,
                                rhow=rhow,
                                mun    = mun,
                                muw    = muw,
                                Ksw=Ksw_types[top_mat_index],
                                psk_model=model,
                                vg_alpha = mvg_alpha_types[top_mat_index],
                                vg_m  = mvg_m_types[top_mat_index],
                                bc_pd  = bc_pd_types[top_mat_index], 
                                bc_lambda = bc_lambda_types[top_mat_index],
                                omega  = thetaS_types[top_mat_index],
                                Sw_max = sw_max_types[top_mat_index],
                                Sw_min = sw_min_types[top_mat_index],
                                density_w_parameters = density_w_exponential,
                                density_n_parameters = density_n_ideal)

if useHet:
    coefficients.setMaterialTypes(Ksw_types = Ksw_types,
                                  omega_types  = thetaS_types,
                                  Sw_max_types = sw_max_types,
                                  Sw_min_types = sw_min_types,
                                  bc_lambda_types = bc_lambda_types,
                                  bc_pd_types = bc_pd_types,
                                  vg_alpha_types = mvg_alpha_types,
                                  vg_m_types = mvg_m_types)

#boundary conditions
#boundary conditions
#air set to 0 pressure head at top
#no flow for water at top
#bottom at equilibrium with water table at some depth

#if try to force psiw_top based on psin_top and pc curve with fixed saturation
psin_top    = 0.0
se_top      = 0.01
sw_top      = se_top*(sw_max_types[top_mat_index]-sw_min_types[top_mat_index]) + sw_min_types[top_mat_index]
psic_top    = pcVGM(se_top,mvg_alpha_types[top_mat_index],mvg_n_types[top_mat_index],mvg_m_types[top_mat_index])
psiw_top    = psin_top - psic_top

#if just try to set psiw to be dry at top and no flow for aqueous phase
set_sw_top  = False#True
#psiw_top    = -5.0

#bottom 
waterTable  = bottom -0.1 #elevation of waterTable
psiwTable   = 0.0
#psi_w at z=0
psiw_bottom = psiwTable + g[2]/gmag*(bottom-waterTable)
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[2]/gmag*(bottom-waterTable)
psic_bottom= psin_bottom-psiw_bottom
se_bottom = seVGM(psic_bottom,mvg_alpha_types[bot_mat_index],mvg_n_types[bot_mat_index],mvg_m_types[bot_mat_index])
sw_bottom = se_bottom*(sw_max_types[bot_mat_index]-sw_min_types[bot_mat_index]) + sw_min_types[bot_mat_index]


print """\t twp_darcy_fc_het_column_8_3d psin_top=%g se_top=%g sw_top=%g psic_top=%g 
\t psiw_top=%g psiw_bottom=%g psic_bottom=%g se_bottom=%g sw_bottom=%g""" % (psin_top,se_top,
                                                                             sw_top,psic_top,
                                                                             psiw_top,
                                                                             psiw_bottom,psic_bottom,
                                                                             se_bottom,sw_bottom)

#for initial conditions
waterTableIC =  top#-elevation of water table 

def getDBC_sw(x,flag):
    #if x[2] >= top-1.0e-10 and set_sw_top:
    #    return lambda x,t: sw_top
    if x[2] <= bottom+1.0e-10:
        return lambda x,t: sw_bottom

def getDBC_psiw(x,flag):
    if x[2] >= top-1.0e-10:
        return (lambda x,t:psin_top,2)
    if x[2] <= bottom+1.0e-10:
     	return (lambda x,t: psiw_bottom,1)
    

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class sw_IC:
    def __init__(self,sw0=sw_bottom):
        self.sw0 = sw0
    def uOfXT(self,x,t):
        #mwf hack desaturate some at top?
        #if x[2] >= 0.0-1.0e-8:
        #    return 0.95
        return self.sw0

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiwTable + g[2]/gmag*(x[2]-waterTableIC)
		
    
initialConditions  = {0:sw_IC(),1:psiw_IC()}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}
analyticalSolution = None
