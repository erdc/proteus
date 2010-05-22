from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
import lens2d 
#default setup
nd = 2
top   = 2.0 #m
right = 4.0 #m
domain = [[0.0,right],[0.0,top]]
slit_left = 1.9
slit_right= 2.1
lensX = [1.5,2.5]; lensY=[1.25,1.5]
lens  = [lensX,lensY]
sourceXs = [slit_left-0.1,slit_left,slit_right,slit_right+0.1]

polyfile = 'lens2d'
L,boundaryTags = lens2d.genPoly(domain=domain,
                                lens=lens,
                                sourceXs=sourceXs,
                                polyfileBase=polyfile)


useHet = False#True
useV2  = False
useForsyth = True
model = 'VGM'#'BCB' #'VGM'
# grab a sand from Forsyth example
#     Type 3: k_i     = 4.898e-12 m^2
#             K_s     = 4.143
#             \theta_s= 0.3250
#             \theta_r= 0.3250*0.2643
#             \alpha  = 3.455     1/m
#             n       = 5.0
#     Type 2: k_i     = 5.55e-12 m^2
#             K_s     = 4.69     m/d
#             \theta_s= 0.3510
#             \theta_r= 0.3510*0.2806
#             \alpha  = 3.63     1/m
#             n       = 1.632

#helmig sanbox example 1 (ch 6)
#     Type 1: k_i     = 6.64e-11 m^2
#             K_s     = 7.0e-4 m/s
#             \theta_s= 0.40
#             \theta_r= 0.40*0.09
#             \lambda = 2.7
#             p_d     = 755 [Pa] -> psi_d= 0.0772 [m]
#             \alpha  = 1/p_d = 13.0 [1/m]
#             n       = lambda+1 = 3.7
#     Type 2: k_i     = 7.15e-12 m^2
#             K_s     = 7.0e-5 m/s
#             \theta_s= 0.39
#             \theta_r= 0.40*0.12
#             \lambda = 2.0
#             p_d     = 2060 [Pa] -> psi_d= 0.211 [m]
#             \alpha  = 1/p_d = 4.75
#             n       = lambda+1 = 3.


nMediaTypes  = 2
alphaVGtypes = numpy.zeros((nMediaTypes,),'d')
nVGtypes     = numpy.zeros((nMediaTypes,),'d')
mVGtypes     = numpy.zeros((nMediaTypes,),'d')
thetaStypes  = numpy.zeros((nMediaTypes,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes,),'d')
KsTypes      = numpy.zeros((nMediaTypes,),'d')
SwMaxTypes   = numpy.ones((nMediaTypes,),'d')
SwMinTypes   = numpy.zeros((nMediaTypes,),'d')

#dummys for Brooks Corey right now
#assign dummy types for Brooks Corey right now
lambdaBCtypes = numpy.ones((nMediaTypes,),'d')
lambdaBCtypes *=3.0
pdBCtypes    = numpy.ones((nMediaTypes,),'d')
pdBCtypes     *= 1.0/3.455

#put everything in [s]
m_per_s_by_m_per_d = 1.1574074e-5

if useForsyth:
    #Forsyth2 Type 3
    it = 0;
    alphaVGtypes[it] = 3.455;  nVGtypes[it]    = 5.0; mVGtypes[it]    = 1.0-1.0/nVGtypes[it];
    thetaStypes[it]  = 0.3250; thetaRtypes[it] = 0.3250*0.2643; SwMinTypes[it] = 0.2643;
    KsTypes[it]      = 4.143*m_per_s_by_m_per_d;
    lambdaBCtypes[it] = nVGtypes[it]-1.0; pdBCtypes[it] = 1.0/alphaVGtypes[it];
    #Forsyth2 Type 2
    it = 1;
    alphaVGtypes[it] = 3.63 ;  nVGtypes[it]    = 1.632; mVGtypes[it]    = 1.0-1.0/nVGtypes[it];
    thetaStypes[it]  = 0.3510; thetaRtypes[it] = 0.3510*0.2806; SwMinTypes[it] = 0.2806
    KsTypes[it]      = 4.69*m_per_s_by_m_per_d;
    lambdaBCtypes[it] = nVGtypes[it]-1.0; pdBCtypes[it] = 1.0/alphaVGtypes[it];
    if not useHet:
        it=1;
        thetaStypes[it]  = thetaStypes[0]; thetaRtypes[it] = thetaRtypes[0];
        KsTypes[it]      = KsTypes[0]; SwMinTypes[it] = SwMinTypes[0]
        nVGtypes[it]     = nVGtypes[0];
        mVGtypes[it]     = 1.0-1.0/nVGtypes[it];
        alphaVGtypes[it] = alphaVGtypes[0]
        lambdaBCtypes[it] = nVGtypes[it]-1.0; pdBCtypes[it] = 1.0/alphaVGtypes[it];
        
else:
    #Helmig example 1 ch 6
    it=0;
    thetaStypes[it]  = 0.40; thetaRtypes[it] = 0.40*0.09; SwMinTypes[it] = 0.09
    KsTypes[it]      = 7.0e-4 #m/s
    lambdaBCtypes[it]= 2.7; pdBCtypes[it] = 0.0772;
    nVGtypes[it] = lambdaBCtypes[it]+1.0; mVGtypes[it]    = 1.0-1.0/nVGtypes[it];
    alphaVGtypes[it] = 1.0/pdBCtypes[it];
    #
    it=1;
    thetaStypes[it]  = 0.39; thetaRtypes[it] = 0.39*0.12; SwMinTypes[it] = 0.12
    KsTypes[it]      = 7.0e-5 #m/s
    lambdaBCtypes[it]= 2.0; pdBCtypes[it] = 0.211;
    nVGtypes[it] = lambdaBCtypes[it]+1.0; mVGtypes[it]    = 1.0-1.0/nVGtypes[it];
    alphaVGtypes[it] = 1.0/pdBCtypes[it];

    if not useHet:
        it=1;
        thetaStypes[it]  = thetaStypes[0]; thetaRtypes[it] = thetaRtypes[0];
        KsTypes[it]      = KsTypes[0] ; SwMinTypes[it] = SwMinTypes[0]
        lambdaBCtypes[it]= lambdaBCtypes[0]; pdBCtypes[0] = pdBCtypes[0];
        nVGtypes[it] = lambdaBCtypes[it]+1.0; mVGtypes[it]    = 1.0-1.0/nVGtypes[it];
        alphaVGtypes[it] = 1.0/pdBCtypes[it];
    
#fluid properties from Helmig, 25 C
#water
viscosity_w     = 1.0e-3  #kg/(m*s)
density_w       = 998.2   #kg/m^3
#TCE
viscosity_n     = 9.0e-4
density_n       = 1460.0

gravity       = 9.8     #m/s^2

g    = [0.0,-gravity,0.0]    #

coefficients = TwophaseDarcy_fc_ff(g=g, 
                                   rhon=density_n,
                                   rhow=density_w,
                                   mun    = viscosity_n,
                                   muw    = viscosity_w,
                                   Ksw=KsTypes[0],
                                   psk_model=model,
                                   vg_alpha = alphaVGtypes[0],
                                   vg_m  = mVGtypes[0],
                                   bc_pd  = pdBCtypes[0], 
                                   bc_lambda = lambdaBCtypes[0],
                                   omega  = thetaStypes[0],
                                   Sw_max = SwMaxTypes[0],
                                   Sw_min = SwMinTypes[0])
if useHet:
    coefficients.setMaterialTypes(Ksw_types = KsTypes,
                                  omega_types  = thetaStypes,
                                  Sw_max_types = SwMaxTypes,
                                  Sw_min_types = SwMinTypes,
                                  bc_lambda_types = lambdaBCtypes,
                                  bc_pd_types = pdBCtypes,
                                  vg_alpha_types = alphaVGtypes,
                                  vg_m_types = mVGtypes)
#     coefficients = TwophaseFFDarcyFCHetV2(Ksw_types=KsTypes,
#                                           rhon=density_n,
#                                           rhow=density_w,
#                                           g=g,
#                                           mvg_alpha_types = alphaVGtypes,
#                                           bc_lambda_types = lambdaBCtypes,
#                                           bc_pd_types = pdBCtypes,
#                                           mvg_n_types = nVGtypes,
#                                           mvg_m_types = mVGtypes,
#                                           thetaS_types = thetaStypes,
#                                           thetaR_types = thetaRtypes,
#                                           mun   = viscosity_n,
#                                           muw   = viscosity_w,
#                                           model = model)

        

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

#napl saturation at top
SnTop = 0.4#0.4
#if not useV2:
#    #convert to effective saturation
#    SeTop = (thetaStypes[0]*(1.0-SnTop)-thetaRtypes[0])/(thetaStypes[0]-thetaRtypes[0])
#    SnTop = 1.0-SeTop
#say napl head at top around ?4 cm
naplHeadTop = 0.3#0.3#0.04# 0.0#m

SeTop = (thetaStypes[0]*(1.0-SnTop)-thetaRtypes[0])/(thetaStypes[0]-thetaRtypes[0])

if model == 'VGM':
    pcTop = pcVGM(SeTop,alphaVGtypes[0],nVGtypes[0],mVGtypes[0])
else:
    pcTop = pcBCB(SeTop,pdBCtypes[0],lambdaBCtypes[0])
    

aqHeadTop = naplHeadTop - pcTop

#Dirichlet boundary conditions
eps = 1.0e-6
def getDBC_sw(x,tag):
    if x[0] <= eps:
        return lambda x,t: 1.0 #fully saturated
    if x[0] >= L[0]-eps:
        return lambda x,t: 1.0
    if (x[1] >= L[1]-eps and
        slit_left <= x[0] and x[0] <= slit_right):
        return lambda x,t: 1.0-SnTop
def getDBC_psiw(x,tag):
    if x[0] <= eps:
        return lambda x,t: L[1]-x[1]
    if x[0] >= L[0]-eps:
        return lambda x,t: L[1]-x[1]
    if (x[1] >= L[1]-eps and
        slit_left <= x[0] and x[0] <= slit_right):
        return lambda x,t: aqHeadTop

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class HydroIC_sw:
    def uOfXT(self,x,t):
        return 1.0
class HydroIC2_sw:
    def uOfXT(self,x,t):
        if (x[1] >= L[1]-eps and
            slit_left <= x[0] and x[0] <= slit_right):
            return 1.0-SnTop
        return 1.0
class HydroIC_psiw:
    def uOfXT(self,x,t):
        return L[1]-x[1]

initialConditions = {0:HydroIC2_sw(),1:HydroIC_psiw()}

#no flow flux boundaries by default
fluxBoundaryConditions = {}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

hour = 3600.0
if useForsyth:
    T = 1.0*hour#24.0*hour
else:
    T = 1.0*hour
