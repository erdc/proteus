from pyadh import *
from pyadh.default_p import *
from twp_darcy_modelParams import *
import lens2d 
name = "re_test"

#default setup
nd = 2
top   = 2.0 #m
right = 4.0 #m
domain = [[0.0,right],[0.0,top]]
slit_left = 1.9
slit_right= 2.1
#lensX = [1.5,2.5]; lensY=[1.25,1.5]
lensX=[1.25,2.25];
lensY=[1.25,1.5];
lens  = [lensX,lensY]
sourceXs = [slit_left-0.1,slit_left,slit_right,slit_right+0.1]

polyfile = 'lens2d'
L,boundaryTags = lens2d.genPoly(domain=domain,
                                lens=lens,
                                sourceXs=sourceXs,
                                polyfileBase=polyfile)


useHet = True#False
model = 'VGM'#'BCB' #'VGM'


nMediaTypes  = 2
alphaVGtypes = numpy.zeros((nMediaTypes,),'d')
nVGtypes     = numpy.zeros((nMediaTypes,),'d')
mVGtypes     = numpy.zeros((nMediaTypes,),'d')
thetaStypes  = numpy.zeros((nMediaTypes,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes,),'d')
KswTypes     = numpy.zeros((nMediaTypes,),'d')
SwMaxTypes   = numpy.ones((nMediaTypes,),'d')
SwMinTypes   = numpy.zeros((nMediaTypes,),'d')
lambdaBCtypes = numpy.ones((nMediaTypes,),'d')
pdBCtypes    = numpy.ones((nMediaTypes,),'d')

#put everything in [s]
m_per_s_by_m_per_d = 1.1574074e-5

if useHet:
    heterogeneityTypes = {0:'sand',1:'clay-3'}
else:
    heterogeneityTypes = {0:PorousMedia,1:PorousMedia}

for index,medium in heterogeneityTypes.iteritems():
    alphaVGtypes[index] = porousMediumDatabase[medium]['mvg_alpha']
    nVGtypes[index]     = porousMediumDatabase[medium]['mvg_n']
    mVGtypes[index]     = porousMediumDatabase[medium]['mvg_m']
    thetaStypes[index]  = porousMediumDatabase[medium]['sw_max']*porousMediumDatabase[medium]['omega']
    thetaRtypes[index]  = porousMediumDatabase[medium]['sw_min']*porousMediumDatabase[medium]['omega']
    KswTypes[index]     = porousMediumDatabase[medium]['Ks']#Ksw is m/s Ks is m/d
    SwMaxTypes[index]   = porousMediumDatabase[medium]['sw_max']
    SwMinTypes[index]   = porousMediumDatabase[medium]['sw_min']
    lambdaBCtypes[index]= porousMediumDatabase[medium]['bc_lambda']
    pdBCtypes[index]    = porousMediumDatabase[medium]['bc_pd']

g    = numpy.array([0.0,-gmag,0.])    #
g_nondim = g/gmag
rhow_nondim = rhow/rhow
thetaSRtypes = thetaStypes-thetaRtypes
coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2withUpwind(nd,
                                                                              KswTypes,
                                                                              nVGtypes,
                                                                              alphaVGtypes,
                                                                              thetaRtypes,
                                                                              thetaSRtypes,
                                                                              g_nondim,
                                                                              rhow_nondim,
                                                                              beta_w,
                                                                              upwindFlag=0,
                                                                              sd = True)
# coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2(KswTypes,
#                                                                     nVGtypes,
#                                                                     alphaVGtypes,
#                                                                     thetaRtypes,
#                                                                     thetaSRtypes,
#                                                                     g_nondim,
#                                                                     rhow_nondim,
#                                                                     beta_w)

psi_top = 0.1
psi_table = 0.0
water_table = 0.0
psi_bottom = psi_table + g[1]/gmag*(0.0-water_table)
hour = 1.0/24 #[d]
T = 1.*hour#24.0e-1*hour #sand 12*hour, clay 24*hour                          # time [s]

#Dirichlet boundary conditions
eps = 1.0e-6
def getDBC_psiw(x,tag):
    if (x[1] >= L[1]-eps and
        slit_left <= x[0] and x[0] <= slit_right and tag==4):
        return lambda x,t: psi_top
    if x[1] <= eps:
        return lambda x,t: psi_bottom

dirichletConditions = {0:getDBC_psiw}

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psi_table + g[1]/gmag*(x[1]-water_table)


initialConditions  = {0:psiw_IC()}
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

