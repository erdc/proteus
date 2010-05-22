from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_modelParams import *
import hole2d#hole2d 
name = "twp_darcy_split_press_hole_sgs"
nd = 2
top = 2.0 #m
right = 4.0 #m
domain = [[0.0,right],[0.0,top]]
slit_left=1.9
slit_right=2.1
holeX=[1.5,2.5]
holeY=[1.5,1.75]
hole = [holeX, holeY]
sourceXs = [slit_left-0.1,slit_left,slit_right,slit_right+0.1]
polyfile='hole2d'
L,boundaryTags=hole2d.genPoly(domain=domain,
                              hole=hole, 
                              sourceXs=sourceXs, 
                              polyfileBase=polyfile)

useHet = False#True#False#True
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
    KswTypes[index]     = porousMediumDatabase[medium]['Ksw']
    SwMaxTypes[index]   = porousMediumDatabase[medium]['sw_max']
    SwMinTypes[index]   = porousMediumDatabase[medium]['sw_min']
    lambdaBCtypes[index]= porousMediumDatabase[medium]['bc_lambda']
    pdBCtypes[index]    = porousMediumDatabase[medium]['bc_pd']

g    = [0.0,-gmag]    #
phase = 'potential'
coefficients = TwophaseDarcy_split_pressure(g=g, 
                                            rhon=rhon,
                                            rhow=rhow,
                                            mun    = mun,
                                            muw    = muw,
                                            Ksw=KswTypes[0],
                                            psk_model=model,
                                            vg_alpha = alphaVGtypes[0],
                                            vg_m  = mVGtypes[0],
                                            bc_pd  = pdBCtypes[0], 
                                            bc_lambda = lambdaBCtypes[0],
                                            omega  = thetaStypes[0],
                                            Sw_max = SwMaxTypes[0],
                                            Sw_min = SwMinTypes[0],
                                            diagonalHet = sd,
                                            sd = sd)
if useHet:
    coefficients.setMaterialTypes(Ksw_types = KswTypes,
                                  omega_types  = thetaStypes,
                                  Sw_max_types = SwMaxTypes,
                                  Sw_min_types = SwMinTypes,
                                  bc_lambda_types = lambdaBCtypes,
                                  bc_pd_types = pdBCtypes,
                                  vg_alpha_types = alphaVGtypes,
                                  vg_m_types = mVGtypes)
        
#convenience functions for background psk relations
se2pc = None
pc2se = None
if model == 'VGM':
    se2pc = lambda se : pcVGM(se,mvg_alpha,mvg_n,mvg_m)
    pc2se = lambda pc : seVGM(pc,mvg_alpha,mvg_n,mvg_m)
elif model == 'BCB':
    se2pc = lambda se : pcBCB(se,bc_pd,bc_lambda)
    pc2se = lambda pc : seBCB(pc,bc_pd,bc_lambda)
#


psin_top   = 0.0
psiw_top   = -1.98
pc_top  = psin_top-psiw_top
Se_top =  pc2se(pc_top)#0.05#1.0e-3                  # effective saturation bottom
Sw_top = Se_top*(sw_max-sw_min) + sw_min

waterTable =  0.02#-0.1#-1.0#elevation of water table 
psiTable= 0.0     #psi_w at water table
#psi_w at z=0
psi_bottom = psiTable + g[1]/gmag*(0.0-waterTable)#waterTable - 0.0
psinTable  = 0.0
psin_bottom= psinTable + rhon/rhow*g[1]/gmag*(0.0-waterTable)#waterTable - 0.0
pc_bottom  = psin_bottom-psi_bottom
Se_bottom =  pc2se(pc_bottom)#0.05#1.0e-3                  # effective saturation bottom
Sw_bottom= Se_bottom*(sw_max-sw_min) + sw_min

hour = 3600.0 #[s]
day  = 24.0*hour
tnList = [0.0,0.1,0.01*day,0.02*day,0.1*day,0.8*day,1.05*day,1.055*day,1.2*day,3.0*day]
T = tnList[-1]                      # time [s]

def getPsiw_top(x,t):
    data = [0.0*day, psiw_top, 0.01*day, 0.02, 0.02*day, 0.02, 0.1*day, -4.0, 0.8*day, -4.0, 1.05*day, 0.02, 1.055*day, 0.02, 1.2*day, -4.0, 3.0*day, -4.0]
    index=2
    while t>data[index] and index < len(data)-2:
        index+=2 
    p0 = data[index-1]
    p1 = data[index+1]
    t0 = data[index-2]
    t1 = data[index]
    li = p0+(t-t0)/(t1-t0)*(p1-p0)
    return li

#Dirichlet boundary conditions
eps = 1.0e-6

def getDBC_psiw(x,tag):
    if x[1] >= L[1]-eps and x[0] > eps and x[0] < L[0]-eps:
       return getPsiw_top
    if x[1] <= eps:
        return lambda x,t: psi_bottom

dirichletConditions = {0:getDBC_psiw}
class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return psiTable + g[1]/gmag*(x[1]-waterTable)


class ConstantIC:
    def __init__(self,val):
        self.val = val
    def uOfXT(self,x,t):
        return self.val

initialConditions  = {0:psiw_IC()}

#no flow flux boundaries by default?
fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

