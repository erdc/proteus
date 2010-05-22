from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_modelParams import *
#import lensAsHole2d as hole2d#hole2d 
import hole2d
name = "twp_darcy_fc_hole"

nd = 2
top = 2.0 #m
right = 4.0 #m
mydomain = [[0.0,right],[0.0,top]]
#slit_left=1.9
#slit_right=2.1
holeX=[1.5,2.5]
holeY=[1.5,1.75]
hole = [holeX, holeY]
#import pdb
#pdb.set_trace()
#sourceXs = [slit_left-0.1,slit_left,slit_right,slit_right+0.1]
polyfile='hole2d'
L,boundaryTags=hole2d.genPoly(domain=mydomain,
                              hole=hole, 
#                              sourceXs=sourceXs, 
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

gmag  = 9.8         # magnitude of gravity             (Earth 9.8 m/s^2) 
g    = [0.0,-gmag]    #
coefficients = TwophaseDarcy_fc(g=g, 
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
                                density_w_parameters = density_w_exponential,
                                density_n_parameters = density_n_ideal,
                                diagonalHet = True,
                                sd = True)
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
psiw_top   = -0.02#-1.98
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
#tnList = [0.0,0.1,0.32*day,0.34*day]
#presList = [-1.98, -1.98, 0.02, 0.02]
tnList = [0.0,0.1,0.32*day,0.34*day,0.54*day,0.63*day,1.0*day,1.22*day,1.26*day,1.54*day,1.63*day,3.0*day]
presList = [-1.98, -1.98, 0.02, 0.02,-2.0, -2.0, -1.98, 0.02, 0.02, -2.0, -2.0, -1.98]
#tnList = [0.0,0.1,0.32*day]
#presList = [psiw_top,psiw_top,psiw_top]
#T = tnList[3]                      # time [s]
T = tnList[-1]                      # time [s]

def getPsiw_top(x,t):
    if t <= 0.0:
        return psiw_top
    #data = [0.0*day, -1.98, 0.1,-1.98, 0.32*day, 0.02]
    #data = [0.0*day, -1.98, 0.32*day, 0.02, 0.34*day, 0.02, 0.54*day, -2.0, 0.63*day, -2.0, 1.0*day, -1.98, 1.22*day, 0.02, 1.26*day, 0.02, 1.54*day, -2.0, 1.63*day, -2.0, 3.0*day, -1.98]
    #import pdb
    #pdb.set_trace()
    index=1
    while index < len(tnList)-1 and t>tnList[index]:
        index+=1
    p0 = presList[index-1]
    p1 = presList[index]
    t0 = tnList[index-1]
    t1 = tnList[index]
    li = p0+(t-t0)/(t1-t0)*(p1-p0)
    #print "getPsiw_top t= %s p0=%s p1=%s t0=%s t1=%s li=%s " % (t,p0,p1,t0,t1,li)
    return li

def getSw_top(x,t):
    pw = getPsiw_top(x,t)
    pc = psin_top-pw
    se = pc2se(pc)
    return se*(sw_max-sw_min) + sw_min



#Dirichlet boundary conditions
eps = 1.0e-6
def getDBC_sw(x,tag):
    if x[1] >= L[1]-eps and x[0] > eps and x[0] < L[0]-eps:
        return getSw_top
    if x[1] <= eps:
        return lambda x,t: Sw_bottom

def getDBC_psiw(x,tag):
    if x[1] >= L[1]-eps and x[0] > eps and x[0] < L[0]-eps:
       return getPsiw_top
    if x[1] <= eps:
        return lambda x,t: psi_bottom

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}
class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= L[1]-eps:
            return getSw_top(x,t)
        psic = pc_bottom + (rhon/rhow - 1.0)*g[1]/gmag*x[1]
        if model == 'VGM':
            se = seVGM(psic,alphaVGtypes[0],nVGtypes[0],mVGtypes[0])
        else:
            se = seBCB(psic,pdBCtypes[0],lambdaBCtypes[0])
        sw = (se*(thetaStypes[0]-thetaRtypes[0]) + thetaRtypes[0])/thetaStypes[0]
        return sw
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

initialConditions  = {0:sw_IC(),1:psiw_IC()}
#initialConditions  = {0:ConstantIC(Sw_bottom),1:ConstantIC(psi_bottom)}

#no flow flux boundaries by default?
fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}
#fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

