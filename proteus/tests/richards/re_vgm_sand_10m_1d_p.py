from proteus import *
from proteus.default_p import *
from proteus.richards import Richards

nd = 1

L=(10.0,1.0,1.0)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
permeability  = (5.04*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
print 'perm',permeability
thetaS        = 0.301   #-
thetaR        = 0.093   #-
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                       0.0,
                                       0.0])
dimensionless_alpha    = mvg_alpha*lengthScale
satRichards = False
optRichards = True
nMediaTypes  = 1
alphaVGtypes = numpy.zeros((nMediaTypes,),'d')
nVGtypes     = numpy.zeros((nMediaTypes,),'d')
thetaStypes  = numpy.zeros((nMediaTypes,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes,),'d')
thetaSRtypes = numpy.zeros((nMediaTypes,),'d')
KsTypes      = numpy.zeros((nMediaTypes,3),'d')

for i in range(nMediaTypes):
    alphaVGtypes[i] = mvg_alpha
    nVGtypes[i]     = mvg_n
    thetaStypes[i]  = thetaS
    thetaRtypes[i]  = thetaR
    thetaSRtypes[i] = thetaStypes[i] - thetaRtypes[i]
    KsTypes[i,:]    = [dimensionless_conductivity,dimensionless_conductivity,dimensionless_conductivity]#m/d?

useSeepageFace = True
if optRichards:
    LevelModelType = Richards.LevelModel
    coefficients = Richards.Coefficients(nd,
                                         KsTypes,
                                         nVGtypes,
                                         alphaVGtypes,
                                         thetaRtypes,
                                         thetaSRtypes,
                                         gravity=dimensionless_gravity,
                                         density=dimensionless_density,
                                         beta=0.0001,
                                         diagonal_conductivity=True)

elif satRichards:
    coefficients = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                             gravity=dimensionless_gravity,
                                                             density=dimensionless_density,
                                                             thetaS=thetaS,
                                                             thetaR=thetaR,
                                                             alpha= dimensionless_alpha,
                                                             n = mvg_n,
                                                             m = mvg_m)
else:
    coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                              gravity=dimensionless_gravity,
                                                              density=dimensionless_density,
                                                              thetaS=thetaS,
                                                              thetaR=thetaR,
                                                              alpha= dimensionless_alpha,
                                                              n = mvg_n,
                                                              m = mvg_m,
                                                              beta=beta)

pondingPressure=0.1#-0.1
bottomPressure = 0.0#0.0
#pondingPressure=-0.1
#bottomPressure = -10.0
pondingSaturation = 0.9
waterTableSaturation = 0.9
initialSaturation = 0.01
#pondingPressure=-0.1
if satRichards:
    def getDBC_Richards_Shock(x,flag):
        if x[0] == L[0]:
            return lambda x,t: pondingSaturation
        if x[0] == 0.0:
            return lambda x,t: waterTableSaturation
else:
    def getDBC_Richards_Shock(x,flag):
        if x[0] == L[0]:
            return lambda x,t: pondingPressure
        if x[0] == 0.0:
            return lambda x,t: bottomPressure
   
dirichletConditions = {0:getDBC_Richards_Shock}

if satRichards:
    class ShockIC_Richards:
        def uOfXT(self,x,t):
            f = getDBC_Richards_Shock(x,0)
            if f:
                return f(x,t)
            return initialSaturation
else:
    class ShockIC_Richards:
        def uOfXT(self,x,t):
            f = getDBC_Richards_Shock(x,0)
            if f:
                return f(x,t)
            return bottomPressure + x[0]*dimensionless_gravity[0]*dimensionless_density

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

def flux(x,flag):
    return None
advectiveFluxBoundaryConditions =  {0:flux}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5/timeScale
