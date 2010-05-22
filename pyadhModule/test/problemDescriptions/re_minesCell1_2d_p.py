from pyadh import *
from pyadh.default_p import *

nd = 2

polyfile = "minesCell1"
top   = 1.8e-1 #m
right = 1.8e-1 #m
L=(right,top,1.0)

analyticalSolution = None

# # # # # # # # media types # # # # # # # # # # 
#     Type 1: fine
#             K_s     = 0.033e-2   m/s  
#             \theta_s= 0.426
#             \theta_r= 0.426*0.081
#             \alpha  = 2.6       1/m
#             m       = 10.67?  m = 1 - 1/n
#             n       = 10.67
#     Type 2: coarse
#             K_s     = 0.209e-2  m/s     
#             \theta_s= 0.411
#             \theta_r= 0.411*0.060
#             \alpha  = 6.4      1/m
#             m       = 11.12?  m = 1 - 1/n
#             n       = 11.12
     

nMediaTypes  = 2
alphaVGtypes = numpy.zeros((nMediaTypes,),'d')
nVGtypes     = numpy.zeros((nMediaTypes,),'d')
thetaStypes  = numpy.zeros((nMediaTypes,),'d')
thetaRtypes  = numpy.zeros((nMediaTypes,),'d')
KsTypes      = numpy.zeros((nMediaTypes,),'d')

#Type ids are base 1

#Type 1
it = 0;
alphaVGtypes[it] = 2.600; nVGtypes[it]    = 10.67;
thetaStypes[it]  = 0.426; thetaRtypes[it] = 0.4260*0.081;
KsTypes[it]      = 0.033e-2;

#Type 2
it = 1;
alphaVGtypes[it] = 6.4 ;  nVGtypes[it]    = 11.12;
thetaStypes[it]  = 0.411; thetaRtypes[it] = 0.411*0.060;
KsTypes[it]      = 0.209e-2;

#set dimensions for types
viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
m_per_s_by_m_per_d = 1.1574074e-5
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)

dimensionless_conductivity  = numpy.zeros(KsTypes.shape,'d')
dimensionless_alpha         = numpy.zeros(alphaVGtypes.shape,'d')
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                      -1.0,
                                      0.0])
for i in range(len(KsTypes.flat)):
    permeability = (KsTypes[i]*viscosity/(gravity*density))  #m^2
    dimensionless_conductivity[i] = KsTypes[i]
    print 'Ks',dimensionless_conductivity
    dimensionless_alpha[i] = alphaVGtypes[i]*lengthScale
    


coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2(KsTypes,
                                                                    nVGtypes,
                                                                    alphaVGtypes,
                                                                    thetaRtypes,
                                                                    thetaStypes-thetaRtypes,
                                                                    dimensionless_gravity,
                                                                    dimensionless_density)

initialBottomPressureHead = L[1]
finalBottomPressureHead   = -0.8 #m
initialTime = 0.0
drainageIncrement = -0.03  #amount to change head
drainageInterval  = 3600.0#how often to drop head
waterTableHeight = L[1]

def drainPressureHead(t):
    nIntervals = math.floor((t-initialTime)/drainageInterval)
    head = nIntervals*drainageIncrement + initialBottomPressureHead
    return head
    #return initialBottomPressureHead

eps = 1.0e-8
def getDBC_drainage(x,tag):
    if tag == 1 and x[1] <= eps:#bottom
        assert abs(x[1]) < eps, "problem tag= %s x=%s L=%s " % (tag,x,L)
        return lambda x,t: drainPressureHead(t)


def getDummyFlux(x,tag):
    pass

dirichletConditions = {0:getDBC_drainage}

class HydroIC_2D_Richards:
    def uOfXT(self,x,t):
        return (x[1]-waterTableHeight)*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:HydroIC_2D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T= 24.0*3600.0/timeScale
