from pyadh import *
from pyadh.default_p import *

nd = 1

L=(10.0,1.0,1.0)
#L=(1.0,1.0,1.0)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
#viscosity     = 1.0 * 10.0#cp -> kg/(m*s)
#density       = 0.433 * 2306.6587   #psi/ft -> kg/m^3
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
#beta          = 0.0
m_per_s_by_m_per_d = 1.1574074e-5
permeability  = 18.0*0.987e-12 #Darcys -> m^2
lengthScale=1.0
timeScale=1.0
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                      0.0,
                                      0.0])
psiD = (2.0*math.sqrt(0.4/18.0)*6894.7573)/(density*gravity)
#psiD = 5.0
print 'psiD',psiD
coefficients = ConservativeHeadRichardsJLeverett(phi_block=numpy.array([0.4]),
                                                 psiD_block=numpy.array([psiD]),
                                                 ns_block=numpy.array([1.0/1.5]),
                                                 #ns_block=numpy.array([1.0]),
                                                 #nk_block=numpy.array([1.0]),
                                                 nk_block=numpy.array([1.5]),
                                                 S_wirr_block=numpy.array([0.0]),
                                                 S_nwr_block=numpy.array([0.0]),
                                                 kr0_block=numpy.array([dimensionless_conductivity]),
                                                 gravity=dimensionless_gravity,
                                                 density=dimensionless_density,
                                                 beta=beta)

pondingPressure= 0.1#-0.01#-1.0
bottomPressure = -psiD#0.0#
#pondingPressure=-0.1
#bottomPressure = -10.0
pondingSaturation = 0.9
waterTableSaturation = 0.9
initialSaturation = 0.01
#pondingPressure=-0.1

def getDBC_Richards_Shock(x,flag):
#     if x[0] == L[0]:
#         return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: bottomPressure

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock(x,0)
        if f:
            return f(x,t)
        return bottomPressure# + x[0]*dimensionless_gravity[0]

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

def waterFlowRate(x,t):
    if t < 2.0:
        return -3.0 * 0.028316847 #ft^3/day -> m^3/day
    else:
        return 0.0

def getAFBC(x,flag):
#    pass
   if x[0] == L[0]:
       return lambda x,t: 0.0
    
def getDFBC(x,flag):
#    pass
   if x[0] == L[0]:
       return waterFlowRate

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}


#T = 0.5/timeScale
#cek testing
#T = 0.0025/timeScale
#mwf testing
#T = 0.25
T = 4.0
