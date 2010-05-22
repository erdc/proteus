from pyadh import *
from pyadh.default_p import *

nd = 2

L=(1.0,5.0,1.0)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
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
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
                                        0.0])
dimensionless_alpha    = mvg_alpha*lengthScale

#useHeadForm = False
useHeadForm = True

if useHeadForm:
    coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                              gravity=dimensionless_gravity,
                                                              density=dimensionless_density,
                                                              thetaS=thetaS,
                                                              thetaR=thetaR,
                                                              alpha= dimensionless_alpha,
                                                              n = mvg_n,
                                                              m = mvg_m)
else:
    coefficients = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                             gravity=dimensionless_gravity,
                                                             density=dimensionless_density,
                                                             thetaS=thetaS,
                                                             thetaR=thetaR,
                                                             alpha= dimensionless_alpha,
                                                             n = mvg_n,
                                                             m = mvg_m)

pondingPressureHead= -0.5 #m
rechargeRate =  -1.0e-2#- means into domain m/d
outFlowSat = 0.9999
outFlowHead = 1.0

eps=1.0e-8
def getDBC_2D_Richards_Shock(x,flag):
    if useHeadForm:
        #if x[1] == L[1]:
        #    return lambda x,t: pondingPressure
        if abs(x[1]-0.0) < eps:
            return lambda x,t: outFlowHead
    else:
        if abs(x[1]-0.0) < eps:
            return lambda x,t: outFlowSat
       
dirichletConditions = {0:getDBC_2D_Richards_Shock}

def getRecharge_2D_Top(x,flag):
    if abs(x[1]-L[1]) < eps:
        return lambda x,t: rechargeRate
    elif abs(x[0]-0.0) < eps or abs(x[0]-L[0]) < eps:
        return lambda x,t: 0.0
def getRecharge_2D_Diff(x,flag):
    return lambda x,t: 0.0
def getDummyFlux(x,flag):
    pass

class HydroIC_2D_Richards:
    def uOfXT(self,x,t):
        return x[1]*dimensionless_gravity[1]*dimensionless_density
        
class ConstIC_2D_Richards:
    def uOfXT(self,x,t):
        if useHeadForm or useTotalHeadForm:
            return outFlowHead
        else:
            return outFlowSat
            
#initialConditions  = {0:HydroIC_2D_Richards()}
initialConditions  = {0:ConstIC_2D_Richards()}

fluxBoundaryConditions = {0:'setFlow'}


advectiveFluxBoundaryConditions =  {0:getRecharge_2D_Top}

diffusiveFluxBoundaryConditions = {0:{0:getDummyFlux}}

T = 0.25/timeScale
