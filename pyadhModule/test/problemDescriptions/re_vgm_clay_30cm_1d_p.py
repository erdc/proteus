from pyadh import *
from pyadh.default_p import *

nd = 1

L=(0.3,1.0,1.0)

analyticalSolutions = None

viscosity     = 8.9e-4  #kg/(m*s) 
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
#mwf check permeability m/d value, can't remember off hand
permeability  = (7.92*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
print 'perm',permeability
thetaS        = 0.368   #-
thetaR        = 0.102   #-
mvg_alpha     = 3.35    #1/m
mvg_n         = 2.0
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

coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                          gravity=dimensionless_gravity,
                                                          density=dimensionless_density,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          alpha= dimensionless_alpha,
                                                          n = mvg_n,
                                                          m = mvg_m,
                                                          beta=beta)

pondingPressure=-0.75
bottomPressure = -10.0
def getDBC_Richards_Shock(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: bottomPressure

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        return -10.0
class HydrostaticIC_Richards:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock(x,0)
        if f:
            return f(x,t)
        return bottomPressure + x[0]*dimensionless_gravity[0]*dimensionless_density

#initialConditions  = {0:ShockIC_Richards()}
initialConditions  = {0:HydrostaticIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

gravityMday2=gravity*(60*60*24)**2
T = 0.25
