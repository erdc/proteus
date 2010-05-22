from pyadh import *
from pyadh.default_p import *

nd = 1

L=(0.3,1.0,1.0)

analyticalSolutions = None

viscosity     = 8.9e-4  #kg/(m*s) 
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
m_per_s_by_m_per_d = 1.1574074e-5
#mwf check permeability m/d value, can't remember off hand
permeability  = (7.92*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
print 'perm',permeability
thetaS        = [0.368]   #-
thetaR        = [0.102]   #-
thetaSR       = [thetaS[i]-thetaR[i] for i in range(1)] 
mvg_alpha     = [3.35]    #1/m
mvg_n         = [2.0]
mvg_m         = [1.0 - 1.0/mvg_n[0]]
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
beta          = 0.0001
#make non-dimensional
dimensionless_conductivity  = [(timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d]
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                       0.0,
                                       0.0])
dimensionless_alpha    = [mvg_alpha[i]*lengthScale for i in range(1)]
upwindFlag = 1
coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2withUpwind(nd,
                                                                              Ks_block=numpy.array(dimensionless_conductivity),
                                                                              vgm_n_block=numpy.array(mvg_n),
                                                                              vgm_alpha_block=numpy.array(mvg_alpha),
                                                                              thetaR_block=numpy.array(thetaR),
                                                                              thetaSR_block=numpy.array(thetaSR),
                                                                              gravity=dimensionless_gravity,
                                                                              density=dimensionless_density,
                                                                              beta=beta,
                                                                              upwindFlag=upwindFlag)


pondingPressure=-0.75

def getDBC_Richards_Shock(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: -10.0

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        return -10.0

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

gravityMday2=gravity*(60*60*24)**2
T = 0.25
