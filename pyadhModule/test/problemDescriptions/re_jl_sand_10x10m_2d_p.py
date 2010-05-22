from pyadh import *
from pyadh.default_p import *

nd = 2

L=(10.0,10.0,1.0)

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
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
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
                                                          beta = beta)

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
dimensionless_gravity  = numpy.array([0.0,
                                      -1.0,
                                      0.0])
psiD = (2.0*math.sqrt(0.4/18.0)*6894.7573)/(density*gravity)
#psiD = 5.0
print 'psiD',psiD
coefficients = ConservativeHeadRichardsJLeverett(phi_block=numpy.array([0.4]),
                                                 psiD_block=numpy.array([psiD]),
                                                 ns_block=numpy.array([1.0/1.5]),
                                                 nk_block=numpy.array([1.5]),
                                                 S_wirr_block=numpy.array([0.02]),
                                                 S_nwr_block=numpy.array([0.0]),
                                                 kr0_block=numpy.array([dimensionless_conductivity]),
                                                 gravity=dimensionless_gravity,
                                                 density=dimensionless_density,
                                                 beta=beta)
# coefficients = ConservativeHeadRichardsJLeverettAni(phi_block=numpy.array([0.4]),
#                                                  psiD_block=numpy.array([psiD]),
#                                                  ns_block=numpy.array([1.0/1.5]),
#                                                  nk_block=numpy.array([1.5]),
#                                                  S_wirr_block=numpy.array([0.02]),
#                                                  S_nwr_block=numpy.array([0.0]),
#                                                  kr0x_block=numpy.array([dimensionless_conductivity]),
#                                                  kr0y_block=numpy.array([0.01*dimensionless_conductivity]),
#                                                  kr0z_block=numpy.array([dimensionless_conductivity]),
#                                                  gravity=dimensionless_gravity,
#                                                  density=dimensionless_density,
#                                                  beta=beta)
thetaS        = 0.4   #-
thetaR        = 0.02*0.4   #-
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
dimensionless_alpha    = mvg_alpha*lengthScale
psiD = 0.0
coefficients = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                          gravity=dimensionless_gravity,
                                                          density=dimensionless_density,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          alpha= dimensionless_alpha,
                                                          n = mvg_n,
                                                          m = mvg_m,
                                                          beta = beta)

bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
dimensionless_pd    = bc_pd/lengthScale
coefficients = ConservativeHeadRichardsBrooksCoreyBurdine(hydraulicConductivity=dimensionless_conductivity,
                                                          gravity=dimensionless_gravity,
                                                          density=dimensionless_density,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          lambdab = bc_lambda,
                                                          pd = dimensionless_pd,
                                                          beta=beta)
coefficients = ConservativeHeadRichardsJLeverettAni(phi_block=numpy.array([0.4]),
                                                 psiD_block=numpy.array([dimensionless_pd]),
                                                 ns_block=numpy.array([bc_lambda]),
                                                 nk_block=numpy.array([-1.0e12]),#not used
                                                 S_wirr_block=numpy.array([0.02]),
                                                 S_nwr_block=numpy.array([0.0]),
                                                 kr0x_block=numpy.array([dimensionless_conductivity]),
                                                 #kr0y_block=numpy.array([0.01*dimensionless_conductivity]),
                                                 kr0y_block=numpy.array([dimensionless_conductivity]),
                                                 kr0z_block=numpy.array([dimensionless_conductivity]),
                                                 gravity=dimensionless_gravity,
                                                 density=dimensionless_density,
                                                 beta=beta)
pondingPressure= 0.1#-0.01#-1.0
bottomPressure = -psiD
#pondingPressure=-0.1
#bottomPressure = -10.0
pondingSaturation = 0.9
waterTableSaturation = 0.9
initialSaturation = 0.01
#pondingPressure=-0.1

def getDBC_Richards_Shock(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: bottomPressure + x[1]*dimensionless_gravity[1]
    else:
        return None
#     if x[1] == L[1]:
#         if x[0] >= L[0]/3.0 and x[0] <= L[0]*2.0/3.0:
#             return lambda x,t: 0.1
#         else:
#             return None

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock(x,0)
        if f:
            return f(x,t)
        return bottomPressure + x[1]*dimensionless_gravity[1]

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

def waterFlowRate(x,t):
    if t < 1.0:
        return -25.0/(L[0]/3.0)
        #return -3.0 * 0.028316847 / (L[0]/3.0) #ft^3/day -> m^3/day -> m/day
    else:
        return 0.0

def getAFBC(x,flag):
    if x[0] in [0.0,L[0]]:
        pass
    if x[1] == L[1]:
        return lambda x,t: 0.0
    if x[1] == 0.0:
        return lambda x,t: 0.0
    
def getDFBC(x,flag):
    if x[0] in [0.0,L[0]]:
        pass
    if x[1] == L[1]:
        if x[0] >= L[0]/3.0 and x[0] <= L[0]*2.0/3.0:
            #pass
            return waterFlowRate
        else:
            return lambda x,t: 0.0
    if x[1] == 0.0:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}


#T = 0.5/timeScale
#cek testing
#T = 0.0025/timeScale
#mwf testing
T = 0.25
T = 3.0
