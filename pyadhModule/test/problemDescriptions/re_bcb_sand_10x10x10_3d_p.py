from pyadh import *
from pyadh.default_p import *

nd = 3

L=(3.0,3.0,3.0)
#L=(10.0,0.1,10.0)

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
bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                        0.0,
                                        -1.0])
dimensionless_alpha    = mvg_alpha*lengthScale
dimensionless_pd    = bc_pd/lengthScale
coefficients = ConservativeHeadRichardsBrooksCoreyBurdine(hydraulicConductivity=dimensionless_conductivity,
                                                          gravity=dimensionless_gravity,
                                                          density=dimensionless_density,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          lambdab = bc_lambda,
                                                          pd = dimensionless_pd,
                                                          beta=beta)

pondingPressure=0.1

def getDBC_3D_Richards_Shock(x,flag):
    if x[2] == L[2]:
        if (x[0] >= L[0]/3.0 and
            x[0] <= 2.0*L[0]/3.0 and
            x[1] >= L[1]/3.0 and
            x[1] <= 2.0*L[1]/3.0):
            return lambda x,t: pondingPressure
    if x[2] == 0.0:
        return lambda x,t: 0.0
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1]):
#         return lambda x,t: x[2]*dimensionless_gravity[2]*dimensionless_density

dirichletConditions = {0:getDBC_3D_Richards_Shock}

class ShockIC_3D_Richards:
    def uOfXT(self,x,t):
        return x[2]*dimensionless_gravity[2]*dimensionless_density

initialConditions  = {0:ShockIC_3D_Richards()}

def getFBC_3D_Richards_Shock(x,flag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0
    if x[2] == L[2]:
        if not (x[0] >= L[0]/3.0 and
                x[0] <= 2.0*L[0]/3.0 and
                x[1] >= L[1]/3.0 and
                x[1] <= 2.0*L[1]/3.0):
            return lambda x,t: 0.0

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {0:getFBC_3D_Richards_Shock}

diffusiveFluxBoundaryConditions = {0:{0:getFBC_3D_Richards_Shock}}

T = 0.25/timeScale
T = 0.3*0.25/timeScale
