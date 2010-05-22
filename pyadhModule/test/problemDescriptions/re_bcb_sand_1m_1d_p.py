from pyadh import *
from pyadh.default_p import *
##\page Tests Test Problems
# \ref re_bcb_sand_1m_1d_p.py "Richards' equation, Brooks-Corey-Burdine constitutive relations, 1m sand column"
#

##\ingroup test
#\file re_bcb_sand_1m_1d_p.py
#\brief Richards' equation with Brooks-Corey-Burdine Coefficients, 1m sand column, 1D
#
#The equation formulation and coefficients are implemented in the ConservativeHeadRichardsBrooksCoreyBurdine
#class. The initial/boundary conditions are
#\f{eqnarray*}
#\psi(x,0) = -x \rho g \\
#\psi(0,t) = 0 \\
#\psi(10,t)= 0.1
#\f}
#

nd = 1

L=(1.0,1.0,1.0)

analyticalSolutions = None

lengthScale   = 1.0     #m
permeability  = 1.0e-7  #m^2
viscosity     = 8.9e-4  #kg/(m*s)
density       = 997.0   #kg/m^3
gravity       = 9.8     #m/s^2
thetaS        = 0.301   #-
thetaR        = 0.093   #-
#thetaR        = 0.21
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
#make non-dimensional
dimensionless_conductivity  = density*gravity*permeability/(viscosity*sqrt(gravity*lengthScale))
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                      0.0,
                                      0.0])
dimensionless_pd    = bc_pd*lengthScale

coefficients = ConservativeHeadRichardsBrooksCoreyBurdine(hydraulicConductivity=dimensionless_conductivity,
                                                          gravity=dimensionless_gravity,
                                                          density=dimensionless_density,
                                                          thetaS=thetaS,
                                                          thetaR=thetaR,
                                                          lambdab = bc_lambda,
                                                          pd = dimensionless_pd)

pondingPressure=0.1

def getDBC_Richards_Shock(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        return x[0]*dimensionless_gravity[0]*dimensionless_density

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.75
