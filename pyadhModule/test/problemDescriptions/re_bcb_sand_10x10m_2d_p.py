from pyadh import *
from pyadh.default_p import *
"""
Richards' equation with Brooks-Corey-Burdine Coefficients, 10x10m sand column, 2D
"""
##\page Tests Test Problems 
# \ref re_bcb_sand_10x10m_p.py "Richards' equation, Brooks-Corey-Burdine constitutive relations, 10x10m sand column"
#

##\ingroup test
#\file re_bcb_sand_10x10m_p.py
#\brief Richards' equation with Brooks-Corey-Burdine Coefficients, 10x10m sand column, 2D
#
#The equation formulation and coefficients are implemented in the ConservativeHeadRichardsBrooksCoreyBurdine
#class. The initial/boundary conditions are
#\f{eqnarray*}
#\psi(x,y,0) = -y \\
#\psi(x,0,t) = 0 \\
#\psi(x,10,t)= 0.1 \quad 1/3 \leq x \leq 2/3 ]]
#\psi(x,10,t)= -10
#\f}
#

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
bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
print 'pd',bc_pd
#make non-dimensional
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                        -1.0,
                                        0.0])
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

def getDBC_2D_Richards_Shock(x,flag):
    if x[1] == L[1]:
        if (x[0] >= L[0]/3.0 and
            x[0] <= 2.0*L[0]/3.0):
            return lambda x,t: pondingPressure
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if (x[0] == 0.0 or
        x[0] == L[0]):
        return lambda x,t: x[1]*dimensionless_gravity[1]*dimensionless_density

dirichletConditions = {0:getDBC_2D_Richards_Shock}

class ShockIC_2D_Richards:
    def uOfXT(self,x,t):
        bc=getDBC_2D_Richards_Shock(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return x[1]*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ShockIC_2D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

def getFBC_2D_Richards_Shock(x,flag):
    if x[1] == L[1]:
        if (x[0] < L[0]/3.0 or
            x[0] > 2.0*L[0]/3.0):
            return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getFBC_2D_Richards_Shock}

diffusiveFluxBoundaryConditions = {0:{0:getFBC_2D_Richards_Shock}}

T = 0.5*0.25/timeScale
