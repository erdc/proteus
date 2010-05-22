from pyadh import *
from pyadh.default_p import *
"""
Richards' equation with Brooks-Corey-Burdine Coefficients, 10m sand column, 1D
"""
##\page Tests Test Problems 
# \ref re_bcb_sand_10m_1d_p.py "Richards' equation, Brooks-Corey-Burdine constitutive relations, 10m sand column"
#

##\ingroup test
#\file re_bcb_sand_10m_1d_p.py
#\brief Richards' equation with Brooks-Corey-Burdine Coefficients, 10m sand column, 1D
#
#The equation formulation and coefficients are implemented in the ConservativeHeadRichardsBrooksCoreyBurdine
#class. The initial/boundary conditions are
#\f{eqnarray*}
#\psi(x,0) = -x \rho g \\
#\psi(0,t) = 0 \\
#\psi(10,t)= 0.1 
#\f}
#
#\todo add important references
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
bc_lambda     = mvg_n-1.0
bc_pd         = 1.0/mvg_alpha
print 'pd',bc_pd
print 'lambda',bc_lambda
(bc_lambda,bc_pd) = VGM_to_BCB_Simple(mvg_alpha,mvg_n)
print 'pd',bc_pd
print 'lambda',bc_lambda
(bc_lambda,bc_pd) = VGM_to_BCB_Johns(mvg_alpha,mvg_n)
print 'pd',bc_pd
print 'lambda',bc_lambda
(bc_lambda,bc_pd) = VGM_to_BCB_MorelSeytoux(mvg_alpha,mvg_n)
print 'pd',bc_pd
print 'lambda',bc_lambda
dimensionless_conductivity  = (timeScale*density*gravity*permeability/(viscosity*lengthScale))/m_per_s_by_m_per_d
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                        0.0,
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
bottomPressure = 0.0
#pondingPressure=-0.1
#bottomPressure = -10.0
#bottomPressure = 0.0
#bottomPressure = -1.0
#Pondingpressure=-bc_pd+1.0e-2
#pondingPressure=-0.0001

def getDBC_Richards_Shock(x,flag):
    if x[0] == L[0]:
        return lambda x,t: pondingPressure
    if x[0] == 0.0:
        return lambda x,t: bottomPressure

dirichletConditions = {0:getDBC_Richards_Shock}

class ShockIC_Richards:
    def uOfXT(self,x,t):
        f = getDBC_Richards_Shock(x,0)
        if f:
            return f(x,t)
        return bottomPressure+x[0]*dimensionless_gravity[0]*dimensionless_density

initialConditions  = {0:ShockIC_Richards()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.25/timeScale
