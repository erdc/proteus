from pyadh import *
from pyadh.default_p import *
import KueperFrindGenPoly3D
##\page Tests Test Problems
# \ref re_kueper_2d_p.py "Richards' equation, Brooks-Corey-Burdine constitutive relations, heterogeneous tank experiment"
#

##\ingroup test
#\file re_kueper_2d_p.py
#\brief Richards' equation with Brooks-Corey-Burdine Coefficients, hterogeneous tank experiment
#
#The equation formulation and coefficients are implemented in the ConservativeHeadRichardsBrooksCoreyBurdineHet
#class. The initial/boundary conditions are
#\f{eqnarray*}
#\psi(x,0) = -x \rho g \\
#\psi(0,t) = 0 \\
#\psi(10,t)= 0.1
#\f}
#

nd = 3

meshfile = "het_column_8"

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = 0.0001#density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
Ks = [1.1808e-1,1.1808e-3,1.1808e-5]
thetaS=[0.5,0.4,0.33]
permeability  = [(K*m_per_s_by_m_per_d)*viscosity/(gravity*density)  for K in Ks] #m^2
thetaR        = [0.105*0.5,0.074*0.4,0.179*0.33]
thetaSR = [(tS-tR) for tS,tR in zip(thetaS,thetaR)]
mvg_alpha     = [4.42*0.3048, 0.487*0.3048, 0.244*0.3048]
mvg_n         = [2.68,1.37,1.09]
lengthScale   = 1.0     #m
timeScale     = 1.0     #d 
#make non-dimensional
dimensionless_conductivity  = [(timeScale*density*gravity*k/(viscosity*lengthScale))/m_per_s_by_m_per_d for k in permeability]
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                      0.0,
                                      -1.0])
dimensionless_alpha    = [mvg*lengthScale for mvg in mvg_alpha]
coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2(Ks_block=numpy.array(dimensionless_conductivity),
                                                                    vgm_n_block=numpy.array(mvg_n),
                                                                    vgm_alpha_block=numpy.array(dimensionless_alpha),
                                                                    thetaR_block=numpy.array(thetaR),
                                                                    thetaSR_block=numpy.array(thetaSR),
                                                                    gravity=dimensionless_gravity,
                                                                    density=dimensionless_density,
                                                                    beta=beta)
pondingpressure= 0.1 #-0.01 #0.1
bottomPressure = 0.0
def getDBC_2D_Richards_KueperFrind_Shock(x,flag):
#     if flag == 1:
#         return lambda x,t: 0.0
#     else:
#         return None
   if x[2] <= -10.0+1.0e-8:
       return lambda x,t: 0.0
   else:
       return None
#     if x[2] >= 0.0-1.0e-8:
#         return lambda x,t: pondingPressure
#    if flag == 1:
#            return lambda x,t: 0.0#bottomPressure
#     if flag == 2:
#             return lambda x,t: pondingPressure

dirichletConditions = {0:getDBC_2D_Richards_KueperFrind_Shock}

class ShockIC_3D_Richards:
    def uOfXT(self,x,t):
        return x[2]*dimensionless_gravity[2]*dimensionless_density

initialConditions  = {0:ShockIC_3D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

def getAFBC_3D_Richards_KueperFrind_Shock(x,flag):
#     if flag in [0,2]:
#         return lambda x,t: 0.0
#     else:
#         return None
   if x[2] <= -10.0+1.0e-8:
       return None
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_3D_Richards_KueperFrind_Shock}

def getDFBC_3D_Richards_KueperFrind_Shock(x,flag):
#     if flag in [0,2]:
#         return lambda x,t: 0.0
#     else:
#         return None
   if x[2] <= -10.0+1.0e-8:
       return None
   return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_3D_Richards_KueperFrind_Shock}}

T = 1.0e20
T = 1000.0
