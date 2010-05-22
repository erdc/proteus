from pyadh import *
from pyadh.default_p import *
name = "re_sand_2d_upwind_ncp1_sc"

nd = 2

L=(10.0,10.0,1.0)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
permeability  = [(5.04*m_per_s_by_m_per_d)*viscosity/(gravity*density)]  #m^2
print 'perm',permeability
thetaS        = [0.301]   #-
thetaR        = [0.093]   #-
thetaSR       = [thetaS[i]-thetaR[i] for i in range(1)] 
mvg_alpha     = [5.47]    #1/m
mvg_n         = [4.264]
mvg_m         = [1.0 - 1.0/mvg_n[0]]
lengthScale   = 1.0     #m
timeScale     = 1.0     #d #1.0/sqrt(g*lengthScale)
#make non-dimensional
dimensionless_conductivity  = [(timeScale*density*gravity*permeability[i]/(viscosity*lengthScale))/m_per_s_by_m_per_d for i in range(1)]
print 'Ks',dimensionless_conductivity
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([0.0,
                                       -1.0,
                                       0.0])
dimensionless_alpha    = [mvg_alpha[i]*lengthScale for i in range(1)]

upwindFlag = 0#1
coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2withUpwind(nd,
                                                                              Ks_block=numpy.array(dimensionless_conductivity),
                                                                              vgm_n_block=numpy.array(mvg_n),
                                                                              vgm_alpha_block=numpy.array(mvg_alpha),
                                                                              thetaR_block=numpy.array(thetaR),
                                                                              thetaSR_block=numpy.array(thetaSR),
                                                                              gravity=dimensionless_gravity,
                                                                              density=dimensionless_density,
                                                                              beta=beta,
                                                                              upwindFlag=upwindFlag,
                                                                              sd=1)

coefficients = ConservativeHeadRichardsMualemVanGenuchtenBlockHetV2(numpy.array(dimensionless_conductivity),
                                                                    numpy.array(mvg_n),
                                                                    numpy.array(mvg_alpha),
                                                                    numpy.array(thetaR),
                                                                    numpy.array(thetaSR),
                                                                    dimensionless_gravity,
                                                                    dimensionless_density,
                                                                    beta)


pondingPressure=0.1#-0.1
bottomPressure = 0.0
#pondingPressure=-0.1
#bottomPressure = -10.0
pondingSaturation = 0.9
waterTableSaturation = 0.9
initialSaturation = 0.01
#pondingPressure=-0.1

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


#T = 0.5/timeScale
#cek testing
#T = 0.0025/timeScale
#mwf testing
T = 0.5*0.25/timeScale
