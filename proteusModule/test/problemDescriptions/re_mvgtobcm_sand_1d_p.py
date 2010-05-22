from pyadh import *
from pyadh.default_p import *

nd = 1

L=(1.0,1.0,1.0)

analyticalSolutions = None

class ConservativeHeadRichardsBCBfromMVG(TransportCoefficients.TC_base):
    """
    Richards' equation with Brooks-Corey-Burdine constitutive relations computed from Mualem-Van Genuchten parameters.
    """
    from pyadh.ctransportCoefficients import  conservativeHeadRichardsBCBfromMVGHomEvaluate
    def __init__(self,
                 hydraulicConductivity,
                 gravity,
                 density,
                 thetaS,
                 thetaR,
                 alpha,
                 n,
                 m):
        self.variableNames=['psi']
        nc=1
        mass={0:{0:'nonlinear'}}
        advection={0:{0:'nonlinear'}}
        diffusion={0:{0:{0:'nonlinear'}}}
        potential={0:{0:'u'}}
        reaction={0:{0:'linear'}}
        hamiltonian={}
        TransportCoefficients.TC_base.__init__(self,
                                               nc,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian)
        self.Ks = hydraulicConductivity
        self.gravity=gravity
        self.rho = density
        self.thetaS = thetaS
        self.thetaR = thetaR
        self.thetaSR = thetaS - thetaR
        self.alpha = alpha
        self.n = n
        self.m = m
    def evaluate(self,t,c):
        self.conservativeHeadRichardsBCBfromMVGHomEvaluate(self.rho,
                                                           self.gravity,
                                                           self.alpha,
                                                           self.n,
                                                           self.m,
                                                           self.thetaR,
                                                           self.thetaSR,
                                                           self.Ks,
                                                           c[('u',0)],
                                                           c[('m',0)],
                                                           c[('dm',0,0)],
                                                           c[('f',0)],
                                                           c[('df',0,0)],
                                                           c[('a',0,0)],
                                                           c[('da',0,0,0)])

lengthScale   = 1.0     #m
permeability  = 1.0e-7  #m^2
viscosity     = 8.9e-4  #kg/(m*s)
density       = 997.0   #kg/m^3
gravity       = 9.8     #m/s^2
thetaS        = 0.301   #-
thetaR        = 0.093   #-
mvg_alpha     = 5.47    #1/m
mvg_n         = 4.264
mvg_m         = 1.0 - 1.0/mvg_n
#make non-dimensional
dimensionless_conductivity  = density*gravity*permeability/(viscosity*sqrt(gravity*lengthScale))
dimensionless_density  = 1.0
dimensionless_gravity  = numpy.array([-1.0,
                                        0.0,
                                        0.0])
dimensionless_alpha    = mvg_alpha*lengthScale
coefficients = ConservativeHeadRichardsBCBfromMVG(hydraulicConductivity=dimensionless_conductivity,
                                                  gravity=dimensionless_gravity,
                                                  density=dimensionless_density,
                                                  thetaS=thetaS,
                                                  thetaR=thetaR,
                                                  alpha= dimensionless_alpha,
                                                  n = mvg_n,
                                                  m = mvg_m)

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

T = 0.25
