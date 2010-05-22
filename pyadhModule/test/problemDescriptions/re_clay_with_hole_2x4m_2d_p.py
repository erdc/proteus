from pyadh import *
from pyadh.default_p import *
import hole2d
import pdb

nd = 2
top = 2.0 #m
right = 4.0 #m
mydomain = [[0.0,right],[0.0,top]]
slit_left=1.9
slit_right=2.1
holeX=[1.5,2.5]
holeY=[1.5,1.75]
hole = [holeX, holeY]
sourceXs = [slit_left-0.1,slit_left,slit_right,slit_right+0.1]
polyfile='hole2d'
L,boundaryTags=hole2d.genPoly(domain=mydomain,
                              hole=hole, 
                              sourceXs=sourceXs, 
                              polyfileBase=polyfile)

analyticalSolution = None

viscosity     = 8.9e-4  #kg/(m*s)
density       = 998.2   #kg/m^3
gravity       = 9.8     #m/s^2
beta          = density*gravity*4.524e-10
m_per_s_by_m_per_d = 1.1574074e-5
permeability  = (7.97*m_per_s_by_m_per_d)*viscosity/(gravity*density)  #m^2
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
PI = atan(1.0)*4.0

def getPressure(x,t):
    data = [0.0, -1.98, 0.01, 0.02, 0.02, 0.02, 0.1, -4.0, 0.8, -4.0, 1.05, 0.02, 1.055, 0.02, 1.2, -4.0, 3.0, -4.0]
    index=2
    while t>data[index]:
        index+=2 
    p0 = data[index-1]
    p1 = data[index+1]
    t0 = data[index-2]
    t1 = data[index]
    li = p0+(t-t0)/(t1-t0)*(p1-p0)
    return li

def getDBC_2D_Richards_Shock(x,flag):
    if x[1] == L[1]:
        return getPressure
    if x[1] == 0.0:
        return lambda x,t: 0.02 


dirichletConditions = {0:getDBC_2D_Richards_Shock}

class ShockIC_2D_Richards:
    def uOfXT(self,x,t):
        bc=getDBC_2D_Richards_Shock(x,0)
        if bc != None:
            return bc(x,t)
        else:
            return (x[1]-0.02)*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ShockIC_2D_Richards()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 3.0
