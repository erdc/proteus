from pyadh import *
from pyadh.default_p import *
from Buckley_Leverett_TransportCoefficients import *
#John Chrispell, Summer 07
L=(1.,1.0,1.0)

nd = 1

a= 1.0e-3
q=[1.0]#[-1.0]

name = "Buckley_Leverett"
coefficients = Nonlinear_VA_linearDR_Coefficients(q=q,a=a,mu_r=0.5)
uLeft = 1.0; uRight = 0.0
#uLeft = 0.0; uRight = 1.0
#now define the Dirichlet boundary conditions
riemFudgeFactor = 0.005
def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: uLeft
def dummyBC(x):
    pass
dirichletConditions = {0:getDBC}#for Riemann problem {0:dummyBC}

class ShockIC:
    def __init__(self,dbc,uRight=0.0):
        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
        self.uRight= uRight
    def uOfXT(self,x,t):
        if x[0] <= riemFudgeFactor:#riemFudgeFactor: #add fudge factor for dg methods to get ic?
            return self.uLeft
        else:
            return self.uRight
class RiemIC:
    def __init__(self,uleft=1.0,uright=0.0):
        self.uLeft=uleft; self.uRight=uright
    def uOfXT(self,x,t):
        if x[0] <= 0.5*L[0]:
            return self.uLeft
        else:
            return self.uRight

#I believe we need something like this right now for dgpk because default projection of step data
#causes oscillations
class RiemSlopeIC:
    def __init__(self,uleft=1.0,uright=0.0,du=0.1):
        self.uLeft=uleft; self.uRight=uright; self.du=du
    def uOfXT(self,x,t):
        if x[0] <= 0.5*L[0]:
            return self.uLeft
        if x[0] > (0.5+self.du)*L[0]:
            return self.uRight
        return self.uLeft + (self.uRight-self.uLeft)/(self.du*L[0])*(x[0]-0.5*L[0])


initialConditions  = {0:ShockIC(getDBC,uRight=uRight)}#{0:RiemIC(uleft=uLeft,uright=uRight)}
#{0:RiemSlopeIC(uleft=uLeft,uright=uRight,du=0.02)}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5#1.0

analyticalSolution = {0:Buckley_Leverett_RiemannSoln(coefficients,uLeft=uLeft,uRight=uRight,
                                                      t0=0.0,x0=riemFudgeFactor,T=T)}
