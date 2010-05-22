from pyadh import *
from pyadh.default_p import *

nd = 2


analyticalSolution = None

polyfile="E6_jun1_50x250_200x300"#"trislicex50to100y200to300"#E6_jun1_50x250_200x300.poly

rho_mm=9.9820e-07
mu_mm = 1.0000e-09
nu_mm = mu_mm/rho_mm
#if want to use weak boundary conditions, need stokesP 
#coefficients = StokesP(g=[0.0,0.0],rho=rho_mm,nu=nu_mm,nd=nd)
#otherwise use Stokes
coefficients = Stokes(g=[0.0,0.0],nd=nd)

#now define the Dirichlet boundary conditions
inflow = 0.1
inflowP= 1.e-4
L = (1.0,1.0,1.0)
EPS = 1.0e-6
def onSolidBoundary(x,tag):
    if tag == 1 and (x[0] > EPS and x[0] < L[0]-EPS and x[1] > EPS and x[1] < L[1]-EPS):
        return True
    return False
def getDBC_pressure(x,tag):
    if x[0] >= L[0]-EPS and tag == 1:
        return lambda x,t: (L[1]-x[1])*coefficients.rho*coefficients.g[1]
    elif x[0] <= EPS and tag == 1:
        return lambda x,t: inflowP+(L[1]-x[1])*coefficients.rho*coefficients.g[1]

def getDBC_u(x,tag):
    if tag == 1 and (x[1] <= EPS or
                     x[1] >= L[1]-EPS):
        return lambda x,t: 0.0
    #elif (x[0] <= EPS):
    #    return lambda x,t: inflow
    #no flow on solid boundary
    elif onSolidBoundary(x,tag):
        return lambda x,t: 0.0
    
def getDBC_v(x,tag):
    if tag == 1 and (x[1] <= EPS or
                     x[1] >= L[1]-EPS): 
        return lambda x,t: 0.0
    #no flow on solid boundary
    elif onSolidBoundary(x,tag):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}


def getAFBC_p(x,tag):
    if tag == 1 and ((x[1] >= L[1]-EPS or
                      x[1] <= EPS) and
                     EPS < x[0] and
                     x[0] < L[0]-EPS): 
        return lambda x,t: 0.0
    if onSolidBoundary(x,tag):
        return lambda x,t: 0.0
    #elif (x[0] <= EPS):
    #    return lambda x,t: -inflow

def dummyBC(x,tag):
    pass
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,1:dummyBC,2:dummyBC}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

class constantIC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        return self.val
    def uOfX(self,x):
        return self.val
    
initialConditions = {0:constantIC(0.),1:constantIC(0.),2:constantIC(0.)}

