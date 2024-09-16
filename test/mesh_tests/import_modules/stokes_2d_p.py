from proteus import *
from proteus import iproteus as ip
from proteus.default_p import *
reload(default_p)
from proteus import Domain
try:
    from . import stokes_2d
except:
    import stokes_2d
from proteus import Context
Context.setFromModule(stokes_2d)
ct=Context.get()
"""
Stokes Poiseuille Equation - this file contains the physics
corresponding to a Stokes' Poiseuille flow.

This model is used primarily as a test to confirm the
accuracy of various numerical schemes.
"""

#######################################################

# variables set from context
name = "poiseulleFlow"
numeric_scheme = "THQuads"
useWeakBoundaryConditions = False

######################################################

# space dimension
nd = 2

# Mesh Creation
# (bottom left corner of mesh)
x0 = [-1.0,-1.0]
# (mesh length in x and y direction)
L = [2,2]
domain=None
polyfile = None
if (numeric_scheme!= "C0Q1C0Q1" and numeric_scheme!="THQuads"):
    rdomain = Domain.RectangularDomain(x=x0[:2],L=L[:2],name="rdomain")
    polyfile="rdomain"
    rdomain.writePoly(polyfile)
    

##################################################

# analytical solution

class pTrue(object):
    def __init__(self):
        self.grad_p = L[0]
        pass
    def uOfX(self,x):
        return self.grad_p*(L[0]-(x[0]-x0[0]))
    def uOfXT(self,x,t):
        return self.uOfX(x)

class uTrue(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return -(x[1]-x0[1])*(x[1]-(x0[1]+L[1]))
    def uOfXT(self,x,t):
        return self.uOfX(x)

class vTrue(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return 0.0
    def uOfXT(self,x,t):
        return self.uOfX(x)

# boundary conditions

def getDBCp(x,flag):
    if (x[0] == x0[0]+L[0]):
        return lambda x,t: pTrue().uOfXT(x,t)

def getDBCu(x,flag):
    if (x[1] == x0[1] or x[1] == x0[1]+L[1]):
        return lambda x,t: uTrue().uOfXT(x,t)

def getDBCv(x,flag):
    if (x[1] == x0[1] or x[1] == x0[1]+L[1]):
        return lambda x,t: vTrue().uOfXT(x,t)

def getAdvFluxBCp(x,flag):
    if (x[0] == x0[0]):
        return lambda x,t: -uTrue().uOfXT(x,t)
    elif (x[1] == x0[1] or x[1] == x0[1]+L[1]):
        return lambda x,t: 0.0
    else:
        pass

def getAdvFluxBCu(x,flag):
    pass

def getAdvFluxBCv(x,flag):
    pass

def getDiffFluxBCu(x,flag):
    if (x[0] == x0[0] + L[0]):
        return lambda x,t: 0.0 

def getDiffFluxBCv(x,flag):
    if (x[0] ==  x0[0] + L[0]):
        return lambda x,t: 0.0

    
# set boundary conditions
dirichletConditions = {0:getDBCp, 1:getDBCu, 2:getDBCv}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBCp, 1:getAdvFluxBCu, 2:getAdvFluxBCv}
diffusiveFluxBoundaryConditions = {0:{}, 1:{1:getDiffFluxBCu} , 2:{2:getDiffFluxBCv}}
# fluxBoundaryConditions = {0:'setFlow'} #options are 'setFlow','noFlow','mixedFlow'

# equation coefficient names
if useWeakBoundaryConditions:
    coefficients = TransportCoefficients.Stokes(rho=1.,nu=1.,g=[0,0],nd=2,steady=True,weakBoundaryConditions=True)
else:
    coefficients = TransportCoefficients.Stokes(rho=1.,nu=1.,g=[0,0],nd=2,steady=True,weakBoundaryConditions=False)

coefficients.variableNames=['p','u','v']

# load analytical solution
analyticalSolution = {0: pTrue(),1: uTrue(), 2:vTrue()}
