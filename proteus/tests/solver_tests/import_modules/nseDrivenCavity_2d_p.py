from proteus import *
from proteus.default_p import *
from proteus import Domain
import nseDrivenCavity_2d
from proteus import Context
Context.setFromModule(nseDrivenCavity_2d)
ct=Context.get()

"""
Stokes Driven Cavity Flow - the file contains the physics
corresponding to a Stokes Driven Cavity Flow.  See notes
for a detailed description.

This model is being used to study the preformance of
Schur complement preconditioners.
"""

#######################################################

# variables set from context
name = ct.name
numeric_scheme = ct.numeric_scheme

######################################################

# space dimension
nd = 2

# Mesh Creation
x0 = [-1.0,-1.0]
L = [2,2]

if (numeric_scheme!= "C0Q1C0Q1" and numeric_scheme!="THQuads"):
    rdomain = Domain.RectangularDomain(x=x0[:2],L=L[:2],name="rdomain")
    polyfile="rdomain"
    rdomain.writePoly(polyfile)

##################################################

class uTrue:
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return 10.0**max(0.0,t-1.0)

class vTrue:
    def __init__(self):
        pass
    def vOfX(self,x):
        return 0.0
    def uOfXT(self,x,t):
        return self.vOfX(x)

# dirichlet boundary condition functions pressure: (x=-1,y),velocity inflow: (x=1,y)

eps = 1.0e-8
he = 0.01
#he *=0.5
#he *=0.5

vent=False

# TODO - *** Add this a physics problem parameter ***
#pressure_has_null_space = True

def top(x,flag):
    if x[1]== x0[1]+L[1]:
        return True
    else:
        return False

def sides(x,flag):
    if (x[0]==x0[0] or x[0]==(x0[0]+L[0])):
        return True
    else:
        return False

def bottom(x,flag):
    if (x[1] == x0[1]):
        return True
    else:
        return False
        
def getDBCp(x,flag):
    pass

def getDBCu(x,flag):
    if sides(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    elif top(x,flag):
        return lambda x,t: uTrue().uOfXT(x,t)

def getDBCv(x,flag):
    if (top(x,flag) or sides(x,flag) or bottom(x,flag)):
        return lambda x,t: 0.0

def getAdvFluxBCp(x,flag):
    pass
    
def getAdvFluxBCu(x,flag):
    pass

def getAdvFluxBCv(x,flag):
    pass

#load dirichlet conditions, flux boundary conditions into the expected variables

dirichletConditions = {0:getDBCp,
                       1:getDBCu,
                       2:getDBCv}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBCp,
                                    1:getAdvFluxBCu,
                                    2:getAdvFluxBCv}
boundaryCreatesNullSpace = True
#diffusiveFluxBoundaryConditions = {0:{},
#                                   1:{1:getDiffFluxBCu} ,
#                                   2:{2:getDiffFluxBCv}}

#fluxBoundaryConditions = {0:'setFlow'} #options are 'setFlow','noFlow','mixedFlow'
# equation coefficient names
#coefficients = TransportCoefficients.Stokes(g=[0.0,0.0,0.0],nd=nd,steady=True,weakBoundaryConditions=False)
coefficients = TransportCoefficients.NavierStokes(rho=1.0,
                                                  nu=1.0,
                                                  g=[0.0,0.0],
                                                  nd=2)
#from proteus.mprans import RANS2P
#LevelModelType = RANS2P.LevelModel
#coefficients = RANS2P.Coefficients(rho_0=1.0,rho_1=1.0, nu_0=1.0, nu_1=1.0,sigma=0.0,g=[0.0,0.0],nd=2, forceStrongDirichlet=False)

# coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
#                                              sigma=0.0,
#                                              rho_0=1.0,
#                                              nu_0=1.0,
#                                              rho_1=1.0,
#                                              nu_1=1.0,
#                                              g=[0.0,0.0],
#                                              nd=2,
#                                              LS_model=None,
#                                              KN_model=None,
#                                              epsFact_density=None,
#                                              stokes=False);

coefficients.variableNames=['p','u','v']
