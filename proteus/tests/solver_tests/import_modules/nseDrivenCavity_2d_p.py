from proteus import *
from proteus.default_p import *
from proteus import defaults
defaults.reset_default_p()
from proteus import Domain
try:
    from . import nseDrivenCavity_2d
    from .TwophaseNavierStokes_ST_LS_SO_VV import TwophaseNavierStokes_ST_LS_SO_VV
except:
    import nseDrivenCavity_2d
    from TwophaseNavierStokes_ST_LS_SO_VV import TwophaseNavierStokes_ST_LS_SO_VV
    
ct = nseDrivenCavity_2d.opts

"""
Navier-Stokes Driven Cavity Flow - this file contains the physics
for a NSE Driven Cavity test simulation.
"""

#######################################################

# variables set from context
name = ct.name
numeric_scheme = ct.numeric_scheme
RE_through_bdy = ct.RE_through_bdy
######################################################

# space dimension
nd = 2

# Mesh Creation
x0 = [-1.0,-1.0]
L = [2,2]
nnx=nny=17

if (numeric_scheme!= "C0Q1C0Q1" and numeric_scheme!="THQuads"):
    rdomain = Domain.RectangularDomain(x=x0[:2],L=L[:2],name="rdomain")
    polyfile="rdomain"
    rdomain.writePoly(polyfile)

##################################################

class uTrue(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return 1. - x[0]**4

class uTrue_RE_through_bdy(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return (10.0**max(0.0,t-1))*(1 - x[0]**4)

class vTrue(object):
    def __init__(self):
        pass
    def vOfX(self,x):
        return 0.0
    def uOfXT(self,x,t):
        return self.vOfX(x)

# dirichlet boundary condition functions pressure: (x=-1,y),velocity
# inflow: (x=1,y)

eps = 1.0e-8
he = 0.01
he *=0.5
he *=0.5
he *=0.5

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

def getDBCu_RE_through_bdy(x,flag):
    if sides(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    elif top(x,flag):
        return lambda x,t: uTrue_RE_through_bdy().uOfXT(x,t)

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

Advectivefluxboundaryconditions =  {0:getAdvFluxBCp,
                                    1:getAdvFluxBCu,
                                    2:getAdvFluxBCv}
boundaryCreatesNullSpace = True
