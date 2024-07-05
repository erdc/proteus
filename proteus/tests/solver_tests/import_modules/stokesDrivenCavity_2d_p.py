from proteus import *
from proteus.default_p import *
from proteus import defaults
defaults.reset_default_p()
from proteus import Domain
import os
try:
    from . import stokesDrivenCavity_2d
except:
    import stokesDrivenCavity_2d
ct = stokesDrivenCavity_2d.opts
"""
Stokes Driven Cavity Flow - this file contains the physics
corresponding to a Stokes Driven Cavity Flow test simulation.
"""

#######################################################

# variables set from context
name = ct.name
numeric_scheme = ct.numeric_scheme
useWeakBoundaryConditions = ct.useWeakBoundaryConditions

######################################################
#                   Mesh Options

# space dimension
nd = 2

# Mesh Properties
# (bottom left part of mesh)
x0 = [-1.,-1.]
# (mesh length in x and y direction)
L = [2,2]

if (numeric_scheme != "C0Q1C0Q1" and numeric_scheme != "THQuads"):
    rdomain = Domain.RectangularDomain(x=x0[:2],L=L[:2],name="rdomain")
    polyfile=os.path.dirname(os.path.abspath(__file__))+"/../"+"rdomain"
    rdomain.polyfile = polyfile
    #rdomain.writePoly(polyfile)
    genMesh=False

##################################################

# TODO - these classes should be renamed as we do not have a true solution for this problem

class uTrue(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return self.uOfX(x)

class vTrue(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return 0.0
    def uOfXT(self,x,t):
        return self.uOfX(x)


#################################################

# boundary conditions

def getDBCp(x,flag):
    pass

def getDBCu(x,flag):
    if x[1] == (x0[1]+L[1]): #inflow/no slip
        return lambda x,t: uTrue().uOfXT(x,t)
    elif (x[0] == x0[0] or x[0] == (x0[0]+L[0]) or x[1] == x0[1] ):
        return lambda x,t: 0.0 # no slip on 'left', 'bottom' and 'right'

def getDBCv(x,flag):
    if (x[0] == x0[0] or x[0] == (x0[0]+L[0]) or x[1] == x0[1] or x[1] == (x0[1]+L[1]) ):
        return lambda x,t: 0.0 # no slip on 'left', 'bottom' and 'right'

def getAdvFluxBCp(x,flag):
    if (x[0] == x0[0] or x[0] == (x0[0]+L[0]) or x[1] == x0[1] or x[1] == (x0[1]+L[1]) ):
        return lambda x,t: 0.0

# def getDiffFluxBCu(x,flag):
#     if (x[0]==1.0):#outflow
#         return lambda x,t: 0.0

# def getDiffFluxBCv(x,flag):
#     if (x[0]==1.0):#outflow
#         return lambda x,t: 0.0


def getAdvFluxBCu(x,flag):
    pass

def getAdvFluxBCv(x,flag):
    pass


#load dirichlet conditions, flux boundary conditions into the expected variables

dirichletConditions = {0:getDBCp, 1:getDBCu, 2:getDBCv}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBCp, 1:getAdvFluxBCu, 2:getAdvFluxBCv}
# diffusiveFluxBoundaryConditions = {0:{}, 1:{1:getDiffFluxBCu} , 2:{2:getDiffFluxBCv}}
# fluxBoundaryConditions = {0:'setFlow'} #options are 'setFlow','noFlow','mixedFlow'

# equation coefficient names
if useWeakBoundaryConditions:
    coefficients = TransportCoefficients.Stokes(rho=1.,nu=1.,g=[0,0],nd=2,steady=True,weakBoundaryConditions=True)
else:
    coefficients = TransportCoefficients.Stokes(rho=1.,nu=1.,g=[0,0],nd=2,steady=True,weakBoundaryConditions=False)

coefficients.variableNames=['p','u','v']
