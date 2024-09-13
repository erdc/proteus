from proteus import *
from proteus.default_p import *
import os
try:
    from .parameters_for_poisson import *
except:
    from parameters_for_poisson import *

name = "poisson"
assert ct.nd==2 or ct.nd==3, "Choose nd=2, or 3"
nd = ct.nd
initialConditions = None

##########
# DOMAIN #
##########
nn=(2**ct.refinement)*10+1
he=1.0/(nn-1.0)
if (nd==2):
    box=Domain.RectangularDomain(L=(1.0,1.0),
                                 x=(0.0,0.0),
                                 name="box");
else:
    nn=5
    he=1.0/(nn-1.0)
    box=Domain.RectangularDomain(L=(1.0,1.0,1.0),
                                 x=(0.0,0.0,0.0),
                                 name="box");
genMesh=ct.genMesh #False
if(genMesh):
    box.writePoly("box")

if ct.unstructured:
    assert ct.useHex==False, "set useHex=False for unstructure meshes"
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box

domain.MeshOptions.nn = domain.MeshOptions.nnx = domain.MeshOptions.nny = domain.MeshOptions.nnz = nn    
 
#the initial test uses triangleFlag=0 from defaults
domain.MeshOptions.triangleFlag=0

nc = 1
##################
# EXACT SOLUTION #
##################
class exact_soln(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if (nd==2):
            return np.sin(2*pi*x[0])*np.sin(2*pi*x[1])
        else:
            return np.sin(2*pi*x[0])*np.sin(2*pi*x[1])*np.sin(2*pi*x[2])

    def duOfXT(self,x,t):
        if (nd==2):
            return [2*pi*np.cos(2*pi*x[0])*np.sin(2*pi*x[1]),
                    2*pi*np.sin(2*pi*x[0])*np.cos(2*pi*x[1])]
        else:
            return [2*pi*np.cos(2*pi*x[0])*np.sin(2*pi*x[1])*np.sin(2*pi*x[2]),
                    2*pi*np.sin(2*pi*x[0])*np.cos(2*pi*x[1])*np.sin(2*pi*x[2]),
                    2*pi*np.sin(2*pi*x[0])*np.sin(2*pi*x[1])*np.cos(2*pi*x[2])]
analyticalSolution = {0:exact_soln()}

#########################
# DIFFUSION COEFFICIENT #
#########################
def A(x):
    if (nd==2):
        return numpy.array([[1.0,0.0],[0.0,1.0]],'d')
    else:
        return numpy.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],'d')
aOfX = {0:A}

##############
# FORCE TERM #
##############
def f(x):
    if (nd==2):
        return 8*pi**2*np.sin(2*pi*x[0])*np.sin(2*pi*x[1])
    else:
        return 12*pi**2*np.sin(2*pi*x[0])*np.sin(2*pi*x[1])*np.sin(2*pi*x[2])
fOfX = {0:f}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    if (nd==2):
        if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
            return lambda x,t: 0.
    else:
        if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
            return lambda x,t: 0.
dirichletConditions = {0:getDBC}

coefficients = TransportCoefficients.PoissonEquationCoefficients({0:A},{0:f},nc,nd,l2proj=[False])

class LevelModelType(OneLevelTransport):
    def getResidual(self,u,r):
        OneLevelTransport.getResidual(self,u,r)
    def calculateElementResidual(self):
        OneLevelTransport.calculateElementResidual(self)
