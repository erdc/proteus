from proteus import *
from proteus.default_p import *
reload(default_p)

##############
# PARAMETERS #
##############
ct=Context.Options([
    ("pDegree",2,"Order of the polynomial approximation"),
    ("refinement",0,"Mesh refinement"),
    ("useHex",True,"Use quads?"),
    ("useBernstein",True,"Use Bernstein polynomials"),
    ("unstructured",False,"Use unstructured triangular mesh")
],mutable=True)

#pDegree=2
#refinement=1
#useHex=True
#useBernstein=True
#xunstructured=False

name = "poisson"
nd = 2
initialConditions = None

##########
# DOMAIN #
##########
nn=(2**ct.refinement)*10+1
he=1.0/(nn-1.0)
box=Domain.RectangularDomain(L=(1.0,1.0),
                             x=(0.0,0.0),
                             name="box");
box.writePoly("box")
if ct.unstructured:
    assert useHex==False, "set useHex=False for unstructure meshes"
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box
    
nc = 1
##################
# EXACT SOLUTION #
##################
class exact_soln:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.sin(2*pi*x[0])*np.sin(2*pi*x[1])

    def duOfXT(self,x,t):
        return [2*pi*np.cos(2*pi*x[0])*np.sin(2*pi*x[1]),
                2*pi*np.sin(2*pi*x[0])*np.cos(2*pi*x[1])]
analyticalSolution = {0:exact_soln()}

#########################
# DIFFUSION COEFFICIENT #
#########################
def A(x):
    return numpy.array([[1.0,0.0],[0.0,1.0]],'d')
aOfX = {0:A}

##############
# FORCE TERM #
##############
def f(x):
    return 8*pi**2*np.sin(2*pi*x[0])*np.sin(2*pi*x[1])
fOfX = {0:f}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: 0.
dirichletConditions = {0:getDBC}

coefficients = TransportCoefficients.PoissonEquationCoefficients({0:A},{0:f},nc,nd,l2proj=[False])

class LevelModelType(OneLevelTransport):
    def getResidual(self,u,r):
        OneLevelTransport.getResidual(self,u,r)
    def calculateElementResidual(self):
        OneLevelTransport.calculateElementResidual(self)
