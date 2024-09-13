import numpy as np
from proteus.iproteus import *
import proteus.default_p as physics
import proteus.default_n as numerics
from proteus.TransportCoefficients import PoissonEquationCoefficients
from proteus  import Profiling
from proteus.Profiling import logEvent
Profiling.logLevel=7
Profiling.verbose=True


import copy

he_max=1.
he_min=0.025
r = 0.25
nSmoothOut = 0
epsTimeStep = 0.1



nd = 2

# use structured mesh
domain = Domain.RectangularDomain(L=[1.0,1.0], x=[0.,0.])
fixedNodeMaterialTypes = np.zeros(10)
he = 0.02

center1 = [0.5, 0.5]

genMesh = True



def my_func(x, t):
    return min(he_max, max(np.abs(np.sqrt((x[0]-0.5)**2+(x[1]-0.5)**2)-r), he_min))
scale = 0.5



from proteus.default_n import *

nn = 33


spaceOrder = 1
useHex = False
nd = domain.nd
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
         basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
         basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
        basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
        basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)



##########################################
# Numerical Options and other parameters #
##########################################

from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# other flags
#----------------------------------------------------
movingDomain = True
checkMass = False
applyCorrection = True
applyRedistancing = True
freezeLevelSet = True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
dt_init = 0.001
T = 0.001
nDTout = 1
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0

#----------------------------------------------------

#  Discretization -- input options
spaceOrder = 1
useHex     = False
useRBLES   = 0.0

# Input checks
if spaceOrder not in [1,2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES)
    sys.exit()

#  Discretization
nd = domain.nd
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
        basis=C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,3)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
        basis=C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
        basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
        basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
