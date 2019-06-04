from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *  
from proteus import Context

'''
flow around a 2D cylinder  benchmark problem.
'''

opts = Context.Options([
    ("T", 4.0, "Time interval [0, T]"),
    ("he",0.01, "maximum size of edges"),
    ("backwardEuler",False,"use backward Euler or not"),
    ("onlySaveFinalSolution",False,"Only save the final solution")
], mutable=True)



nd = 2
spaceOrder=1
useHex=False
DX = opts.he
usePETSc = False#True


parallelPartitioningType = MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1

L = 2.2
H = 0.41
Um = 1.5
radius = 0.05
fl_H = H

# Input checks
if spaceOrder not in [1,2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit() 
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

nLevels = 1
#from cylinder2dDomain import *
try:
    from .symmetricDomain_john import *
except:
    from symmetricDomain_john import *
    
domain = symmetric2D(box=(2.2,0.41),
                     L= 0.2,
                     H = 0.2,
                     r = 0.05,
                     C = (0.2,0.2),
                     DX = DX,
                     refinement_length=0.5,
                     DX_coarse = DX)
boundaryTags=domain.boundaryFlags

# Time stepping
T= opts.T
runCFL = 0.9
dt_fixed = 0.005
dt_init = 0.005
nDTout = int(old_div(T,dt_fixed))
dt_init = min(dt_init,0.5*dt_fixed)
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)] 
if opts.onlySaveFinalSolution == True:
    tnList = [0.0,dt_init,opts.T]


useBackwardEuler = opts.backwardEuler
# Numerical parameters
ns_shockCapturingFactor  = 0.0
ns_lag_shockCapturing = True#False
ns_lag_subgridError = True

epsFact_density    = 1.5
epsFact_viscosity  = 1.5
epsFact_redistance = 0.33
epsFact_consrv_heaviside = 1.5
epsFact_consrv_dirac     = 1.5
epsFact_consrv_diffusion = 10.0

# Fluid
rho = 1.0
#mu  =rho*0.2
nu = 1.0e-3


# Gravity
g = [0.0,0.0]

triangleOptions="pAq30ena"#D=Delaunay gives bad results for this composite meshing approach
genMesh=True
domain.writePLY('cylinder2D')
domain.writePoly('cylinder2D')
