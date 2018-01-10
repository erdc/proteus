from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *  
'''
flow around a 2D cylinder  benchmark problem.
'''

nd = 2
spaceOrder=1
Refinement=2
useHex=False
points_on_grain = 21
DX = 0.03
DX *=0.5
#DX *=0.5
usePETSc = False#True

parallelPartitioningType = MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1

L = 2.2
H = 0.41
Um = 1.0e-1
radius = 0.05
fl_H = H

# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
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
from cylinder2dDomain import *
domain = cylinder2D(L=L,
                    H=H,
                    C=(0.5*H,0.5*H),#C=(0.2,0.2),
                    r=radius,
                    points_on_grain=4*int(ceil(2.0*pi*radius/DX))+1)
boundaryTags=domain.boundaryFlags

# Time stepping
T= 8.0
runCFL = 0.33
dt_fixed = 2.5e-1
dt_init = 1.0e-3
nDTout = int(T/dt_fixed)
dt_init = min(dt_init,0.5*dt_fixed)
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)] 


useBackwardEuler = False

# Numerical parameters
ns_shockCapturingFactor  = 0.2
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


triangleOptions="pAq30.0Dena%f" % (.5*DX**2)  #% (0.5*(DX)**2,)
print triangleOptions
genMesh=True
domain.writePLY('cylinder2D')
domain.writePoly('cylinder2D')
