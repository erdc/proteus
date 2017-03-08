from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
'''
flow around a 2D cylinder  benchmark problem.
'''
x0 = [-1.0,-1.0]
L  = [2,2]

# Fluid
rho = 1.0e0
##mu  =rho*0.2
nu = 1.0e-3
he = 0.0005
nd = 2
spaceOrder=1
Refinement=2
useHex=False
usePETSc = True

rdomain = Domain.RectangularDomain(x=x0[:2],L=L[:2],name="rdomain")
polyfile="rdomain"
rdomain.writePoly(polyfile)

# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         triangleOptions = "VApq30Dena%8.8f" % he
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
        triangleOptions = "VApq30Dena%8.8f" % he

nLevels = 1
# Time stepping
T=  100.00
#runCFL = 0.33
dt_fixed = 4.5
dt_init = 4.5
#dt_fixed = 1.0e-2 #2.5e-1
#dt_init = 1.0e-2
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

# Gravity
g = [0.0,0.0]
