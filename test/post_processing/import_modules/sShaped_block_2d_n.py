from proteus import *
from proteus.default_n import *
from .sShaped_block_2d_p import *

######################################################

# context variables
nLevels = 1
finiteElement = 'c0p1'
conservativeFluxSchemeFlag = 'pwl-bdm'
numericalFluxTypeFlag = 'none'
parallel = 'false'

#####################################################

# Define finite element space

if finiteElement=="c0p1":
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif finiteElement=="ncp1":
    femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}
elif finiteElement=="dgp1":
    femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
elif finiteElement=="c0p2":
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
else:
    print('INVALID FINITE ELEMENT SELECTED')

#######################################################

# Define conservative flux scheme

if conservativeFluxSchemeFlag=='pwl-bdm':
    conservativeFlux = {0:'pwl-bdm'}
elif conservativeFluxSchemeFlag=='pwl-bdm2':
    conservativeFlux = {0:'pwl-bdm2'}
elif conservativeFluxSchemeFlag=='pwc':
    conservativeFlux = {0:'pwc'}
elif conservativeFluxSchemeFlag=='p1-nc':
    conservativeFlux = {0:'p1-nc'}
elif conservativeFluxSchemeFlag=='dg-point-eval':
    conservativeFlux = {0:'dg-point-eval'}#{0:'pwl'}
elif conservativeFluxSchemeFlag=='none':
    pass
else:
    print('INVALID CONSERVATIVE FLUX SCHEME')


######################################################

# Define numerical flux type

if numericalFluxTypeFlag=='none':
    pass
elif numericalFluxTypeFlag=='sipg':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
elif numericalFluxTypeFlag=='iipg':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
elif numericalFluxTypeFlag=='ldg':
    numericalFluxType = Diffusion_LDG
elif numericalFluxTypeFlag=='nipg':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_NIPG
else:
    print('INVALID NUMERICAL FLUX TYPE')

######################################################

timeIntegration = NoIntegration
nDTout = 1

quad_order = 4
#mwf everything gets SimplexGaussQuadrature
elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

triangleOptions = "q30Dena0.005A"

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc

levelLinearSolver = LU
#levelLinearSolver = PETSc

linTolFac = 1.0e-10
