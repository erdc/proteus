#number of space dimensions
nd=1
#time integration
#
#quadrature
#timeIntegration_vof = "RK"
#timeIntegration_ls = "RK"
#
#timeIntegration_vof = "FLCBDF"
#timeIntegration_ls = "FLCBDF"
rtol_u = {0:1.0e-2}
atol_u = {0:1.0e-2}
rtol_res = {0:1.0e-2}
atol_res = {0:1.0e-2}
#
timeIntegration_vof = "BE"
timeIntegration_ls = "BE"
#try special purpose non-conservative transport operator or not
useNCLS = True
useRDLS = True
cDegree_ls=0
pDegree_ls=1
cDegree_vof=0
pDegree_vof=1
#runCFL=0.1
runCFL = 0.3#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
vortex_quad_order = 2*max(pDegree_vof,pDegree_ls)+1
#number nodes in each direction
nn=41#26#17
#nn=81
#number of meshes in multilevel mesh
nLevels = 1
#end time of simulation
T = 2.0
#number of output time steps, ignored for adaptive/CFL based runs
nDTout = 1
#mass correction
applyCorrection=False#True
applyRedistancing=True
#eps
epsFactHeaviside=3.0
oepsFactDirac=3.0
epsFactDiffusion=10.0
epsFactRedistance=0.0
#
shockCapturingFactor_vof=0.99
shockCapturingFactor_ls=0.99
shockCapturingFactor_rd=0.99
#atol
atolRedistance = 0.01/float(nn-1)
atolConservation = 1.0e-8
#redist solver
fmmFlag=0
checkMass=False#True
