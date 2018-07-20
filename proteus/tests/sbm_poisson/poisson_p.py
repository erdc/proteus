from proteus import *
from proteus.default_p import *
import math
from poisson2D import *
import Poisson_M as LV

#===============================================================================
# Level model type
#===============================================================================
LevelModelType = LV.LevelModel
logEvent = Profiling.logEvent
name = soname + "_ls"

nd = 2

#===============================================================================
# My Coefficient
#===============================================================================
class My_coefficient(LV.Coefficients):#######

    from proteus.ctransportCoefficients import ncLevelSetCoefficientsEvaluate

    def __init__(self,):

        mass = {0: {0: 'linear'},}
        advection = {0: {0: 'linear'},}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {0: {0: 'linear'},}
       
        LV.Coefficients.__init__(self)#######

        self.outputQuantDOFs = True
#         self.LUMPED_MASS_MATRIX = True
        self.USE_SBM = USE_SBM

    def save_dirichlet_dofs(self):
        mesh = self.model.mesh
        fes = self.model.u[0].femSpace ##P1
        self.dirichlet_bc_dofs = {'dof': [], 'xyz': [], 'label': [], 'value': []}

        for eN in range(mesh.nElements_global):
            for k in range(fes.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                i = fes.referenceFiniteElement.interpolationConditions.quadrature2DOF_element(k)
                dofN = fes.dofMap.l2g[eN, i]
                x = fes.interpolationPoints[eN, k]

                for ebN_element in range(mesh.nElementBoundaries_element):
                    if fes.referenceFiniteElement.interpolationConditions.definedOnLocalElementBoundary(k, ebN_element) == True:
                        ebN = mesh.elementBoundariesArray[eN, ebN_element]
                        materialFlag = mesh.elementBoundaryMaterialTypes[ebN]
                        if materialFlag in [1, 2, 3, 4, 7]:###############7='obstacle'
                            self.dirichlet_bc_dofs['dof'].append(dofN)
                            self.dirichlet_bc_dofs['xyz'].append(x)
                            self.dirichlet_bc_dofs['label'].append(materialFlag)
                            self.dirichlet_bc_dofs['value'].append(analyticalSolution[0].uOfXT(x,0.0))###########boundary values

        self.dirichlet_bc_dofs['dof'] = numpy.asarray(
            self.dirichlet_bc_dofs['dof'], 'i')
        self.dirichlet_bc_dofs['xyz'] = numpy.asarray(
            self.dirichlet_bc_dofs['xyz'], 'd')
        self.dirichlet_bc_dofs['label'] = numpy.asarray(
            self.dirichlet_bc_dofs['label'], 'd')
        self.dirichlet_bc_dofs['value'] = numpy.asarray(
            self.dirichlet_bc_dofs['value'], 'd')


        
    def attachModels(self, modelList):
        self.model = modelList[0]
        self.mesh = self.model.mesh

        self.save_dirichlet_dofs()
        
        self.ML_new = np.zeros_like(self.model.u[0])
        self.q_v = numpy.zeros((self.model.u[0].dof.shape[0], 3), 'd')
        self.ebqe_v = numpy.zeros_like(self.model.ebqe[('dH', 0, 0)])
        
    def preStep(self, t, firstStep=False):
        # Serious error: Veclocity depends on time, since the mesh is moving.
        #         self.q_v[..., 0] = -2.0 * math.pi * \
        #             (self.model.q['x'][..., 1] - rotation_center[0])
        #         self.q_v[..., 1] = 2.0 * math.pi * \
        #             (self.model.q['x'][..., 0] - rotation_center[1])
        #         self.q_v[..., 0] = 1.0
        #         self.q_v[..., 1] = 0.0
        #end-fixed_mesh
        copyInstructions = {}
        return copyInstructions
#

    def postStep(self, t, firstStep=False):
        
#         self.model.q[('u',0)] = true_solution().uOfXT(self.model.q['x'],0)
#         self.model.q[('grad(u)',0)][:,:,0] = true_solution().duOfXT(self.model.q['x'],0)[0]
#         self.model.q[('grad(u)',0)][:,:,1] = true_solution().duOfXT(self.model.q['x'],0)[1]
        self.model.calculateSolutionAtQuadrature()
        copyInstructions = {}
        return copyInstructions

    def evaluate(self, t, c):
        pass

#     def calculateResidual(self, *args):#This is in Coefficients
# 
#         import burgers2D_GP_ALE as M
# #         import burgers2D_GP_ALE_M2 as M
# 
#         M.getResidual(
#             *args,
#             ad_function=analyticalSolution[0].advection_function,
#             moving_function=analyticalSolution[0].mesh_velocity)  # t is already updated in prestep
#         # Dirichlete BC
#         args[49][self.dirichlet_bc_dofs['dof']
#                  ] = self.dirichlet_bc_dofs['value']
                 
    def apply_dirichlet_bc(self,u):
        u[self.dirichlet_bc_dofs['dof']] = self.dirichlet_bc_dofs['value']


coefficients = My_coefficient()
#===============================================================================
# BC
#===============================================================================

def getDBC_u(x, flag):
    if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'], boundaryTags['obstacle']]:
        return 0.0

dirichletConditions = {0: getDBC_u,}


advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {}

#===============================================================================
# Initial condition
#===============================================================================
class initial_u:
    def uOfXT(self, x, t=0):
        return 0.0

initialConditions = {
                    0: initial_u(),
                    }
#===============================================================================
# Analytic solution
#===============================================================================
class true_solution(object):
    def uOfXT(self, x, t):
        R2 = 0.2**2
        r2 = x[0]*x[0]+x[1]*x[1]
        d = 0.0
        if r2<R2:
            d = 0.5*(R2-r2)
        return d

    def duOfXT(self, x, t):
        R2 = 0.2**2
        r2 = x[0]*x[0]+x[1]*x[1]
        d = [0.0, 0.0]
        if r2<R2:
            d = [-x[0], -x[1]]
        return d

analyticalSolution = {0: true_solution()}