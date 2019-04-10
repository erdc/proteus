from __future__ import division
from past.utils import old_div
from builtins import object
import numpy as np
from proteus.Profiling import logEvent
from proteus.MeshTools import MeshOptions
from proteus.defaults import (Physics_base,
                              Numerics_base,
                              System_base)
# models
from proteus.mprans import (RANS2P,
                            RANS3PF,
                            VOF,
                            RDLS,
                            NCLS,
                            MCorr,
                            CLSVOF,
                            AddedMass,
                            MoveMesh,
                            MoveMeshMonitor,
                            Pres,
                            PresInit,
                            PresInc)
# numerical options
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools,
                     NumericalFlux)

# default values for several models
epsFact = 1.5
sc_uref = 1.
sc_beta = 1.5
shockCapturingFactor = 0.5
minTol = 1e-8


class ParametersHolder:
    """
    """
    def __init__(self, ProblemInstance=None):
        # gain access to problem class if necessary
        self._Problem = ProblemInstance
        # default options
        self.nModels = 0
        self.models_list = []
        self.Models = ParametersModelsHolder(self._Problem)
        self.physical = ParametersPhysical()
        self.mesh = self._Problem.domain.MeshOptions

    def initializeParameters(self):
        all_models = [self.Models.rans2p,
                      self.Models.rans3p,
                      self.Models.clsvof,
                      self.Models.vof,
                      self.Models.ncls,
                      self.Models.rdls,
                      self.Models.moveMeshElastic,
                      self.Models.moveMeshMonitor,
                      self.Models.addedMass,
                      self.Models.pressureInitial,
                      self.Models.pressure,
                      self.Models.pressureIncrement,
                      self.Models.mcorr]
        logEvent('----------')
        logEvent('Mesh Options')
        for key, value in self.mesh.__dict__.items():
            if key[0] != '_':  # do not print hidden attributes
                logEvent('{key}: {value}'. format(key=key, value=value))
        logEvent('----------')
        logEvent('Physical Parameters')
        for key, value in self.physical.__dict__.items():
            if key[0] != '_':  # do not print hidden attributes
                logEvent('{key}: {value}'. format(key=key, value=value))
        logEvent('----------')
        self.nModels = 0
        self.models_list = []
        n_base = Numerics_base()
        p_base = Physics_base()
        for i in range(len(all_models)):
            model = all_models[i]
            if model['index'] is not None:
                model.initializePhysics()
                model.initializeNumerics()
                self.nModels += 1
                self.models_list += [model]
                logEvent('TwoPhaseFlow parameters for model: {name}'.format(name=model['name']))
                logEvent('-----')
                logEvent('{name} PHYSICS'.format(name=model.name))
                logEvent('-----')
                logEvent('COEFFICIENTS OPTIONS')
                for key, value in model.p.CoefficientsOptions.__dict__.items():
                    if key[0] != '_':  # do not print hidden attributes
                        logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('END OF COEFFICIENTS OPTIONS')
                for key, value in model.p.__dict__.items():
                    if key[0] != '_':  # do not print hidden attributes
                        if key in p_base.__dict__.keys():
                            if value != p_base.__dict__[key]:
                                logEvent('(!) {key}: {value}'. format(key=key, value=value))
                            else:
                                logEvent('{key}: {value}'. format(key=key, value=value))
                        else:
                            logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('-----')
                logEvent('{name} NUMERICS'.format(name=model.name))
                logEvent('-----')
                for key, value in model.n.__dict__.items():
                    if key[0] != '_':  # do not print hidden attributes
                        if key in n_base.__dict__.keys():
                            if value != n_base.__dict__[key]:
                                logEvent('(!) {key}: {value}'. format(key=key, value=value))
                            else:
                                logEvent('{key}: {value}'. format(key=key, value=value))
                        else:
                            logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('----------')

        logEvent('-----')
        logEvent('-----')
        logEvent('-----')
        logEvent('-----')
        logEvent('-----')
        n = Numerics_base()
        p = Physics_base()
        for i in range(len(all_models)):
            model = all_models[i]
            if model['index'] is not None:
                logEvent('TwoPhaseFlow parameters for model: {name}'.format(name=model['name']))
                logEvent('-----')
                for key, value in p.__dict__.items():
                    if key[0] != '_' and value != model.p.__dict__[key]:  # do not print hidden attributes
                        logEvent('{key}: {value}'. format(key=key, value=model.p[key]))
                logEvent('-----n')
                for key, value in n.__dict__.items():
                    if key[0] != '_' and value != model.n.__dict__[key]:  # do not print hidden attributes
                        logEvent('{key}: {value}'. format(key=key, value=model.n[key]))


class ParametersModelsHolder:
    """
    """
    def __init__(self, ProblemInstance=None):
        self._Problem = ProblemInstance
        self.rans2p = ParametersModelRANS2P(Problem=self._Problem)
        self.vof = ParametersModelVOF(Problem=self._Problem)
        self.ncls = ParametersModelNCLS(Problem=self._Problem)
        self.rdls = ParametersModelRDLS(Problem=self._Problem)
        self.addedMass = ParametersModelAddedMass(Problem=self._Problem)
        self.moveMeshMonitor = ParametersModelMoveMeshMonitor(Problem=self._Problem)
        self.moveMeshElastic = ParametersModelMoveMeshElastic(Problem=self._Problem)
        self.clsvof = ParametersModelCLSVOF(Problem=self._Problem)
        self.rans3p = ParametersModelRANS3P(Problem=self._Problem)
        self.pressureInitial = ParametersModelPressureInitial(Problem=self._Problem)
        self.pressure = ParametersModelPressure(Problem=self._Problem)
        self.pressureIncrement = ParametersModelPressureIncrement(Problem=self._Problem)
        self.mcorr = ParametersModelMCorr(Problem=self._Problem)


class FreezableClass(object):
    """Base class for all parameters class, enforces attribute freezing
    """
    __frozen = False

    def __init__(self):
        pass

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, val):
        self.__setattr__(key, val)

    def __setattr__(self, key, val):
        if self.__frozen and not hasattr(self, key):
            raise TypeError("{key} is not an option for model {name}".format(key=key, name=self.name))
        object.__setattr__(self, key, val)

    def _freeze(self):
        self.__frozen = True

    def addOption(self, name, value):
        self.__frozen = False
        self.__setattr__(name, value)
        self._freeze()

class ParametersModelBase(FreezableClass):
    """
    """
    def __init__(self,
                 name=None,
                 index=None,
                 Problem=None):
        super(ParametersModelBase, self).__init__()
        self.name = name
        self.index = index
        self.auxiliaryVariables = []
        self._Problem = Problem
        self.p = Physics_base(nd=self._Problem.domain.nd)
        self.p.myTpFlowProblem = self._Problem
        self.p.name = name
        self.p.CoefficientsOptions = FreezableClass()
        self.p._freeze()
        self.n = Numerics_base()
        self.n.name = name
        self.n.ShockCapturingOptions = FreezableClass()
        self.n.SubgridErrorOptions = FreezableClass()
        self.n._freeze()
        # NON LINEAR SOLVERS
        # self.n.fullNewtonFlag = True
        # self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # self.n.levelNonlinearSolver = NonlinearSolvers.Newton
        # self.n.nonlinearSmoother = None
        # self.n.levelNonlinearSolverConvergenceTest = 'r'
        # self.n.nonlinearSolverConvergenceTest = 'r'
        # LINEAR ALGEBRA
        # self.n.matrix = LinearAlgebraTools.SparseMatrix
        # self.n.linearSmoother = None
        # NUMERICAL FLUX
        # self.n.massLumping = False
        # self.n.conservativeFlux = None
        # TOLERANCES
        self.n.nl_atol_res = None
        self.n.l_atol_res = None
        # self.n.linTolFac = 0.001
        # self.n.useEisenstatWalker = False
        # self.n.tolFac = 0.
        # self.n.maxNonlinearIts = 50
        # self.n.maxLineSearches = 0

    def initializePhysics(self):
        self.p.domain = self._Problem.domain
        self.p.nd = self._Problem.nd
        self.p.movingDomain = self._Problem.movingDomain
        self.p.genMesh = self._Problem.Parameters.mesh.genMesh
        # initialize extra parameters
        self._initializePhysics()
        self.p._unfreeze()

    def _initializePhysics(self):
        # to overwrite for each models
        pass

    def initializeNumerics(self):
        self.n.runCFL = self._Problem.cfl
        # MESH
        mesh = self._Problem.Parameters.mesh
        self.n.triangleFlag = mesh.triangleFlag
        self.n.nnx = mesh.nnx
        self.n.nny = mesh.nny
        self.n.nnz = mesh.nnz
        self.n.triangleOptions = mesh.triangleOptions
        self.n.parallelPartitioningType = mesh.parallelPartitioningType
        self.n.nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
        self.n.restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
        # TIME INTEGRATION
        self.n.runCFL = self._Problem.cfl
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.elementQuadrature = FESpace['elementQuadrature']
        self.n.elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
        # SUPERLU
        if self._Problem.useSuperlu:
            self.n.multilevelLinearSolver = LinearSolvers.LU
            self.n.levelLinearSolver = LinearSolvers.LU
        # AUXILIARY VARIABLES
        self.n.auxiliaryVariables = self.auxiliaryVariables
        # initialize extra parameters
        self._initializeNumerics()
        self.n._unfreeze()

    def _initializeNumerics(self):
        # to overwrite for each models
        pass


class ParametersModelRANS2P(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelRANS2P, self).__init__(name='rans2p', index=None,
                                                    Problem=Problem)
        self.timeDiscretization = 'be'
        copts = self.p.CoefficientsOptions
        copts.useMetrics = 1.
        copts.epsFact_viscosity = epsFact
        copts.epsFact_density = epsFact
        copts.forceStrongDirichlet = False
        copts.weak_bc_penalty_constant = 1e6
        copts.useRBLES = 0
        copts.useRANS = 0
        copts.closure = 0
        copts.useVF = 1
        copts.timeDiscretization = None
        copts.timeOrder = 2
        copts.stokes = False
        copts.eb_adjoint_sigma = 1.
        copts.Closure_0_model = None
        copts.Closure_1_model = None
        copts._freeze()
        scopts = self.n.ShockCapturingOptions
        scopts.shockCapturingFactor = shockCapturingFactor
        scopts.lag = True
        scopts._freeze()
        seopts = self.n.SubgridErrorOptions
        seopts.lag = True
        seopts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = RANS2P.LevelModel
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # NUMERICAL FLUX
        self.n.numericalFluxType = RANS2P.NumericalFlux
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'rans2p_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.linTolFac = 0.01
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 50
        self.n.maxLineSearches = 0
        self._freeze()

    def _initializePhysics(self):
        mparams = self._Problem.Parameters.Models
        pparams = self._Problem.Parameters.physical # physical parameters
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEX
        ME_model = mparams.rans2p.index
        assert ME_model is not None, 'rans2p model index was not set!'
        CLSVOF_model = mparams.clsvof.index
        VF_model = mparams.vof.index
        LS_model = mparams.ncls.index
        # POROSITY / RELAXATION
        if hasattr(domain, 'porosityTypes'):
            porosityTypes = domain.porosityTypes
            dragAlphaTypes = domain.dragAlphaTypes
            dragBetaTypes = domain.dragBetaTypes
            epsFact_solid = domain.epsFact_solid
        else:
            porosityTypes = None
            dragAlphaTypes = None
            dragBetaTypes = None
            epsFact_solid = None
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = RANS2P.Coefficients(epsFact=copts.epsFact_viscosity,
                                                  sigma=pparams.surf_tension_coeff,
                                                  rho_0=pparams.densityA,
                                                  nu_0=pparams.kinematicViscosityA,
                                                  rho_1=pparams.densityB,
                                                  nu_1=pparams.kinematicViscosityB,
                                                  g=pparams.gravity,
                                                  nd=nd,
                                                  ME_model=ME_model,
                                                  CLSVOF_model=CLSVOF_model,
                                                  VF_model=VF_model,
                                                  LS_model=LS_model,
                                                  Closure_0_model=copts.Closure_0_model,
                                                  Closure_1_model=copts.Closure_1_model,
                                                  epsFact_density=copts.epsFact_density,
                                                  stokes=copts.stokes,
                                                  useVF=copts.useVF,
                                                  useRBLES=copts.useRBLES,
                                                  useMetrics=copts.useMetrics,
                                                  eb_adjoint_sigma=copts.eb_adjoint_sigma,
                                                  eb_penalty_constant=copts.weak_bc_penalty_constant,
                                                  forceStrongDirichlet=copts.forceStrongDirichlet,
                                                  turbulenceClosureModel=copts.closure,
                                                  movingDomain=self.p.movingDomain,
                                                  porosityTypes=porosityTypes,
                                                  dragAlphaTypes=dragAlphaTypes,
                                                  dragBetaTypes=dragBetaTypes,
                                                  epsFact_solid=epsFact_solid,
                                                  barycenters=domain.barycenters)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['pressure'],
                                    1: IC['vel_u'],
                                    2: IC['vel_v']}
        if nd == 3:
            self.p.initialConditions[3] = IC['vel_w']
        # BOUNDARY CONDITIONS
        boundaryConditions = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: boundaryConditions['pressure_DBC'],
                                          1: boundaryConditions['vel_u_DBC'],
                                          2: boundaryConditions['vel_v_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC'],
                                                      1: boundaryConditions['vel_u_AFBC'],
                                                      2: boundaryConditions['vel_v_AFBC']}
            self.p.diffusiveFluxBoundaryConditions = {0: {},
                                                      1: {1: boundaryConditions['vel_u_DFBC']},
                                                      2: {2: boundaryConditions['vel_v_DFBC']}}
            if nd == 3:
                self.p.dirichletConditions[3] = boundaryConditions['vel_w_DBC']
                self.p.advectiveFluxBoundaryConditions[3] = boundaryConditions['vel_w_AFBC']
                self.p.diffusiveFluxBoundaryConditions[3] = {3: boundaryConditions['vel_w_DFBC']}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
                                          1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                                          2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
                                                      1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                                      2: lambda x, flag: domain.bc[flag].v_advective.init_cython()}
            self.p.diffusiveFluxBoundaryConditions = {0: {},
                                                      1: {1:lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                                      2: {2:lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}
            if nd == 3:
                self.p.dirichletConditions[3] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
                self.p.advectiveFluxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
                self.p.diffusiveFluxBoundaryConditions[3] = {3: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}

    def _initializeNumerics(self):
        nd = self._Problem.domain.nd
        # TIME
        if self.timeDiscretization=='vbdf':
            self.n.timeIntegration = TimeIntegration.VBDF
            self.n.timeOrder = self.p.CoefficientsOptions.timeOrder
        else: #backward euler
            self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['pBasis'],
                            1: FESpace['velBasis'],
                            2: FESpace['velBasis']}
        if nd == 3:
            self.n.femSpaces[3] = FESpace['velBasis']
        # NUMERICAL FLUX
        seopts = self.n.SubgridErrorOptions
        self.n.subgridError = RANS2P.SubgridError(coefficients=self.p.coefficients,
                                                  nd=nd,
                                                  lag=seopts.lag,
                                                  hFactor=FESpace['hFactor'])
        scopts = self.n.ShockCapturingOptions
        self.n.shockCapturing = RANS2P.ShockCapturing(coefficients=self.p.coefficients,
                                                      nd=nd,
                                                      shockCapturingFactor=scopts.shockCapturingFactor,
                                                      lag=scopts.lag)
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.01*self.n.nl_atol_res


class ParametersModelRANS3P(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelRANS3P, self).__init__(name='rans3p', index=None,
                                                    Problem=Problem)
        self.timeOrder = 2
        self.timeDiscretization = 'vbdf'
        copts = self.p.CoefficientsOptions
        copts.useMetrics = 1.
        copts.epsFact_viscosity = epsFact
        copts.epsFact_density = epsFact
        copts.forceStrongDirichlet = False
        copts.ns_sed_forceStrongDirichlet = False
        copts.weak_bc_penalty_constant = 1e6
        copts.useRBLES = 0
        copts.useRANS = 0
        copts.closure = 0
        copts.useVF = 1
        copts.PSTAB = 0
        copts.ARTIFICIAL_VISCOSITY = 3
        copts.INT_BY_PARTS_PRESSURE = 1
        copts.cE = 1.
        copts.cMax = 1.
        copts.forceTerms = None
        copts.stokes = False
        copts.eb_adjoint_sigma = 1.
        copts.MULTIPLY_EXTERNAL_FORCE_BY_DENSITY = 0
        copts.USE_SUPG = False
        copts._freeze()
        scopts = self.n.ShockCapturingOptions
        scopts.shockCapturingFactor = shockCapturingFactor
        scopts.lag = True
        scopts._freeze()
        seopts = self.n.SubgridErrorOptions
        seopts.lag = True
        seopts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = RANS3PF.LevelModel
        # NUMERICAL FLUX
        self.n.numericalFluxType = RANS3PF.NumericalFlux
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'rans3p_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        self.n.nonlinearSolverConvergenceTest = 'rits'
        self.n.levelNonlinearSolverConvergenceTest = 'rits'
        # TOLERANCES
        self.n.linTolFac = 0.
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 1 # This is a linear problem
        self.n.maxLineSearches = 0
        # OTHER
        self.n.addOption('forceTerms', None)
        self.p.addOption('forceTerms', None)
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        mparams = self._Problem.Parameters.Models
        pparams = self._Problem.Parameters.physical # physical parameters
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEX
        VOF_model=mparams.vof.index
        LS_model=mparams.ncls.index
        RD_model=mparams.rdls.index
        MCORR_model=mparams.mcorr.index
        SED_model=None
        VOS_model=None
        CLSVOF_model = mparams.clsvof.index
        V_model = mparams.rans3p.index
        PINC_model = mparams.pressureIncrement.index
        PRESSURE_model = mparams.pressure.index
        # COEFFICIENTS
        Closure_0_model = None
        Closure_1_model = None
        copts = self.p.CoefficientsOptions
        if copts.forceTerms is not None:
            self.p.forceTerms = copts.forceTerms
            copts.MULTIPLY_EXTERNAL_FORCE_BY_DENSITY = 1
        self.p.coefficients = RANS3PF.Coefficients(epsFact=copts.epsFact_viscosity,
                                                   sigma=pparams.surf_tension_coeff,
                                                   rho_0=pparams.densityA,
                                                   nu_0=pparams.kinematicViscosityA,
                                                   rho_1=pparams.densityB,
                                                   nu_1=pparams.kinematicViscosityB,
                                                   g=pparams.gravity,
                                                   nd=nd,
                                                   ME_model=V_model,
                                                   PRESSURE_model=PRESSURE_model,
                                                   SED_model=SED_model,
                                                   CLSVOF_model=CLSVOF_model,
                                                   VOF_model=VOF_model,
                                                   VOS_model=VOS_model,
                                                   LS_model=LS_model,
                                                   Closure_0_model=Closure_0_model,
                                                   Closure_1_model=Closure_1_model,
                                                   epsFact_density=copts.epsFact_density,
                                                   stokes=copts.stokes,
                                                   useVF=copts.useVF,
                                                   useRBLES=copts.useRBLES,
                                                   useMetrics=copts.useMetrics,
                                                   eb_adjoint_sigma=copts.eb_adjoint_sigma,
                                                   eb_penalty_constant=copts.weak_bc_penalty_constant,
                                                   forceStrongDirichlet=copts.forceStrongDirichlet,
                                                   turbulenceClosureModel=copts.closure,
                                                   movingDomain=self.p.movingDomain,
                                                   PSTAB=copts.PSTAB,
                                                   USE_SUPG=copts.USE_SUPG,
                                                   ARTIFICIAL_VISCOSITY=copts.ARTIFICIAL_VISCOSITY,
                                                   INT_BY_PARTS_PRESSURE=copts.INT_BY_PARTS_PRESSURE,
                                                   cE=copts.cE,
                                                   cMax=copts.cMax,
                                                   MULTIPLY_EXTERNAL_FORCE_BY_DENSITY=copts.MULTIPLY_EXTERNAL_FORCE_BY_DENSITY)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['vel_u'],
                                    1: IC['vel_v']}
        if nd == 3:
            self.p.initialConditions[2] = IC['vel_w']
        # BOUNDARY CONDITIONS
        boundaryConditions = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: boundaryConditions['vel_u_DBC'],
                                          1: boundaryConditions['vel_v_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: boundaryConditions['vel_u_AFBC'],
                                                      1: boundaryConditions['vel_v_AFBC']}
            self.p.diffusiveFluxBoundaryConditions = {0: {0: boundaryConditions['vel_u_DFBC']},
                                                      1: {1: boundaryConditions['vel_v_DFBC']}}
            if nd == 3:
                self.p.dirichletConditions[2] = boundaryConditions['vel_w_DBC']
                self.p.advectiveFluxBoundaryConditions[2] = boundaryConditions['vel_w_AFBC']
                self.p.diffusiveFluxBoundaryConditions[2] = {2: boundaryConditions['vel_w_DFBC']}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                                          1: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                                      1: lambda x, flag: domain.bc[flag].v_advective.init_cython()}
            self.p.diffusiveFluxBoundaryConditions = {0: {0:lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                                      1: {1:lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}
            if nd == 3:
                self.p.dirichletConditions[2] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
                self.p.advectiveFLuxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
                self.p.diffusiveFluxBoundaryConditions[2] = {2: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}

    def _initializeNumerics(self):
        self.n.forceTerms = self.p.forceTerms
        nd = self._Problem.domain.nd
        # TIME
        if self.timeDiscretization=='vbdf':
            self.n.timeIntegration = TimeIntegration.VBDF
            self.n.timeOrder = 2
        else: #backward euler
            self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['velBasis'],
                            1: FESpace['velBasis']}
        if nd == 3:
            self.n.femSpaces[2] = FESpace['velBasis']
        # NUMERICAL FLUX
        seopts = self.n.SubgridErrorOptions
        self.n.subgridError = RANS3PF.SubgridError(coefficients=self.p.coefficients,
                                                  nd=nd,
                                                  lag=seopts.lag,
                                                  hFactor=FESpace['hFactor'])
        scopts = self.n.ShockCapturingOptions
        self.n.shockCapturing = RANS3PF.ShockCapturing(coefficients=self.p.coefficients,
                                                      nd=nd,
                                                      shockCapturingFactor=scopts.shockCapturingFactor,
                                                      lag=scopts.lag)
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.01*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.1*self.n.nl_atol_res


class ParametersModelPressure(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelPressure, self).__init__(name='pressure', index=None,
                                                      Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.useRotationalForm = False
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = Pres.LevelModel
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'pressure_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # NUMERICAL FLUX 
        self.n.numericalFluxType = NumericalFlux.ConstantAdvection_exterior
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.linTolFac = 0.
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        PRESSURE_model = mparams.pressure.index
        V_model = mparams.rans3p.index
        PINC_model = mparams.pressureIncrement.index
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = Pres.Coefficients(modelIndex=PRESSURE_model,
                                                fluidModelIndex=V_model,
                                                pressureIncrementModelIndex=PINC_model,
                                                useRotationalForm=copts.useRotationalForm)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['pressure']}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['pressure_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: BC['pressure_AFBC']}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython()}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.stepController = StepControl.FixedStep
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['pBasis']}
        # TOLERANCE
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.01*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.1*self.n.nl_atol_res

class ParametersModelPressureInitial(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelPressureInitial, self).__init__(name='pressureInitial', index=None,
                                                             Problem=Problem)
        # NUMERICAL FLUX 
        self.n.numericalFluxType = NumericalFlux.ConstantAdvection_exterior
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linearSmoother = LinearSolvers.NavierStokesPressureCorrection # pure neumann laplacian solver
        self.n.linear_solver_options_prefix = 'pinit_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.linTolFac = 0.
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        PRESSURE_model = mparams.pressure.index
        V_model = mparams.rans3p.index
        PINIT_model = mparams.pressureInitial.index
        # COEFFICIENTS
        self.p.coefficients=PresInit.Coefficients(nd=nd,
                                                  modelIndex=PINIT_model,
                                                  fluidModelIndex=V_model,
                                                  pressureModelIndex=PRESSURE_model)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['pressure']}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['pressure_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: BC['pressure_AFBC']}
            self.p.diffusiveFluxBoundaryConditions = {0:{0: BC['pressure_increment_DFBC']}}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython()}
            self.p.diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython()}}
        # freeze attributes
        self._freeze()

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.stepController = StepControl.FixedStep
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['pBasis']}
        # LINEAR ALGEBRA
        if self._Problem.useSuperlu:
            self.n.linearSmoother = None
        # TOLERANCE
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.01*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.1*self.n.nl_atol_res

class ParametersModelPressureIncrement(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelPressureIncrement, self).__init__(name='pressureIncrement', index=None,
                                                               Problem=Problem)
        # LEVEL MODEL
        self.p.LevelModelType = PresInc.LevelModel
        # NUMERICAL FLUX 
        self.n.numericalFluxType = PresInc.NumericalFlux
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'phi_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.linTolFac = 0.
        self.n.maxNonlinearIts = 50
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        pparams = self._Problem.Parameters.physical
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        V_model = mparams.rans3p.index
        PINC_model = mparams.pressureIncrement.index
        # COEFFICIENTS
        self.p.coefficients = PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*pparams.densityB,
                                                   rho_s_min = (1.0-1.0e-8)*pparams.densityA,
                                                   nd = nd,
                                                   modelIndex=PINC_model,
                                                   fluidModelIndex=V_model,
                                                   fixNullSpace=False)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['pressure_increment']}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['pressure_increment_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: BC['pressure_increment_AFBC']}
            self.p.diffusiveFluxBoundaryConditions = {0:{0: BC['pressure_increment_DFBC']}}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].pInt_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pInc_advective.init_cython()}
            self.p.diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython()}}
        # freeze attributes
        self._freeze()

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.stepController = StepControl.FixedStep
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['pBasis']}
        # LINEAR ALGEBRA
        if self._Problem.useSuperlu:
            self.n.linearSmoother = None
        # TOLERANCE
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.01*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.01*self.n.nl_atol_res


class ParametersModelCLSVOF(ParametersModelBase):
    def __init__(self, Problem):
        super(ParametersModelCLSVOF, self).__init__(name='clsvof', index=None,
                                                    Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.useMetrics = 1.
        copts.epsFactHeaviside = epsFact
        copts.epsFactDirac = epsFact
        copts.epsFactRedist = 0.33
        copts.lambdaFact = 10.
        copts.outputQuantDOFs = True
        copts.computeMetrics = 1
        copts.computeMetricsForBubble = False
        copts.eps_tolerance_clsvof = False
        copts.disc_ICs = False
        copts._freeze()
        # LEVEL MODEL 
        self.p.LevelModelType = CLSVOF.LevelModel
        # NUMERICAL FLUX
        self.n.numericalFluxType = CLSVOF.NumericalFlux
        # NON LINEAR SOLVER
        self.n.levelNonlinearSolver = NonlinearSolvers.CLSVOFNewton
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        self.n.addOption('updateJacobian', True)
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'clsvof_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 50
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        CLSVOF_model = mparams.clsvof.index
        V_model = mparams.rans2p.index
        if V_model is None:
            V_model = mparams.rans3p.index
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = CLSVOF.Coefficients(V_model=V_model,
                                                  ME_model=CLSVOF_model,
                                                  useMetrics=copts.useMetrics,
                                                  epsFactHeaviside=copts.epsFactHeaviside,
                                                  epsFactDirac=copts.epsFactDirac,
                                                  epsFactRedist=copts.epsFactRedist,
                                                  lambdaFact=copts.lambdaFact,
                                                  outputQuantDOFs=copts.outputQuantDOFs,
                                                  computeMetrics=copts.computeMetrics,
                                                  disc_ICs=copts.disc_ICs)
        self.p.coefficients.variableNames = ['phi']
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['clsvof']}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['clsvof_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: BC['clsvof_AFBC']}
            self.p.diffusiveFluxBoundaryConditions = {0:{0: BC['clsvof_DFBC']}}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vof_advective.init_cython()}
            self.p.diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].clsvof_diffusive.init_cython()}}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # TOLERANCE
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            if self.p.CoefficientsOptions.eps_tolerance_clsvof:
                self.nl_atol_res = 1e-12
            else:
                self.n.nl_atol_res = max(minTol, 0.01*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.1*self.n.nl_atol_res


class ParametersModelVOF(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelVOF, self).__init__(name='vof', index=None,
                                                 Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.useMetrics = True
        copts.checkMass = True
        copts.sc_uref = sc_uref
        copts.sc_beta = sc_beta
        copts.epsFact = epsFact
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = VOF.LevelModel
        # TIME INTEGRATION
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController  = StepControl.Min_dt_cfl_controller
        # NUMERICS
        self.n.ShockCapturingOptions.shockCapturingFactor = shockCapturingFactor
        self.n.ShockCapturingOptions.lag = True
        # NUMERICAL FLUX
        self.n.numericalFluxType = VOF.NumericalFlux
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'vof_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.maxNonlinearIts = 50
        self.n.maxLineSearches = 0
        self.n.tolFac = 0.
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        ME_model = mparams.vof.index
        assert ME_model is not None, 'ls model index was not set!'
        if mparams.rans2p.index is not None:
            V_model = mparams.rans2p.index
        elif mparams.rans3p.index is not None:
            V_model = mparams.rans3p.index
        else:
            assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'
        RD_model = mparams.rdls.index
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = VOF.Coefficients(V_model=V_model,
                                               RD_model=RD_model,
                                               ME_model=ME_model,
                                               checkMass=copts.checkMass,
                                               useMetrics=copts.useMetrics,
                                               epsFact=copts.epsFact,
                                               sc_uref=copts.sc_uref,
                                               sc_beta=copts.sc_beta,
                                               movingDomain=self.p.movingDomain)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['vof']}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['vof_DBC']}
            self.p.advectiveFluxBoundaryConditions = {0: BC['vof_AFBC']}
        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython()}
            self.p.advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vof_advective.init_cython()}
        self.p.diffusiveFluxBoundaryConditions = {0: {}}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # NUMERICAL FLUX
        self.n.subgridError = VOF.SubgridError(coefficients=self.p.coefficients,
                                               nd=nd)
        scopts = self.n.ShockCapturingOptions
        self.n.shockCapturing = VOF.ShockCapturing(coefficients=self.p.coefficients,
                                                   nd=nd,
                                                   shockCapturingFactor=scopts.shockCapturingFactor,
                                                   lag=scopts.lag)
        # TOLERANCE
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res


class ParametersModelNCLS(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelNCLS, self).__init__(name='ncls', index=None,
                                                  Problem=Problem)
        # PHYSICS
        copts = self.p.CoefficientsOptions
        copts.useMetrics = True
        copts.checkMass = False
        copts.sc_uref = sc_uref
        copts.sc_beta = sc_beta
        copts.epsFact = epsFact
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = NCLS.LevelModel
        # TIME INTEGRATION
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController  = StepControl.Min_dt_cfl_controller
        # NUMERICAL FLUX
        self.n.numericalFluxType = NCLS.NumericalFlux
        self.n.ShockCapturingOptions.shockCapturingFactor = shockCapturingFactor
        self.n.ShockCapturingOptions.lag = True
        # NUMERICAL FLUX
        self.n.numericalFluxType = NCLS.NumericalFlux
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'ncls_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.maxNonlinearIts = 50
        self.n.maxLineSearches = 0
        self.n.tolFac = 0.
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        ME_model = mparams.ncls.index
        assert ME_model is not None, 'ls model index was not set!'
        if mparams.rans2p.index is not None:
            V_model = mparams.rans2p.index
        elif mparams.rans3p.index is not None:
            V_model = mparams.rans3p.index
        else:
            assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'
        RD_model = mparams.rdls.index
        copts = self.p.CoefficientsOptions
        self.p.coefficients = NCLS.Coefficients(V_model=V_model,
                                                RD_model=RD_model,
                                                ME_model=ME_model,
                                                checkMass=copts.checkMass,
                                                useMetrics=copts.useMetrics,
                                                epsFact=copts.epsFact,
                                                sc_uref=copts.sc_uref,
                                                sc_beta=copts.sc_beta,
                                                movingDomain=self.p.movingDomain)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['ncls']}
        # BOUNDARY CONDITIONS
        self.p.dirichletConditions = {0: lambda x, flag: None}
        self.p.advectiveFluxBoundaryConditions = {}
        self.p.diffusiveFluxBoundaryConditions = {0: {}}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # NUMERICAL FLUX
        self.n.subgridError = NCLS.SubgridError(coefficients=self.p.coefficients,
                                                nd=nd)
        scopts = self.n.ShockCapturingOptions
        self.n.shockCapturing = NCLS.ShockCapturing(coefficients=self.p.coefficients,
                                                    nd=nd,
                                                    shockCapturingFactor=scopts.shockCapturingFactor,
                                                    lag=scopts.lag)
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res


class ParametersModelRDLS(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelRDLS, self).__init__(name='rdls', index=None,
                                                  Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.useMetrics = True
        copts.applyRedistancing = True
        copts.backgroundDiffusionFactor = 0.01
        copts.epsFact = 0.33
        copts.ELLIPTIC_REDISTANCING = 0
        copts._freeze()
        scopts = self.n.ShockCapturingOptions
        scopts.shockCapturingFactor = 0.9
        scopts.lag = False
        scopts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = RDLS.LevelModel
        # TIME INTEGRATION
        self.n.timeIntegration = TimeIntegration.NoIntegration
        self.n.stepController = StepControl.Newton_controller
        # NONLINEAR SOLVERS
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        self.n.nonlinearSmoother = NonlinearSolvers.NLGaussSeidel
        # NUMERICAL FLUX
        self.n.numericalFluxType = NumericalFlux.DoNothing
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'rdls_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 25
        self.n.maxLineSearches = 25
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        nModelId = mparams.ncls.index
        assert nModelId is not None, 'ncls model index was not set!'
        rdModelId = mparams.rdls.index
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = RDLS.Coefficients(applyRedistancing=copts.applyRedistancing,
                                                epsFact=copts.epsFact,
                                                nModelId=nModelId,
                                                rdModelId=rdModelId,
                                                useMetrics=copts.useMetrics,
                                                backgroundDiffusionFactor=copts.backgroundDiffusionFactor,
                                                ELLIPTIC_REDISTANCING=copts.ELLIPTIC_REDISTANCING)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        self.p.initialConditions = {0: IC['rdls']}
        # BOUNDARY CONDITIONS
        self.p.dirichletConditions = {0: lambda x, flag: None}
        self.p.weakDirichletConditions = {0: RDLS.setZeroLSweakDirichletBCsSimple}
        self.p.advectiveFluxBoundaryConditions = {}
        self.p.diffusiveFluxBoundaryConditions = {0: {}}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # TIME
        self.n.timeIntegration = TimeIntegration.NoIntegration
        self.n.stepController = StepControl.Newton_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # NON LINEAR SOLVER
        self.n.nonlinearSmoother = NonlinearSolvers.NLGaussSeidel
        # NUMERICAL FLUX
        self.n.subgridError = RDLS.SubgridError(coefficients=self.p.coefficients,
                                                nd=nd)
        scopts = self.n.ShockCapturingOptions
        self.n.shockCapturing = RDLS.ShockCapturing(coefficients=self.p.coefficients,
                                                    nd=nd,
                                                    shockCapturingFactor=scopts.shockCapturingFactor,
                                                    lag=scopts.lag)
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.01*mesh.he)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res

class ParametersModelMCorr(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelMCorr, self).__init__(name='mcorr', index=None,
                                                   Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.useMetrics = True
        copts.checkMass = False
        copts.applyCorrection = True
        copts.epsFactHeaviside = epsFact
        copts.epsFactDirac = epsFact
        copts.epsFactDiffusion = 10.
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = MCorr.LevelModel
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController  = StepControl.Min_dt_cfl_controller
        # NUMERICAL FLUX
        self.n.numericalFluxType = NumericalFlux.DoNothing
        # NON LINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'mcorr_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.linTolFac = 0.
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 50
        self.n.maxLineSearches = 0
        self.n.useEisenstatWalker = True
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        ME_model = mparams.mcorr.index
        assert ME_model is not None, 'vof model index was not set!'
        LS_model = mparams.ncls.index
        VOF_model = mparams.vof.index
        if mparams.rans2p.index is not None:
            V_model = mparams.rans2p.index
        elif mparams.rans3p.index is not None:
            V_model = mparams.rans3p.index
        else:
            assert mparams.rans2p.index is not None or params.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = MCorr.Coefficients(LSModel_index=LS_model,
                                                 V_model=V_model,
                                                 me_model=ME_model,
                                                 VOFModel_index=VOF_model,
                                                 applyCorrection=copts.applyCorrection,
                                                 nd=nd,
                                                 checkMass=copts.checkMass,
                                                 useMetrics=copts.useMetrics,
                                                 epsFactHeaviside=copts.epsFactHeaviside,
                                                 epsFactDirac=copts.epsFactDirac,
                                                 epsFactDiffusion=copts.epsFactDiffusion)
        # INITIAL CONDITIONS
        class zero_phi:
            def __init__(self):
                pass
            def uOfX(self,X):
                return 0.0
            def uOfXT(self,X,t):
                return 0.0
        self.p.initialConditions  = {0:zero_phi()}
        # BOUNDARY CONDITIONS
        # N/A

    def _initializeNumerics(self):
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.0001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res


class ParametersModelAddedMass(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelAddedMass, self).__init__(name='addedMass', index=None,
                                                       Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.flags_rigidbody = None
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = AddedMass.LevelModel
        # TIME
        self.n.timeIntegration = TimeIntegration.BackwardEuler_cfl
        self.n.stepController = StepControl.Min_dt_cfl_controller
        # NONLINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.AddedMassNewton
        self.n.levelNonlinearSolver = NonlinearSolvers.AddedMassNewton
        self.n.nonlinearSmoother = NonlinearSolvers.AddedMassNewton
        # NUMERICAL FLUX
        self.n.numericalFluxType = AddedMass.NumericalFlux
        # LINEAR ALGEBRA
        self.n.linearSmoother = LinearSolvers.NavierStokesPressureCorrection
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'am_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.linTolFac = 0.
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 1
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        if mparams.rans2p.index is not None:
            V_model = mparams.rans2p.index
        elif mparams.rans3p.index is not None:
            V_model = mparams.rans3p.index
        else:
            assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with addedMass'
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = AddedMass.Coefficients(nd=nd,
                                                     V_model=V_model,
                                                     barycenters=domain.barycenters,
                                                     flags_rigidbody=copts.flags_rigidbody)
        # INITIAL CONDITIONS
        IC = self._Problem.initialConditions
        class dp_IC:
            def uOfXT(self, x, t):
                return 0.0
        self.p.initialConditions = {0: dp_IC()}
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].pAddedMass_dirichlet.init_cython()}
        self.p.advectiveFluxBoundaryConditions = {}
        def getFlux_am(x, flag):
            #the unit rigid motions will applied internally
            #leave this set to zero
            return lambda x,t: 0.0
        self.p.diffusiveFluxBoundaryConditions = {0: {0: getFlux_am}}

    def _initializeNumerics(self):
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.0001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = self.n.nl_atol_res

class ParametersModelMoveMeshMonitor(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelMoveMeshMonitor, self).__init__(name='moveMeshMonitor', index=None,
                                                             Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.func = lambda x, t: 1000.
        copts.he_min = 0.
        copts.he_max = 1000.
        copts.epsFact = epsFact
        copts.epsTimeStep = 0.1
        copts.nSmoothOut = 0.
        copts.nSmoothIn = 0.
        copts.grading = 1.1
        copts.grading_type = 2
        copts.resetNodeVelocityArray = None
        copts.useLS = True
        copts.fixedNodeMaterialTypes = None
        copts.fixedElementMaterialTypes = None
        copts.noNodeVelocityNodeMaterialTypes = None
        copts.scale_with_nd = False
        copts.do_firstStep = False
        copts.ntimes_solved = 1
        copts._freeze()
        # TIME INTEGRATION
        self.n.timeIntegration = TimeIntegration.NoIntegration
        # NONLINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.MoveMeshMonitorNewton
        self.n.levelNonlinearSolver = NonlinearSolvers.MoveMeshMonitorNewton
        self.n.nonlinearSmoother = NonlinearSolvers.MoveMeshMonitorNewton
        # NUMERICAL FLUX
        self.n.numericalFluxType = NumericalFlux.Diffusion_SIPG_exterior
        # LINEAR ALGEBRA
        self.n.linearSmoother = LinearSolvers.NavierStokesPressureCorrection
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'mesh2_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.linTolFac = 0.
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 1
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        ME_MODEL = mparams.moveMeshMonitor.index
        assert ME_MODEL is not None, 'moveMeshMonitor model index was not set!'
        if myparams.useLS is True:
            LS_MODEL = mparams.ncls.index
        else:
            LS_MODEL = None
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = MoveMeshMonitor.Coefficients(func=copts.func,
                                                           nd=nd,
                                                           he_max=copts.he_max,
                                                           he_min=copts.he_min,
                                                           ntimes_solved=copts.ntimes_solved,
                                                           LS_MODEL=LS_MODEL,
                                                           ME_MODEL=ME_MODEL,
                                                           fixedNodeMaterialTypes=copts.fixedNodeMaterialTypes,
                                                           fixedElementMaterialTypes=copts.fixedElementMaterialTypes,
                                                           noNodeVelocityNodeMaterialTypes=copts.noNodeVelocityNodeMaterialTypes,
                                                           nSmoothOut=copts.nSmoothOut,
                                                           nSmoothIn=copts.nSmoothIn,
                                                           epsTimeStep=copts.epsTimeStep,
                                                           epsFact_density=copts.epsFact,
                                                           grading=copts.grading,
                                                           grading_type=copts.grading_type,
                                                           scale_with_nd=copts.scale_with_nd,
                                                           do_firstStep=copts.do_firstStep)
        # INITIAL CONDITIONS
        self.p.initialConditions = None
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        self.p.dirichletConditions = {0: lambda x, flag: None}
        # self.p.advectiveFluxBoundaryConditions = {}
        self.p.diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: None}}

    def _initializeNumerics(self):
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['lsBasis']}
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.0001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res


class ParametersModelMoveMeshElastic(ParametersModelBase):
    """
    """
    def __init__(self, Problem):
        super(ParametersModelMoveMeshElastic, self).__init__(name='moveMeshElastic', index=None,
                                                             Problem=Problem)
        copts = self.p.CoefficientsOptions
        copts.E = 1.
        copts.nu = 0.3
        copts._freeze()
        # LEVEL MODEL
        self.p.LevelModelType = MoveMesh.LevelModel
        # TIME INTEGRATION
        self.n.timeIntegration = TimeIntegration.NoIntegration
        # NONLINEAR SOLVER
        self.n.multilevelNonlinearSolver = NonlinearSolvers.Newton
        # NUMERICAL FLUX
        self.n.numericalFluxType = NumericalFlux.Stress_IIPG_exterior
        # LINEAR ALGEBRA
        self.n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.levelLinearSolver = LinearSolvers.KSP_petsc4py
        self.n.linear_solver_options_prefix = 'mesh_'
        self.n.linearSolverConvergenceTest = 'r-true'
        # TOLERANCES
        self.n.tolFac = 0.
        self.n.maxNonlinearIts = 4
        self.n.maxLineSearches = 0
        # freeze attributes
        self._freeze()

    def _initializePhysics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # NUM PARAMS
        nMediaTypes = len(domain.regionFlags)  # (!) should be region flags
        smTypes = np.zeros((nMediaTypes+1, 2), 'd')
        smFlags = np.zeros((nMediaTypes+1,), 'i')
        copts = self.p.CoefficientsOptions
        smTypes[:, 0] = copts.E
        smTypes[:, 1] = copts.nu
        # MODEL INDEXING
        mparams = self._Problem.Parameters.Models
        ME_model = mparams.moveMeshElastic.index
        assert ME_model is not None, 'vof model index was not set!'
        if mparams.rans2p.index is not None:
            V_model = mparams.rans2p.index
        elif mparams.rans3p.index is not None:
            V_model = mparams.rans3p.index
        else:
            assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'
        # COEFFICIENTS
        copts = self.p.CoefficientsOptions
        self.p.coefficients = MoveMesh.Coefficients(nd=nd,
                                                    V_model=V_model,
                                                    modelType_block=smFlags,
                                                    modelParams_block=smTypes,
                                                    meIndex=ME_model)
        # INITIAL CONDITIONS
        self.p.initialConditions = None
        # BOUNDARY CONDITIONS
        BC = self._Problem.boundaryConditions
        if domain.useSpatialTools is False or self._Problem.useBoundaryConditionsModule is False:
            self.p.dirichletConditions = {0: BC['hx'],
                                          1: BC['hy']}
            self.p.stressFluxBoundaryConditions = {0: BC['u_stress'],
                                                   1: BC['v_stress']}
            if nd == 3:
                self.p.dirichletConditions[2] = BC['hz']
                self.p.stressFluxBoundaryConditions[2] = BC['w_stress']

        else:
            self.p.dirichletConditions = {0: lambda x, flag: domain.bc[flag].hx_dirichlet.init_cython(),
                                          1: lambda x, flag: domain.bc[flag].hy_dirichlet.init_cython()}
            self.p.stressFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_stress.init_cython(),
                                                   1: lambda x, flag: domain.bc[flag].v_stress.init_cython()}
            if nd == 3:
                self.p.dirichletConditions[2] = lambda x, flag: domain.bc[flag].hz_dirichlet.init_cython()
                self.p.stressFluxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_stress.init_cython()
        self.p.fluxBoundaryConditions = {0: 'noFlow',
                                         1: 'noFlow'}
        self.p.advectiveFluxBoundaryConditions = {}

        self.p.diffusiveFluxBoundaryConditions = {0: {},
                                                  1: {}}
        if nd == 3:
            self.p.fluxBoundaryConditions[2] = 'noFlow'
            self.p.diffusiveFluxBoundaryConditions[2] = {}

    def _initializeNumerics(self):
        domain = self._Problem.domain
        nd = domain.nd
        # FINITE ELEMENT SPACES
        FESpace = self._Problem.FESpace
        self.n.femSpaces = {0: FESpace['velBasis'],
                            1: FESpace['velBasis']}
        if nd == 3:
            self.n.femSpaces[2] = FESpace['velBasis']
        # TOLERANCES
        mesh = self._Problem.Parameters.mesh
        if self.n.nl_atol_res is None:
            self.n.nl_atol_res = max(minTol, 0.0001*mesh.he**2)
        if self.n.l_atol_res is None:
            self.n.l_atol_res = 0.001*self.n.nl_atol_res


class ParametersPhysical(FreezableClass):
    """
    """
    def __init__(self):
        super(ParametersPhysical, self).__init__()
        self.name = 'physical'
        self.densityA = 998.2
        self.densityB = 1.205
        self.kinematicViscosityA = 1.004e-6
        self.kinematicViscosityB = 1.500e-5
        self.surf_tension_coeff = 72.8e-3
        self.gravity = [0., -9.81, 0.]
        # freeze attributes
        self._freeze()

