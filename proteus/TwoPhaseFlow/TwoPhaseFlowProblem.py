from __future__ import division
from past.utils import old_div
from proteus import FemTools as ft
from proteus import MeshTools as mt
from proteus.Profiling import logEvent

class TwoPhaseFlowProblem:
    """ TwoPhaseFlowProblem """

    def __init__(self,
                 ns_model=0, #0: rans2p, 1: rans3p
                 nd=2,
                 # TIME STEPPING #
                 cfl=0.33,
                 outputStepping=None,
                 # DOMAIN AND MESH #
                 structured=False,
                 he=None,
                 nnx=None,
                 nny=None,
                 nnz=None,
                 domain=None,
                 triangleFlag=1,
                 # INITIAL CONDITIONS #
                 initialConditions=None,
                 # BOUNDARY CONDITIONS #
                 boundaryConditions=None,
                 # OTHERS #
                 useSuperlu=None):
        """ Constructor for structured meshes  """
        # ***** SET OF ASSERTS ***** #
        assert ns_model in [0,1], "ns_model={0,1} for rans2p or rans3p respectively"
        assert nd in [2,3], "nd={2,3}"
        assert cfl <= 1, "Choose cfl <= 1"
        assert isinstance (outputStepping,OutputStepping), "Provide an object from the OutputStepping class"
        assert type(he)==float , "Provide (float) he (characteristic mesh size)"
        if structured:
            assert type(nnx)==int and type(nny)==int, "Provide (int) nnx and nny"
            if nd==3:
                assert type(nnz)==int, "Provide (int) nnz"
        assert domain is not None, "Provide a domain"
        assert triangleFlag in [0,1,2], "triangleFlag must be 1, 2 or 3"
        if initialConditions is not None:
            assert type(initialConditions)==dict, "Provide dict of initial conditions"
            self.assert_initialConditions(ns_model,nd,initialConditions)
        if boundaryConditions is not None:
            assert type(boundaryConditions)==dict, "Provide dict of boundary conditions"
            self.assert_boundaryConditions(ns_model,nd,boundaryConditions)

        # ***** SAVE PARAMETERS ***** #
        self.ns_model=ns_model
        self.nd=nd
        self.cfl=cfl
        self.outputStepping=outputStepping.getOutputStepping()
        self.he=he
        self.nnx=nnx
        self.nny=nny
        self.nnz=nnz
        self.domain=domain
        self.triangleFlag=triangleFlag
        self.initialConditions=initialConditions
        self.boundaryConditions=boundaryConditions
        self.useSuperlu = useSuperlu
        self.Parameters = Parameters()

        # ***** CHOOSE SOME DEFAULT OPTIONS FOR PARALLEL RUNS ***** #
        self.parallelPartitioningType = mt.MeshParallelPartitioningTypes.node
        #parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
        self.nLayersOfOverlapForParallel = 0

        # ***** CREATE FINITE ELEMENT SPACES ***** #
        self.FESpace = FESpace(self.ns_model,self.nd).getFESpace()

        # ***** DEFINE PHYSICAL AND NUMERICAL PARAMETERS ***** #
        self.physical_parameters = default_physical_parameters
        self.rans2p_parameters = default_rans2p_parameters
        self.rans3p_parameters = default_rans3p_parameters
        self.clsvof_parameters = default_clsvof_parameters

    def assert_initialConditions(self,ns_model,nd,initialConditions):
        assert 'pressure' in initialConditions, 'Provide pressure in ICs'
        assert 'vel_u' in initialConditions, 'Provide vel_u in ICs'
        assert 'vel_v' in initialConditions, 'Provide vel_v in ICs'
        if nd==3:
            assert 'vel_w' in initialConditions, 'Provide vel_w in ICs'
        if ns_model == 1:
            assert 'clsvof' in initialConditions, 'Provide clsvof in ICs'
            if ns_model==1: #rans3p
                assert 'pressure_increment' in initialConditions, 'Provide pressure_increment in ICs'
    #
    def assert_boundaryConditions(self,ns_model,nd,boundaryConditions):
        # check dirichlet BCs
        assert 'pressure_DBC' in boundaryConditions, "Provide pressure_DBC"
        assert 'vel_u_DBC' in boundaryConditions, "Provide vel_u_DBC"
        assert 'vel_v_DBC' in boundaryConditions, "Provide vel_v_DBC"
        if nd==3:
            assert 'vel_w_DBC' in boundaryConditions, "Provide vel_w_DBC"
        assert 'clsvof_DBC' in boundaryConditions, "Provide clsvof_DBC"
        # check advective flux BCs
        assert 'pressure_AFBC' in boundaryConditions, "Provide pressure_AFBC"
        assert 'vel_u_AFBC' in boundaryConditions, "Provide vel_u_AFBC"
        assert 'vel_v_AFBC' in boundaryConditions, "Provide vel_v_AFBC"
        if nd==3:
            assert 'vel_w_AFBC' in boundaryConditions, "Provide vel_w_AFBC"
        assert 'clsvof_AFBC' in boundaryConditions, "Provide clsvof_AFBC"
        # check diffusive flux BCs
        assert 'vel_u_DFBC' in boundaryConditions, "Provide vel_u_DFBC"
        assert 'vel_v_DFBC' in boundaryConditions, "Provide vel_v_DFBC"
        if nd==3:
            assert 'vel_w_DFBC' in boundaryConditions, "Provide vel_w_DFBC"
        assert 'clsvof_DFBC' in boundaryConditions, "Provide clsvof_DFBC"
        if ns_model==1: #rans3p
            # check dirichlet BCs
            assert 'pressure_increment_DBC' in boundaryConditions, "Provide pressure_increment_DBC"
            # check advective flux BCs
            assert 'pressure_increment_AFBC' in boundaryConditions,"Provide pressure_increment_AFBC"
            # check diffusive flux BCs
            assert 'pressure_increment_DFBC' in boundaryConditions,"Provide pressure_increment_DFBC"

class OutputStepping:
    """
    OutputStepping handles how often the solution is outputted.
    """
    def __init__(self,
                 final_time,
                 dt_init=0.001,
                 dt_output=None,
                 nDTout=None,
                 dt_fixed=None):
        self.final_time=final_time
        self.dt_init=dt_init
        assert not (dt_output is None and nDTout is None), "Provide dt_output or nDTout"
        self.dt_output=dt_output
        self.nDTout = nDTout
        self.dt_fixed = dt_fixed
    #
    def getOutputStepping(self):
        # COMPUTE dt_init #
        dt_init = min(0.1 * self.dt_output, self.dt_init)
        if self.nDTout is None:
            self.nDTout = int(round(old_div(self.final_time, self.dt_output)))
        else:
            self.dt_output = float(self.final_time)/float(self.nDTout)
        #
        return {'final_time': self.final_time,
                'dt_init': dt_init,
                'dt_output': self.dt_output,
                'nDTout': self.nDTout,
                'dt_fixed': self.dt_fixed}
#

class FESpace:
    """
    Create FE Spaces.
    """
    def __init__(self,ns_model,nd):
        assert ns_model == 0 or ns_model == 1, "ns_model must be 0 (rans2p) or 1 (rans3p)"
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        self.ns_model=ns_model
        self.nd=nd
        # For now we just support rans2p with: p1-p1 and rans3p with: p2-p1
        if ns_model == 0: #rans2p
            self.velSpaceOrder=1
            self.pSpaceOrder=1
        else: #rans3p
            self.velSpaceOrder=2
            self.pSpaceOrder=1

    def getFESpace(self):
        ##################
        # VELOCITY SPACE #
        ##################
        if self.velSpaceOrder == 1: # p1 space
            hFactor = 1.0
            velBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            hFactor = 0.5
            velBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ##################
        # PRESSURE SPACE #
        ##################
        if self.pSpaceOrder == 1: # p1 space
            pBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            pBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ###################
        # LEVEL SET SPACE #
        ###################
        lsBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis # p1 space
        ###################
        # QUADRATURE RULE #
        ###################
        if max(self.velSpaceOrder,self.pSpaceOrder)==1:
            elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 3)
            elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 3)
        else:
            elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 5)
            elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 5)
        #
        return {'lsBasis': lsBasis,
                'hFactor': hFactor,
                'velBasis': velBasis,
                'pBasis': pBasis,
                'elementQuadrature': elementQuadrature,
                'elementBoundaryQuadrature': elementBoundaryQuadrature}

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
default_physical_parameters ={'densityA': 998.2,
                              'viscosityA': 1.004e-6,
                              'densityB': 1.205,
                              'viscosityB': 1.500e-5,
                              'surf_tension_coeff': 72.8E-3,
                              'gravity': [0.0, -9.8, 0.0]}

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
default_rans2p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 3.0,
                             'epsFact_density': 3.0,
                             'ns_forceStrongDirichlet': False,
                             'weak_bc_penalty_constant': 1.0E6,
                             'useRBLES': 0.0,
                             'useRANS': 0.0,
                             'ns_closure': 0,
                             'useVF': 1.0,
                             'ns_shockCapturingFactor': 0.25,
                             'ns_lag_shockCapturing': True,
                             'ns_lag_subgridError': True,
                             'timeDiscretization': 'vbdf',
                             'timeOrder': 2}
default_rans3p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 3.0,
                             'epsFact_density': 3.0,
                             'ns_forceStrongDirichlet': False,
                             'ns_sed_forceStrongDirichlet': False,
                             'weak_bc_penalty_constant': 1.0E6,
                             'useRBLES': 0.0,
                             'useRANS': 0.0,
                             'ns_closure': 0,
                             'useVF': 1.0,
                             'ns_shockCapturingFactor': 0.5,
                             'ns_lag_shockCapturing': True,
                             'ns_lag_subgridError': True,
                             'timeDiscretization': 'vbdf',
                             'timeOrder': 2,
                             'PSTAB': 0,
                             'USE_SUPG': True,
                             'ARTIFICIAL_VISCOSITY': 2,
                             'cE': 1.0,
                             'cMax': 1.0}
default_clsvof_parameters = {'useMetrics': 1.0,
                             'epsFactHeaviside': 1.5,
                             'epsFactDirac': 1.5,
                             'epsFactRedist': 0.33,
                             'lambdaFact': 10.0,
                             'outputQuantDOFs': True,
                             'computeMetrics': 1,
                             'eps_tolerance_clsvof': False}

class Parameters:
    def __init__(self):
        self.minTol = 1e-8
        self.nModels = 0
        self.models = []
        # ***************************************** #
        # ********** PHYSICAL PARAMETERS ********** #
        # ***************************************** #
        self.physical ={'densityA': 998.2,
                        'viscosityA': 1.004e-6,
                        'densityB': 1.205,
                        'viscosityB': 1.500e-5,
                        'surf_tension_coeff': 72.8E-3,
                        'gravity': [0.0, -9.8, 0.0]}

        # ****************************************** #
        # ********** NUMERICAL PARAMETERS ********** #
        # ****************************************** #
        # shared variables between different models
        epsFact = 1.5
        sc_uref = 1.
        sc_beta = 1.5
        shockCapturingFactor = 0.5
        # RANS2P: full doc in proteus.mprans.RANS2P.py
        self.rans2p = {'index': None,
                       'tolFac': 0.001,
                       'name': 'rans2p',
                       'useMetrics': 1.0,
                       'epsFact_viscosity': epsFact,
                       'epsFact_density': epsFact,
                       'ns_forceStrongDirichlet': False,
                       'weak_bc_penalty_constant': 1.0E6,
                       'useRBLES': 0.0,
                       'useRANS': 0.0,
                       'ns_closure': 0,
                       'useVF': 1.0,
                       'ns_shockCapturingFactor': shockCapturingFactor,
                       'ns_lag_shockCapturing': True,
                       'ns_lag_subgridError': True,
                       'timeDiscretization': 'vbdf',
                       'timeOrder': 2}
        # RANS3P: full doc in proteus.mprans.RANS3P.py
        self.rans3p = {'index': None,
                       'name': 'rans3p',
                       'useMetrics': 1.0,
                       'epsFact_viscosity': epsFact,
                       'epsFact_density': epsFact,
                       'ns_forceStrongDirichlet': False,
                       'ns_sed_forceStrongDirichlet': False,
                       'weak_bc_penalty_constant': 1.0E6,
                       'useRBLES': 0.0,
                       'useRANS': 0.0,
                       'ns_closure': 0,
                       'useVF': 1.0,
                       'ns_shockCapturingFactor': shockCapturingFactor,
                       'ns_lag_shockCapturing': True,
                       'ns_lag_subgridError': True,
                       'timeDiscretization': 'vbdf',
                       'timeOrder': 2,
                       'PSTAB': 0,
                       'USE_SUPG': True,
                       'ARTIFICIAL_VISCOSITY': 2,
                       'cE': 1.0,
                       'cMax': 1.0}
        self.clsvof = {'index': None,
                       'name': 'clsvof',
                       'useMetrics': 1.0,
                       'epsFactHeaviside': epsFact,
                       'epsFactDirac': epsFact,
                       'epsFactRedist': 0.33,
                       'lambdaFact': 10.0,
                       'outputQuantDOFs': True,
                       'computeMetrics': 1,
                       'eps_tolerance_clsvof': False}
        self.pressure = {
            # index of model
            'index': None,
            # name of model
            'name': 'pressure',
        }
        self.pressureIncrement = {
            # index of model
            'index': None,
            # name of model
            'name': 'pressureincrement',
        }
        self.pressureInitial = {
            # index of model
            'index': None,
            # name of model
            'name': 'pressureInitial',
        }
        # VOF: full doc in proteus.mprans.VOF.py
        self.vof = {
            # index of model
            'index': None,
            # name of model
            'name': 'vof',
            # tolerance factor
            'tolFac': 0.001,
            'useMetrics': True,
            'checkMass': True,
            'sc_uref': sc_uref,
            'sc_beta': sc_beta,
            'epsFact': epsFact,
            'shockCapturingFactor': shockCapturingFactor,
            'lag_shockCapturing': True,
        }
        # NCLS: full doc in proteus.mprans.NCLS.py
        self.ncls = {
            # index of model
            'index': None,
            # name of model
            'name': 'ncls',
            # tolerance factor
            'tolFac': 0.001,
            'useMetrics': True,
            'checkMass': False,
            'sc_uref': sc_uref,
            'sc_beta': sc_beta,
            'epsFact': epsFact,
            'shockCapturingFactor': shockCapturingFactor,
            'lag_shockCapturing': True,
        }
        # RDLS: full doc in proteus.mprans.RDLS.py
        self.rdls = {
            # index of model
            'index': None,
            # name of model
            'name': 'rdls',
            # tolerance factor
            'tolFac': 0.01,
            'useMetrics': True,
            'applyRedistancing': True,
            'backgroundDiffusionFactor': 0.01,
            'epsFact': 0.33,
            'shockCapturingFactor': shockCapturingFactor,
            'lag_shockCapturing': False,
        }
        # AddedMass: full doc in proteus.mprans.AddedMass.py
        self.addedmass = {
            # index of model
            'index': None,
            # name of model
            'name': 'addedmass',
            # tolerance factor
            'tolFac': 0.001,
            # flags of rigid body (1: rigid body)
            'flags_rigidbody': None,
        }
        # MoveMesh: full doc in proteus.mprans.MoveMesh.py
        self.movemeshelastic = {
            # index of model
            'index': None,
            # name of model
            'name': 'movemeshelastic',
            # tolerance factor
            'tolFac': 0.001,
            # Young's modulus
            'E': 1.,
            # Poisson ratio
            'nu': 0.3,
        }
        # MoveMeshMonitor: full doc in proteus.mprans.MoveMeshMonitor.py
        self.movemeshmonitor = {
            # index of model
            'index': None,
            # name of model
            'name': 'movemeshmonitor',
            # tolerance factor
            'tolFac': 0.001,
            # tolerance factor
            'tolFac': 0.001,
            # monitor function f(x, t)
            'func': lambda x, t: 1000.,
            # min element size
            'he_min': 0.,
            # max element size
            'he_max': 1.,
            # epsFact around f=0 where he=he_min
            'epsFact': epsFact,
            # time step size from 0 to 1
            'epsFact': 0.1,
            # number of Laplace smoothing steps
            'nSmoothOut': 0.,
            # grading factor of mesh (increase in volume)
            'grading': 1.,
            # grading type (0: no grading, 1: hyperbolic, 2: log)
            'grading_type': 0,
            # reset node velocity at beginning of step
            'resetNodeVelocityArray': True,
            # use LS distance to smooth mesh on top of func
            'useLS': True,
            # list of nodeMaterialTypes to fix (1 to fix)
            'fixedNodeMaterialTypes': None,
            # list of elementMaterialTypes to fix (1 to fix)
            'fixedElementMaterialTypes': None,
            # list of nodeMaterialTypes where nodeVelocity should be 0.
            'noNodeVelocityNodeMaterialTypes': None,
            # scale distance function with nd
            'scale_with_nd': False,
            # do the first step
            'do_firstStep': False,
        }

    def initializeParameters(self):
        all_models = [self.rans2p,
                      self.rans3p,
                      self.clsvof,
                      self.vof,
                      self.ncls,
                      self.rdls,
                      self.movemeshelastic,
                      self.movemeshmonitor,
                      self.addedmass,
                      self.pressureInitial,
                      self.pressure,
                      self.pressureIncrement]
        self.nModels = 0
        self.models = []
        for i in range(len(all_models)):
            model = all_models[i]
            if model['index'] >= 0:
                self.nModels += 1
                self.models += [model]
                logEvent('TwoPhaseFlow parameters for model: {name}'.format(name=model['name']))
                for key, value in model.items():
                    logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('----------')
