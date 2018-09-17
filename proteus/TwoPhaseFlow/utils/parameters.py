from proteus.Profiling import logEvent

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
physical ={'densityA': 998.2,
           'viscosityA': 1.004e-6,
           'densityB': 1.205,
           'viscosityB': 1.500e-5,
           'surf_tension_coeff': 72.8E-3,
           'gravity': [0.0, -9.8, 0.0]}

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
rans2p = {'useMetrics': 1.0,
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
rans3p = {'useMetrics': 1.0,
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
clsvof = {'useMetrics': 1.0,
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
        # VOF: full doc in proteus.mprans.VOF.py
        self.vof = {
            # index of model
            'index': None,
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
                      self.addedmass]
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

