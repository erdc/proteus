from __future__ import division
from past.utils import old_div
from builtins import object
from proteus.Profiling import logEvent
from proteus.MeshTools import MeshOptions

# default values for several models
epsFact = 1.5
sc_uref = 1.
sc_beta = 1.5
shockCapturingFactor = 0.5


class ParametersHolder:
    """
    """
    def __init__(self, ProblemInstance=None):
        # gain access to problem class if necessary
        self.__Problem = ProblemInstance
        # default options
        self.minTol = 1e-8
        self.nModels = 0
        self.models_list = []
        self.Models = ParametersModelsHolder()
        self.physical = ParametersPhysical()
        self.mesh = MeshOptions(nd=self.__Problem.domain.nd)

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
        for i in range(len(all_models)):
            if i == 0:
                logEvent('----------')
            model = all_models[i]
            if model['index'] >= 0:
                self.nModels += 1
                self.models_list += [model]
                logEvent('TwoPhaseFlow parameters for model: {name}'.format(name=model['name']))
                for key, value in model.__dict__.items():
                    if key[0] != '_':  # do not print hidden attributes
                        logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('----------')


class ParametersModelsHolder:
    """
    """
    def __init__(self):
        self.rans2p = ParametersModelRANS2P()
        self.vof = ParametersModelVOF()
        self.ncls = ParametersModelNCLS()
        self.rdls = ParametersModelRDLS()
        self.addedMass = ParametersModelAddedMass()
        self.moveMeshMonitor = ParametersModelMoveMeshMonitor()
        self.moveMeshElastic = ParametersModelMoveMeshElastic()
        self.clsvof = ParametersModelCLSVOF()
        self.rans3p = ParametersModelRANS3P()
        self.pressureInitial = ParametersModelPressureInitial()
        self.pressure = ParametersModelPressure()
        self.pressureIncrement = ParametersModelPressureIncrement()
        self.mcorr = ParametersModelMCorr()


class ParametersBase(object):
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


class ParametersModelBase(ParametersBase):
    """
    """
    def __init__(self,
                 name=None,
                 index=None):
        super(ParametersModelBase, self).__init__()
        self.name = name
        self.index = index
        self.tolFac = 0.001
        self.minTol = 1e-8
        self.auxiliaryVariables = []


class ParametersModelRANS2P(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelRANS2P, self).__init__(name='rans2p', index=None)
        self.useMetrics = 1.
        self.epsFact_viscosity = epsFact
        self.epsFact_density = epsFact
        self.ns_forceStrongDirichlet = False
        self.weak_bc_penalty_constant = 1e6
        self.useRBLES = 0
        self.useRANS = 0
        self.ns_closure = 0
        self.useVF = 1
        self.ns_shockCapturingFactor = shockCapturingFactor
        self.ns_lag_shockCapturing = True
        self.ns_lag_subgridError = True
        self.timeDiscretization = None
        self.timeOrder = 2
        self.stokes = False
        self.eb_adjoint_sigma = 1.
        # freeze attributes
        self._freeze()


class ParametersModelRANS3P(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelRANS3P, self).__init__(name='rans3p', index=None)
        self.useMetrics = 1.
        self.epsFact_viscosity = epsFact
        self.epsFact_density = epsFact
        self.ns_forceStrongDirichlet = False
        self.ns_sed_forceStrongDirichlet = False
        self.weak_bc_penalty_constant = 1e6
        self.useRBLES = 0
        self.useRANS = 0
        self.ns_closure = 0
        self.useVF = 1
        self.ns_shockCapturingFactor = shockCapturingFactor
        self.ns_lag_shockCapturing = True
        self.ns_lag_subgridError = True
        self.timeDiscretization = 'vbdf'
        self.timeOrder = 2
        self.PSTAB = 0
        self.USE_SUPG = True
        self.ARTIFICIAL_VISCOSITY = 2
        self.cE = 1.
        self.cMax = 1.
        # freeze attributes
        self._freeze()


class ParametersModelPressure(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressure, self).__init__(name='pressure', index=None)
        # freeze attributes
        self._freeze()


class ParametersModelPressureInitial(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressureInitial, self).__init__(name='pressureInitial', index=None)
        # freeze attributes
        self._freeze()


class ParametersModelPressureIncrement(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressureIncrement, self).__init__(name='pressureIncrement', index=None)
        # freeze attributes
        self._freeze()


class ParametersModelCLSVOF(ParametersModelBase):
    def __init__(self):
        super(ParametersModelCLSVOF, self).__init__(name='clsvof', index=None)
        self.useMetrics = 1.
        self.epsFactHeaviside = epsFact
        self.epsFactDirac = epsFact
        self.epsFactRedist = 0.33
        self.lambdaFact = 10.
        self.outputQuantDOFs = True
        self.computeMetrics = 1
        self.eps_tolerance_clsvof = False
        # freeze attributes
        self._freeze()


class ParametersModelVOF(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelVOF, self).__init__(name='vof', index=None)
        self.useMetrics = True
        self.checkMass = True
        self.sc_uref = sc_uref
        self.sc_beta = sc_beta
        self.epsFact = epsFact
        self.shockCapturingFactor = shockCapturingFactor
        self.lag_shockCapturing = True
        # freeze attributes
        self._freeze()


class ParametersModelNCLS(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelNCLS, self).__init__(name='ncls', index=None)
        self.useMetrics = True
        self.checkMass = False
        self.sc_uref = sc_uref
        self.sc_beta = sc_beta
        self.epsFact = epsFact
        self.shockCapturingFactor = shockCapturingFactor
        self.lag_shockCapturing = True
        # freeze attributes
        self._freeze()


class ParametersModelRDLS(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelRDLS, self).__init__(name='rdls', index=None)
        self.tolFac = 0.01
        self.useMetrics = True
        self.applyRedistancing = True
        self.backgroundDiffusionFactor = 0.01
        self.epsFact = 0.33
        self.shockCapturingFactor = shockCapturingFactor
        self.lag_shockCapturing = False
        # freeze attributes
        self._freeze()


class ParametersModelMCorr(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelMCorr, self).__init__(name='mcorr', index=None)
        self.tolFac = 0.0001
        self.useMetrics = True
        self.checkMass = False
        self.applyCorrection = True
        self.epsFactHeaviside = epsFact
        self.epsFactDirac = epsFact
        self.epsFactDiffusion = 10.
        # freeze attributes
        self._freeze()


class ParametersModelAddedMass(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelAddedMass, self).__init__(name='addedMass', index=None)
        self.flags_rigidbody = None
        # freeze attributes
        self._freeze()


class ParametersModelMoveMeshMonitor(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelMoveMeshMonitor, self).__init__(name='moveMeshMonitor', index=None)
        self.func = lambda x, t: 1000.
        self.he_min = 0.
        self.he_max = 1000.
        self.epsFact = epsFact
        self.epsTimeStep = 0.1
        self.nSmoothOut = 0.
        self.nSmoothIn = 0.
        self.grading = 1.1
        self.grading_type = 2
        self.resetNodeVelocityArray = None
        self.useLS = True
        self.fixedNodeMaterialTypes = None
        self.fixedElementMaterialTypes = None
        self.noNodeVelocityNodeMaterialTypes = None
        self.scale_with_nd = False
        self.do_firstStep = False
        self.ntimes_solved = 1
        # freeze attributes
        self._freeze()


class ParametersModelMoveMeshElastic(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelMoveMeshElastic, self).__init__(name='moveMeshElastic', index=None)
        self.E = 1.
        self.nu = 0.3
        # freeze attributes
        self._freeze()


class ParametersPhysical(ParametersBase):
    """
    """
    def __init__(self):
        super(ParametersPhysical, self).__init__()
        self.densityA = 998.2
        self.densityB = 1.205
        self.viscosityA = 1.004e-6
        self.viscosityB = 1.500e-5
        self.surf_tension_coeff = 72.8e-3
        self.gravity = [0., -9.81, 0.]
        # freeze attributes
        self._freeze()

