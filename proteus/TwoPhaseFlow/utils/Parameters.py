from __future__ import division
from past.utils import old_div
from builtins import object
from proteus.Profiling import logEvent

# default values for several models
epsFact = 1.5
sc_uref = 1.
sc_beta = 1.5
shockCapturingFactor = 0.5

class ParametersHolder:
    """
    """
    def __init__(self):
        self.minTol = 1e-8
        self.nModels = 0
        self.models_list = []
        self.Models = ParametersModelsHolder()
        self.physical = ParametersPhysical()

    def initializeParameters(self):
        all_models = [self.Models.rans2p,
                      self.Models.rans3p,
                      self.Models.clsvof,
                      self.Models.vof,
                      self.Models.ncls,
                      self.Models.rdls,
                      self.Models.movemeshelastic,
                      self.Models.movemeshmonitor,
                      self.Models.addedmass,
                      self.Models.pressureInitial,
                      self.Models.pressure,
                      self.Models.pressureIncrement]
        self.nModels = 0
        self.models_list = []
        for i in range(len(all_models)):
            model = all_models[i]
            if model['index'] >= 0:
                self.nModels += 1
                self.models_list += [model]
                logEvent('TwoPhaseFlow parameters for model: {name}'.format(name=model['name']))
                for key, value in model.__dict__.items():
                    logEvent('{key}: {value}'. format(key=key, value=value))
                logEvent('----------')

class ParametersModelBase(object):
    """
    """
    def __init__(self,
                 name=None,
                 index=None):
        self.name = name
        self.index = index
        self.tolFac = 0.001
        self.minTol = 1e-8

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, val):
        self.__dict__[key] = val

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
        self.useRBLES = 0.
        self.useRANS = 0
        self.ns_closure = 0.
        self.useVF = 1
        self.ns_shockCapturingFactor = shockCapturingFactor
        self.ns_lag_shockCapturing = True
        self.ns_lag_subgridError = True
        self.timeDiscretization = None
        self.timeOrder = 2

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
        self.useRBLES = 0.
        self.useRANS = 0
        self.ns_closure = 0.
        self.useVF = 1
        self.ns_shockCapturingFactor = shockCapturingFactor
        self.ns_lag_shockCapturing = True
        self.ns_lag_subgridError = True
        self.timeDiscretization = 'vdbf'
        self.timeOrder = 2
        self.PSTAB = 0
        self.USE_SUPG = True
        self.ARTIFICIAL_VISCOSITY = 2
        self.cE = 1.
        self.cMax = 1.

class ParametersModelPressure(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressure, self).__init__(name='pressure', index=None)

class ParametersModelPressureInitial(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressureInitial, self).__init__(name='pressureInitial', index=None)

class ParametersModelPressureIncrement(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelPressureIncrement, self).__init__(name='pressureincrement', index=None)

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

class ParametersModelAddedMass(ParametersModelBase):
    """
    """
    def __init__(self):
        super(ParametersModelAddedMass, self).__init__(name='addedmass', index=None)
        self.flags_rigidbody = None

class ParametersModelMoveMeshMonitor(ParametersModelBase):
    def __init__(self):
        super(ParametersModelMoveMeshMonitor, self).__init__(name='movemeshmonitor', index=None)
        self.fund = lambda x, t: 1000.
        self.he_min = 0.
        self.he_max = 1000.
        self.epsFact = epsFact
        self.epsTimeStep = 0.1
        self.nSmoothOut = 0.
        self.grading = 1.1
        self.grading_type = 2
        self.resetNodeVelocityArray = None
        self.useLS = True
        self.fixedNodeMaterialTypes = None
        self.fixedElementMaterialTypes = None
        self.noNodeVelocityNodeMaterialTypes = None
        self.scale_with_nd = None
        self.do_firstStep = False

class ParametersModelMoveMeshElastic(ParametersModelBase):
    def __init__(self):
        super(ParametersModelMoveMeshElastic, self).__init__(name='movemeshelastic', index=None)
        self.E = 1.
        self.nu = 0.3

class ParametersPhysical:
    def __init__(self):
        self.densityA = 998.2
        self.densityB = 1.205
        self.viscosityA = 1.004e-6
        self.viscosityB = 1.500e-5
        self.surf_tension_coeff = 72.8e-3
        self.gravity = [0., -9.81, 0.]

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, val):
        self.__dict__[key] = val

class ParametersModelsHolder:
    def __init__(self):
        self.rans2p = ParametersModelRANS2P()
        self.vof = ParametersModelVOF()
        self.ncls = ParametersModelNCLS()
        self.rdls = ParametersModelRDLS()
        self.addedmass = ParametersModelAddedMass()
        self.movemeshmonitor = ParametersModelMoveMeshMonitor()
        self.movemeshelastic = ParametersModelMoveMeshElastic()
        self.clsvof = ParametersModelCLSVOF()
        self.rans3p = ParametersModelRANS3P()
        self.pressureInitial = ParametersModelPressureInitial()
        self.pressure = ParametersModelPressure()
        self.pressureIncrement = ParametersModelPressureIncrement()
