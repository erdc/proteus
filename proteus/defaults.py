import sys, recordtype, os, imp, copy, inspect
from . import (TransportCoefficients,
               Transport,
               default_p,
               TimeIntegration,
               Quadrature,
               FemTools,
               SubgridError,
               ShockCapturing,
               NumericalFlux,
               NonlinearSolvers,
               LinearAlgebraTools,
               LinearSolvers,
               clapack,
               StepControl,
               AuxiliaryVariables,
               MeshTools,
               default_n,
               SplitOperator,
               default_so)

from .Archiver import ArchiveFlags
from .Profiling import logEvent

import sys

if sys.version_info.major < 3:  # Python 2?
    # Using exec avoids a SyntaxError in Python 3.
    exec("""def reraise(exc_type, exc_value, exc_traceback=None):
                raise exc_type, exc_value, exc_traceback""")
else:
    def reraise(exc_type, exc_value, exc_traceback=None):
        if exc_value is None:
            exc_value = exc_type()
        if exc_value.__traceback__ is not exc_traceback:
            raise exc_value.with_traceback(exc_traceback)
        raise exc_value

physics_default_keys = []
physics_excluded_keys = []

for k in (set(dir(default_p)) -
          (set(dir(FemTools))|
           set(dir(MeshTools))|
           set(dir(TransportCoefficients))|
           set(dir(Transport)))):
    if (k[:2] != '__' and k not in ['FemTools',
                                    'MeshTools',
                                    'TransportCoefficients',
                                    'Transport']):
        physics_default_keys.append(k)
    else:
        physics_excluded_keys.append(k)

_Physics_base = recordtype.recordtype('Physics_base',
                                     [(k,default_p.__dict__[k])
                                      for k in physics_default_keys],
                                     use_slots=False)
class Physics_base(_Physics_base):
    __frozen = False

    def __init__(self, **args):
        super(Physics_base,self).__init__(**args)
        for k in set(physics_default_keys) - set(args.keys()):
            v = default_p.__dict__[k]
            if not inspect.isclass(v):
                try:
                    self.__dict__[k] = copy.deepcopy(v)
                except:
                    pass

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, val):
        self.__setattr__(key, val)

    def __setattr__(self, key, val):
        if self.__frozen and not hasattr(self, key):
            raise TypeError("{key} is not an option".format(key=key))
        object.__setattr__(self, key, val)

    def _freeze(self):
        self.__frozen = True

    def _unfreeze(self):
        self.__frozen = False

    def addOption(self, name, value):
        if self.__frozen is True:
            frozen = True
            self.__frozen = False
        self.__setattr__(name, value)
        if frozen is True:
            self._freeze()
                
def reset_default_p():
    for k,v in Physics_base().__dict__.items():
        default_p.__dict__[k] = v

def load_physics(pModule, path='.'):
    reset_default_p()
    sys.path.append(path)
    p = imp.load_source(pModule, os.path.join(path, pModule+".py"))
    sys.path.remove(path)
    physics_object = Physics_base()
    for k,v in p.__dict__.items():
        if k not in physics_excluded_keys:
            physics_object.__dict__[k] = v
    return physics_object

numerics_default_keys = []
numerics_excluded_keys = []

for k in (set(dir(default_n)) -
          (set(dir(TimeIntegration))|
           set(dir(Quadrature))|
           set(dir(FemTools))|
           set(dir(SubgridError))|
           set(dir(ShockCapturing))|
           set(dir(NumericalFlux))|
           set(dir(NonlinearSolvers))|
           set(dir(LinearAlgebraTools))|
           set(dir(LinearSolvers))|
           set(dir(clapack))|
           set(dir(StepControl))|
           set(dir(AuxiliaryVariables))|
           set(dir(MeshTools)))):
    if (k[:2] != '__' and k not in ['TimeIntegration',
                                    'Quadrature',
                                    'FemTools',
                                    'SubgridError',
                                    'ShockCapturing',
                                    'NumericalFlux',
                                    'NonlinearSolvers',
                                    'LinearAlgebraTools',
                                    'LinearSolvers',
                                    'clapack',
                                    'StepControl',
                                    'AuxiliaryVariables',
                                    'MeshTools']):
        numerics_default_keys.append(k)
    else:
        numerics_excluded_keys.append(k)

_Numerics_base = recordtype.recordtype('Numerics_base',
                                      [(k,default_n.__dict__[k])
                                       for k in numerics_default_keys],
                                      use_slots=False)
class Numerics_base(_Numerics_base):
    __frozen = False

    def __init__(self, **args):
        super(Numerics_base,self).__init__(**args)
        for k in set(numerics_default_keys) - set(args.keys()):
            v = default_n.__dict__[k]
            if not inspect.isclass(v):
                try:
                    self.__dict__[k] = copy.deepcopy(default_n.__dict__[k])
                except:
                    pass

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, val):
        self.__setattr__(key, val)

    def __setattr__(self, key, val):
        if self.__frozen and not hasattr(self, key):
            raise TypeError("{key} is not an option".format(key=key))
        object.__setattr__(self, key, val)

    def _freeze(self):
        self.__frozen = True

    def _unfreeze(self):
        self.__frozen = False

    def addOption(self, name, value):
        if self.__frozen is True:
            frozen = True
            self.__frozen = False
        self.__setattr__(name, value)
        if frozen is True:
            self._freeze()

def reset_default_n():
    for k,v in Numerics_base().__dict__.items():
        default_n.__dict__[k] = v

def load_numerics(nModule, path='.'):
    reset_default_n()
    sys.path.append(path)
    n = imp.load_source(nModule, os.path.join(path, nModule+".py"))
    sys.path.remove(path)
    numerics_object = Numerics_base()
    for k,v in n.__dict__.items():
        if k not in numerics_excluded_keys:
            numerics_object.__dict__[k] = v
    return numerics_object

system_default_keys = []
system_excluded_keys = []

for k in (set(dir(default_so)) -
          (set(dir(SplitOperator))|
           set(dir(ArchiveFlags)))):
    if (k[:2] != '__' and k not in ['SplitOperator',
                                    'ArchiveFlags']):
        system_default_keys.append(k)
    else:
        system_excluded_keys.append(k)

_System_base = recordtype.recordtype('System_base',
                                    [(k,default_so.__dict__[k])
                                     for k in system_default_keys],
                                    use_slots=False)

class System_base(_System_base):
    def __init__(self, **args):
        super(System_base,self).__init__(**args)
        for k in set(system_default_keys) - set(args.keys()):
            v = default_so.__dict__[k]
            if not inspect.isclass(v):
                try:
                    self.__dict__[k] = copy.deepcopy(default_so.__dict__[k])
                except:
                    pass
def reset_default_so():
    for k,v in System_base().__dict__.items():
        default_so.__dict__[k] = v

def load_system(soModule, path='.'):
    reset_default_so()
    sys.path.append(path)
    so = imp.load_source(soModule, os.path.join(path, soModule+".py"))
    sys.path.remove(path)
    system_object = System_base()
    for k,v in so.__dict__.items():
        if k not in system_excluded_keys:
            system_object.__dict__[k] = v
    return system_object

