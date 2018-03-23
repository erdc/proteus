import sys, recordtype, os, imp
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
               lapackWrappers,
               StepControl,
               AuxiliaryVariables,
               MeshTools,
               default_n,
               SplitOperator,
               default_so)

from .Archiver import ArchiveFlags

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

Physics_base = recordtype.recordtype('Physics_base',
                                     [(k,default_p.__dict__[k])
                                      for k in physics_default_keys],
                                     use_slots=False)

def reset_default_p():
    for k,v in Physics_base().__dict__.iteritems():
        default_p.__dict__[k] = v

def load_physics(pModule, path='.'):
    reset_default_p()
    sys.path.append(path)
    p = imp.load_source(pModule, os.path.join(path, pModule+".py"))
    sys.path.remove(path)
    physics_object = Physics_base()
    for k,v in p.__dict__.iteritems():
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
           set(dir(lapackWrappers))|
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
                                    'lapackWrappers',
                                    'StepControl',
                                    'AuxiliaryVariables',
                                    'MeshTools']):
        numerics_default_keys.append(k)
    else:
        numerics_excluded_keys.append(k)

Numerics_base = recordtype.recordtype('Numerics_base',
                                      [(k,default_n.__dict__[k])
                                       for k in numerics_default_keys],
                                      use_slots=False)

def reset_default_n():
    for k,v in Numerics_base().__dict__.iteritems():
        default_n.__dict__[k] = v

def load_numerics(nModule, path='.'):
    reset_default_n()
    sys.path.append(path)
    n = imp.load_source(nModule, os.path.join(path, nModule+".py"))
    sys.path.remove(path)
    numerics_object = Numerics_base()
    for k,v in n.__dict__.iteritems():
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

System_base = recordtype.recordtype('System_base',
                                    [(k,default_so.__dict__[k])
                                     for k in system_default_keys],
                                    use_slots=False)

def reset_default_so():
    for k,v in System_base().__dict__.iteritems():
        default_so.__dict__[k] = v

def load_system(soModule, path='.'):
    reset_default_so()
    sys.path.append(path)
    so = imp.load_source(soModule, os.path.join(path, soModule+".py"))
    sys.path.remove(path)
    modulefile.close()
    system_object = System_base()
    for k,v in so.__dict__.iteritems():
        if k not in system_excluded_keys:
            system_object.__dict__[k] = v
    return system_object

