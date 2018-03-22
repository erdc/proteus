import recordtype
from . import MeshTools, FemTools, TransportCoefficients, Transport, default_p

physics_default_keys = []
excluded_keys = []

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
        excluded_keys.append(k)

Physics_base = recordtype.recordtype('Physics',
                                     [(k,default_p.__dict__[k])
                                      for k in physics_default_keys],
                                     use_slots=False)

def reset_default_p():
    for k,v in Physics_base().__dict__.iteritems():
        default_p.__dict__[k] = v

def load_from_module(pModule):
    reset_default_p()
    p = __import__(pModule)
    physics_object = Physics_base()
    for k,v in p.__dict__.iteritems():
        if k not in excluded_keys:
            physics_object.__dict__[k] = v
    return physics_object
