import os

if 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('garnet'):
    from garnet import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('topaz'):
    from topaz import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('onyx'):
    from onyx import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('spirit'):
    from spirit import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('stampede'):
    from stampede import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('copper'):
    from copper import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('lightning'):
    from lightning import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('excalibur'):
    from excalibur import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('viutill'):
    from viutill import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('tamucluster'):
    from tamucluster import *
elif 'PROTEUS_ARCH' in os.environ and os.environ['PROTEUS_ARCH'].startswith('centos'):
    from centos import *
else:
    from default import *
