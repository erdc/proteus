import os

if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('garnet')):
    from garnet import *
elif 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('spirit')):
    from spirit import *
elif 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('copper')):
    from copper import *
elif 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('lightning')):
    from lightning import *
elif 'HOSTNAME' in os.environ and os.environ['HOSTNAME'].startswith('viutill'):
    from viutill import *
else:
    from default import *
