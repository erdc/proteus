import os

if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('garnet')):
    from garnet import *
if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('copper')):
    from copper import *
if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('lightning')):
    from lightning import *
else:
    from default import *
