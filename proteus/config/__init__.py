import os

if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('garnet')
                                 or
                                 os.environ['HOSTNAME'].startswith('copper')
                                 or
                                 os.environ['HOSTNAME'].startswith('lightning')):
    from garnet import *
else:
    from default import *
