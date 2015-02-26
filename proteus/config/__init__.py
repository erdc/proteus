import platform

if platform.node().startswith('garnet') or platform.node().startswith('copper'):
    from garnet import *
else:
    from default import *