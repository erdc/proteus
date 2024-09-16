#!/usr/bin/env python
"""
Fine-scale viscous Burgers equation solver.

The equation is

.. math: \frac{\partial u}{\partial t} + \nabla \cdot (\frac{1}{2} u^2 - \epsilon \nabla u) = 0
"""
from burgers_init import *

physics.name = "burgers_{0}d".format(physics.nd)
so.name = physics.name
Profiling.logLevel=5
Profiling.verbose=True
    
ns = NumericalSolution.NS_base(so,[physics],[numerics],so.sList,opts,simFlagsList=simFlagsList)

import time
start = time.time()
failed=ns.calculateSolution('run1')
assert(not failed)
end = time.time() # we measure time required to obtain the fully resolved solution
print('Time required was %s seconds' % (end - start))

if physics.nd == 3:
    plot_burgers3d()
#saving mass and stiffness matrices below

Model = ns.modelList[0].levelModelList[-1]
mj = Model.initializeMassJacobian()
Model.getMassJacobian(mj)
kj = Model.initializeSpatialJacobian()
Model.getSpatialJacobian(kj)

rowptr,colind,nzval = mj.getCSRrepresentation()
np.savetxt('iam',rowptr,'%d', ' ')
np.savetxt('jam',colind,'%d', ' ')
np.savetxt('am',nzval,delimiter=' ')
rowptr_s,colind_s,nzval_s = kj.getCSRrepresentation()
np.savetxt('ias',rowptr_s,'%d',' ')
np.savetxt('jas',colind_s,'%d',' ')
np.savetxt('as',nzval_s,delimiter=' ')
