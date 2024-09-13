#!/usr/bin/env python
"""
Fine-scale viscous Burgers equation solver.

The equation is

.. math: \frac{\partial u}{\partial t} + \nabla \cdot (\frac{1}{2} u^2 - \epsilon \nabla u) = 0
"""
from burgers_init import *

physics.name = "pod_burgers_{0}d".format(physics.nd)
if use_deim:
    physics.name = "pod_burgers_DEIM_{0}d".format(physics.nd)
so.name = physics.name
Profiling.logLevel=5
Profiling.verbose=True
if use_deim:
    numerics.multilevelNonlinearSolver = NonlinearSolvers.POD_DEIM_Newton
    numerics.levelNonlinearSolver = NonlinearSolvers.POD_DEIM_Newton
else:
    numerics.multilevelNonlinearSolver = NonlinearSolvers.POD_Newton
    numerics.levelNonlinearSolver = NonlinearSolvers.POD_Newton
numerics.tolFac = 0.0#relative tolerance
numerics.nl_atol_res = 1.0e-4 #nonlinear solver rtolerance
ns = NumericalSolution.NS_base(so,[physics],[numerics],so.sList,opts)

import time
start = time.time()
failed=ns.calculateSolution('run1')
assert(not failed)
end = time.time() # we measure time required to obtain the fully resolved solution
print('Time required was %s seconds' % (end - start))

if physics.nd == 3:
    burgers_plot3d()
