#!/usr/bin/env python
"""
Fine-scale heat equation solver
The equation is du/du - Laplace u + u + f(x,y,z,t) = 0
"""
from heat_init import *

physics.name = "heat_3d"
so.name = physics.name
ns = NumericalSolution.NS_base(so,[physics],[numerics],so.sList,opts)

import time
start = time.time()
failed=ns.calculateSolution('run1')
assert(not failed)
end = time.time() # we measure time required to obtain the fully resolved solution
print('Time required was %s seconds' % (end - start))

#arrays for using matplotlib's unstructured plotting interface
x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,0]
y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,1]
z = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,2]
u = ns.modelList[0].levelModelList[-1].u[0].dof

#we want to build a 3d plot of f(x,y,z0), when z0 = 0.5
u_range = u[4410:4851]
x_range = x[4410:4851]
y_range = y[4410:4851]
u_range = u_range.reshape(21,21)
x_range = x_range.reshape(21,21)
y_range = y_range.reshape(21,21)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
surf=ax.plot_surface(x_range, y_range, u_range, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('x'); plt.ylabel('y')
plt.title('approximate solution at t = 1');

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.05f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig("solution.png")
plt.show()

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
