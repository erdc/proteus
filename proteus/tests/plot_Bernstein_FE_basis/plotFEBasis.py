import proteus
from proteus.FemTools import QuadraticOnSimplexWithNodalBasis
from proteus.FemTools import LagrangeOnCubeWithNodalBasis
from proteus.FemTools import LinearOnCubeWithNodalBasis
from proteus.FemTools import BernsteinOnCube

#psi = LagrangeOnCubeWithNodalBasis(nd=2)
psi = BernsteinOnCube(nd=2)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

npts = 40
x = np.arange(-1, 1., 1./(npts-1))
y = np.arange(-1, 1., 1./(npts-1))
xi,yi = np.meshgrid(x,y)
zi = xi.copy()
zi[:]=0.0
for i in range(9):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for j,xj in enumerate(x):
        for k,yk in enumerate(y):
            zi[j,k]=psi.basis[i]([xj,yk])
    ax.plot_surface(xi,yi,zi, rstride=1, cstride=1, cmap=cm.RdBu, linewidth=0, antialiased=False)
    plt.axis('tight')
    plt.savefig('basis_0'+str(i)+'.png')
    

#ax.plot_wireframe(xi,yi,0*zi, cmap=cm.coolwarm,zorder=10)
#ax.view_init(30, 0)

#plt.show()



