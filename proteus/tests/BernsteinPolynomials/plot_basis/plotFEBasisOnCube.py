import proteus
from proteus.FemTools import QuadraticOnSimplexWithNodalBasis
from proteus.FemTools import LagrangeOnCubeWithNodalBasis
from proteus.FemTools import LinearOnCubeWithNodalBasis
from proteus.FemTools import BernsteinOnCube

#psi = LinearOnCubeWithNodalBasis(nd=2)
#psi = LagrangeOnCubeWithNodalBasis(nd=2)
psi = BernsteinOnCube(nd=2,order=4)
nDOFs_per_element = len(psi.basis)

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

#2D plot
for i in range(nDOFs_per_element):
    fig = plt.figure()
    for j,xj in enumerate(x):
        for k,yk in enumerate(y):
            zi[j,k]=psi.basis[i]([xj,yk])
    plt.pcolor(xi, yi, zi, cmap='RdBu',vmin=zi.min(),vmax=1.0)
    plt.colorbar()
    plt.axis('tight')
    plt.savefig('2DView_basis_0'+str(i)+'.png')
    plt.close()
    
#3D plot
for i in range(nDOFs_per_element):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for j,xj in enumerate(x):
        for k,yk in enumerate(y):
            zi[j,k]=psi.basis[i]([xj,yk])
    ax.plot_surface(xi,yi,zi, rstride=1, cstride=1, cmap=cm.RdBu, linewidth=0, antialiased=False)
    plt.axis('tight')
    plt.savefig('3DView_basis_0'+str(i)+'.png')    
    plt.close()
