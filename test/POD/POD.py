#!/usr/bin/env python
"""
Proper orthogonal decomposition for the heat equation solver
"""
from heat_init import *
from read_hdf5 import *

#we construct this NS object again ONLY to have access to the load vector at each time step
physics.name = "heat_3d_reduction"
so.name = physics.name
ns = NumericalSolution.NS_base(so,[physics],[numerics],so.sList,opts)
save_projected_soln = True #do we save the projection of the reduced order solution on the fine grid in xmf?

#doing reduction in the remaining part of the code

import petsc4py
from petsc4py import PETSc
from petsc4py.PETSc import Mat
import numpy as np

DB = 4 #the number of basis vectors! Increase or decrease as you wish...

U = np.loadtxt('SVD_basis')
U = U[:,0:DB]
U_transpose = U.conj().T
#load mass matrix and reduce it
 
rp = np.loadtxt('iam',int)
rp = rp.astype(petsc4py.PETSc.IntType)
ci = np.loadtxt('jam',int)
ci = ci.astype(petsc4py.PETSc.IntType)
nz = np.loadtxt('am')

nr = rp.shape[0] - 1
nc = nr

M = Mat().createAIJ(size=(nr,nc),csr=(rp,ci,nz))

M_intermediate = np.zeros((nr,DB),'d')

for i in range(0,nr):
    res = M.getRow(i)
    M_intermediate[i,:] = np.dot(res[1],U[res[0],:])

M_pod = np.dot(U_transpose, M_intermediate)
np.savetxt('M_pod', M_pod, delimiter=' ')

#load stiffness matrix and reduce it

rp = np.loadtxt('ias',int)
rp = rp.astype(petsc4py.PETSc.IntType)
ci = np.loadtxt('jas',int)
ci = ci.astype(petsc4py.PETSc.IntType)
nz = np.loadtxt('as')

S = Mat().createAIJ(size=(nr,nc),csr=(rp,ci,nz))

S_intermediate = np.zeros((nr,DB),'d')

for i in range(0,nr):
    res = S.getRow(i)
    S_intermediate[i,:] = np.dot(res[1],U[res[0],:])

S_pod = np.dot(U_transpose, S_intermediate)
np.savetxt('S_pod', M_pod, delimiter=' ')

del M_intermediate
del S_intermediate

K = M_pod + DT*S_pod #generalized mass matrix
K = np.linalg.inv(K) #its inverse

archive = Archiver.XdmfArchive(".","heat_3d",readOnly=True)

#reading initial condition
u = read_from_hdf5(archive.hdfFile,'/u0')
u_pod = np.dot(U_transpose,u)


f = np.zeros((nr,),'d') #rhs vector
Model = ns.modelList[0].levelModelList[-1]
ns.tCount=0
if save_projected_soln:
    Model.u[0].dof[:] = u
    for index,model in enumerate(ns.modelList):
        ns.archiveInitialSolution(model,index)
import time
start = time.time()

#Backward Euler time loop
for i in range(1,nDTout+1):
    t = i*DT
    Model.timeIntegration.t = t
    Model.getLoadVector(f)
    f_pod = np.dot(U_transpose,f)
    rhs = np.dot(M_pod,u_pod)-DT*f_pod
    u_pod = np.dot(K,rhs)

    u_approx = np.dot(U,u_pod)

    ns.tCount += 1
    ns.tn_last = t
    if save_projected_soln:
        Model.u[0].dof[:] = u_approx
        for index,model in enumerate(ns.modelList):
            ns.archiveSolution(model,index,t)

    if (i == 50):
	U_approx = u_approx
    label="/%s%d" % ('u',i)
    print('trying to read from %s ' % label)
    u = read_from_hdf5(archive.hdfFile,label)

    err = u-u_approx
    err *= err
    err *= 1.0/9261.0
    L2approx = np.sqrt(err.sum())
    print('error is %s ' % L2approx)

end = time.time() # we measure time required to obtain the reduced solution
elapsed = end - start
print('time required was %s seconds' % (end-start))

#arrays for using matplotlib's unstructured plotting interface
u_range = U_approx[4410:4851]
u_range = u_range.reshape(21,21)
#x and y are needed for plotting only, z = 0.5
label="/%s%d" % ('nodesSpatial_Domain',0)
print('trying to read from %s ' % label)
coord = read_from_hdf5(archive.hdfFile,label)
x = coord[4410:4851,0]
y = coord[4410:4851,1]
x = x.reshape(21,21)
y = y.reshape(21,21)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
surf=ax.plot_surface(x, y, u_range, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('x'); plt.ylabel('y')
plt.title('reduced solution at t = 0.5');

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.05f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig("solution_reduced_t0_5.png")
plt.show()

#also look at final time step
u_range = u_approx[4410:4851]
u_range = u_range.reshape(21,21)

fig.clf()
ax = fig.gca(projection='3d')
surf=ax.plot_surface(x, y, u_range, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('x'); plt.ylabel('y')
plt.title('reduced solution at t = %s' % (nDTout*DT));

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.05f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig("solution_reduced_final.png")
plt.show()

#cleanup numerical solution files since calculateSolution wasn't called
for index,model in enumerate(ns.modelList):
    ns.closeArchive(model,index)

