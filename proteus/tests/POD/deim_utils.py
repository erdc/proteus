#!/usr/bin/env python
"""
utility module for generating deim interpolants
"""
from read_hdf5 import *
import numpy as np

def read_snapshots(archive,nsnap,val_name):
    """
    assumes nsnap values of array in val_name are stored in h5file as
    /val_name'i' for i=0,nspap-1

    loads these into a matrix and returns
    """
    label_base="/%s%d"
    u = read_from_hdf5(archive.hdfFile,label_base % (val_name,0))
    S = np.reshape(u,(u.shape[0],1))
    for i in range(1,nsnap):
        label=label_base % (val_name,i)
        u = read_from_hdf5(archive.hdfFile,label)
        u = np.reshape(u,(u.shape[0],1))
        S = np.append(S,u,axis=1)
    #
    return S


def generate_svd_decomposition(archive,nsnap,val_name,outbase):
    """
    assumes nsnap values of array in val_name are stored in h5file as
    /val_name'i' for i=0,nspap-1

    loads these into a matrix, performs an SVD, and stores the output in outbase_SVD_basis, 
      outbase_singular_values in numpy's binary format

    returns U,s,V svd decomposition of snapshots
    """
    S = read_snapshots(archive,nsnap,val_name)

    U, s, V= np.linalg.svd(S,full_matrices=False)
    
    np.savetxt(outbase+'_SVD_basis',U,delimiter=' ')
    np.savetxt(outbase+'_SVD_singular_values',s,delimiter=' ')

    return U,s,V

def calculate_deim_indices(Uin):
    """
    input: Uin n x m array of basis vectors for nonlinear function snapshots
    output: rho, m vector of indices \rho_i for extracting $\vec F$ values

    """
    n,m=Uin.shape
    rind = np.argmax(np.absolute(Uin[:,0]))
    U=np.array(Uin[:,0])
    rho=np.array([rind],'i')
    #Pt = np.zeros((1,n),'d')
    #P[0,rind]=1.0
    for j in range(1,m):
        u = Uin[:,j]        
        Up=U[rho]#Up= np.dot(Pt,U)
        up=u[rho]#up= np.dot(Pt,u)
        if j==1:
            c=up/Up
            r=u-U*c
        else:
            c =np.linalg.solve(Up,up)
            r=u-np.dot(U,c)          
        rind=np.argmax(np.absolute(r))
        rho_new = np.zeros(j+1,'i'); 
        rho_new[:-1]=rho; rho_new[-1]=rind; rho = rho_new
        U_new=np.zeros((n,j+1),'d')
        U_new[:,:-1]=U.reshape(n,j); U_new[:,-1]=u
        U=U_new
    #
    return rho

def deim_alg(Uin,m):
    """
    Basic procedure

    - given $m$, dimension for $F$ reduced basis $\mathbf{U}_m$
    - call DEIM algorithm to determine $\vec \rho$. 
    - build $\mathbf{P}$ from $\rho$ as 
      $$
      \mathbf{P} = [\vec e_{\rho_1},\vec e_{\rho_2},\dots,\vec e_{\rho_m}]
      $$
    - invert $\mathbf{P}^T\mathbf{U}_m$
    - return \rho and $\mathbf{P}_F=\mathbf{U}_m(\mathbf{P}^T\mathbf{U}_m)^{-1}$

    
    """
    assert m <= Uin.shape[1]
    Um = Uin[:,0:m]
    rho = calculate_deim_indices(Um)
    PtUm = Um[rho]
    assert PtUm.shape == (m,m)
    PtUmInv = np.linalg.inv(PtUm)
    PF= np.dot(Um,PtUmInv)
    return rho,PF

def calculate_gpod_indices(Uin):
    """
    input: Uin n x m array of basis vectors for nonlinear function snapshots
    output: rho, vector of indices \rho_i for extracting $\vec F$ values, and nind, vector of index quantity rate per basis vector

    """
    n,m=Uin.shape
    indices = set()
    rind = np.argmax(np.absolute(Uin[:,0]))
    indices.add(rind)
    neigs = set([max(0,rind-1)]+[min(n-1,rind+1)])
    indices |= neigs
    rho = np.array(list(indices),dtype='i')
    l = len(indices)
    nind = np.array([l])
    for j in range(1,m):
	U = Uin[:,0:j]
	u = Uin[:,j]
	c,l2res,rank,svals = np.linalg.lstsq(U[rho],u[rho])
	r = u-np.dot(U,c)
	r[rho] = 0.0
	rind = np.argmax(np.absolute(r))
	if not rind in indices:
	    indices.add(rind)
	else:
	    break
	neigs = set([max(0,rind-1)]+[min(n-1,rind+1)])
	indices |= neigs
	rho = np.array(list(indices),dtype='i')
	nind = np.append(nind, len(indices)-l)
	l = len(indices)
    #
    return rho, nind

def gpod_alg(Uin,m):
    """
    Basic procedure

    - given $m$, dimension for $F$ reduced basis $\mathbf{U}_m$
    - call Gappy POD algorithm to determine $\vec \rho$.
    - build $\mathbf{P}$ from $\rho$ as
      $$
      \mathbf{P} = [\vec e_{\rho_1},\vec e_{\rho_2},\dots,\vec e_{\rho_m}]
      $$
    - invert $\mathbf{P}^T\mathbf{U}_m$ in a presudo sense
    - return \rho and $\mathbf{P}_F=\mathbf{U}_m(\mathbf{P}^T\mathbf{U}_m)^{-1}$


    """
    assert m <= Uin.shape[1]
    Um = Uin[:,0:m]
    rho, nind = calculate_gpod_indices(Um)
    PtUm = Um[rho]
    assert PtUm.shape == (np.sum(nind,dtype='i'),m)
    PtUmInv = np.linalg.pinv(PtUm)
    PF= np.dot(Um,PtUmInv)
    return rho,nind,PF

def visualize_zslice(variable,nnx,nny,iz,x=None,y=None,name=None):
    """
    convenience function for plotting a slice
    """
    istart = nnx*nny*iz
    iend   = nnx*nny*(iz+1)
    v_slice= variable[istart:iend]
    v_slice= v_slice.reshape(nnx,nny)
    if x == None:
        x = np.outer(np.arange(nnx),np.arange(nnx))
    if y == None:
        y = np.outer(np.arange(nny),np.arange(nny))
    assert x.shape == v_slice.shape
    assert y.shape == v_slice.shape

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf=ax.plot_surface(x,y,v_slice,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
    plt.xlabel('x'); plt.ylabel('y')
    if name == None:
        name = 'deim_slice_z={0}.png'.format(iz)
    plt.savefig(name)


    return surf
