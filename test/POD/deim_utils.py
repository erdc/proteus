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
