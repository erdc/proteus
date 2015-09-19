#!/usr/bin/env python

from read_hdf5 import *
import numpy as np

"""
assumes nsnap values of array in val_name are stored in h5file as
/val_name'i' for i=0,nspap-1

loads these into a matrix and returns
"""
val_name = 'u'
nspap = 101
archive = Archiver.XdmfArchive(".","pod_burgers_1d",readOnly=True)
label_base="/%s%d"
u = read_from_hdf5(archive.hdfFile,label_base % (val_name,0))
S = np.reshape(u,(u.shape[0],1))
for i in range(1,nspap):
    label=label_base % (val_name,i)
    u = read_from_hdf5(archive.hdfFile,label)
    u = np.reshape(u,(u.shape[0],1))
    S = np.append(S,u,axis=1)
    #
np.savetxt('pod_Solutions', S, delimiter=' ')
