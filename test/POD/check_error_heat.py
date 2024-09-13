#!/usr/bin/env python
from read_hdf5 import *

T = 1.0
nDTout = 100
DT = T/float(nDTout)

def uex0(x,t):
    """
    Exact solution
    """
    return 128.0*(1.0-t)*x[...,1]*(1.0-x[...,1])*x[...,2]*(1.0-x[...,2])*x[...,0]*(1.0-x[...,0])

archive = Archiver.XdmfArchive(".","heat_3d",readOnly=True)

label="/%s%d" % ('nodesSpatial_Domain',0)
print('trying to read from %s ' % label)
coord = read_from_hdf5(archive.hdfFile,label)

import numpy as np
u = read_from_hdf5(archive.hdfFile,'/u0')
uex_vals = np.zeros(u.shape,'d')

for i in range(0,nDTout+1):
    time_level_to_read=i
    label="/%s%d" % ('u',time_level_to_read)
    print('trying to read from %s ' % label)
    u = read_from_hdf5(archive.hdfFile,label)
    uex_vals = uex0(coord,i*DT)
    err = u-uex_vals
    err *= err
    err *= 1.0/9261.0 #9261 = 21^3
    L2approx = np.sqrt(err.sum())
    print("Trapezoidal approximation for error at dofs for nx=21 ny=21 nz=21 is %s " % L2approx)

