#!/usr/bin/env python

from read_hdf5 import *
import burgers_init,deim_utils


T = 1.0
nDTout = 100
DT = T/float(nDTout)

archive = Archiver.XdmfArchive(".",burgers_init.physics.name,readOnly=True)

import numpy as np

#generate snapshots for solution
S = deim_utils.read_snapshots(archive,nDTout+1,'u',)
U, s, V = np.linalg.svd(S, full_matrices=False)
print 'SVD for solution done!'
np.savetxt('Snapshots', S, delimiter=' ')
np.savetxt('SVD_basis', U, delimiter=' ')
np.savetxt('Singular_values', s, delimiter=' ')

if burgers_init.use_hyper:
    Sf = deim_utils.read_snapshots(archive,nDTout+1,'spatial_residual0')
    Uf,sf,Vf = np.linalg.svd(Sf,full_matrices=False)
    print 'SVD for spatial residual done!'
    np.savetxt('Fs_Snapshots', Sf, delimiter=' ')
    np.savetxt('Fs_SVD_basis', Uf, delimiter=' ')
    np.savetxt('Fs_Singular_values', sf, delimiter=' ')
    
