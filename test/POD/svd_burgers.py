#!/usr/bin/env python
from burgers_init import use_deim
from proteus import deim_utils,Archiver
from proteus.deim_utils import read_snapshots

T = 1.0
nDTout = 100
DT = T/float(nDTout)

archive = Archiver.XdmfArchive(".",burgers_init.physics.name,readOnly=True)

import numpy as np

#generate snapshots for solution
S = read_snapshots(archive,nDTout+1,'u',)
U, s, V = np.linalg.svd(S, full_matrices=False)
print('SVD for solution done!')
np.savetxt('SVD_basis', U, delimiter=' ')
np.savetxt('Singular_values', s, delimiter=' ')

if use_deim:
    Sf = read_snapshots(archive,nDTout+1,'spatial_residual0')
    Uf,sf,Vf = np.linalg.svd(Sf,full_matrices=False)
    print('SVD for spatial residual done!')
    np.savetxt('Fs_SVD_basis', Uf, delimiter=' ')
    np.savetxt('Fs_Singular_values', sf, delimiter=' ')
    Sm = read_snapshots(archive,nDTout+1,'mass_residual0')
    Um,sm,Vm = np.linalg.svd(Sm,full_matrices=False)
    print('SVD for mass residual done!')
    np.savetxt('Fm_SVD_basis', Um, delimiter=' ')
    np.savetxt('Fm_Singular_values', sm, delimiter=' ')
    
