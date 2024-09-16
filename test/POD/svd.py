#!/usr/bin/env python
from read_hdf5 import *

T = 1.0
nDTout = 100
DT = T/float(nDTout)

archive = Archiver.XdmfArchive(".","heat_3d",readOnly=True)

import numpy as np

u = read_from_hdf5(archive.hdfFile,'/u0')
S = np.reshape(u,(np.shape(u)[0],1))
for i in range(1,nDTout+1):
    time_level_to_read=i
    label="/%s%d" % ('u',time_level_to_read)
    print('trying to read from %s ' % label)
    u = read_from_hdf5(archive.hdfFile,label)
    u = np.reshape(u, (np.shape(u)[0],1))
    S = np.append(S,u,axis=1)

U, s, V = np.linalg.svd(S, full_matrices=False)
print('SVD done!')
np.savetxt('SVD_basis', U, delimiter=' ')
np.savetxt('Singular_values', s, delimiter=' ')

