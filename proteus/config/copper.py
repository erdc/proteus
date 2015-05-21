from default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=['-L/cm/shared/apps/openblas/dynamic/0.2.6/lib','-lopenblas'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/cm/shared/apps/openblas/dynamic/0.2.6/lib','-lopenblas']
PROTEUS_BLAS_LIB_DIR = '/cm/shared/apps/openblas/dynamic/0.2.6/lib'
PROTEUS_BLAS_LIB   = 'openblas'
PROTEUS_LAPACK_LIB_DIR = '/cm/shared/apps/openblas/dynamic/0.2.6/lib'
PROTEUS_LAPACK_LIB = 'openblas'
