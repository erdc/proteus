from default import *

PROTEUS_EXTRA_LINK_ARGS=['-L/opt/acml/5.3.1/gfortran64/lib','-lacml','-Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/opt/acml/5.3.1/gfortran64/lib','-lacml']
PROTEUS_BLAS_LIB_DIR = '/opt/acml/5.3.1/gfortran64/lib'
PROTEUS_BLAS_LIB   = 'acml'
PROTEUS_LAPACK_LIB_DIR = '/opt/acml/5.3.1/gfortran64/lib'
PROTEUS_LAPACK_LIB = 'acml'
