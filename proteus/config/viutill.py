from .default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=['-L/app/acml/acml-4.4.0/gfortran64/lib','-lacml'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/app/acml/acml-4.4.0/gfortran64/lib','-lacml']
PROTEUS_BLAS_LIB_DIR = '/app/acml/acml-4.4.0/gfortran64/lib'
PROTEUS_BLAS_LIB   = 'acml'
PROTEUS_LAPACK_LIB_DIR = '/app/acml/acml-4.4.0/gfortran64/lib'
PROTEUS_LAPACK_LIB = 'acml'
