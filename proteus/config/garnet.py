from default import *

PROTEUS_EXTRA_COMPILE_ARGS=['-Wall',
                            '-DF77_POST_UNDERSCORE',
                            '-DUSE_BLAS',
                            '-DCMRVEC_BOUNDS_CHECK',
                            '-DMV_VECTOR_BOUNDS_CHECK',
                            '-I/opt/cray/hdf5-parallel/1.8.13/GNU/49/include']+platform_extra_compile_args
PROTEUS_EXTRA_LINK_ARGS=['-L/opt/cray/hdf5-parallel/1.8.13/GNU/49/lib',
                         '-lhdf5',
                         '-Wl,-rpath=/opt/cray/hdf5-parallel/1.8.13/GNU/49/lib',
                         '-L/opt/acml/5.3.1/gfortran64/lib','-lacml',
                         '-Wl,-rpath=/opt/acml/5.3.1/gfortran64/lib'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/opt/acml/5.3.1/gfortran64/lib','-lacml']
PROTEUS_BLAS_LIB_DIR = '/opt/acml/5.3.1/gfortran64/lib'
PROTEUS_BLAS_LIB   = 'acml'
PROTEUS_LAPACK_LIB_DIR = '/opt/acml/5.3.1/gfortran64/lib'
PROTEUS_LAPACK_LIB = 'acml'
