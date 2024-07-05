from .default import *
    
platform_extra_compile_args = ['-DPETSC_INCLUDE_AS_C']
    
PROTEUS_PRELOAD_LIBS=['/p/home/apps/sgi/mpt-2.12-sgi712r26/lib/libmpi.so',
'/p/home/apps/sgi/mpt-2.12-sgi712r26/lib/libmpi_lustre.so',
'/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/libmkl_rt.so']
#,
#'/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/libmkl_sequential.so',
#'/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/libmkl_core.so',
#'/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/libmkl_sequential.so',
#'/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/libmkl_avx2.so']
#PROTEUS_PRELOAD_LIBS=['/p/home/apps/sgi/mpt-2.12-sgi712r26/lib/libmpi.so','/p/home/apps/sgi/mpt-2.12-sgi712r26/lib/libmpi_lustre.so']
#PROTEUS_PRELOAD_LIBS=['/p/home/apps/sgi/mpt-2.12-sgi712r26/lib/libmpi_lustre.so']

PROTEUS_EXTRA_LINK_ARGS=['-L'+pjoin(prefix,'lib'),'-llapacke','-lrefblas','-ltmglib','-lpthread','-lm','-lgfortran'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=[]
PROTEUS_BLAS_LIB_DIR = '/app/unsupported/COST/gotoblas2/1.13/gnu/lib'
PROTEUS_BLAS_LIB   = 'goto2'
PROTEUS_LAPACK_LIB_DIR = '/app/unsupported/COST/lapack/3.5.0/gnu/lib'
PROTEUS_LAPACK_LIB = 'lapack'
