#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 1:00:00
#SBATCH -p workq
#SBATCH -A hpc_proteus02o
#SBATCH -J proteus_build
#load proteus module and ensure proteus's python is in path
module load proteus/fct
ll="================================================================================"
N=64
info () {
    echo $ll
    echo "$*"
    echo $ll
    }
platform="$(uname -s)"
case "${platform}" in
    Linux*)     LIB_SUFFIX="so";;
    Darwin*)    LIB_SUFFIX="dyld";;
    *)          LIB_SUFFIX="so";;
esac
info "building on " ${platform}
export PROTEUS_PREFIX=${PWD}/proteus_env
export PROTEUS_ARCH='mike'
PYVER=$(python -c 'from distutils import sysconfig; print(sysconfig.get_python_version())')
info "creating virtual env"
python -m venv ${PROTEUS_PREFIX}
source ${PWD}/proteus_env/bin/activate
info "installing petsc python dependencies"
python -m pip install --upgrade pip
pip install cython numpy mpi4py
info "installing petsc"
git clone https://gitlab.com/petsc/petsc.git
cd petsc
./configure --prefix=${PROTEUS_PREFIX} --with-mpi-compilers -CC=${MPICC} -CXX=${MPICXX} -F90=${MPIF90} -F77=${MPIF77} --download-hdf5 --download-tetgen --download-triangle  --download-petsc4py --download-superlu --download-parmetis --download-metis --download-superlu_dist --download-zlib --download-hypre --download-zoltan --download-eigen --download-cmake
make PETSC_DIR=${PWD} all
make install
cd ..
info "Installing post-petsc python packages"
export PYTHONPATH=${PROTEUS_PREFIX}/lib
CC=${MPICC} HDF5_MPI="ON" HDF5_DIR=${PROTEUS_PREFIX} pip install h5py scipy pybind11 swig future
info "installing xtensor stack: xtl, xtensor, xtensor-python"
info "installing xtl"
git clone https://github.com/xtensor-stack/xtl.git
cd xtl
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} ..
make install
cd ../..
info "installing xtensor"
git clone https://github.com/xtensor-stack/xtensor.git
cd xtensor
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} ..
make install
cd ../..
info "installing xtensor-python"
git clone https://github.com/xtensor-stack/xtensor-python.git
cd xtensor-python
mkdir build
cd build
pybind11_DIR=${PROTEUS_PREFIX}/lib/python${PYVER}/site-packages/pybind11 cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} ..
make install
cd ../..
info "installing chrono"
git clone https://github.com/projectchrono/chrono.git
cd chrono
git checkout 8.0.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} -DBUILD_DEMOS:BOOL=OFF -DENABLE_MODULE_CASCADE:BOOL=OFF -DENABLE_UNIT_CASCADE:BOOL=ON -DENABLE_MODULE_IRRLICHT:BOOL=OFF -DENABLE_MODULE_POSTPROCESS:BOOL=OFF -DENABLE_MODULE_VEHICLE:BOOL=OFF -DENABLE_MODULE_FSI:BOOL=OFF -DENABLE_OPENMP:BOOL=ON -DENABLE_MODULE_PYTHON:BOOL=ON -DENABLE_MODULE_COSIMULATION:BOOL=OFF -DENABLE_MODULE_MATLAB:BOOL=OFF -DENABLE_MODULE_MKL:BOOL=OFF -DENABLE_MODULE_PARALLEL:BOOL=OFF -DENABLE_MODULE_OPENGL:BOOL=OFF -DENABLE_MODULE_OGRE:BOOL=OFF -DCMAKE_BUILD_TYPE:STRING=Release -DPYTHON_EXECUTABLE:PATH=${PROTEUS_PREFIX}/bin/python -DPYTHON_INCLUDE_DIR:PATH=$(python -c 'from distutils import sysconfig; print(sysconfig.get_python_inc())') -DPYTHON_LIBRARY:PATH=$(python -c 'from distutils import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')/libpython${PYVER}.${LIB_SUFFIX} -DSWIG_EXECUTABLE=${PROTEUS_PREFIX}/bin/swig ..
make -j ${N} all
make install
cd ../..
info "installing scorec"
git clone https://github.com/SCOREC/core.git
cd core
git checkout v2.2.7
mkdir build
rm -rf pumi-meshes
git clone https://github.com/SCOREC/pumi-meshes.git pumi-meshes
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_EXES:BOOL=ON -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -DCMAKE_C_COMPILER:STRING=${MPICC} -DCMAKE_CXX_COMPILER:STRING=${MPICXX} -DSCOREC_CXX_WARNINGS:BOOL=OFF -DENABLE_ZOLTAN:BOOL=ON -DSCOREC_CXX_OPTIMIZE:BOOL=ON -DMDS_ID_TYPE:STRING=int -DPCU_COMPRESS=ON -DENABLE_FPP=ON -DENABLE_SIMMETRIX:BOOL=OFF -DZOLTAN_PREFIX:PATH=${PROTEUS_PREFIX} -DPARMETIS_PREFIX:PATH=${PROTEUS_PREFIX} -DCMAKE_BUILD_TYPE:STRING=Release ..
make -j ${N} all
make install
cd ../..
info "installing gmsh"
git clone https://gitlab.onelab.info/gmsh/gmsh.git
cd gmsh
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${PROTEUS_PREFIX} -DCMAKE_CXX_COMPILER:FILEPATH=${MPICXX} -DCMAKE_C_COMPILER:FILEPATH=${MPICC} -DENABLE_BUILD_SHARED:BOOL=ON -DENABLE_WRAP_PYTHON:BOOL=ON ..
make -j ${N} all
make install
cd ../..
info "intalling triangle executable"
mkdir triangle
cd triangle
curl https://netlib.org/voronoi/triangle.zip --output triangle.zip
unzip triangle.zip
sed -e s/-DLINUX//g makefile > makefile.darwin
make -f makefile.darwin triangle
cp triangle ${PROTEUS_PREFIX}/bin
cd ..
info "installing tetgen executable"
curl https://wias-berlin.de/software/tetgen/1.5/src/tetgen1.6.0.tar.gz --output tetgen1.6.0.tar.gz
tar xzvf tetgen1.6.0.tar.gz
cd tetgen1.6.0
make tetgen
cp tetgen ${PROTEUS_PREFIX}/bin
cd ..
info "install proteus"
git clone https://github.com/cekees/proteus
cd proteus
N=${N} python setup.py build_ext
pip install -v .
