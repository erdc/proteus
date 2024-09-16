# Proteus: Computational Methods and Simulation Toolkit [![Build Status](https://travis-ci.com/cekees/proteus.svg?branch=main)](https://app.travis-ci.com/github/cekees/proteus) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/erdc/proteus_tutorial/master?filepath=index.ipynb)  [![DOI](https://zenodo.org/badge/2212385.svg)](https://zenodo.org/badge/latestdoi/2212385)


Proteus (http://proteustoolkit.org) is a Python package for
rapidly developing computer models and numerical methods.

# Installation



```bash
mamba install proteus -c conda-forge
```

For a development installation, you want to install Proteus's dependencies and compile Proteus from source:

```bash
mamba env create -f environment-dev.yml #environment-dev-up.yml to try unpinned dependencies
mamba activate proteus-dev
pip install -v -e .
```

# HPC Installation

For installation on high performance environments you may want to install Proteus's dependencies from source as well. We recommend using the PETSc build system to install most of the dependencies. The following is general outline:

Create a basic python environment you can install into:

```
mamba env create -f petsc-dev.yml
mamba activate petsc-dev
```

or

```
python -m venv petsc-dev
pip install setuptools make cmake cython swig pybind11 numpy
```

Next, build petsc from source

```
bash scripts/petsc_config_linux_conda_seq.sh #see https://petsc.org/release/install/
```

Next, build additional C++ dependencies and install into environment prefix (e.g $CONDA_PREFIX or $VIRTUAL_ENV)

https://github.com/projectchrono/chrono #>=9.0.1, with python bindings
https://github.com/scorec/core #>=2.2.8, shared, python bindinds not needed
https://github.com/xtensor-stack/xtl
https://github.com/xtensor-stack/xtensor
https://github.com/xtensor-stack/xtensor-python

Finally, locally build and install remaining python dependencies into environment

```
CC=mpicc CXX=mpicxx MPI_DIR=$MPI_ROOT pip install -v mpi4py==3.1.6 --no-build-isolation --no-binary=:all:
PETSC_DIR=$CONDA_PREFIX PETSC_ARCH="" pip install -v ../petsc/src/binding/petsc4py --no-build-isolation --no-binary=:all:
HDF5_MPI=ON HDF5_DIR=${CONDA_PREFIX} CC=mpicc CXX=mpicxx pip install -v h5py --no-build-isolation --no-binary=:all: #note h5py depends on mpi4py for this config
CC=mpicc CXX=mpicxx pip install -v . --no-build-isolation --no-binary=:all:
```

Some optional packages can be installed with pip/mamba: py2gmsh, memory_profiler, scipy, pytest.

See https://github.com/erdc/proteus/wiki/How-to-Build-Proteus for old information on building the entire stack.

# Developer Information

The source code, wiki, and issue tracker are on GitHub at

https://github.com/erdc/proteus.
