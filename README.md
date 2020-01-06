# Proteus: Computational Methods and Simulation Toolkit [![Build Status](https://travis-ci.com/erdc/proteus.svg?branch=master)](https://travis-ci.com/erdc/proteus) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/erdc/training_proteus/master?filepath=index.ipynb)  [![DOI](https://zenodo.org/badge/2212385.svg)](https://zenodo.org/badge/latestdoi/2212385)


Proteus (http://proteustoolkit.org) is a Python package for
rapidly developing computer models and numerical methods.

# Installation

The recommended way is the following:

```bash
make develop
make test
```
https://github.com/erdc/proteus/wiki/How-to-Build-Proteus-Using-HashDist

However, we are transitioning to a conda-based environment. Proteus can be installed with:

```bash
conda install proteus -c conda-forge
```

For a development installation, you want to install Proteus's dependencies and compile Proteus from source:

```bash
conda env create -f environment-dev.yml
conda activate proteus-dev
python setup.py install
```

# Developer Information

The source code, wiki, and issue tracker are on GitHub at
https://github.com/erdc/proteus.
