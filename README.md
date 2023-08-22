# Proteus: Computational Methods and Simulation Toolkit [![Build Status](https://travis-ci.com/cekees/proteus.svg?branch=main)](https://travis-ci.com/cekees/proteus) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/erdc/proteus_tutorial/master?filepath=index.ipynb)  [![DOI](https://zenodo.org/badge/2212385.svg)](https://zenodo.org/badge/latestdoi/2212385)


Proteus (http://proteustoolkit.org) is a Python package for
rapidly developing computer models and numerical methods.

# Installation



```bash
conda install proteus -c conda-forge
```

For a development installation, you want to install Proteus's dependencies and compile Proteus from source:

```bash
conda env create -f environment-dev.yml
conda activate proteus-dev
make develop-conda # or pip install -v -e .
```

You can also build proteus and dependencies from source (without conda) with 

```bash
make develop
make test
```

See https://github.com/erdc/proteus/wiki/How-to-Build-Proteus for more information on building the entire stack.

# Developer Information

The source code, wiki, and issue tracker are on GitHub at
https://github.com/erdc/proteus.
