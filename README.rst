Proteus: Computational Methods and Simulation Toolkit
======================================================

Proteus (http://proteus.usace.army.mil) is a Python package for
rapidly developing computer models and numerical methods.

Installation
=============

Set the following environment variables and type 'make' in the root directory.

::

  setenv PROTEUS ${HOME}/proteus #or wherever you put it
  setenv PROTEUS_ARCH linux
  setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
  setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
  setenv PATH .:${PROTEUS_PREFIX}/bin:${PATH}
  setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
  setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib


Developer Information
======================

The source code, wiki, and issue tracker are on github at
https://github.com/erdc-cm/proteus. The developers' mailing list is
http://groups.google.com/group/proteus-dev. Both require approval at
this time.
