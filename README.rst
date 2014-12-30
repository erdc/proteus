Proteus: Computational Methods and Simulation Toolkit [![Build Status](https://api.shippable.com/projects/545531c244927f89db3e7d70/badge?branchName=master)](https://app.shippable.com/projects/545531c244927f89db3e7d70/builds/latest) [![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/erdc-cm/proteus-public)
======================================================

Proteus (http://proteus.usace.army.mil) is a Python package for
rapidly developing computer models and numerical methods.

Installation
=============

Make sure the following environment variables ``PROTEUS,
PROTEUS_PREFIX, PROTEUS_ARCH,`` and  ``PROTEUS_PYTHON`` are
undefined.

Then, change to the root directory and type ``make; make check``.

This will use the HashDist build system to automatically build all
Proteus dependencies for you.

After installation, you should receive instructions on how to update your 
``PATH`` environment variable so that the Proteus-installed python and
other packages are available.  You may want to put this command in your 
bash startup profile.

Developer Information
======================

The source code, wiki, and issue tracker are on GitHub at
https://github.come/erdc-cm/proteus. The developers' mailing list is
http://groups.google.com/group/proteus-dev.
