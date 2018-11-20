BoundaryConditions
******************

:doc:`BoundaryConditions`

Usage
=====

This module offers a framework to deal with boundary conditions that can either
be set manually or by using predefined functions that have been optimised
through cython/C++ code.

Import
------

The base `BoundaryConditions` module can be imported as follows:

.. code-block:: python

   from proteus import BoundaryConditions as bc

The `mprans` version of `BoundaryConditions`, adding functionality specific to
multi-phase flows, can be imported as follows:

.. code-block:: python

   from proteus.mprans import BoundaryConditions as bc

.. note::

   The base `proteus.BoundaryConditions` only initialises empty
   `BoundaryCondition` instances. In practice,
   `proteus.mprans.BoundaryConditions` is always the module imported for
   multi-phase flow applications. It is also usually used through
   :doc:`SpatialTools <./spatial_tools>` without being directly imported.
