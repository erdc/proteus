.. _boundary_conditions:

BoundaryConditions
******************

source: :py:mod:`proteus.BoundaryConditions` and
:py:mod:`proteus.mprans.BoundaryConditions`

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

   The base :py:mod:`proteus.BoundaryConditions` only initialises empty
   :py:class:`proteus.BoundaryConditions.BC_Base` instances. In practice,
   :py:mod:`proteus.mprans.BoundaryConditions` is always the module imported
   for multi-phase flow applications, with pre-populated
   :py:mod:`proteus.mprans.BoundaryConditions.BC_RANS` that include the most
   common BCs used in Proteus. It is also usually used through
   :ref:`spatial_tools` without being directly imported.
