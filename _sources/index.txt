.. Proteus documentation master file. This file generates the Proteus homepage (index.html) and serves as the root of the toctree

.. _intro-sec:

Introduction
============

Proteus is a Python package for rapidly developing computer models and
numerical methods. It is focused on models of continuum mechanical
processes described by partial differential equations and on
discretizations and solvers for computing approximate solutions to
these equations. Proteus consists of a collection of Python modules
and scripts. Proteus also uses several C, C++, and Fortran libraries,
which are either external open source packages or part of Proteus, and
several open source Python packages.

The design of Proteus is organized around two goals:

* Make it easy to solve new model equations with existing numerical methods
* Make it easy to solve existing model equations with new numerical methods

We want to improve the development process for models *and*
methods. Proteus is not intended to be an expert system for solving
partial differential equations. In fact, effective numerical methods
are often physics-based. Nevertheless many physical models are
mathematically represented by the same small set of differential
operators, and effective numerical methods can be developed with minor
adjustments to existing methods. The problem with much existing
software is that the physics and numerics are completely intertwined,
which makes it difficult to extend (and maintain). In Proteus the
description of the physical model and initial-boundary value problems
are nearly "method agnostic".  This approach has been used in the
developement of a variety of mathematical models and numerical
methods, both of which are described in more detail below
(:ref:`capabilities-sec`).

.. _obtaining-sec:

Obtaining and Installing Proteus
================================

For now, Proteus is only available as source from our public `GitHub
<https://github.com/erdc-cm/proteus>`_ repository  because we haven't
had the resources to devote to producing regular binary releases.  We
use the `HashDist <https://hashdist.github.io/>`_ software environment
installation tool to bootstrap any needed dependencies for you.  If
you already have compilers and Git installed on your system, you can
install Proteus with the following commands:

More information is available on our `Wiki <https://github.com/erdc-cm/proteus/wiki>`_, and you can ask for help on
the `Developers' Mailing List <https://groups.google.com/forum/#!forum/proteus-dev>`_.

.. _getting_start-sec:

Getting Started
===============

The following IPython Notebooks provide a gentle introduction to
getting started with Proteus, from setting up a model domain to
solving your own problems with a comprehensive example.  The notebooks
can all be downloaded as an `archive <https://github.com/erdc-cm/proteus-notebooks/archive/master.zip>`_.

* `Defining the Model Equations <http://nbviewer.ipython.org/github/erdc-cm/proteus-notebooks/blob/master/Equation_Tutorial.ipynb>`_
* `Defining a Model Domain <http://nbviewer.ipython.org/github/erdc-cm/proteus-notebooks/blob/master/Domain_Tutorial.ipynb>`_

.. _running-sec:

Running
=======

If you have successfully compiled Proteus then you should be able to do::

   %cd $PROTEUS/test
   %$PROTEUS_PREFIX/bin/parun poisson_3d_p.py poisson_3d_c0p1_n.py

The solution will be saved in a file ending in .xmf, which can be
opened with ParaView or Ensight.

.. _capabilities-sec:

Capabilities
============

Test problems and some analytical solutions have been implemented for

* Poisson's equation
* The heat equation
* Linear advection-diffusion-reaction equations
* Singly degenerate nonlinear advection-diffusion-reaction equations (including various forms of Burger's equation)
* Doubly degenerate nonlinear advection-diffusion-reaction equations
* The eikonal (signed distance) equation
* The diffusive wave equations for overland flow
* 1D and 2D Shallow Water Equations
* Richards' equation (mass conservative head- and saturation-based)
* Two-phase flow in porous media with diffuse interface (fully coupled  and IMPES formulations)
* Two-phase flow in porous media with a sharp interface (level set formulation)
* Stokes equations
* Navier-Stokes equations
* Reynolds-Averged Navier-Stokes equations
* Two-phase Stokes/Navier-Stokes/RANS flow with a sharp interface (level set/VOF formulation) 
* Linear elasticity

These problems are solved on unstructured simplicial meshes. Simple
meshes can be generated with tools included with Proteus, and more
complex meshes can by imported from other mesh generators. The finite
elements implemented are

Classical and vartiational multiscale methods

* :math:`C_0 P_1`
* :math:`C_0 P_2`
* :math:`C_0 Q_1`
* :math:`C_0 Q_2`

Discontinuous Galerkin methods

* :math:`C_{-1} P_0`
* :math:`C_{-1} P_1`  (Lagrange Basis)
* :math:`C_{-1} P_2`  (Lagrange Basis)
* :math:`C_{-1} P_k`  (Monomial Basis)

Non-conforming and mixed methods

* :math:`P_1` non-conforming 
* :math:`C_0 P_1 C_0 P_2` Taylor-Hood

The time integration methods are

* Backward Euler
* Forward Euler
* :math:`\Theta` Methods
* Strong Stability Preserving Runge-Kutta Methods
* Adaptive BDF Methods
* Pseudo-transient continuation

The linear solvers are

* Jacobi
* Gauss-Seidel
* Alternating Schwarz
* Full Multigrid
* Wrappers to LAPACK, SuperLU, and PETSc

The nonlinear solvers are

* Jacobi
* Gauss-Seidel
* Alternating Schwarz
* Newton's method
* Nonlinear Multigrid (Full Approximation Scheme)
* Fast Marching and Fast Sweeping

Additional tools are included for pre- and post-processings meshes and
solutions files generated by Proteus and other models, including methods for
obtaining locally-conservative velocity fields from :math:`C_0` finite
elements.

.. _release-sec:

Release Policy
==============

The releases are numbered major.minor.revision

* A revision increment indicates a bug fix or extension that shouldn't break any models working with the same major.minor number.
* A minor increment introduces significant new functionality but is mostly backward compatible
* A major increment may require changes to input files and significantly change the underlying Proteus implementation.

These are not hard and fast rules, and there is no time table for releases.

References
==========

* `Numerical modeling of drag for flow through vegetated domains and
  porous structures
  <http://dx.doi.org/10.1016/j.advwatres.2012.01.002>`_ (2012)
  S.A. Mattis, C. N. Dawson, C. E. Kees, and M. W. Farthing, *Advances
  in Water Resources*, 39, pp44-59
* `Parallel Computational Methods and Simulation for Coastal and
  Hydraulic Applications Using the Proteus Toolkit
  <http://www.dlr.de/sc/en/Portaldata/15/Resources/dokumente/pyhpc2011/submissions/pyhpc2011_submission_11.pdf>`_
  (2011) C. E. Kees and M. W. Farthing (2011) *Supercomputing11:
  Proceedings of the PyHPC11 Workshop*
* `A Conservative Level Set Method for Variable-Order Approximations
  and Unstructured Meshes
  <http://dx.doi.org/10.1016/j.jcp.2011.02.030>`_ (2011) C.E. Kees,
  I. Akkerman, Y. Bazilevs, and M. W. Farthing *Journal of
  Computational Physics* 230(12), pp4536–4558
* `Locally Conservative, Stabilized Finite Element Methods For
  Variably Saturated Flow
  <http://dx.doi.org/10.1016/j.cma.2008.06.005>`_ (2008) Kees, C.E.,
  M. W. Farthing, and C. N. Dawson, *Computer Methods in Applied
  Mechanics and Engineering*, 197, pp4610-4625
* `Locally Conservative, Stabilized Finite Element Methods for a Class
  of Variable Coefficient Navier-Stokes Equations
  <http://chl.erdc.usace.army.mil/ERDC-CHL-TR-09-12>`_ (2009)
  C. E. Kees, M. W. Farthing, and M. T. Fong, *ERDC/CHL TR-09-12*
* `Evaluating Finite Element Methods for the Level Set Equation
  <http://chl.erdc.usace.army.mil/ERDC-CHL-TR-09-11>`_ (2009)
  M. W. Farthing and C. E. Kees, *ERDC/CHL TR-09-11*
* `A Review of Methods for Moving Boundary Problems
  <http://chl.erdc.usace.army.mil/ERDC-CHL-TR-09-10>`_ (2009)
  C. E. Kees, M. W. Farthing, T. C. Lackey, and R. C. Berger,
  *ERDC/CHL TR-09-10*
* `Implementation of Discontinuous Galerkin Methods for the Level Set
  Equation on Unstructured Meshes
  <http://chl.erdc.usace.army.mil/library/publications/chetn/pdf/chetn-xiii-2.pdf>`_
  (2008) M. W. Farthing and C. E. Kees, *ERDC/CHL CHETN-XIII-2*

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Source Code Documentation
=========================

.. toctree::
   :maxdepth: 1

   api/modules
   cfiles
   

