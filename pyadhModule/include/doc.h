/*!
\file doc.h

\brief This file is for general documentation on PyADH and for placing
documentation outside of source files.

The primary purpose of this file is for the PyADH main page. If for
some reason we want to place documentation for a source code file
outside of the source code, this documentation can be placed after the
main page definition in this file.

\mainpage 

\section intro_sec Introduction

PyADH is a set of modules for rapidly developing computer models and
numerical methods for partial differential equations. The modules are
implemented as C/C++/FORTRAN libraries, Python modules, and Python
\ref scripts "scripts". PyADH has been used to solve three general
classes of partial differential equations:

- Systems of linear and nonlinear advection-diffusion-reaction equations 
- Linear and nonlinear Hamilton-Jacobi equations 
- Stokes and Navier-Stokes equations
- Combinations of the above using split operator formulations

\section obtaining_sec Obtaining PyADH

The best way to get started is to download a pre-built package for you
platform below. PyADH currently depends on many other open source projects,
so that obtaining and building all the depencies can be time
consuming. The pre-built packages contain all the binaries and source
code necessary to develop in PyADH, including the svn repository
information and the pyadhGraphical module for visualization and
graphics post-processing. With one of these packages you can update
the source and rebuild portions of the code if necessary. If one of these doesn't work or if you wish to start from scratch or update your source, see \ref svn below.

Mac OSX 10.5 platforms:

- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-darwin_i386_macports.tar.gz </tt>
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-darwin_ppc_macports.tar.gz </tt>
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-darwin_x86_64.tar.gz </tt>
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-darwin_ppc64.tar.gz </tt>

GNU/Linux platforms:

- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-linux.tar.gz </tt> (RHEL 5)
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-jade.tar.gz </tt> (Cray XT4 CNL at ERDC MSRC)
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-lonestar.tar.gz </tt> (Dell Cluster at TACC)
- <tt> https://adh.usace.army.mil/pyadhReleases/pyadh-0.5.0-ranger.tar.gz </tt> (Sun Constellation at TACC)

\section installing_sec Installing and Running

This example is for the linux64 package and csh.

<tt> \%tar xzf pyadh-0.5.0-linux64.tar.gz </tt> \n
<Tt> \%setenv PYADH $PWD/pyadh </tt> \n
<tt> \%setenv PYADH_ARCH linux64 </tt> \n
<tt> \%setenv PATH $PYADH/$PYADH_ARCH/bin:$PATH </tt> \n
<tt> \%setenv LD_LIBRARY_PATH $PYADH/$PYADH_ARCH/lib:$LD_LIBRARY_PATH </tt> \n

Replace linux64 with darwin_i386_macports, jade, etc. for one of the other architechtures. On MacOS also you also need to add

<tt> \%setenv PATH $PYADH/$PYADH_ARCH/Python.framework/Versions/Current/bin:$PATH </tt>

To use runtime visualization you may also need to add

<tt> \%setenv LD_LIBRARY_PATH $PYADH/$PYADH_ARCH/lib/vtk-5.3:$LD_LIBRARY_PATH </tt>

To run an example do

<tt> \%cd $PYADH/pyadhModule/test/problemDescriptions </tt> \n
<tt> \%parun navier_stokes_cylinder_2d_p.py navier_stokes_cylinder_2d_c0p1c0p1_n.py -l 3 -v</tt> \n

The solution will be saved in a file ending in .xmf, which can be opened with ParaView3 or Ensight. ParaView3 is included in the binary package in $PYADH/$PYADH_ARCH.

\section building_sec Compiling

To get updated source and recompile the pyadh module you can do

<tt> \%cd $PYADH/pyadhModule </tt> \n
<Tt> \%svn update </tt> \n
<tt> \%python setup.py install </tt> \n

\section release Release Policy

The releases are numbered major.minor.revision

- A revision increment indicates a bug fix or extension that shouldn't break any models working with the same major.minor number.
- A minor increment introduces significant new functionality but is mostly backward compatible
- A major increment may require changes to input files and significantly change the underlying PyADH implementation.

These are not hard and fast rules, and there is no time table for releases.

\section svn SVN Repository Access

To obtain or update the source use an svn client to connect to our SVN repository. You can request an account at 

<tt> https://adh.usace.army.mil/svnaccess  </tt>

Select PYADH_R if you only need read access or PYADH_RW if you might
want to make changes. You can then run 'svn update' from your pyadh or
pyadhGraphical subdirectory. To check out a clean copy of the source
you can do the following:

<tt> \%svn co https://adh.usace.army.mil/svn/pyadh/trunk pyadh </tt>

You may also want pyadhGraphical:

<tt> \%svn co https://adh.usace.army.mil/svn/pyadhGraphical/trunk pyadhGraphical </tt>

Follow the directions in the source to build the packages.

\section testing_sec Testing PyADH

Fully coupled models are specified by defining a p-file ("problem
file") and a n-file ("numerics file").  By convention, p-files are
given an ending <tt> _p.py</tt>, while numerics file names consist of
the problem name and a trailing identifier <tt> _meth_n.py </tt> where
<tt> meth </tt> is a tag identifying the spatial approximation used
for the test problem. For multi-physics models an sso-file ("split
operator file") describes a list of p- and n-files as well as
additional information on how the solutions are to be coordinated.

See the \ref Tests "Test problems" section for a list of existing
problems.

To run a problem consisting of a single mathematical model and
visualize the output with vtk one can use the <tt> parun
</tt> script 

<tt> %parun problem_p.py method_n.py -V vtk </tt>

For multiple PDE models one can use <tt> parun </tt> with the sso-file

<tt> %parun splitting.sso -V vtk </tt>

where the file <tt> splitting.sso </tt> defines the <tt> problem_p.py
method_n.py </tt> pairs for each model in the simulation.

Running in parallel depends on the mpi environment, but generally looks  something like

<tt> %mpiexec -n 16 python parun splitting.sso </tt>

Use the <tt> parun --help </tt> to get a list of options.

\section capabalities_sec Capabilities

\ref Tests "Test problems" and some \ref analyticalSolutions "analytical solutions" have been
implemented for

- Poisson's equation
- The heat equation
- Linear advection-diffusion-reaction equations
- Singly degenerate nonlinear advection-diffusion-reaction equations (including various forms of Burger's equation)
- Doubly degenerate nonlinear advection-diffusion-reaction equations
- The eikonal (signed distance) equation
- The travel time equation
- Richards' equation (mass conservative head- and saturation-based)
- Two-phase flow in porous media with diffuse interface (fully coupled  and IMPES formulations)
- Two-phase flow in porous media with a sharp interface (level set formulation)
- Stokes equations
- Navier-Stokes equations
- Two-phase Stokes flow with a sharp interface (level set/VOF formulation) 
- Two-phase Navier-Stokes flow with a sharp interface (level set/VOF formulation)

These problems are solved on unstructured simplicial meshes. Simple
meshes can be generated with tools included with PyADH, and more
complex meshes can by imported from other mesh generators. The finite
elements implemented are

Classical and vartiational multiscale methods:
- \f$ C_0 P_1 \f$
- \f$ C_0 P_2  \f$

Discontinuous Galerkin methods
- \f$ C_{-1} P_0 \f$
- \f$ C_{-1} P_1\f$  (Lagrange Basis)
- \f$ C_{-1} P_2\f$  (Lagrange Basis)
- \f$ C_{-1} P_k\f$  (Monomial Basis)

Non-conforming and mixed methods
- \f$ P_1 \f$ non-conforming 
- \f$ C_0 P_1 C_0 P_2 \f$ Taylor-Hood
- Mini element

The time integration methods are

- Backward Euler
- Forward Euler
- \f$ \Theta \f$ Methods
- Strong Stability Preserving Runge-Kutta Methods
- Adaptive BDF Methods
- Pseudo-transient continuation

The linear solvers are

- Jacobi
- Gauss-Seidel
- Alternating Schwarz
- Full Multigrid
- Wrappers to LAPACK, SuperLU, and PETSc

The nonlinear solvers are

- Jacobi
- Gauss-Seidel
- Alternating Schwarz
- Newton's method
- Nonlinear Multigrid (Full Approximation Scheme)
- Fast Marching and Fast Sweeping

Additional tools are included for pre- and post-processings meshes and
solutions files generated by ADH and PyADH, including methods for
obtaining locally-conservative velocity fields from \f$C_0\f$ finite
elements.

\section references_sec References

- <A HREF="http://dx.doi.org/10.1016/S0309-1708(01)00022-7"> Mixed finite element methods and higher-order temporal approximations</A> (2002) M. W. Farthing, C. E. Kees, and C. T. Miller, Advances in Water Resources 25(1), 85-101.
- <A HREF="http://dx.doi.org/10.1016/S0309-1708(02)00187-2"> Mixed finite element methods and higher order temporal approximations for variably saturated groundwater flow</A> (2003) M. W. Farthing, C. E. Kees, and C. T. Miller, Advances in Water Resources 26(4), 373-394.
- <A HREF="http://dx.doi.org/10.1016/S0309-1708(03)00076-9"> Efficient steaady-state solution techniques for variably saturated groundwater flow</A> (2003) M. W. Farthing, C. E. Kees, T. S. Coffrey, C. T. Kelley, and C. T. Miller, Advances in Water Resources 26(8), 833-849.
- <A HREF="http://dx.doi.org/10.1016/S0309-1708(01)00054-9"> Higher order time integration methods for two-phase flow</A> (2002) C. E. Kees and C. T. Miller, Advances in Water Resources 25(2), 159-177.
- <A HREF="http://doi.acm.org/10.1145/332242.334001"> C++ implementations of numerical methods for solving differential-algebraic equations: design and optimization considerations</A> (1999) C. E. Kees and C. T. Miller, ACM Transactions on Mathematical Software 25(4), 377-403.

*/
/*
  hidden
- <A HREF="http://www.springerlink.com/content/w8j2057163458656"> Inexact Newton methods and the method of lines for solving Richards' equation in two space dimensions</A> (1998) M. D. Tocci, C. T. Kelley, C. T. Miller, and C. E. Kees, Computational Geosciences 2(4), 291-309.
- <A HREF="http://dx.doi.org/10.1016/j.advwatres.2008.01.010"> Comparison of derivative-free optimization methods for groundwater supply and hydraulic capture community problems</A> (2008) K. R. Fowler, J. P. Reese, C. E. Kees, J. E. Dennis Jr., C. T. Kelley, C. T. Miller, C. Audet, A. J. Booker, G. Couture, R. W. Darwin, M. W. Farthing, D. E. Finkel, J. M. Gablonsky, G. Gray, and T. G. Kolda, Advances in Water Resources 31(5), 743-757.
- <A HREF="http://www.springerlink.com/content/m6676t73436567r3"> Multiphase Flow and Transport Modeling in Heterogeneous Porous Media</A> (2006) R. Helmig, C. T. Miller, H. Jakobs, H. Class, M. Hilpert, C. E. Kees, and J. Niessner, Mathematics in Industry 8(1), 449-488. 
- <A HREF="http://www.springerlink.com/content/q042015437386q73"> Versatile Two-Level Schwarz Preconditioners for Multiphase Flow</A> (2003) C. E. Kees, C. T. Miller, E. W. Jenkins, and C. T. Kelley, Computational Geosciences 7(2), 1573-1499.
- <A HREF="http://dx.doi.org/10.1016/j.advwatres.2005.07.001"> An ELLAM approximation for advective-dispersive transport with nonlinear sorption</A> (2006) M. W. Farthing, C. E. Kees, T. F. Russell, and C. T. Miller, Advances in Water Resources 29(5), 657-675.
- <A HREF="http://www.springerlink.com/content/t5j1732160xq3238"> Solution of a Well-Field Design Problem with Implicit Filtering</A> (2004) K. R. Fowler, C. T. Kelley, C. T. Miller, C. E. Kees, R. W. Darwin, J. P. Reese, M. W. Farthing, and M. S. C. Reed, Optimization and Engineering 5(2), 207-234.
  
 */
/*! \defgroup scripts Scripts 
  
  \brief Python scripts for various tasks. 
*/

/*! \defgroup test test
  \brief A test set for 
verifying the PyADH implementation
*/

/* Put external documentation below this line */
