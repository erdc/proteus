Dambreak flow - Collagrosi and Landrini (2003)
==============================================

Description
-----------
The problem comprises a 0.60m x 1.20m (height x width) column of
water in a 1.8 m high and 3.22 m wide container. The column collapses under the action of gravity
and impacts to a wall. The top of the domain is
left open, when the rest of the boundary patches act as free slip walls.
In the following figure, a sketch of the dambreak initial conditions
is shown.

.. figure:: ./dambreakColagrossi.bmp
   :width: 100%
   :align: center

This case tests the ability of PROTEUS to simulate the free-surface
evolution and forces / pressures on structures, according to data that
are available in the following references.  For more details, see
runfiles or references.

Tests
-------
The python test file named ``test_dambreak_Colagrossi.py`` is made up of 
two tests:

* The first test is to check that the run is successfully completed.
* The second test is to validate the results comparing them to reference values. For this case we will compare the numerical and reference maximum pressure for a given point.
One can run this test file typing ``py.test --boxed ../../../Tests/1st_set/test_dambreak_Colagrossi.py``.

Execution
---------

The case is run using the command:

parun dambreak_Colagrossi_so.py -l 5 -O ../../../inputTemplates/petsc.options.superlu_dist

* Executing parun -h in the command line will give a list of all avaliable options as well as a description of the option.

References
----------

- Colagrossi A and Landrini M (2003) Numerical simulation of
  interfacial flows by smoothed particle hydrodynamics, Journal of
  Computational Physics,191,448-475.

- Martin, J. C. & Moyce, W. J., (1952) Part IV. An Experimental Study
  of the Collapse of Liquid Columns on a Rigid Horizontal Plane
  Phil. Trans. R. Soc. Lond. A 244 (882) 312-324.

- Zhou, Z. Q., De Kat, J. O. and Buchner, B. (1999) A nonlinear 3-D
  approach to simulate green water dynamics on deck in: J. Piquet
  (Ed.), Proc. 7th Int. Conf. Num. Ship Hydrod., Nantes, 5.11, 15.
