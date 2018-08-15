Nonlinear wave generation/absorption
====================================

Description
-----------

Plane regular nonlinear waves are mild and high steepness waves that 
propagate in a single direction, in uniform wave fronts.  The wave 
profile deviates from the sinusoidal shape, and it typically exhibits 
high and sharp wave crests and low and flat wave troughs.  It is not 
always accurate to calculate the wave celerity from the linear 
dispersion theory, especially for highly nonlinear waves. 
Fenton (1988) proposes a method for calculating the nonlinear wave 
properties and profile, which is adopted for the generation of 
nonlinear waves within Proteus. 

In terms of classification, linear waves are in the right top area of the following diagram (Lé Méhauté 1976), where the vertical axis corresponds to the non dimensional wave height and the horizontal to the non dimensional water depth, with gT\ :sup:`2` being proportional to the wavelength.


.. figure:: ./Mehaute_nonlinear_waves_01.png
   :width: 100%
   :align: center

The numerical wave flume represents the geometry used in Higuera et al 2013 for their numerical tests. The file nonlinearTest is a batch script that runs all these tests using context options. 

This case tests demonstrates the ability of PROTEUS to simulate the 
generation, propagation and absorption of regular non-linear waves. 

Tests
------
The python test file named ``test_nonlinearWaves.py`` is made up of three tests:

* The first test checks that the run is completed successfully.
* The second test is to validate the results comparing them to the theory. For this case we will compare the numerical and theoretical wave height in the middle of the tank.
* The third test evaluates wave reflection and compares to a threshold. The calculation of reflection is performed by applying Isaacson's 3rd method (Isaacson 1991) to the primary harmonic of the signal.

One can run this test file typing ``py.test --boxed test_nonlinearWaves.py``.

References
----------

- Fenton JD (1988) The numerical solution of steady water wave 
  problems, Comp and Geosc, 14(3), 357-368
  
- Lé Méhauté, B., (1976). “Introduction to Hydrodynamics and water waves”, Springer-Verlag, New York.

- Isaacson (1991), Measurement of regular wave reflection, Journal of Waterway Port Coastal and Ocean Engineering 117(6), 553-569





