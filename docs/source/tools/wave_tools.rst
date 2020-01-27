.. _wave_tools:

WaveTools
*********

source: :py:mod:`proteus.WaveTools`

Usage
=====

This module offers a framework for calculating the free-surface elevation and wave velocities based on various wave theories. Wave theories are organised in classes within the module. The module is written in Python / Cython and C++ for optimising calculation speed.  

Import in command line
----------------------

Once installed Proteus, you can open a python or ipython command line and type 

.. code-block:: python

   from proteus import WaveTools as wt

You can see information on the module (including available classes) by typing:

.. code-block:: python

   help(wt)

and you can see a list of classes and functions by typing 

.. code-block:: python

   wt.__all__
   
Each function or class has documentation info which is accessible by typing

.. code-block:: python

   help(wt.function)
   help(wt.class)
   help(wt.class.function)
  
List of wave theories
---------------------
Available classes that correspond to wave theories are

``SteadyCurrent`` - Introduce steady currents with ramp time

``SolitaryWave`` - Generate 1st order solitary waves

``MonochromaticWaves`` - Generate linear and nonlinear monochromatic waves. Nonlinear wave theory according to `Fenton's Fourier transform <http://johndfenton.com/Steady-waves/Fourier.html>`_

``NewWave`` - Generate waves according to NewWave theory `Tromans et al 1991 <https://www.onepetro.org/conference-paper/ISOPE-I-91-154>`_

``RandomWaves`` - Generate plane random waves from JONSWAP or Pierson Moskovitch 
        
``MultiSpectraRandomWaves`` - Generate random waves by overlaying multiple frequency spectra. Wave spectra can meet in different angles
 
``DirectionalWaves`` - Generate random waves using JOSNWAP / PM spectra for frequencies and cos-2s/Mitsuyashu spectra for directions
 
``TimeSeries`` - Generate waves from a given free-surface time series. Time series can be reconstructed using direct or windowed methods 
 
``RandomWavesFast`` - Same as ``RandomWaves`` only much more computationally efficient, see `Dimakopoulos et al 2019 <https://www.icevirtuallibrary.com/doi/abs/10.1680/jencm.17.00016>`_
 
``RandomNLWaves``  - Generate plane random waves from JONSWAP or RM, using 2nd order theory, following `Dalzell's formulae <https://www.sciencedirect.com/science/article/abs/pii/S0141118799000085>`_
 
``RandomNLWavesFast`` - Same as  ``RandomNLWaves`` only much more computationally efficient, by using the approach of `Dimakopoulos et al 2019 <https://www.icevirtuallibrary.com/doi/abs/10.1680/jencm.17.00016>`_
 
``CombineWaves`` - Generate waves by combining any of the wave theories above

How to use in Proteus
---------------------

The wave tools module is loaded at the preample as in the case of the command line:

.. code-block:: python

   from proteus import WaveTools as wt
   
 
Then the target wave theory is set, by initializing the relevant class as follows

.. code-block:: python

   wave = wt.MonochromaticWaves(period=1.,
                                 waveHeight=0.1,
                                 mwl=0.5,
                                 depth=0.5,
                                 g=np.array([0,-9.81,0]),
                                 waveDir=np.array([1,0,0]),
			         waveType="Fenton", 
			         autoFenton=True,
			         Nf=8)
            
The wave theory is passed through the ``setUnsteadyTwoPhaseVelocityInlet`` boundary condition as follows:

.. code-block:: python

   tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

If the relaxation zone method is used, then the class should be passed through the relevant ``setGenerationZones`` function

.. code-block:: python

   tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing = smoothing)
   
   
Guidance for using the ``setUnsteadyTwoPhaseVelocityInlet`` and ``setGenerationZones`` functions are given in the :doc:`BoundaryConditions <./boundary_conditions>` and :doc:`Spatial Tools <./spatial_tools>` sections of the documentation

Simple examples of usage within the context of a 2D numerical tank can be found in `air-water-vv <https://github.com/erdc/air-water-vv/tree/master/2d/numericalTanks>`_

How to use as stand-alone tool
------------------------------

After importing the tool in a python interface (command line, editor) following the examples above, you can load a class that corresponds to a wave theory, as follows:

.. code-block:: python
   
   wave = wt.RandomWavesFast(Tstart=0.,
                         Tend=5000.,
                         x0=np.array([0.,0.,0.])
                         Tp=2.5,
                         Hs=0.1,
                         mwl=0.5,
                         depth=0.5,
                         waveDir=np.array([1,0,0]),
                         g=np.array([0,-9.81,0]),
                         N=2000,
                         bandFactor=2.,
                         spectName="JONSWAP",
                         Lgen=1.,
                         Nwaves=16,
                         Nfreq=32,
                         checkAcc=True,
                         fast=True)

Then the free surface and velocity for a point in space and time can be calculated as follows:

.. code-block:: python

   x0 = [1.,0.,0.]
   t0 = 0.
   U = wave.u(x0,t0)

Full time series can be calculated and plotted by appropriately manipulating the calculations and storing in arrays, e.g.:

.. code-block:: python

   x0 = [1.,0.,0.]
   time_array = np.linspace(0,10,1000)
   eta = np.zeros(len(time_array),)
   for i,t in enumerate(time_array):
   	eta[i] = wave.eta(x0,t)
   import matplotlib.pyplot as plt
   plt.plot(time_array,eta)
   plt.xlabel("Time (s)")
   plt.ylabel("Free-surface elevation (m)")
   plt.savefig("Free-surface.pdf")
   plt.show()
