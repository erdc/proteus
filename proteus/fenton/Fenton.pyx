"""
Scripts for creating Fenton waves.
Modified from johndfenton.com/Steady-waves/Fourier.html to work with python and Proteus
Used in proteus.WaveTools
"""
import os
from subprocess import check_call
from shutil import copy
import numpy as np
cimport numpy as np
from proteus import Profiling

cdef extern from "Fourier.cpp":
    cdef void runfourier()

def writeInput(waveheight,
               depth,
               period=None,
               wavelength=None,
               mode='Period',
               current_criterion=1,
               current_magnitude=0,
               ncoeffs=8,
               height_steps=1,
               g=9.81,
               niter=40,
               conv_crit=1.e-05,
               points_freesurface=50,
               points_velocity=16,
               points_vertical=20):
    '''
    Creates input files for Fourier script

    Parameters
    ----------
    waveheight: double
        Height of wave
    depth: double
        Water depth
    period: double
        Wave period
    mode: string
        'Period' or 'Wavelength'
    current_criterion: int
        1: Euler, 2: Stokes
    current_magnitude: double
        Magnitude of current
    ncoeffs: int
        Number of Fourier coefficients
    height_steps: int
        Number of height steps to reach H/d
    g: double
        Gravity
    niter: int
        Max number of iterations 
    conv_crit: double
        Criterion for convergence
    points_freesurface: int
        Number of points on free surface
    points_velocity: int
        Number of velocity/acceleration profiles to print out
    points_vertical: int
        Number of vertical points in each profile
    '''
    # Data input file
    assert period is not None or wavelength is not None, 'Period or wavelength must be set for Fenton wave'
    if period is None and wavelength is not None:
        mode = 'Wavelength'
    if period is not None and wavelength is None:
        mode = 'Period'
    filename = 'Data.dat'
    with open(filename, 'w') as f:
        waveheight_dimless = waveheight/depth
        mode = mode
        if mode == 'Period':
            length_dimless = period*np.sqrt(g/depth)
        elif mode == 'Wavelength':
            length_dimless = wavelength/depth
        current_magnitude_dimless = current_magnitude/np.sqrt(g*depth)
        f.write('''Wave
{waveheight}
{mode}
{length}
{current_criterion}
{current_magnitude}
{ncoeffs}
{height_steps}
'''.format(waveheight=waveheight_dimless, mode=mode, length=length_dimless, current_criterion=current_criterion,
current_magnitude=current_magnitude_dimless, ncoeffs=ncoeffs, height_steps=height_steps))
    # Convergence options file
    filename = 'Convergence.dat'
    with open(filename, 'w') as f:
        f.write('''Control file to control convergence and output of results
{niter}		Maximum number of iterations for each height step; 10 OK for ordinary waves, 40 for highest
{conv_crit}	Criterion for convergence, typically 1.e-4, or 1.e-5 for highest waves
'''.format(niter=niter, conv_crit=conv_crit))
    # Points number file
    filename = 'Points.dat'
    with open(filename, 'w') as f:
        f.write('''Control file to control convergence and output of results
{niter}		Maximum number of iterations for each height step; 10 OK for ordinary waves, 40 for highest
{conv_crit}	Criterion for convergence, typically 1.e-4, or 1.e-5 for highest waves
'''.format(niter=niter, conv_crit=conv_crit))

def runFourier():
    '''
    Runs Fourier.cpp script to get Fenton wave solution
    (!) must be called after writeInput(...)
    '''
    runfourier()

def getBYCoeffs():
    '''
    Get B and Y coeffs of solution
    (!) must be called after runFourier()

    Returns
    -------
    BCoeffs: array_like
        B coeffs of solution
    YCoeffs: array_like
        Y coeffs of solution
    '''
    FFT = np.genfromtxt('Solution.res', delimiter=None, comments='#')
    BCoeffs = FFT[:,1]
    YCoeffs = FFT[:,2]
    return BCoeffs, YCoeffs

def getWavelength():
    '''
    Get wavelength/depth (dimensionless)
    (!) must be called after runFourier()

    Parameters
    ----------
    depth: double
        water depth used to compute solution

    Returns
    -------
    wavelength: double
        wavelength
    '''
    wavelength = None
    with open('Solution.res', 'r') as f:
        for line in f:
            l = line.split()
            if 'Wave' in l and 'length' in l:
                wavelength = l[5]
    return float(wavelength)


def copyFiles():
    '''
    Copy the files to log directory
    '''
    Profiling.logDir = './log'
    copy('Data.dat',Profiling.logDir)
    copy('Convergence.dat',Profiling.logDir)
    copy('Points.dat',Profiling.logDir)
    copy('Solution.res',Profiling.logDir)
    copy('Surface.res',Profiling.logDir)
    copy('Flowfield.res',Profiling.logDir)

def __get_dir():
    return os.path.dirname(__file__)
