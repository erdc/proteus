# -*- coding: cp1252 -*-
# A type of -*- python -*- file
#cython: embedsignature=True
"""Tools for working with water waves.

The primary objective of this module is to provide solutions (exact and
approximate) for the free surface deformation and subsurface velocity
components of water waves. These can be used as boundary conditions, wave
generation sources, and validation solutions for numerical wave codes.

.. inheritance-diagram:: ShockCapturing
    :parts: 2
"""
from math import pi, tanh, sqrt, exp, log, sin, cos, cosh, sinh
import numpy as np
import cmath as cmath
from Profiling import logEvent
import time as tt
import sys as sys


def loadExistingFunction(funcName, validFunctions):
    """ Checks if a function name is known function and returns it

    Checks if a function name is present in a list of functions.
    If True, the function is returned. If False, raises SystemExit. 
    
    Parameters
    ----------    
    funcName : string
            Function name 
    validFunctions : List[function]
            List of valid functions (list of objects) 

    Returns
    --------
    function

    Raises
    ---------
    SystemExit
    

    """
    funcNames = []
    for func  in validFunctions:
            funcNames.append(func.__name__)
            if func.__name__ == funcName:
                func_ret = func
    if funcName not in funcNames:
        logEvent("ERROR! Wavetools.py: Wrong function type (%s) given: Valid wavetypes are %s" %(funcName,funcNames), level=0)
        sys.exit(1)
    return func_ret
       


def setVertDir(g):
    """ Returns the unit vector for the vertical direction 

    The vertical direction is opposite to the gravity direction

    Parameters
    ----------    
    g : numpy.ndarray 
        Gravitational acceleration vector (must have 3 components)

    Returns
    --------
    numpy.ndarray
        
    """
    return -np.array(g/(sqrt(g[0]**2 + g[1]**2 + g[2]**2)))


def setDirVector(vector):
    """ Returns the direction of a vector in the form of a unit vector

    Parameters
    ----------    
    vector : numpy.ndarray
           1D numpy array with three components

    Returns
    --------
    numpy.ndarray
   
    """
    return vector/(sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2))

def dirCheck(v1, v2):
    """ Checks if two vectors are vertical raises SystemError if True

    Parameters
    ----------    
    v1 : numpy.ndarray
        1st vector with three components
           
    v2 : numpy.ndarray 
        2nd vector with three components
        
    Returns
    --------
    None

    Raises
    ---------
    SystemExit
    
    """
    dircheck = abs(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
        #print self.dircheck
    if dircheck > 1e-10:
        logEvent("Wave direction is not perpendicular to gravity vector. Check input",level=0)
        return sys.exit(1)
    else:
        return None
def reduceToIntervals(fi,df):
    """ Prepares the x-axis array with N elements for numerical integration

    Integration is performed along he x- axis. If
    fi = [a1, a2, a3,...,a_N-1 a_N] then it returns the array 
    [a1, 0.5(a1+a2), 0.5(a2+a3),...0.5(a_N-1+a_N), a_N].
    Input array must have constant step

    Parameters
    ----------    
    fi : numpy.ndarray
        x- array with N elements
    df : float
        Constant step of array
    Returns
    --------
    numpy.ndarray
            
    """
    fim_tmp = (0.5*(fi[1:]+fi[:-1])).tolist()
    return np.array([fim_tmp[0]-0.5*df]+fim_tmp+[fim_tmp[-1]+0.5*df])
def returnRectangles(a,x):
    """ Returns 2D discrete integral array using the rectangle method

    The calculation for each array element is
    :math:`(\Delta y_i = 0.5(a_{n-1}+a_{n})*(x_{n-1}-x_{n})`

    Parameters
    ----------    
    a : numpy.ndarray
        Description: Array of y(x) function with N+1 elements
    x : numpy.ndarray
        Description: x- coordinate array with N elements

    Returns
    --------
    numpy.ndarray
            
    
    """
    return 0.5*(a[1:]+a[:-1])*(x[1:]-x[:-1])
def returnRectangles3D(a,x,y):
    """ Returns 3D discrete integrals using the rectangle method

    The calculation for each array element is
    :math: `(\Delta y = 0.25*(a_(n-1,m-1)+a_(n,m-1)+a_(n-1,m) ...
    ...+a_(n,m))*(x_n-1-x_n) *(z_m-1-z_m))`
    
    Parameters
    ----------    
    a : numpy.ndarray
        2D Array of y(x,y) function with (N+1)x(M+1)elements
            Type: Numpy array or list
    x : numpy.ndarray
        Description: x- coordinate array with N+1 elements
    y : numpy.ndarray
            Description: x- coordinate array with N+1 elements
            
    Returns
    --------
    numpy.ndarray
    """
    ai = 0.5*(a[1:,:]+a[:-1,:])
    ai = 0.5*(ai[:,1:]+ai[:,:-1])
    for ii in range(len(x)-1):
        ai[ii,:] *= (y[1:]-y[:-1]) 
    for jj in range(len(y) - 1):
        ai[:,jj] *= (x[1:]-x[:-1])    
    return ai
def normIntegral(f,dom):
    """Returns a normalised 2D function

    The calculation is  :math: `(\int_\Omega f d\Omega =1)`

    Parameters
    ----------    
    f : numpy.ndarray
        Discrete 2D function
        Numpy array or list
    
    dom : float
        Discrete function step

    Returns
    --------
    numpy.ndarray
    """
    G0 = 1./sum(returnRectangles(f,dom))
    return G0*f



def eta_mode(x, t, kDir, omega, phi, amplitude):
    """Calculates the free surface elevation for a single frequency mode

    Parameters
    ----------    
    x : numpy.ndarray
        Position vector
    t : float
        Time variable
    kDir : numpy.ndarray
        Wave number vector
    omega : float
        Angular frequency
    phi : float
        Description: Wave phase 
    amp : float
        Description: Wave amplitude
        
    Returns
    --------
    float
        The free surface elevation at x,t

    """
    phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi
    return amplitude*cos(phase)


def vel_mode(x, t, kDir, kAbs, omega, phi, amplitude, mwl, depth, g, vDir):
    """Calculates the wave velocity components for a single frequency mode

    Parameters
    ----------    
    x : numpy.ndarray
        Position vector
    t : float
        Time variable
    kDir : numpy.ndarray
        Wave number vector
    kAbs : floatkAbs
        Wave number magnitude
    omega : float
        Angular frequency
    phi : float
        Description: Wave phase 
    amplidute : float
        Description: Wave amplitude
    mwl : float
        Mean water level
    depth : float
        Water depth
    g : numpy.ndarray
        Gravitational accelaration (dummy parameter, to be removed)
    vDir : numpy.ndarray
        Unit vector aligned with vertical direction
        

    Returns
    --------
    numpy.ndarray
        1D Numpy array of the velocity vector at x,t
    """

    phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi
    Z =  (vDir[0]*x[0] + vDir[1]*x[1]+ vDir[2]*x[2]) - mwl
    UH = 0.
    UV=0.
    ii=0.
    UH=amplitude*omega*cosh(kAbs*(Z + depth))*cos( phase )/sinh(kAbs*depth)
    UV=amplitude*omega*sinh(kAbs*(Z + depth))*sin( phase )/sinh(kAbs*depth)
    waveDir = kDir/kAbs
#waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.waveDir,wavelength=self.wi[ii], phi0 = self.phi[ii]).u(x,y,z,t)
    V = np.array([UH*waveDir[0]+UV*vDir[0],
                  UH*waveDir[1]+UV*vDir[1],
                  UH*waveDir[2]+UV*vDir[2]])
    return V



def sigma(omega,omega0):
    """Calculates sigma function for JONSWAP spectrum

       See http://www.wikiwaves.org/Ocean-Wave_Sectra
    Parameters
    ----------    
    omega : numpy.ndarray
        Angular frequency array
    omega0 : numpy.ndarray
        Peak angular frequency
    Returns
    --------
    numpy.ndarray
        1D Numpy array of simga function with respect to f
  
    """
    sigmaReturn = np.where(omega > omega0,0.09,0.07)
    return sigmaReturn


def JONSWAP(f,f0,Hs,gamma=3.3,TMA=False, depth = None):
    """Calculates the JONSWAP frequency spectrum (Goda 2009)

    The calculation includes the TMA modification, if TMA =True
    
    Parameters
    ----------    
    f : numpy.ndarray
        Frequency array
    f0 : float
        Peak frequency
    Hs : float
        Significant wave height
    gamma : Optional[float]
        Peak enhancement factor
    TMA : bool
            Description: TMA switch
        Type: boolean
    depth : Optional[float]
        Water depth

    Returns
    --------
    numpy.ndarray
        1D Numpy array of the spectrum in frequency domain

    """
    Tp = 1./f0
    bj = 0.0624*(1.094-0.01915*log(gamma))/(0.23+0.0336*gamma-0.185/(1.9+gamma))
    r = np.exp(-(Tp*f-1.)**2/(2.*sigma(f,f0)**2))
    tma = 1.
    if TMA:
        if (depth == None):
            logEvent("Wavetools:py. Provide valid depth definition definition for TMA spectrum")
            logEvent("Wavetools:py. Stopping simulation")
            sys.exit(1)
        k = dispersion(2*pi*f,depth)
        tma = np.tanh(k*depth)*np.tanh(k*depth)/(1.+ 2.*k*depth/np.sinh(2.*k*depth))

    return tma * bj*(Hs**2)*(1./((Tp**4) *(f**5)))*np.exp(-1.25*(1./(Tp*f)**(4.)))*(gamma**r)

def PM_mod(f,f0,Hs):
    """Calculates the Pierson-Moskovitz spectrum (or Bretschneider or ISSC)

    Reference:
    http://www.orcina.com/SoftwareProducts/OrcaFlex/Documentation/Help/Content/html/Waves,WaveSpectra.htm
    And then to Tucker M J, 1991. Waves in Ocean Engineering. Ellis Horwood Ltd. (Chichester).
    
    f : numpy.ndarray
        Frequency array
    f0 : float
        Peak frequency
    Hs : float
        Significant wave height

    Returns
    --------
    numpy.ndarray
        1D Numpy array of the spectrum in frequency domain

    """
    return (5.0/16.0)*Hs**2*(f0**4/f**5)*np.exp((-5.0/4.0)*(f0/f)**4)

def cos2s(theta,f,s=10):
    """Calculates the cos-2s directional spreading function
    see USACE - CETN-I-28 http://chl.erdc.usace.army.mil/library/publications/chetn/pdf/cetn-i-28.pdf

    Parameters
    ----------    
    theta : numpy.ndarray
        Wave angle array
    f : numpy.ndarray
        Frequency array
    s : Optional[float]
        Spreading parameter

    Returns
    --------
    numpy.ndarray
        2D Numpy array of cos-2s spectrum
    """
    fun = np.zeros((len(theta),len(f)),)
    for ii in range(len(fun[0,:])):
        fun[:,ii] = np.cos(theta/2)**(2*s)
    return fun
def mitsuyasu(theta,fi,f0,smax=10):
    """The cos2s wave directional spread with wave frequency dependency

    Equation from "Random Seas and Design of Maritime Structures" - Y. Goda - 2010
    (3rd ed) eq. 2.22 - 2.25

    Parameters
    ----------    
    theta : numpy.ndarray
        Wave angle array
    fi : numpy.ndarray
        Frequency array
    f0 : float
        Peak frequency
    smax : Optional[float]
        Spreading parameter
    Returns
    --------
    numpy.ndarray
        2D Numpy array of Mitsuyashu-type spectrum
    """

    s = smax * (fi/f0)**(5)
    ii = np.where(fi>f0)[0][0]
    s[ii:] = smax * (fi[ii:]/f0)**(-2.5)
    fun = np.zeros((len(theta),len(fi)),)
    for ii in range(len(fun[0,:])):
        fun[:,ii] = np.cos(theta/2)**(2.*s[ii])
    return fun





def dispersion(w,d, g = 9.81,niter = 1000):
    """Calculates the wave number for single or multiple frequencies using linear dispersion relation. 

    Parameters
    ----------    
    w : float or np.ndarray
        Angular frequency 
    d : float
        Water depth
    g : Optional[float]
        Gravitational acceleration
    niter : Optional[int]
        Solution iterations

    Returns
    --------
    float or numpy.ndarray
        Wavenumber as a float or 1D array for multiple frequencies
    """
    w_aux = np.array(w)
    K = w_aux**2/g
    for jj in range(niter):
        K =  w_aux**2/(g*np.tanh(K*d))
    if type(K) is float:
        return K[0]
    else:
        return K


def tophat(l,cutoff):
    """ Calculates and returns a top hat filter array

    Parameters
    ----------    
    l : int
        Length of array
    cutoff : float
         Cut off fraction at both the leading and tailing part of the array

    Returns
    --------
    numpy.ndarray

    """
    a = np.zeros(l,)
    cut = int(cutoff*l)
    a[cut:-cut] = 1.
    return a

def costap(l,cutoff=0.1):
    """ Calculates and returns a top hat filter array

    Parameters
    ----------    
    l : int
        Length of array
    cutoff : float
         Cut off fraction at both the leading and tailing part of the array

    Returns
    --------
    numpy.ndarray
    """
    npoints = int(cutoff*l)
    wind = np.ones(l)
    for k in range(l): # (k,np) = (n,N) normally used
        if k < npoints:
            wind[k] = 0.5*(1.-cos(pi*float(k)/float(npoints)))
        if k > l - npoints -1:
            wind[k] = 0.5*(1.-cos(pi*float(l-k-1)/float(npoints)))
    return wind

def decompose_tseries(time,eta,dt):
    """ Performs spectral analysis and calculates angular frequency components, amplitude, phase and mean level power
    of a time series with constant sampling.

    Parameters
    ----------    
    time : numpy.ndarray
        Time array
    eta :numpy.ndarray 
        Signal array
    dt : float
        Sampling interval
        
    Returns
    --------
    List
        A list with results with four entries:
         0 -> numpy array with angular frequency components 
         1 -> numpy array with amplitude of each component aa
         2 -> numpy array with phase of each component pp
         3 -> float of the 0th fourier mode (wave setup) 

         """
    nfft = len(time) 
    results = []
    fft_x = np.fft.fft(eta,nfft)
    freq = np.fft.fftfreq(nfft,dt)                              #%complex spectrum
    iend = np.where(freq<0)[0][0]
    setup = np.real(fft_x[0])/nfft
    fft_x = fft_x[1:iend]
    freq = freq[1:iend]
                              #%retaining only first half of the spectrum
    aa = 2.*abs(fft_x)/nfft                                 #%amplitudes (only the ones related to positive frequencies)
    ww = 2*pi*freq
    

    pp = np.zeros(len(aa),complex)
    for k in range(len(aa)):
        pp[k]=cmath.phase(fft_x[k])                       #% Calculating phases phases
    pp = np.real(pp)                                         # Append results to list
    results.append(ww)
    results.append(aa)
    results.append(pp)
    results.append(setup)
    return results

    



class MonochromaticWaves:
    """
    This class is used for generating regular waves in both linear and nonlinear regimes
    
    Parameters
    ----------    
    param : period
            Description: Regular wave period
            Type: Float
    param : waveHeight
            Description: Regular wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : g
            Description: Gravitational acceleration vector 
            Type: Numpy array
    param : wavelength
            Description: Regular wave lenght, calculated from linear dispersion if set to None
            Type: float
    param : waveType
            Description: Set to "Linear" for linear wave theory and "Fenton" for using Fenton Fourier appoximation
            Type: string
    param : Ycoeff
            Description: Fenton Fourier coefficients for free-surface elevation (set to None for linear wave theory)
            Type: Numpy array
    param : Bcoeff
            Description: Fenton Fourier coefficients for velocity (set to None for linear wave theory)
            Type: Numpy array
    param : meanVelocity
            Description: Mean velocity for Fenton Fourier approximation
            Type: Numpy array
    param : phi0
            Description: Regular wave phase (0 by default)
            Type: float

    """
    def __init__(self,
                 period,
                 waveHeight,
                 mwl,
                 depth,
                 g,
                 waveDir,
                 wavelength=None,
                 waveType="Linear",
                 Ycoeff = None, 
                 Bcoeff =None, meanVelocity = np.array([0.,0,0.]),
                 phi0 = 0.):

        self.knownWaveTypes = ["Linear","Fenton"]
        self.waveType = waveType
        if self.waveType not in self.knownWaveTypes:
            logEvent("Wrong wavetype given: Valid wavetypes are %s" %(self.knownWaveTypes), level=0)
            sys.exit(1)
        self.g = np.array(g)
        self.waveDir =  setDirVector(np.array(waveDir))
        self.vDir = setVertDir(g)
        self.gAbs = sqrt(self.g[0]*self.g[0]+self.g[1]*self.g[1]+self.g[2]*self.g[2])

#Checking if g and waveDir are perpendicular
        dirCheck(self.waveDir,self.vDir)
        self.phi0=phi0
        self.period = period
        self.waveHeight = waveHeight
        self.mwl = mwl
        self.depth = depth
        self.omega = 2.0*pi/period

#Calculating / checking wavelength data
        if  self.waveType== "Linear":
            self.k = dispersion(w=self.omega,d=self.depth,g=self.gAbs)
            self.wavelength = 2.0*pi/self.k
        else:
            try:
                self.k = 2.0*pi/wavelength
                self.wavelength=wavelength
            except:
                logEvent("ERROR! Wavetools.py: Wavelenght is not defined for nonlinear waves. Enter wavelength in class arguments",level=0)
                sys.exit(1)
        self.kDir = self.k * self.waveDir
        self.amplitude = 0.5*self.waveHeight
        self.meanVelocity = np.array(meanVelocity)
#Checking that meanvelocity is a vector

        if(len(meanVelocity) != 3):
            logEvent("ERROR! Wavetools.py: meanVelocity should be a vector with 3 components. ",level=0)
            sys.exit(1)

        self.Ycoeff = Ycoeff
        self.Bcoeff = Bcoeff

# Checking for
        if (Ycoeff==None) or (Bcoeff==None):
            if self.waveType!= "Linear":
                logEvent("ERROR! Wavetools.py: Need to define Fenton Fourier coefficients Ycoeff and Bcoeff (free-surface and velocity) for nonlinear waves",level=0)
                sys.exit(1)
    def eta(self, x, t):
        """Calculates free surface elevation (MonochromaticWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """

        if self.waveType == "Linear":
            return eta_mode(x,t,self.kDir,self.omega,self.phi0,self.amplitude)
        elif self.waveType == "Fenton":
            HH = 0.
            ii =0.
            for Y in self.Ycoeff:
                ii+=1
                HH+=eta_mode(x,t,ii*self.kDir,ii*self.omega,ii*self.phi0,Y)
            return HH/self.k

    def u(self, x, t):
        """Calculates wave velocity vector (MonochromaticWaves class).
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """
    
        if self.waveType == "Linear":
            return vel_mode(x, t, self.kDir,self.k,self.omega,self.phi0,self.amplitude,self.mwl,self.depth,self.g,self.vDir)
        elif self.waveType == "Fenton":
            Ufenton = self.meanVelocity.copy()
            ii = 0
            for B in self.Bcoeff:
                ii+=1
                wmode = ii*self.omega
                kmode = ii*self.k
                kdir = self.waveDir*kmode
                amp = tanh(kmode*self.depth)*sqrt(self.gAbs/self.k)*B/self.omega
                Ufenton+= vel_mode(x,t,kdir,kmode,wmode,ii*self.phi0,amp,self.mwl,self.depth,self.g,self.vDir)
            return Ufenton # + self.meanVelocity[comp]


class RandomWaves:
    """
    This class is used for generating plane random waves using linear reconstruction of components from a
    wave spectrum
    
    Parameters
    ----------    
    param : Tp
            Description: Peak wave period
            Type: Float
    param : Hs
            Description: Significant wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir
            Description: Wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : N
            Description: Number of frequency components
            Type: integer
    param : bandfactor
            Description: Spectral band factor. fmax = bandFactor/Tp, fmin = 1/(bandFactor*Tp)
            Type: float
    param : spectName
            Description: Name of spectral distribution
            Type: string
    param : spectral_params
            Description: Dictionary of arguments specific to the spectral distribution 
            Example for JONSWAP = {"gamma": 3.3, "TMA":True,"depth": depth}
            TMA=True activates the TMA modification, which in turn needs the depth as a parameter
            Type: Dictionary 
    param : phi
            Description: Component phases (if set to None, phases are picked at random)
            Type: numpy array

    """
    def __init__(self,
                 Tp,
                 Hs,
                 mwl,#m significant wave height
                 depth ,           #m depth
                 waveDir,
                 g,      #peak  frequency
                 N,
                 bandFactor,         #accelerationof gravity
                 spectName ,# random words will result in error and return the available spectra 
                 spectral_params =  None, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 phi=None
                 ):
        validSpectra = [JONSWAP,PM_mod]
        spec_fun =loadExistingFunction(spectName, validSpectra)                 
        self.g = np.array(g)
        self.waveDir =  setDirVector(np.array(waveDir))
        self.vDir = setVertDir(g)
        dirCheck(self.waveDir,self.vDir)
        self.gAbs = sqrt(self.g[0]*self.g[0]+self.g[1]*self.g[1]+self.g[2]*self.g[2])
        self.Hs = Hs
        self.depth = depth
        self.Tp = Tp
        self.fp = 1./Tp
        self.bandFactor = bandFactor
        self.N = N
        self.mwl = mwl
        self.fmax = self.bandFactor*self.fp
        self.fmin = self.fp/self.bandFactor
        self.df = (self.fmax-self.fmin)/float(self.N-1)
        self.fi = np.linspace(self.fmin,self.fmax,self.N)
        self.omega = 2.*pi*self.fi
        self.ki = dispersion(self.omega,self.depth,g=self.gAbs)
        if phi == None:
            self.phi = 2.0*pi*np.random.random(self.fi.shape[0])
            logEvent('ERROR! Wavetools.py: No phase array is given. Assigning random phases. Outputing the phasing of the random waves')
        else:
            try: 
                self.phi = np.array(phi)
                if self.phi.shape[0] != self.fi.shape[0]:
                    logEvent('ERROR! Wavetools.py: Phase array must have N elements')
                    sys.exit(1)
                    
            except:
                logEvent('ERROR! Wavetools.py: phi argument must be an array with N elements')
                sys.exit(1)

        #ai = np.sqrt((Si_J[1:]+Si_J[:-1])*(fi[1:]-fi[:-1]))
        self.fim = reduceToIntervals(self.fi,self.df)
        if (spectral_params == None):
            self.Si_Jm = spec_fun(self.fim,self.fp,self.Hs)
        else:
            try:
                self.Si_Jm = spec_fun(self.fim,self.fp,self.Hs,**spectral_params)
            except:
                logEvent('ERROR! Wavetools.py: Additional spectral parameters are not valid for the %s spectrum' %spectName)
                sys.exit(1)
        

        self.ai = np.sqrt(2.*returnRectangles(self.Si_Jm,self.fim))
        self.kDir = np.zeros((len(self.ki),3),)
        for ii in range(3):
             self.kDir[:,ii] = self.ki[:] * self.waveDir[ii] 
    def eta(self, x, t):
        """Calculates free surface elevation (RandomWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        Eta=0.
        for ii in range(self.N):
            Eta+= eta_mode(x, t,self.kDir[ii],self.omega[ii],self.phi[ii],self.ai[ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """Calculates wave velocity vector (RandomWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """

        U=0.
        for ii in range(self.N):
            U+= vel_mode(x, t, self.kDir[ii], self.ki[ii],self.omega[ii],self.phi[ii],self.ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U
    def writeEtaSeries(self,Tstart,Tend,x0,fname,Vgen= np.array([0.,0,0])):
        """Writes a timeseries of the free-surface elevation

        It also returns the free surface elevation as a time-eta array.
        If Vgen !=[0.,0.,0.,] then Tstart is modified to account for the
        wave transformation at the most remote point of the relaxation zone.
    
        Parameters
        ----------    
        Tstart : float
            Start time
        Tend : float
            End time
        x0 : numpy.ndarray
            Position vector of the time series
        fname : string
            Filename for timeseries file
        Vgen : Optional[numpy.ndarray]
            Length vector of relaxation zone 
            

        Returns
        ----------
        numpy.ndarray
            2D numpy array Nx2 containing free-surface elevation in time. 
        """
        if sum(Vgen[:]*self.waveDir[:])< 0 :
                logEvent('ERROR! Wavetools.py: Location vector of generation zone should not be opposite to the wave direction')
                sys.exit(1)        
        dt = self.Tp/50.
        Tlag = np.zeros(len(self.omega),)
        for j in range(len(self.omega)):
            Tlag[j] = sum(self.kDir[j,:]*Vgen[:])/self.omega[j]
        Tlag = max(Tlag)
        Tstart = Tstart - Tlag
        Np = int((Tend - Tstart)/dt)
        time = np.linspace(Tstart,Tend,Np )
        etaR  = np.zeros(len(time), )
        for jj in range(len(time)):
            etaR[jj] = self.eta(x0,time[jj])
        np.savetxt(fname,zip(time,etaR))
        series = np.zeros((len(time),2),)
        series[:,0] = time
        series[:,1] = etaR
        return series



class MultiSpectraRandomWaves(RandomWaves):
    """This class is used for generating random waves
    by combining multiple spectra with different distributions and directions
    
    Parameters
    ----------
    param : Nspectra
            Description: Total number of spectra
            Type: integer
    param : Tp
            Description: List of peak wave periods
            Type: list
    param : Hs
            Description: List of significant wave heights
            Type: list
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir
            Description: Lsit of wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: Numpy array
    param : N
            Description: List of numbers of frequency components
            Type: list
    param : bandfactor
            Description: List of spectral band factors
            Type: list
    param : spectName
            Description: List of names of spectral distribution
            Type: list
    param : spectral_params
            Description: List of names of spectral distribution (see RandomWaves class)
            Type: list
    param : phi
            Description: List of component phases 
            Type: list

    """
    def __init__(self,
                 Nspectra,
                 Tp, # np array with 
                 Hs,
                 mwl,#m significant wave height
                 depth ,           #m depth
                 waveDir,
                 g,      #peak  frequency
                 N,
                 bandFactor,         #accelerationof gravity
                 spectName ,# random words will result in error and return the available spectra 
                 spectral_params, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 phi
                 ):
# Checking length of arrays / lists to be equal to NSpectra
        try:
            if (len(Tp) != Nspectra) or (len(Hs) != Nspectra) or (len(waveDir) != Nspectra) or \
               (len(N) != Nspectra) or (len(bandFactor) != Nspectra) or \
               (len(spectName) != Nspectra) or (len(spectral_params) != Nspectra) or(len(phi) != Nspectra):

                logEvent('ERROR! Wavetools.py: Parameters passed in MultiSpectraRandomWaves must be in array or list form with length Nspectra  ')
                sys.exit(1)
               
        except:
            logEvent('ERROR! Wavetools.py: Parameters passed in MultiSpectraRandomWaves must be in array or list form with length Nspectra  ')
            sys.exit(1)
        # Initialize numpy arrays for complete reconstruction
        self.Nall = 0 
        for nn in N:
            self.Nall+=nn
        

        self.omegaM = np.zeros(self.Nall,float)
        self.kiM = np.zeros(self.Nall,float)
        self.aiM = np.zeros(self.Nall,float)
        self.kDirM = np.zeros((self.Nall,3),float)
        self.phiM= np.zeros(self.Nall,float)


        NN = 0
        for kk in range(Nspectra):
            logEvent("ERROR! Wavetools.py: Reading spectra No %s" %kk)
            NN1 = NN
            NN +=N[kk]
            RandomWaves.__init__(self,
                                 Tp[kk], # np array with 
                                 Hs[kk],
                                 mwl,#m significant wave height
                                 depth,           #m depth
                                 waveDir[kk],
                                 g,      #peak  frequency
                                 N[kk],
                                 bandFactor[kk],         #accelerationof gravity
                                 spectName[kk],# random words will result in error and return the available spectra 
                                 spectral_params[kk], #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                                 phi[kk]
                             )
            self.omegaM[NN1:NN] = self.omega
            self.kiM[NN1:NN] = self.ki
            self.aiM[NN1:NN] = self.ai
            self.kDirM[NN1:NN,:] =self.kDir[:,:]
            self.phiM[NN1:NN] = self.phi
        

    def eta(self, x, t):
        """Calculates free surface elevation(MultiSpectraRandomWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """

        Eta=0.
        for ii in range(self.Nall):
            Eta+= eta_mode(x, t, self.kDirM[ii],self.omegaM[ii],self.phiM[ii],self.aiM[ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """Calculates wave velocity vector (MultiSpectraRandomWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """

        U=0.
        for ii in range(self.Nall):
            U+= vel_mode(x,t,self.kDirM[ii], self.kiM[ii],self.omegaM[ii],self.phiM[ii],self.aiM[ii],self.mwl,self.depth,self.g,self.vDir)
        return U



class DirectionalWaves(RandomWaves):
    """
    This class is used for generating directional random waves using linear reconstruction of components from a
    wave spectrum
    
    Parameters
    ----------    
    param : M
            Description: Number of directional components
            Type: integer
    param : Tp
            Description: Peak wave period
            Type: Float
    param : Hs
            Description: Significant wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir0
            Description: Leading wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : N
            Description: Number of frequency components
            Type: integer
    param : bandfactor
            Description: Spectral band factor. fmax = bandFactor/Tp, fmin = 1/(bandFactor*Tp)
            Type: float
    param : spectName
            Description: Name of spectral distribution
            Type: string
    param : spreadName
            Description: Name of spreading distribution
            Type: string
    param : spectral_params
            Description: Dictionary of arguments specific to the spectral distribution (see RandomWaves class)
            Type: Dictionary 
    param : spread_params
            Description: Dictionary of arguments specific to the spreading distribution
            Example for Cos-2s = {"s": 10}
            Example for Mitsuyashu-type = {"fp": 1/Tp, "smax":10}
            Type: Dictionary 
    param : phi
            Description: Component phases (if set to None, phases are picked at random)
            Type: numpy array
    param : phiSymm
            Description: Switch for enabling a symmetric phase allocation across directional components 
            Type: boolean

    """
    def __init__(self,
                 M,  #half bin of frequencies
                 Tp, # np array with 
                 Hs, # 
                 mwl,#m significant wave height
                 depth ,           #m depth
                 waveDir0,  # Lead direction
                 g,      #peak  frequency
                 N,    # Number of frequencies
                 bandFactor,         #accelerationof gravity
                 spectName ,# random words will result in error and return the available spectra 
                 spreadName ,# random words will result in error and return the available spectra 
                 spectral_params = None, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 spread_params = None,  
                 phi=None, # phi must be an (2*M+1)*N numpy array
                 phiSymm = False # When true, phi[-pi/2,0] is symmetric to phi[0,pi/2]
                 ):   
        validSpread = [cos2s,mitsuyasu]
        spread_fun =  loadExistingFunction(spreadName, validSpread)
        self.M = M
        self.Mtot = 2*M+1
        self.waveDir0 = setDirVector(waveDir0)
        self.vDir = setVertDir(g) 


 # Loading Random waves to get the frequency array the wavelegnths and the frequency spectrum
        RandomWaves.__init__(self,
                             Tp, # np array with 
                             Hs,
                             mwl,#m significant wave height
                             depth,           #m depth
                             self.waveDir0,
                             g,      #peak  frequency
                             N,
                             bandFactor,         #accelerationof gravity
                             spectName,# random words will result in error and return the available spectra 
                             spectral_params, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                             phi = None 
        )

       
        
        # Directional waves propagate usually in a plane -90 to 90 deg with respect to the direction vector, normal to the gavity direction. Rotating the waveDir0 vector around the g vector to produce the directional space
        from SpatialTools import rotation3D
        self.thetas = np.linspace(-pi/2,pi/2,2*M+1)
        self.dth = (self.thetas[1] - self.thetas[0])
        self.waveDirs = np.zeros((2*M+1,3),)
        self.phiDirs = np.zeros((2*M+1,N),)
        self.aiDirs = np.zeros((2*M+1,N),)
        

        temp_array = np.zeros((1,3),)
        temp_array[0,:] = waveDir0
        directions = range(0,self.Mtot)

# initialising wave directions
        for rr in directions: 
            theta = self.thetas[rr]            
            self.waveDirs[rr,:] = rotation3D(temp_array,theta,self.vDir)[0,:]
            self.waveDirs[rr,:]=setDirVector( self.waveDirs[rr,:])


# Initialising phasing
        if phi == None:
            self.phiDirs = 2.0*pi*np.random.rand(self.Mtot,self.fi.shape[0])
        elif np.shape(phi) == (2*M+1,self.fi.shape[0]):
            self.phiDirs = phi
        else:
            logEvent("ERROR! Wavetools.py: phi in DirectionalWaves class must be given either as None or as a list with 2*M + 1 numpy arrays with length N")
            sys.exit(1)
            
        if (phiSymm):
            for i in range(0,M):
                self.phiDirs[M+1+i,:] = self.phiDirs[self.M - 1 - i,:]
            
            


        self.theta_m = reduceToIntervals(self.thetas,self.dth)        
        if (spread_params == None):
            self.Si_Sp = spread_fun(self.theta_m,self.fim)
        else:
            try:
                self.Si_Sp = spread_fun(self.theta_m,self.fim, **spread_params)
            except:
                logEvent('ERROR! Wavetools.py: Additional spread parameters are not valid for the %s spectrum' %spectName)
                sys.exit(1)

        # Setting amplitudes 
        #Normalising the spreading function
        freq = range(0,self.N)
    # Normalising integral over all frequencies
        for ii in freq:            
            self.Si_Sp[:,ii] = normIntegral(self.Si_Sp[:,ii],self.theta_m)
            self.Si_Sp[:,ii]*= self.Si_Jm[ii] 
    # Creating amplitudes spectrum
        self.aiDirs[:] = np.sqrt(2.*returnRectangles3D(self.Si_Sp,self.theta_m,self.fim))
    def eta(self, x, t):
        """Calculates free surface elevation (DirectionalWaves class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        Eta=0.
        for jj in range(self.Mtot):
            for ii in range(self.N):
                kDiri = self.waveDirs[jj]*self.ki[ii]
                Eta+= eta_mode(x,t,kDiri,self.omega[ii],self.phiDirs[jj,ii],self.aiDirs[jj,ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """Calculates wave velocity vector (Directional class).
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """
        U=0.
        for jj in range(self.Mtot):
            for ii in range(self.N):
                kDiri = self.waveDirs[jj]*self.ki[ii]
                U+= vel_mode(x,t,kDiri, self.ki[ii],self.omega[ii],self.phiDirs[jj,ii],self.aiDirs[jj,ii],self.mwl,self.depth,self.g,self.vDir)
        return U
     


            
                

        

class TimeSeries:
    """This class is used for generating waves from an arbirtrary free-surface elevation time series
    
    Parameters
    ----------    
    param : timeSeriesFile
            Description: Time series file name (csv or txt)
            Type: string
    param : skiprows
            Description: Number of header rows
            Type: float
    param : timeSeriesPosition
            Description: Coordinates of the gauge / signal location
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : N
            Description: Number of frequency components
            Type: integer
    param : mwl
            Description: Still water level
            Type: float
    param : waveDir
            Description: Leading wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : cutoffTotal
            Description: Cut off fraction, applied both at the leading and tailing parts of the series
            Type: float
    param : rec_direct
            Description: Switch for activating direct decomposition
            Type: string
    param : window_params
            Description: Dictionary of parameters for window method
            e.g.  window_params = {"Nwaves":15, "Tm": Tp/1.1, "Window":"costap"} (minimum parameters required)
            e.g.  window_params = {"Nwaves":15, "Tm": Tp/1.1, "Window":"costap", "Overlap":0.5, "Cutoff":0.2} (full range of parameters)
            Type: Dictionary 
    param : arrayData
            Description: Switch for passing the time series as an array (False by default)
            Type: boolean
    param : seriesArray
            Description: Free surface elevation time series given in an array format (None by default) 
            Type: numpy array

    """
    def __init__(self,
                 timeSeriesFile, # e.g.= "Timeseries.txt",
                 skiprows,
                 timeSeriesPosition,
                 depth  ,
                 N ,          #number of frequency bins
                 mwl ,        #mean water level
                 waveDir, 
                 g,
                 cutoffTotal = 0.01,
                 rec_direct = True,
                 window_params = None, #If rec_direct = False then wind_params = {"Nwaves":Nwaves,"Tm":Tm,"Window":wind_filt,"Overlap":overlap,"Cutoff":cutoff}
                 arrayData = False,
                 seriesArray = None
                 ):

        # Setting the depth
        self.depth = depth
        self.rec_direct = rec_direct
        # Number of wave components
        self.N = N
        self.Nwaves = None
        # Position of timeSeriesFile
        if(len(timeSeriesPosition)==3):
            self.x0 = timeSeriesPosition[0]
            self.y0 = timeSeriesPosition[1]
            self.z0 = timeSeriesPosition[2]
        else:
            logEvent("ERROR! Wavetools.py: Location vector for timeSeries must have three-components",level=0)
            sys.exit(1)
            

        # Mean water level
        self.mwl = mwl
        # Wave direction
        self.waveDir = setDirVector(waveDir)
        # Gravity
        self.g = np.array(g)
        # Derived variables
        # Gravity magnitude
        self.gAbs = sqrt(sum(g * g))
        # Definition of gravity direction
        self.vDir = setVertDir(g)
        dirCheck(self.waveDir,self.vDir)
        #Reading time series

        
        self.arrayData = arrayData
        if(self.arrayData):
            tdata = seriesArray
        else:
            filetype = timeSeriesFile[-4:]
            fid = open(timeSeriesFile,"r")
            if (filetype !=".txt") and (filetype != ".csv"):
                logEvent("ERROR! Wavetools.py: File %s must be given in .txt or .csv format" % (timeSeriesFile),level=0)
                sys.exit(1)
            elif (filetype == ".csv"):
                tdata = np.loadtxt(fid,skiprows=skiprows,delimiter=",")
            else:
                tdata = np.loadtxt(fid,skiprows=skiprows)
            fid.close()
        #Checks for tseries file
        # Only 2 columns: time & eta
        ncols = len(tdata[0,:])
        if ncols != 2:
            logEvent("ERROR! Wavetools.py: Timeseries file (%s) must have only two columns [time, eta]" % (timeSeriesFile),level=0)
            sys.exit(1)
        time_temp = tdata[:,0]
        self.dt = (time_temp[-1]-time_temp[0])/(len(time_temp)-1)



        # If necessary, perform interpolation
        doInterp = False
        for i in range(1,len(time_temp)):
            dt_temp = time_temp[i]-time_temp[i-1]
        #check if time is at first column
            if time_temp[i]<=time_temp[i-1]:
                logEvent("ERROR! WaveTools.py:  Found not consistent time entry between %s and %s row in %s file. Time variable must be always at the first column of the file and increasing monotonically" %(i-1,i,timeSeriesFile) )
                sys.exit(1)
        #check if sampling rate is constant
            if dt_temp!=self.dt:
                doInterp = True
        if(doInterp):
            logEvent("ERROR! WaveTools.py: Not constant sampling rate found, proceeding to signal interpolation to a constant sampling rate",level=0)
            self.time = np.linspace(time_temp[0],time_temp[-1],len(time_temp))
            self.eta = np.interp(self.time,time_temp,tdata[:,1])
        else:
            self.time = time_temp
            self.eta = tdata[:,1]

        self.t0  = self.time[0]        
        # Remove mean level from raw data
        self.eta -= np.mean(self.eta)
        # Filter out first 2.5 % and last 2.5% to make the signal periodic
        self.eta *= costap(len(self.time),cutoff=cutoffTotal)
        # clear tdata from memory
        del tdata
        # Calculate time lenght
        self.tlength = (self.time[-1]-self.time[0])
        # Matrix initialisation
        self.windows_handover = []
        self.windows_rec = []

    # Direct decomposition of the time series for using at reconstruct_direct
        if (self.rec_direct):
            Nf = self.N
            self.nfft=len(self.time)
            logEvent("ERROR! WaveTools.py: performing a direct series decomposition")
            self.decomp = decompose_tseries(self.time,self.eta,self.dt)
            self.ai = self.decomp[1]
            ipeak = np.where(self.ai == max(self.ai))[0][0]
            imax = min(ipeak + Nf/2,len(self.ai))
            imin = max(0,ipeak - Nf/2)
            self.ai = self.ai[imin:imax]
            self.omega = self.decomp[0][imin:imax]
            self.phi = - self.decomp[2][imin:imax]
            self.ki = dispersion(self.omega,self.depth,g=self.gAbs)
            self.Nf = imax - imin
            self.setup = self.decomp[3]
            self.kDir = np.zeros((len(self.ki),3),"d")
            for ii in range(len(self.ki)):
                self.kDir[ii,:] = self.ki[ii]*self.waveDir[:]


                # Spectral windowing
        else:
            if (window_params==None):
                logEvent("ERROR! WaveTools.py: Set parameters for spectral windowing. Argument window_params must be a dictionary")
                sys.exit(1)
            try:
                self.Nwaves = window_params["Nwaves"]
            except:
                logEvent("ERROR! WaveTools.py: Dictionary key 'Nwaves' (waves per window) not found in window_params dictionary")
                sys.exit(1)

            try:           
                self.Tm = window_params["Tm"]
            except:
                logEvent("ERROR! WaveTools.py: Dictionary key 'Tm' (mean or characteristic wave period) not found in window_params dictionary")
                sys.exit(1)

            try:           
                self.windowName = window_params["Window"]
            except:
                logEvent("ERROR! WaveTools.py: Dictionary key 'Window' (window function type) not found in window_params dictionary")
                sys.exit(1)

            if(self.Nwaves > 0.5*self.tlength / self.Tm):
                logEvent("ERROR! WaveTools.py: Reconstruction is expected to have two windows or more. Plese reduce the number of waves per window or switch to direct decomposition )")
                sys.exit(1)



            validWindows = [costap, tophat]
            wind_filt =  loadExistingFunction(self.windowName, validWindows) 
            logEvent("ERROR! WaveTools.py: performing series decomposition with spectral windows")
            # Portion of overlap, compared to window time
            try:
                self.overlap = window_params["Overlap"]            
            except:
                self.overlap = 0.25
                logEvent("ERROR! WaveTools.py: Overlap entry in window_params dictionary not found. Setting default value of 0.25 (1/4 of the window length)")

            try:
                self.cutoff = window_params["Cutoff"]            
            except:
                self.cutoff= 0.1
                logEvent("ERROR! WaveTools.py: Cutoff entry in window_params dictionary not found. Setting default value of 0.1 (1/10 of the window length)")
                
                

            # Portion of window filtered with the Costap filter
            # Setting the handover time, either at the middle of the overlap or just after the filter
            self.handover = max(1.1 *self.cutoff,  self.overlap / 2.)
            if (self.handover > 0.9 * self.overlap):
                logEvent("ERROR! Wavetools.py: Window handover is not optimal as the cutoff is too close to the overlap. Decrease cutoff or increase overlap")
                sys.exit(1)
            self.Twindow =  self.Tm * self.Nwaves            # setting the window duration (approx.). Twindow = Tmean * Nwaves
            self.Toverlap = self.overlap * self.Twindow             
            self.Nwindows = int( (self.tlength -   self.Twindow ) / (self.Twindow - self.Toverlap) ) + 1             #Getting the actual number of windows  (N-1) * (Twindow - Toverlap) + Twindow = total time
            self.Twindow = self.tlength/(1. + (1. - self.overlap)*(self.Nwindows-1))            # Correct Twindow and Toverlap for duration and integer number of windows
            self.Toverlap = self.overlap*self.Twindow
            logEvent("ERROR! Wavetools.py: Correcting window duration for matching the exact time range of the series. Window duration correspond to %s waves approx." %(self.Twindow / self.Tm) )
            diff = (self.Nwindows-1.)*(self.Twindow -self.Toverlap)+self.Twindow - self.tlength
            logEvent("ERROR! Wavetools.py: Checking duration of windowed time series: %s per cent difference from original duration" %(100*diff) )
            logEvent("ERROR! Wavetools.py: Using %s windows for reconstruction with %s sec duration and %s per cent overlap" %(self.Nwindows, self.Twindow,100*self.overlap ))
# Setting where each window starts and ends
            for jj in range(self.Nwindows):
                span = np.zeros(2,"d")
                tfirst = self.time[0] + self.Twindow
                tlast = self.time[-1] - self.Twindow
                if jj == 0:
                    ispan1 = 0
                    ispan2 = np.where(self.time> tfirst)[0][0]
                elif jj == self.Nwindows-1:
                    ispan1 = np.where(self.time > tlast)[0][0]
                    ispan2 = len(self.time)-1
                else:
                    tstart = self.time[ispan2] - self.Toverlap
                    ispan1 = np.where(self.time > tstart)[0][0]
                    ispan2 = np.where(self.time > tstart + self.Twindow )[0][0]
                span[0] = ispan1
                span[1] = ispan2
# Storing time series in windows and handover times
                self.windows_handover.append( self.time[ispan2] - self.handover*self.Twindow )
                self.windows_rec.append(np.array(zip(self.time[ispan1:ispan2],self.eta[ispan1:ispan2])))
# Decomposing windows to frequency domain
            self.decompose_window = []
#            style = "k-"
#            ii = 0
            
            for wind in self.windows_rec:
                self.nfft=len(wind[:,0])
                wind[:,1] *=wind_filt(self.nfft,cutoff = self.cutoff)
                decomp = decompose_tseries(wind[:,0],wind[:,1],self.dt)
                self.N = min(self.N, len(decomp[0]))
                Nftemp = self.N
                ipeak =  np.where(decomp[1] == max(decomp[1]))[0][0]
                imax = min(ipeak + Nftemp/2,len(decomp[1]))
                imin = max(0,ipeak - Nftemp/2)
                self.Nf = imax-imin
                if (self.Nf < self.N):
                    if imin == 0:
                        imax = imax + (self.N - self.Nf)
                    else:
                        imin = imin - (self.N - self.Nf)
                    self.Nf = self.N

                decomp[1] = decomp[1][imin:imax]
                decomp[0] = decomp[0][imin:imax]
                decomp[2] = -decomp[2][imin:imax]
                ki = dispersion(decomp[0],self.depth,g=self.gAbs)
                kDir = np.zeros((len(ki),3),"d")
                for ii in range(len(ki)):
                    kDir[ii,:] = ki[ii]*self.waveDir[:]
                decomp.append(kDir)
                decomp.append(ki)

                self.decompose_window.append(decomp)
                
            
        if(self.rec_direct):
            self.eta = self.etaDirect
            self.u = self.uDirect
        else:
            self.eta =  self.etaWindow
            self.u = self.uWindow
     
#                if style == "k-":
#                    style = "kx"
#                else:
#                    style ="k-"
#                plt.plot(wind[:,0],wind[:,1],style)
#                plt.plot(self.time,self.eta,"bo",markersize=2)
#                plt.plot([self.windows_handover[ii],self.windows_handover[ii]] , [-1000,1000],"b--")
#                ii+=1
#            plt.ylim(-1,2)
#            plt.grid()
#            plt.savefig("rec.pdf")
#            self.Twindow = self.Npw*self.dt
#            self.Noverlap = int(self.Npw *0.25)

    def etaDirect(self, x, t):
        """Calculates free surface elevation(Timeseries class-direct method
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        Eta=0.        
        for ii in range(0,self.Nf):
            x1 = np.array(x)-[self.x0, self.y0, self.z0]
            Eta+= eta_mode(x1,t-self.t0,self.kDir[ii],self.omega[ii],self.phi[ii],self.ai[ii])
        return Eta

    def uDirect(self, x, t):
        """Calculates wave velocity vector (Timeseries class-direct method)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """
        U=0.
        for ii in range(0,self.Nf):
            x1 = x-[self.x0, self.y0, self.z0]
            U+= vel_mode(x1, t-self.t0, self.kDir[ii],self.ki[ii], self.omega[ii],self.phi[ii],self.ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U

    def findWindow(self,t):
        """Returns the current spectral window in TimeSeries class."

        Parameters
        ----------

        t : float
                Time variable

        Returns
        -------
        int
            Index of window as an integer
        
        """
        term = 1. - self.handover
        if t-self.time[0] >= term*self.Twindow:
            Nw = min(int((t-self.time[0] - term*self.Twindow)/(self.Twindow - 2. * self.handover * self.Twindow)) + 1, self.Nwindows-1)
            if t-self.time[0] < self.windows_handover[Nw-1] - self.time[0]:
                Nw-=1
        else:
            Nw = 0
        return Nw
        
    def etaWindow(self, x, t):
        """Calculates free surface elevation(Timeseries class-window method
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        Nw = self.findWindow(t)
        ai =  self.decompose_window[Nw][1]
        omega = self.decompose_window[Nw][0]
        phi = self.decompose_window[Nw][2]
        kDir = self.decompose_window[Nw][4]
        t0 = self.windows_rec[Nw][0,0]
        Eta=0.        
        for ii in range(0,self.Nf):
            x1 = np.array(x)-[self.x0, self.y0, self.z0]
            Eta+= eta_mode(x1, t-t0, kDir[ii], omega[ii], phi[ii], ai[ii])
        return Eta

    def uWindow(self, x, t):
        """Calculates wave velocity vector (Timeseries class-window method)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """
        Nw = self.findWindow(t)
        ai =  self.decompose_window[Nw][1]
        omega = self.decompose_window[Nw][0]
        phi = self.decompose_window[Nw][2]
        kDir = self.decompose_window[Nw][4]
        ki = self.decompose_window[Nw][5]
        t0 = self.windows_rec[Nw][0,0]
        U=0.
        for ii in range(0,self.Nf):
            x1 =  np.array(x)-[self.x0, self.y0, self.z0]
            U+= vel_mode(x1, t-t0, kDir[ii],ki[ii],omega[ii],phi[ii],ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U



class RandomWavesFast(RandomWaves):
    """
    This class is used for generating plane random waves in an optimised manner
    using linear reconstruction of components from a wave spectrum
    
    Parameters
    ----------    
    param : Tstart
            Description: Start time
            Type: Float    
    param : Tend
            Description: End time
            Type: Float    
    param : x0
            Description: Position vector for the time series
            Type: Numpy array    
    param : Tp
            Description: Peak wave period
            Type: Float
    param : Hs
            Description: Significant wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir
            Description: Wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : N
            Description: Number of frequency components
            Type: integer
    param : bandfactor
            Description: Spectral band factor. fmax = bandFactor/Tp, fmin = 1/(bandFactor*Tp)
            Type: float
    param : spectName
            Description: Name of spectral distribution
            Type: string
    param : spectral_params
            Description: Dictionary of arguments specific to the spectral distribution 
            Example for JONSWAP = {"gamma": 3.3, "TMA":True,"depth": depth}
            TMA=True activates the TMA modification, which in turn needs the depth as a parameter
            Type: Dictionary 
    param : phi
            Description: Component phases (if set to None, phases are picked at random)
            Type: numpy array
    param : Lgen
            Description: Length of the generation zone (np.array([0., 0., 0.]) by default
            Type: numpy array
    """

    def __init__(self,
                 Tstart,
                 Tend,
                 x0,
                 Tp,
                 Hs,
                 mwl,#m significant wave height
                 depth ,           #m depth
                 waveDir,
                 g,      #peak  frequency
                 N,
                 bandFactor,         #accelerationof gravity
                 spectName ,# random words will result in error and return the available spectra 
                 spectral_params =  None, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 phi=None,
                 Lgen = np.array([0., 0. ,0. ])
                 ):
            RandomWaves.__init__(self,
                                 Tp, # np array with 
                                 Hs,
                                 mwl,#m significant wave height
                                 depth,           #m depth
                                 waveDir,
                                 g,      #peak  frequency
                                 N,
                                 bandFactor,         #accelerationof gravity
                                 spectName,# random words will result in error and return the available spectra 
                                 spectral_params, #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                                 phi
                             )
            fname = "RandomSeries"+"_Hs_"+str(self.Hs)+"_Tp_"+str(self.Tp)+"_depth_"+str(self.depth)
            series = self.writeEtaSeries(Tstart,Tend,x0,fname,Lgen)
            cutoff = 0.2*self.Tp/(series[-1,0]-series[0,0])
            TS = TimeSeries(
                 fname, # e.g.= "Timeseries.txt",
                 0,
                 x0,
                 self.depth ,
                 32 ,          #number of frequency bins
                 self.mwl ,        #mean water level
                 self.waveDir, 
                 self.g,
                 cutoffTotal = cutoff,
                 rec_direct = False,
                 window_params = {"Nwaves":15 ,"Tm":self.Tp/1.1,"Window":"costap"},
                 arrayData = True,
                 seriesArray = series
                 )
            self.eta = TS.eta
            self.u = TS.u





class RandomNLWaves(RandomWaves):
    """
    This class is contains functions for calculating random waves with 2nd order corrections
    
    Parameters
    ----------    
    param : Tstart
            Description: Start time
            Type: Float    
    param : Tend
            Description: End time
            Type: Float    
    param : Tp
            Description: Peak wave period
            Type: Float
    param : Hs
            Description: Significant wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir
            Description: Wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : N
            Description: Number of frequency components
            Type: integer
    param : bandfactor
            Description: Spectral band factor. fmax = bandFactor/Tp, fmin = 1/(bandFactor*Tp)
            Type: float
    param : spectName
            Description: Name of spectral distribution
            Type: string
    param : spectral_params
            Description: Dictionary of arguments specific to the spectral distribution 
            Example for JONSWAP = {"gamma": 3.3, "TMA":True,"depth": depth}
            TMA=True activates the TMA modification, which in turn needs the depth as a parameter
            Type: Dictionary 
    param : phi
            Description: Component phases (if set to None, phases are picked at random)
            Type: numpy array
    """
    def __init__(self,
                 Tstart,
                 Tend,
                 Tp,                      #wave period
                 Hs,                      #significant wave height
                 mwl,                     #mean water level
                 depth,                   #water depth          
                 waveDir,                 #wave direction vector with three components
                 g,                       #gravitational accelaration vector with three components      
                 N,                       #number of frequency bins
                 bandFactor,              #width factor for band around peak frequency fp       
                 spectName,               #random words will result in error and return the available spectra 
                 spectral_params=None,    #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 phi=None                 #array of component phases
                 ):
        RandomWaves.__init__(self,Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params,phi)
        
        
        self.eta_linear = self.eta
        self.eta = self.wtError 
        self.u = self.wtError 
    
    def eta_2ndOrder(self,x,t):
        """Calculates the free surface elevation for 2nd-order terms

        Uses 2nd order random wave theory

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """

        Eta2nd = 0.
        for i in range(0,self.N):
            ai_2nd = (self.ai[i]**2*self.ki[i]*(2+3/sinh(self.ki[i]*self.depth)**2))/(4*tanh(self.ki[i]*self.depth))
            wwi_2ndOrder = eta_mode(x,t,2*self.kDir[i],2*self.omega[i],2*self.phi[i],ai_2nd)
            Eta2nd += wwi_2ndOrder  
        return Eta2nd



    #higher harmonics
    def eta_short(self,x,t):
        """Calculates the free surface elevation for higher-order terms

        Uses 2nd order random wave theory

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """

        Etashort = 0.
        for i in range(0,self.N-1):
            for j in range(i+1,self.N):
                Dp = (self.omega[i]+self.omega[j])**2 - self.gAbs*(self.ki[i]+self.ki[j])*tanh((self.ki[i]+self.ki[j])*self.depth)
                Bp = (self.omega[i]**2+self.omega[j]**2)/(2*self.gAbs) - ((self.omega[i]*self.omega[j])/(2*self.gAbs)) *(1-1./(tanh(self.ki[i]*self.depth)*tanh(self.ki[j]*self.depth))) *(((self.omega[i]+self.omega[j])**2 + self.gAbs*(self.ki[i]+self.ki[j])*tanh((self.ki[i]+self.ki[j])*self.depth))/Dp) + ((self.omega[i]+self.omega[j])/(2*self.gAbs*Dp))*((self.omega[i]**3/sinh(self.ki[i]*self.depth)**2) + (self.omega[j]**3/sinh(self.ki[j]*self.depth)**2))	
                ai_short = self.ai[i]*self.ai[j]*Bp
                wwi_short = eta_mode(x,t,self.kDir[i]+self.kDir[j],self.omega[i]+self.omega[j],self.phi[i]+self.phi[j],ai_short)
                Etashort += wwi_short
        return Etashort



    #lower harmonics
    def eta_long(self,x,t):
        """Calculates the free surface elevation for lower-order terms

        Uses 2nd order random wave theory

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """


        Etalong = 0.
        for i in range(0,self.N-1):
            for j in range(i+1,self.N):
                Dm = (self.omega[i]-self.omega[j])**2 - self.gAbs*(self.ki[i]-self.ki[j])*tanh((self.ki[i]-self.ki[j])*self.depth)	
                Bm = (self.omega[i]**2+self.omega[j]**2)/(2*self.gAbs) + ((self.omega[i]*self.omega[j])/(2*self.gAbs))*(1+1./(tanh(self.ki[i]*self.depth)*tanh(self.ki[j]*self.depth)))*(((self.omega[i]-self.omega[j])**2 + self.gAbs*(self.ki[i]-self.ki[j])*tanh((self.ki[i]-self.ki[j])*self.depth))/Dm) + ((self.omega[i]-self.omega[j])/(2*self.gAbs*Dm))*((self.omega[i]**3/sinh(self.ki[i]*self.depth)**2) - (self.omega[j]**3/sinh(self.ki[j]*self.depth)**2))
                ai_long = self.ai[i]*self.ai[j]*Bm 
                wwi_long = eta_mode(x,t,self.kDir[i]-self.kDir[j],self.omega[i]-self.omega[j],self.phi[i]-self.phi[j],ai_long)
                Etalong += wwi_long
        return Etalong


    #set-up calculation
    def eta_setUp(self,x,t):
        """Calculates the free surface elevation set up

        Uses 2nd order random wave theory

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """


        EtasetUp = 0.
        for i in range(0,self.N):
            wwi_setUp = (self.ai[i]**2*self.ki[i])/(2*sinh(2*self.ki[i]*self.depth))
            EtasetUp += wwi_setUp
        return EtasetUp



    #overall free surface elevation
    def eta_overall(self,x,t,setUp=False):
        """Calculates the free surface elevation with 2nd order corrections

        Uses 2nd order random wave theory

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        Etaoverall =  self.eta_linear(x,t) + self.eta_2ndOrder(x,t) + self.eta_short(x,t) + self.eta_long(x,t)
        if setUp:   
            Etaoverall -= self.eta_setUp(x,t)
        return Etaoverall



    def writeEtaSeries(self,Tstart,Tend,dt,x0,fname, mode="all",setUp=False,Vgen=np.array([0.,0.,0.])):
        """Writes a timeseries of the free-surface elevation

        It also returns the free surface elevation as a time-eta array.
        If Vgen !=[0.,0.,0.,] then Tstart is modified to account for the
        wave transformation at the most remote point of the relaxation zone.
    
        Parameters
        ----------    
        Tstart : float
            Start time
        Tend : float
            End time
        dt : float
            Sampling interval
        x0 : numpy.ndarray
            Position vector of the time series
        fname : string
            Filename for timeseries file
        mode: Optional[string]
            Mode of set up calculations (all, long, short, setup)
        setUp: Optional[bool]
            Switch for activating setup calculation
        Vgen : Optional[numpy.ndarray]
            Length vector of relaxation zone 
            

        Returns
        ----------
        numpy.ndarray
            2D numpy array Nx2 containing free-surface elevation in time. 
        """
        if sum(Vgen[:]*self.waveDir[:])< 0 :
            logEvent('ERROR! Wavetools.py: Location vector of generation zone should not be opposite to the wave direction')
            sys.exit(1)        

        Tlag = np.zeros(len(self.omega),)
        for j in range(len(self.omega)):
            Tlag[j] = sum(self.kDir[j,:]*Vgen[:])/self.omega[j]
        Tlag = max(Tlag)
        Tstart = Tstart - Tlag

        Nseries = int(Tend - Tstart)/dt + 1
        timelst=np.linspace(Tstart, Tend, Nseries)
        series = np.zeros((Nseries,2),)
        series[:,0] = timelst
        for i in range(len(timelst)):
            time = series[i,0]
            if mode == "all":
                series[i,1] = self.eta_overall(x0,time,setUp)
            elif mode == "setup":
                series[i,1] = self.eta_setUp(x0,time)
            elif mode == "short":
                series[i,1] = self.eta_short(x0,time) + self.eta_2ndOrder(x0,time)
            elif mode == "long":
                series[i,1] = self.eta_long(x0,time) 
            elif mode == "linear":
                series[i,1] = self.eta_linear(x0,time)
            else:
                logEvent('ERROR! Wavetools.pyx: Argument mode in RandomNLWaves.writeEtaSeries should be "all", "setup", "short", "long" or "linear"')
                sys.exit(1)        
        delimiter =" "
        if fname[-4:]==".csv":
            delimiter = ","        
        np.savetxt(fname,series,delimiter=delimiter)
        return series

    def wtError(self,x,t):
        """Raises error for using RandomNLWavesFast class instead 

        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable


        Returns
        --------
        None

        Raises
        --------
        SystemExit

        """

        logEvent("ERROR! Wavetools.py: eta and u functions not available for this class. Please use RandomNLWavesFast for generating random waves with nonlinear correction",0)
        sys.exit(1)
    
    

class RandomNLWavesFast:
    """
    This class is used for generating plane random waves with 2ns order correction in an optimised manner
    using linear reconstruction of components from a wave spectrum
    
    Parameters
    ----------    
    param : Tstart
            Description: Start time
            Type: Float    
    param : Tend
            Description: End time
            Type: Float    
    param : x0
            Description: Position vector for the time series
            Type: Numpy array    
    param : Tp
            Description: Peak wave period
            Type: Float
    param : Hs
            Description: Significant wave height
            Type: float
    param : mwl
            Description: Still water level
            Type: float
    param : depth
            Description: Water depth
            Type: float
    param : waveDir
            Description: Wave direction vector
            Type: numpy array
    param : g
            Description: Gravitational acceleration vector 
            Type: numpy array
    param : N
            Description: Number of frequency components
            Type: integer
    param : bandfactor
            Description: Spectral band factor. fmax = bandFactor/Tp, fmin = 1/(bandFactor*Tp)
            Type: float
    param : spectName
            Description: Name of spectral distribution
            Type: string
    param : spectral_params
            Description: Dictionary of arguments specific to the spectral distribution 
            Example for JONSWAP = {"gamma": 3.3, "TMA":True,"depth": depth}
            TMA=True activates the TMA modification, which in turn needs the depth as a parameter
            Type: Dictionary 
    param : phi
            Description: Component phases (if set to None, phases are picked at random)
            Type: numpy array
    param : Vgen
            Description: Length of the generation zone (np.array([0., 0., 0.]) by default
            Type: numpy array
    """
    def __init__(self,
                 Tstart,
                 Tend,
                 x0,
                 Tp,                      #wave period
                 Hs,                      #significant wave height
                 mwl,                     #mean water level
                 depth,                   #water depth          
                 waveDir,                 #wave direction vector with three components
                 g,                       #gravitational accelaration vector with three components      
                 N,                       #number of frequency bins
                 bandFactor,              #width factor for band around peak frequency fp       
                 spectName,               #random words will result in error and return the available spectra 
                 spectral_params=None,    #JONPARAMS = {"gamma": 3.3, "TMA":True,"depth": depth} 
                 phi=None,
                 Vgen = np.array([0.,0.,0.])    #array of component phases
                 ):
        aR = RandomWaves(Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params,phi)
        aRN = RandomNLWaves(Tstart,Tend,Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params,phi)
        self.omega = aR.omega

        Tmax = 4.*pi/(max(self.omega)-min(self.omega))
        modes = ["short","linear","long"]
        periods = [Tp/2./1.1,Tp/1.1, Tmax]
        self.TS= []
        ii = -1
        for mode in modes:
            ii+=1
            fname = "randomNLWaves_"+mode+".csv"
            dt = periods[ii]/50.
            series = aRN.writeEtaSeries(Tstart,Tend,dt,x0,fname,mode,False,Vgen)
            Tstart_temp = series[0,0]
            cutoff = 0.2*Tp/(Tend-Tstart_temp)
            self.TS.append(TimeSeries(
                    fname,
                    0,
                    x0,
                    depth,
                    32,
                    mwl,
                    waveDir,
                    g,
                    cutoffTotal = cutoff,
                    rec_direct = False,
                    window_params = {"Nwaves":15 ,"Tm":periods[ii],"Window":"costap"},
                    arrayData = True,
                    seriesArray = series)
                           )
        self.series =  [ Tstart,Tend,
                         fname,
                    0,
                    x0,
                    depth,
                    32,
                    mwl,
                    waveDir,
                    g,
                    cutoff,
                     False,
                     {"Nwaves":15 ,"Tm":periods[ii],"Window":"costap"},
                     True,
                     series]

#        self.series = ii


    def eta(self,x,t):
        """Calculates free surface elevation (RandomNLWavesFast class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        float
            Free-surface elevation as a float

        """
        etaR = self.TS[0].eta(x,t)+ self.TS[1].eta(x,t)+self.TS[2].eta(x,t)
        return etaR


    def u(self,x,t):
        """Calculates wave velocity vector (RandomNLWavesFast class)
        Parameters
        ----------    
        x : numpy.ndarray
            Position vector
        t : float
            Time variable

        Returns
        --------
        numpy.ndarray
            Velocity vector as 1D array

        """
        uR = self.TS[0].u(x,t)+ self.TS[1].u(x,t)+self.TS[2].u(x,t)
        return uR
        



