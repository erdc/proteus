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
    """ Checks if a function name  is present in a list of known functions, returns system exit if not present
    param: funcName : function name in form of string under consideration
    param: validFunctions: list of valid functions objects (not names in strings)

    """
    funcNames = []
    for func  in validFunctions:
            funcNames.append(func.__name__)
            if func.__name__ == funcName:
                func_ret = func
    if funcName not in funcNames:
        logEvent("WaveTools.py: Wrong function type (%s) given: Valid wavetypes are %s" %(funcName,funcNames), level=0)
        sys.exit(1)
    return func_ret
       


def setVertDir(g):
    """ Sets the unit vector for the vertical direction, opposite to the gravity vector
    param: g : gravitational acceleration vector [L/T^2] (must have 3 components)
    """
    return -g/(sqrt(g[0]**2 + g[1]**2 + g[2]**2))

def setDirVector(vector):
    """ Returns the direction of a vector in the form of a unit vector
    param: vector : Any vector with three components
    """
    return vector/(sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2))

def dirCheck(v1, v2):
    """ Checks if to vectors are vertical and returns system exit if not
    param: v1 : 1st vector  [-]  with three components
    param: v2 : 2nd vector  [-]  with three components
    """
    dircheck = abs(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
        #print self.dircheck
    if dircheck > 1e-10:
        logEvent("Wave direction is not perpendicular to gravity vector. Check input",level=0)
        return sys.exit(1)
    else:
        return None
def reduceToIntervals(fi,df):
    """ Prepares the x- axis array of size N for numerical integration along he x- axis. 
    If fi = [a1, a2, a3,...,a_N-1 a_N] then it returns the array 
    [a1, 0.5(a1+a2), 0.5(a2+a3),...0.5(a_N-1+a_N), a_N]. Input array must have constant step
    param: fi : x- array  [-]
    param: df : dx constant step of array  [-]
    """
    fim_tmp = (0.5*(fi[1:]+fi[:-1])).tolist()
    return np.array([fim_tmp[0]-0.5*df]+fim_tmp+[fim_tmp[-1]+0.5*df])
def returnRectangles(a,x):
    """ Returns \delta y of y(x) using the rectangle method (\delta y = 0.5*(a_n-1+a_n)*(x_n-1-x_n) 
    param: a : y(x) function   [-]
    param: x : x- coordinate  [-]
    """
    return 0.5*(a[1:]+a[:-1])*(x[1:]-x[:-1])
def returnRectangles3D(a,x,y):
    """ Returns \delta y of  y(x,z) using the rectangle method 
    \delta y = 0.25*(a_(n-1,m-1)+a_(n,m-1)+a_(n-1,m)+a_(n,m))*(x_n-1-x_n) *(z_m-1-z_m)
    param: a : a(x,y) function   [-]
    param: x : x- coordinate  [-]
    param: y : y- coordinate  [-]
    """
    ai = 0.5*(a[1:,:]+a[:-1,:])
    ai = 0.5*(ai[:,1:]+ai[:,:-1])
    for ii in range(len(x)-1):
        ai[ii,:] *= (y[1:]-y[:-1]) 
    for jj in range(len(y) - 1):
        ai[:,jj] *= (x[1:]-x[:-1])    
    return ai
def normIntegral(Sint,th):
    """Given an Sint(th) function, it returns Sint_n, such as \int (Sint_n dth = 1)
    param: Sint : Sint(th) function   [-]
    param: th : th- coordinate  [-]
    """
    G0 = 1./sum(returnRectangles(Sint,th))
    return G0*Sint



def eta_mode(x, t, kDir, omega, phi, amplitude):
    """Returns a single frequency mode for free-surface elevation at point x,y,z,t
    :param kDir: wave number vector [1/L] with three components
    :param omega: angular frequency [1/T]
    :param phi: phase [0,2*pi]
    :param amplitude: wave amplitude [L/T^2]
    """
    phase = x[0]*kDir[0]+x[1]*kDir[1]+x[2]*kDir[2] - omega*t  + phi
    return amplitude*cos(phase)


def vel_mode(x, t, kDir, kAbs, omega, phi, amplitude, mwl, depth, g, vDir):
    """Returns a single frequency mode for velocity at point x,y,z,t
    :param kDir: wave number vector [1/L] with three components
    :param omega: angular frequency [1/T]
    :param phi: phase [0,2*pi]
    :param amplitude: wave amplitude [L/T^2]
    :param mwl: mean water level [L]
    :param depth: water depth [L]
    :param g: gravity vector
    :param vDir (vertical direction - opposite to the gravity vector)
    :param comp: component "x", "y" or "z"
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
    """sigma function for JONSWAP spectrum
       http://www.wikiwaves.org/Ocean-Wave_Sectra
    """
    sigmaReturn = np.where(omega > omega0,0.09,0.07)
    return sigmaReturn


def JONSWAP(f,f0,Hs,gamma=3.3,TMA=False, depth = None):
    """The wave spectrum from Joint North Sea Wave Observation Project
    Jonswap equation from "Random Seas and Design of Maritime Structures" - Y. Goda - 2010 (3rd ed) eq. 2.12 - 2.15
    TMA modification from "Random Seas and Design of Maritime Structures" - Y. Goda - 2010 (3rd ed) eq. 2.19
    :param f: wave frequency [1/T] (not angular frequency)
    :param f0: direcpeak frequency [1/T] (not angular frequency)
    :param Hs: significant wave height [L]
    :param g: gravity [L/T^2]
    :param gamma: peak enhancement factor [-]
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
    """modified Pierson-Moskovitz spectrum (or Bretschneider or ISSC)
    Reference http://www.orcina.com/SoftwareProducts/OrcaFlex/Documentation/Help/Content/html/Waves,WaveSpectra.htm
    And then to Tucker M J, 1991. Waves in Ocean Engineering. Ellis Horwood Ltd. (Chichester).
    :param f: frequency [1/T]
    :param f0: peak frequency [1/T]
    :param alpha: alpha fitting parameter [-]
    :param beta: beta fitting parameter [-]
    :param g: graivty [L/T^2]
    """
    return (5.0/16.0)*Hs**2*(f0**4/f**5)*np.exp((-5.0/4.0)*(f0/f)**4)

def cos2s(theta,f,s=10):
    """The cos2s wave directional Spread 
    see USACE - CETN-I-28 http://chl.erdc.usace.army.mil/library/publications/chetn/pdf/cetn-i-28.pdf
    :param theta: ange of wave direction, with respect to the peak direction
    :param f: wave frequency [1/T] (not angular frequency). Dummy variable in this one
    :param s: directional peak parameter. as s ->oo the distribution converges to 
    """
    fun = np.zeros((len(theta),len(f)),)
    for ii in range(len(fun[0,:])):
        fun[:,ii] = np.cos(theta/2)**(2*s)
    return fun
def mitsuyasu(theta,fi,f0,smax=10):
    """The cos2s wave directional spread with wave frequency dependency (mitsuyasu spread) 
    Equation from "Random Seas and Design of Maritime Structures" - Y. Goda - 2010 (3rd ed) eq. 2.22 - 2.25
    :param theta: ange of wave direction, with respect to the peak direction
    :param f: wave frequency [1/T] (not angular frequency). Dummy variable in this one
    :param s: directional peak parameter. as s ->oo the distribution converges to 
    """

    s = smax * (fi/f0)**(5)
    ii = np.where(fi>f0)[0][0]
    s[ii:] = smax * (fi[ii:]/f0)**(-2.5)
    fun = np.zeros((len(theta),len(fi)),)
    for ii in range(len(fun[0,:])):
        fun[:,ii] = np.cos(theta/2)**(2.*s[ii])
    return fun





def dispersion(w,d, g = 9.81,niter = 1000):
    """Calculates wave number magnitude as a scallar or an arry of modes linear dispersion relation

    :param w: cyclical frequency (can be scalar or an aray of frequency modes)
    :param d: depth [L]
    :param niter: number  of solution iterations
    :param g: gravity [L/T^2
    """
#    print("Initiating dispersion")
    w_aux = np.array(w)
    K = w_aux**2/g
#    print("Initial dispersion value = %s" %str(Kd/d))
    for jj in range(niter):
       #Kdn_1 = Kd
        K =  w_aux**2/(g*np.tanh(K*d))
        #Kdn_1 /=0.01*Kd
        #Kdn_1 -= 100.
        #Kdn_1 = abs(Kdn_1)
        #try: Kdn_1 = mean(Kdn_1)
        #except: continue
    #print "Solution convergence for dispersion relation %s percent" % Kdn_1
#    print("Final k value = %s" %str(Kd/d))
#    print("Wavelength= %s" %str(2.*pi*d/Kd))
    if type(K) is float:
        return K[0]
    else:
        return K


def tophat(l,cutoff):
    """ returns a top hat filter 
    :param l: array length
    :param cutoff: cut off fraction at either side of the array zero values will be imposed at the first and last cutoff*l array elements

    """
    a = np.zeros(l,)
    cut = int(cutoff*l)
    a[cut:-cut] = 1.
    return a

def costap(l,cutoff=0.1):
    """ Cosine taper filter Goda (2010), Random Seas and Design of Maritime Structures equation 11.40   
    :param l: array length
    :param cutoff: cut off fraction at either side of the array zero values will be imposed at the first and last cutoff*l array elements"""
    npoints = int(cutoff*l)
    wind = np.ones(l)
    for k in range(l): # (k,np) = (n,N) normally used
        if k < npoints:
            wind[k] = 0.5*(1.-cos(pi*float(k)/float(npoints)))
        if k > l - npoints -1:
            wind[k] = 0.5*(1.-cos(pi*float(l-k-1)/float(npoints)))
    return wind

def decompose_tseries(time,eta,dt):
    """ This function does a spectral decomposition of a time series with constant sampling.
     It returns a list with results with four entries:
         0 -> numpy array with frequency components ww
         1 -> numpy array with amplitude of each component aa
         2 -> numpy array with phase of each component pp
         3 -> float of the 0th fourier mode (wave setup) 
         :param time: time array [T]
         :param eta: signal array
         :param dt: sampling frequency [1/T] 
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
    """Generate a monochromatic wave train in the linear regime
    :param period: Monochromatic wave period
    :param waveHeight: Monochromatic wave height
    :param mwl: Mean water level
    :param depth: Mean water depth
    :param g: gravitational acceleration
    :param waveDir: wave direction vector (all 3 components needed)
    :param wavelength: wavelength for nonlinear (Fenton) waves. Can assume None if waves are linear, need to declare if waveType is Fenton
    :param waveType: can be Linear or Fenton (nonlinear). Linear by default
    :param Ycoeff: Y coefficient array for Fenton waves (see JD Fenton (1988) THE NUMERICAL SOLUTION OF STEADY WATER WAVE PROBLEMS, Computer and Geosciences, 14(3), 357-368, 
                   http://johndfenton.com/Papers/Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf    
    :param BCoeff: B coefficient array for Fenton waves (see reference above)
    :meanVelocity: Current velocity. Recommended use with Fenton waves
    :phi0: Phase of the wave                 
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
        if  self.waveType is "Linear":
            self.k = dispersion(w=self.omega,d=self.depth,g=self.gAbs)
            self.wavelength = 2.0*pi/self.k
        else:
            try:
                self.k = 2.0*pi/wavelength
                self.wavelength=wavelength
            except:
                logEvent("WaveTools.py: Wavelenght is not defined for nonlinear waves. Enter wavelength in class arguments",level=0)
                sys.exit(1)
        self.kDir = self.k * self.waveDir
        self.amplitude = 0.5*self.waveHeight
        self.meanVelocity = np.array(meanVelocity)
#Checking that meanvelocity is a vector

        if(len(meanVelocity) != 3):
            logEvent("WaveTools.py: meanVelocity should be a vector with 3 components. ",level=0)
            sys.exit(1)

        self.Ycoeff = Ycoeff
        self.Bcoeff = Bcoeff

# Checking for
        if (Ycoeff is None) or (Bcoeff is None):
            if self.waveType is not "Linear":
                logEvent("WaveTools.py: Need to define Ycoeff and Bcoeff (free-surface and velocity) for nonlinear waves",level=0)
                sys.exit(1)
    def eta(self, x, t):
        if self.waveType is "Linear":
            return eta_mode(x,t,self.kDir,self.omega,self.phi0,self.amplitude)
        elif self.waveType is "Fenton":
            HH = 0.
            ii =0.
            for Y in self.Ycoeff:
                ii+=1
                HH+=eta_mode(x,t,ii*self.kDir,ii*self.omega,self.phi0,Y)
            return HH/self.k

    def u(self, x, t):
        if self.waveType is "Linear":
            return vel_mode(x, t, self.kDir,self.k,self.omega,self.phi0,self.amplitude,self.mwl,self.depth,self.g,self.vDir)
        elif self.waveType is "Fenton":
            Ufenton = self.meanVelocity
            ii = 0
            for B in self.Bcoeff:
                ii+=1
                wmode = ii*self.omega
                kmode = ii*self.k
                kdir = self.waveDir*kmode
                amp = tanh(kmode*self.depth)*sqrt(self.gAbs/self.k)*B/self.omega
                Ufenton+= vel_mode(x,t,kdir,kmode,wmode,self.phi0,amp,self.mwl,self.depth,self.g,self.vDir)
            return Ufenton # + self.meanVelocity[comp]


class RandomWaves:
    """Generate approximate random wave solutions
    :param Tp: frequency [1/T]
    :param Hs: significant wave height [L]
    :param mwl: mean water level [L]
    :param  depth: depth [L]
    :param waveDir:wave Direction vector with three components [-]
    :param g: Gravitational acceleration vector with three components [L/T^2]
    :param N: number of frequency bins [-]
    :param bandFactor: width factor for band  around fp [-]
    :param spectName: Name of spectral function in string format. Use a random word and run the code to obtain the vaild spectra names
    :param spectral_params: Dictionary of additional arguments for spectral function, specific to each spectral function, except from Hs and Tp e.g. {"gamma": 3.3, "TMA" = True, "depth" = 1} for Jonswap. Check spectral function arguments 
    :param phi: Array of component phases - if set to none, random phases are assigned
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
            logEvent('WaveTools.py: No phase array is given. Assigning random phases. Outputing the phasing of the random waves')
        else:
            try: 
                self.phi = np.array(phi)
                if self.phi.shape[0] != self.fi.shape[0]:
                    logEvent('WaveTools.py: Phase array must have N elements')
                    sys.exit(1)
                    
            except:
                logEvent('WaveTools.py: phi argument must be an array with N elements')
                sys.exit(1)

        #ai = np.sqrt((Si_J[1:]+Si_J[:-1])*(fi[1:]-fi[:-1]))
        self.fim = reduceToIntervals(self.fi,self.df)
        if (spectral_params == None):
            self.Si_Jm = spec_fun(self.fim,self.fp,self.Hs)
        else:
            try:
                self.Si_Jm = spec_fun(self.fim,self.fp,self.Hs,**spectral_params)
            except:
                logEvent('WaveTools.py: Additional spectral parameters are not valid for the %s spectrum' %spectName)
                sys.exit(1)
        

        self.ai = np.sqrt(2.*returnRectangles(self.Si_Jm,self.fim))
        self.kDir = np.zeros((len(self.ki),3),)
        for ii in range(3):
             self.kDir[:,ii] = self.ki[:] * self.waveDir[ii] 
    def eta(self, x, t):
        """Free surface displacement

        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.
        for ii in range(self.N):
            Eta+= eta_mode(x, t,self.kDir[ii],self.omega[ii],self.phi[ii],self.ai[ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        U=0.
        for ii in range(self.N):
            U+= vel_mode(x, t, self.kDir[ii], self.ki[ii],self.omega[ii],self.phi[ii],self.ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U

class MultiSpectraRandomWaves(RandomWaves):
    """Generate a random wave timeseries from multiple spectra. 
    Same input parameters as RandomWaves class but they have to be all in lists with the same lenght as the spectra (except from g!)
    :param Nspectra, number of spectra
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

                logEvent('WaveTools.py: Parameters passed in MultiSpectraRandomWaves must be in array or list form with length Nspectra  ')
                sys.exit(1)
               
        except:
            logEvent('WaveTools.py: Parameters passed in MultiSpectraRandomWaves must be in array or list form with length Nspectra  ')
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
            logEvent("WaveTools.py: Reading spectra No %s" %kk)
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
        """Free surface displacement

        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.
        for ii in range(self.Nall):
            Eta+= eta_mode(x, t, self.kDirM[ii],self.omegaM[ii],self.phiM[ii],self.aiM[ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        U=0.
        for ii in range(self.Nall):
            U+= vel_mode(x,t,self.kDirM[ii], self.kiM[ii],self.omegaM[ii],self.phiM[ii],self.aiM[ii],self.mwl,self.depth,self.g,self.vDir)
        return U



class DirectionalWaves(RandomWaves):
    """Generate a random wave timeseries from directional waves
    Same input parameters as RandomWaves with the addition of:
    :param M: number of discrete directions
    :param waveDir0: lead direction in vector form (3 components required)
    :spreadName: Spreading function name (can be cos2s or mitsuyashu), given in string format
    :spread_params: Parameters specific to each spread functions, e.g. {"s":15} or {fi0: 1, smax=20}, except from f and theta. Check spread functions for more info
    :phiSymm: Logical variable, by default False, when set to True it generated same phase for symmetrically arranged directions, with respect to the lead direction
    
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
            logEvent("WaveTools.py: phi in DirectionalWaves class must be given either as None or as a list with 2*M + 1 numpy arrays with length N")
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
                logEvent('WaveTools.py: Additional spread parameters are not valid for the %s spectrum' %spectName)
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
        """Free surface displacement

        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.
        for jj in range(self.Mtot):
            for ii in range(self.N):
                kDiri = self.waveDirs[jj]*self.ki[ii]
                Eta+= eta_mode(x,t,kDiri,self.omega[ii],self.phiDirs[jj,ii],self.aiDirs[jj,ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()

    def u(self, x, t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        U=0.
        for jj in range(self.Mtot):
            for ii in range(self.N):
                kDiri = self.waveDirs[jj]*self.ki[ii]
                U+= vel_mode(x,t,kDiri, self.ki[ii],self.omega[ii],self.phiDirs[jj,ii],self.aiDirs[jj,ii],self.mwl,self.depth,self.g,self.vDir)
        return U
     


            
                

        

class TimeSeries:
    """Generate a time series by using spectral windowing method.
    :param timeSeriesFile: Time series file name
    :param skiprows: How many rows to skip while reading time series
    :param: timeSeriesPosition: 3D vector showing the position of the surface elevation sampling
    :param  depth: depth [L]
    :param N: number of frequency bins [-]
    :param mwl: mean water level [L]
    :param waveDir: wave Direction vector
    :param g: Gravitational acceleration vector (3 components required)
    :param rec_direct: Logical variable, True for direct reconstruction, False for windowed reconstrunction
    :window_params: dictionary for window reconstruction parameters. Mandatory definition for Nwaves (how many waves per window) Tm (mean wave period), wind_filt (window filter name in string form). Optional: Overlap (window overlap as a percentage of window lenght), Cutoff (length of domain wher filter is applied, as a percentage of the 1/2 of window length)
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
                 rec_direct = True,
                 window_params = None #If rec_direct = False then wind_params = {"Nwaves":Nwaves,"Tm":Tm,"Window":wind_filt,"Overlap":overlap,"Cutoff":cutoff}
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
            logEvent("WaveTools.py: Location vector for timeSeries must have three-components",level=0)
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
        filetype = timeSeriesFile[-4:]
        logEvent("WaveTools.py: Reading timeseries from %s file: %s" % (filetype,timeSeriesFile),level=0)
        fid = open(timeSeriesFile,"r")
        if (filetype !=".txt") and (filetype != ".csv"):
                logEvent("WaveTools.py: File %s must be given in .txt or .csv format" % (timeSeriesFile),level=0)
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
            logEvent("WaveTools.py: Timeseries file (%s) must have only two columns [time, eta]" % (timeSeriesFile),level=0)
            sys.exit(1)
        time_temp = tdata[:,0]
        self.dt = (time_temp[-1]-time_temp[0])/(len(time_temp)-1)



        # If necessary, perform interpolation
        doInterp = False
        for i in range(1,len(time_temp)):
            dt_temp = time_temp[i]-time_temp[i-1]
        #check if time is at first column
            if time_temp[i]<=time_temp[i-1]:
                logEvent("WaveTools.py:  Found not consistent time entry between %s and %s row in %s file. Time variable must be always at the first column of the file and increasing monotonically" %(i-1,i,timeSeriesFile) )
                sys.exit(1)
        #check if sampling rate is constant
            if dt_temp!=self.dt:
                doInterp = True
        if(doInterp):
            logEvent("WaveTools.py: Not constant sampling rate found, proceeding to signal interpolation to a constant sampling rate",level=0)
            self.time = np.linspace(time_temp[0],time_temp[-1],len(time_temp))
            self.eta = np.interp(self.time,time_temp,tdata[:,1])
        else:
            self.time = time_temp
            self.eta = tdata[:,1]

        self.t0  = self.time[0]        
        # Remove mean level from raw data
        self.eta -= np.mean(self.eta)
        # Filter out first 2.5 % and last 2.5% to make the signal periodic
        self.eta *= costap(len(self.time),cutoff=0.025)
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
            logEvent("WaveTools.py: performing a direct series decomposition")
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
                logEvent("WaveTools.py: Set parameters for spectral windowing. Argument window_params must be a dictionary")
                sys.exit(1)
            try:
                self.Nwaves = window_params["Nwaves"]
            except:
                logEvent("WaveTools.py: Dictionary key 'Nwaves' (waves per window) not found in window_params dictionary")
                sys.exit(1)

            try:           
                self.Tm = window_params["Tm"]
            except:
                logEvent("WaveTools.py: Dictionary key 'Tm' (mean or characteristic wave period) not found in window_params dictionary")
                sys.exit(1)

            try:           
                self.windowName = window_params["Window"]
            except:
                logEvent("WaveTools.py: Dictionary key 'Window' (window function type) not found in window_params dictionary")
                sys.exit(1)

            if(self.Nwaves > 0.5*self.tlength / self.Tm):
                logEvent("WaveTools.py: Reconstruction is expected to have two windows or less. Plese reduce the number of waves per window or switch to direct decomposition )")
                sys.exit(1)



            validWindows = [costap, tophat]
            wind_filt =  loadExistingFunction(self.windowName, validWindows) 
            logEvent("WaveTools.py: performing series decomposition with spectral windows")
            # Portion of overlap, compared to window time
            try:
                self.overlap = window_params["Overlap"]            
            except:
                self.overlap = 0.25
                logEvent("WaveTools.py: Overlap entry in window_params dictionary not found. Setting default value of 0.25 (1/4 of the window length)")

            try:
                self.cutoff = window_params["Cutoff"]            
            except:
                self.cutoff= 0.1
                logEvent("WaveTools.py: Cutoff entry in window_params dictionary not found. Setting default value of 0.1 (1/10 of the window length)")
                
                

            # Portion of window filtered with the Costap filter
            # Setting the handover time, either at the middle of the overlap or just after the filter
            self.handover = max(1.1 *self.cutoff,  self.overlap / 2.)
            if (self.handover > 0.9 * self.overlap):
                logEvent("WaveTools.py: Window handover is not optimal as the cutoff is too close to the overlap. Decrease cutoff or increase overlap")
                sys.exit(1)
            self.Twindow =  self.Tm * self.Nwaves            # setting the window duration (approx.). Twindow = Tmean * Nwaves
            self.Toverlap = self.overlap * self.Twindow             
            self.Nwindows = int( (self.tlength -   self.Twindow ) / (self.Twindow - self.Toverlap) ) + 1             #Getting the actual number of windows  (N-1) * (Twindow - Toverlap) + Twindow = total time
            self.Twindow = self.tlength/(1. + (1. - self.overlap)*(self.Nwindows-1))            # Correct Twindow and Toverlap for duration and integer number of windows
            self.Toverlap = self.overlap*self.Twindow
            logEvent("WaveTools.py: Correcting window duration for matching the exact time range of the series. Window duration correspond to %s waves approx." %(self.Twindow / self.Tm) )
            diff = (self.Nwindows-1.)*(self.Twindow -self.Toverlap)+self.Twindow - self.tlength
            logEvent("WaveTools.py: Checking duration of windowed time series: %s per cent difference from original duration" %(100*diff) )
            logEvent("WaveTools.py: Using %s windows for reconstruction with %s sec duration and %s per cent overlap" %(self.Nwindows, self.Twindow,100*self.overlap ))
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
        """Free surface displacement
        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.        
        for ii in range(0,self.Nf):
            x1 = np.array(x)-[self.x0, self.y0, self.z0]
            Eta+= eta_mode(x1,t-self.t0,self.kDir[ii],self.omega[ii],self.phi[ii],self.ai[ii])
        return Eta

    def uDirect(self, x, t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        U=0.
        for ii in range(0,self.Nf):
            x1 = x-[self.x0, self.y0, self.z0]
            U+= vel_mode(x1, t-self.t0, self.kDir[ii],self.ki[ii], self.omega[ii],self.phi[ii],self.ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U

    def findWindow(self,t):
        term = 1. - self.handover
        if t-self.time[0] >= term*self.Twindow:
            Nw = min(int((t-self.time[0] - term*self.Twindow)/(self.Twindow - 2. * self.handover * self.Twindow)) + 1, self.Nwindows-1)
            if t-self.time[0] < self.windows_handover[Nw-1] - self.time[0]:
                Nw-=1
        else:
            Nw = 0
        return Nw
        
    def etaWindow(self, x, t):
        """Free surface displacement
        :param x: floating point x coordinate
        :param t: time"""
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
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
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
            x1 = x-[self.x0, self.y0, self.z0]
            U+= vel_mode(x1, t-t0, kDir[ii],ki[ii],omega[ii],phi[ii],ai[ii],self.mwl,self.depth,self.g,self.vDir)
        return U





