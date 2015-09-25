"""Tools for working with water waves
The primary objective of this module is to provide solutions (exact and

approximate) for the free surface deformation and subsurface velocity
components of water waves. These can be used as boundary conditions,
wave generation sources, and validation solutions for numerical wave
codes.

"""
from math import pi, tanh, sqrt, exp, log, sin, cos, cosh, sinh
import numpy as np
from Profiling import logEvent
import time as tt
import sys as sys


def sigma(omega,omega0):
    """sigma function for JONSWAP spectrum
    """
    sigmaReturn = np.where(omega > omega0,0.09,0.07)
    return sigmaReturn

def JONSWAP(f,f0,Hs,g,gamma):
    """The wave spectrum from Joint North Sea Wave Observation Project

    :param f: wave frequency [1/T]
    :param f0: peak frequency [1/T]
    :param Hs: significant wave height [L]
    :param g: gravity [L/T^2]
    :param gamma: peak enhancement factor [-]
    """
    omega = 2.0*pi*f
    omega0 = 2.0*pi*f0
    alpha = 2.0*pi*0.0624*(1.094-0.01915*log(gamma))/(0.23+0.0336*gamma-0.0185/(1.9+gamma))
    r = np.exp(- (omega-omega0)**2/(2*sigma(omega,omega0)**2*omega0**2))
    return (alpha*Hs**2*omega0**4/omega**5)*np.exp(-(5.0/4.0)*(omega0/omega)**4)*gamma**r

def piersonMoskovitz(f,f0,alpha=8.1e-3,beta=0.74,g=9.8):
    """Pierson-Moskovitz spectrum

    :param f: frequency [1/T]
    :param f0: peak frequency [1/T]
    :param alpha: alpha fitting parameter [-]
    :param beta: beta fitting parameter [-]
    :param g: graivty [L/T^2]
    """
    return (5.0/16.0)*Hs**2*(f0**4/f**5)*np.exp((-5.0/4.0)*(f0/f)**4)

def cos2s(theta,s):
    return cos(theta/2)**(2*s)

def normInt(thetas,dir_fun,s,N):
    G0 = 0.
    theta = 0.
    for ii in range(N):
        G0+= dir_fun(theta,s)*dth
        theta+=dth
    return 1./G0
    

def dispersion(w,d, g = 9.81,niter = 1000):
    """Calculates wave vector k from linear dispersion relation

    :param w: cyclical frequency
    :param d: depth [L]
    :param niter: number  of solution iterations
    :param g: gravity [L/T^2
    """
#    print("Initiating dispersion")
    Kd = w*sqrt(d/g)
#    print("Initial dispersion value = %s" %str(Kd/d))
    for jj in range(niter):
       #Kdn_1 = Kd
        Kd = w*w*d/g/np.tanh(Kd)
        #Kdn_1 /=0.01*Kd
        #Kdn_1 -= 100.
        #Kdn_1 = abs(Kdn_1)
        #try: Kdn_1 = mean(Kdn_1)
        #except: continue
    #print "Solution convergence for dispersion relation %s percent" % Kdn_1
#    print("Final k value = %s" %str(Kd/d))
#    print("Wavelength= %s" %str(2.*pi*d/Kd))
    return(Kd/d)

class MonochromaticWaves:
    """Generate a monochromatic wave train in the linear regime
    """
    def __init__(self,period,waveHeight,mwl,depth,g,waveDir,wavelength=None,waveType="Linear",Ycoeff = None, Bcoeff =None, meanVelocity = 0.,phi0 = 0.):
        self.knownWaveTypes = ["Linear","Fenton","userDefined"]
        self.waveType = waveType
        self.g = g
        self.gAbs = sqrt(sum(g * g))
        self.waveDir = waveDir/sqrt(sum(waveDir * waveDir))
        self.phi0=phi0
        if self.waveType not in self.knownWaveTypes:
            logEvent("Wrong wavetype given: Valid wavetypes are %s" %(self.knownWaveTypes), level=0)
            sys.exit(1)
        self.dircheck = abs(sum(g * waveDir))
        #print self.dircheck
        if self.dircheck > 1e-6:
            logEvent("Wave direction is not perpendicular to gravity vector. Check input",level=0)
            sys.exit(1)
        self.period = period
        self.waveHeight = waveHeight
        self.mwl = mwl
        self.depth = depth
        self.omega = 2.0*pi/period
        if  self.waveType is "Linear":
            self.k = dispersion(w=self.omega,d=self.depth,g=self.gAbs)
            self.wavelength = 2.0*pi/self.k
        else:
            try:
                self.k = 2.0*pi/wavelength
                self.wavelength=wavelength
            except:
                logEvent("Wavelenght is not defined for nonlinear waves. Enter wavelength in class arguments",level=0)  
                sys.exit(1)
        self.kDir = self.k * self.waveDir 
        self.amplitude = 0.5*self.waveHeight
        self.meanVelocity = meanVelocity
        self.vDir = self.g/self.gAbs
        self.Ycoeff = Ycoeff
        self.Bcoeff = Bcoeff
        if (Ycoeff is None) or (Bcoeff is None):
            if self.waveType is not "Linear":
                pr.logEvent("Need to define Ycoeff and Bcoeff (free-surface and velocity) for nonlinear waves",level=0)                            
                sys.exit(1)
    def phase(self,x,y,z,t):        
#        return y*self.kDir[1] - self.omega*t + self.phi0
        return x*self.kDir[0]+y*self.kDir[1]+z*self.kDir[2] - self.omega*t + self.phi0
    
    def eta(self,x,y,z,t):
        if self.waveType is "Linear":
            global VTIME
            a = tt.time()
            aa = self.amplitude*cos(self.phase(x,y,z,t))
            VTIME = (tt.time() - a)
            return aa
        else:
            HH = 0.
            ii =0.
            for Y in self.Ycoeff:
                ii+=1
                HH+=Y*cos(ii*self.phase(x,y,z,t))
            if self.waveType is "Fenton": return HH/self.k
            else: return HH
    def Z(self,x,y,z):
        return   -(self.vDir[0]*x + self.vDir[1]*y+ self.vDir[2]*z) - self.mwl
                    
    def u(self,x,y,z,t,ss = "x"):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        UH=0.
        UV=0.
        ii=0.
        if self.waveType is "Linear":
            UH+=self.amplitude*self.omega*cosh(self.k*(self.Z(x,y,z)+self.depth))*cos(self.phase(x,y,z,t))/sinh(self.k*self.depth)
            UV+=self.omega*self.amplitude*sinh(self.k*(self.Z(x,y,z)+self.depth))*sin(self.phase(x,y,z,t))/sinh(self.k*self.depth)
#waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.waveDir,wavelength=self.wi[ii], phi0 = self.phi[ii]).u(x,y,z,t)
            Vcomp = {
                "x":UH*self.waveDir[0] + UV*self.vDir[0],
                "y":UH*self.waveDir[1] + UV*self.vDir[1],
                "z":UH*self.waveDir[2] + UV*self.vDir[2],
                }

        elif self.waveType is "Fenton":
            for B in self.Bcoeff:
                ii+=1
                UV+=ii*B*sinh(self.k*(self.Z(x,y,z)+self.depth))*sin(self.phase(x,y,z,t))/cosh(self.k*self.depth)
                UH+=ii*B*cosh(self.k*(self.Z(x,y,z)+self.depth))*cos(self.phase(x,y,z,t))/cosh(self.k*self.depth)
#waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.waveDir,wavelength=self.wi[ii], phi0 = self.phi[ii]).u(x,y,z,t)
                Vcomp = {
                    "x":UH*self.waveDir[0] + UV*self.vDir[0],
                    "y":UH*self.waveDir[1] + UV*self.vDir[1],
                    "z":UH*self.waveDir[2] + UV*self.vDir[2],

                    }
        else:
            logEvent("Check Wave types. Available wave types are %s" % waveType,level=0)
            exit(1)
        return Vcomp[ss]

class RandomWaves:
    """Generate approximate random wave solutions

    :param Tp: peak period [T]
    :param Hs: significant wave height [L]
    :param  d: depth [L]
    :param fp: frequency [1/T]
    :param bandFactor: width factor for band  around fp [-]
    :param N: number of frequency bins [-]
    :param mwl: mean water level [L]"""
    
    def __init__(self,
                 Tp = 5.0,         #s peak period
                 Hs = 2.0,         #m significant wave height
                 d = 2.0,           #m depth
                 fp = 1.0/5.0,      #peak  frequency
                 bandFactor = 2.0, #controls width of band  around fp
                 N = 101,          #number of frequency bins
                 mwl = 0.0,        #mean water level
                 waveDir = np.array([1,0,0]),
                 g = np.array([0, -9.81, 0]),         #accelerationof gravity
                 spec_fun = JONSWAP,
                 gamma=3.3
                 ): 
        self.gamma=gamma
        self.waveDir = waveDir/sqrt(sum(waveDir * waveDir))
        self.g = np.array(g)
        self.gAbs = sqrt(sum(g * g))
        self.vDir = self.g/self.gAbs
        self.Tp = Tp
        self.Hs = Hs
        self.depth = d
        self.fp = fp
        self.bandFactor = bandFactor
        self.N = N
        self.mwl = mwl
        self.fmax = self.bandFactor*self.fp
        self.fmin = self.fp/self.bandFactor
        self.df = (self.fmax-self.fmin)/float(self.N-1)
        self.fi=np.zeros(self.N,'d')
        for i in range(self.N):
            self.fi[i] = self.fmin+self.df*i
        self.omega = 2.*pi*self.fi
        self.ki = dispersion(2.0*pi*self.fi,self.depth,g=self.gAbs)
        self.wi = 2.*pi/self.ki
        self.phi = 2.0*pi*np.random.random(self.fi.shape[0])
        #ai = np.sqrt((Si_J[1:]+Si_J[:-1])*(fi[1:]-fi[:-1]))
        fim_tmp = (0.5*(self.fi[1:]+self.fi[:-1])).tolist()
        self.fim = np.array([fim_tmp[0]-0.5*self.df]+fim_tmp+[fim_tmp[-1]+0.5*self.df])
        self.Si_Jm = spec_fun(self.fim,f0=self.fp,Hs=self.Hs,g=self.g,gamma=self.gamma)
        self.ai = np.sqrt((self.Si_Jm[1:]+self.Si_Jm[:-1])*(self.fim[1:]-self.fim[:-1]))
        self.waves = MonochromaticWaves
        self.kDir = np.zeros((self.N, 3) , "d")
        for k in range(N):
            self.kDir[k,:] = self.ki[k]*self.waveDir[:]
    def Z(self,x,y,z):
        return   -(self.vDir[0]*x + self.vDir[1]*y+ self.vDir[2]*z) - self.mwl
    def eta(self,x,y,z,t):
        """Free surface displacement
        
        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.
        for ii in range(self.N):
            Eta+=self.ai[ii]*cos(x*self.kDir[ii,0]+y*self.kDir[ii,1]+z*self.kDir[ii,2] - self.omega[ii]*t + self.phi[ii])
        return Eta
#        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()
    
    def u(self,x,y,z,t,ss = "x"):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        UH=0.
        UV=0.
        for ii in range(self.N):
            UH+=self.ai[ii]*self.omega[ii]*cosh(self.ki[ii]*(self.Z(x,y,z)+self.depth))*cos(x*self.kDir[ii,0]+y*self.kDir[ii,1]+z*self.kDir[ii,2] - self.omega[ii]*t + self.phi[ii])/sinh(self.ki[ii]*self.depth)
            UV+=self.ai[ii]*self.omega[ii]*sinh(self.ki[ii]*(self.Z(x,y,z)+self.depth))*sin(x*self.kDir[ii,0]+y*self.kDir[ii,1]+z*self.kDir[ii,2] - self.omega[ii]*t + self.phi[ii])/sinh(self.ki[ii]*self.depth)
#waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.waveDir,wavelength=self.wi[ii], phi0 = self.phi[ii]).u(x,y,z,t)
        Vcomp = {
            "x":UH*self.waveDir[0] + UV*self.vDir[0],
            "y":UH*self.waveDir[1] + UV*self.vDir[1],
            "z":UH*self.waveDir[2] + UV*self.vDir[2],
            }
        return Vcomp[ss]
#        Z = z - self.mwl
#        return (2.0*pi*self.fi*self.ai*np.cos(2.0*pi*self.fi*t-self.ki*x+self.phi)*
#                np.cosh(self.ki*(self.d+Z))/np.sinh(self.ki*self.d)).sum()


class timeSeries:
    """Generate a time series by using spectral windowing method.

    :param ts: time-series array [T]
    :param  d: depth [L]
    :param Npeaks: Number of spectral peaks
    :param bandFactor[Npeaks]: width factor for band  around spectral peaks [-]
    :param peakFrequencies[Npeaks]: expected peak frequencies
    :param N: number of frequency bins [-]
    :param Nwaves: Number of waves per window (Approx)
    :param mwl: mean water level [L]"""
    
    def __init__(self,
                 timeSeriesFile= "Timeseries.txt",
                 skiprows = 0,
                 d = 2.0,
                 Npeaks = 1, #m depth
                 bandFactor = [2.0], #controls width of band  around fp
                 peakFrequencies = [1.0],
                 N = 32,          #number of frequency bins
                 Nwaves = 20,
                 mwl = 0.0,        #mean water level
                 waveDir = np.array([1,0,0]),
                 g = np.array([0, -9.81, 0])         #accelerationof gravity
                 ): 
        self.depth = d
        self.Npeaks = Npeaks
        self.bandFactor = np.array(bandFactor)        
        self.peakFrequencies = np.array(peakFrequencies)
        self.N = N
        self.Nwaves = Nwaves
        self.mwl = mwl
        self.waveDir = waveDir/sqrt(sum(waveDir * waveDir))
        self.g = np.array(g)

#derived variables
        self.gAbs = sqrt(sum(g * g))
        self.vDir = self.g/self.gAbs
        self.fmax = self.bandFactor*self.peakFrequencies
        self.fmin = self.peakFrequencies/self.bandFactor

        self.df = (self.fmax-self.fmin)/float(self.N-1)
        self.fi=np.zeros((self.N,self.Npeaks),'d')

        for i,j in np.ndenumerate(fi):
            self.fi[i,j] = self.fmin[i,j]+self.df[i,j]*i

        self.omegai = 2.*pi*self.fi
        self.ki = dispersion(self.omegai,self.depth,g=self.gAbs)
        self.wi = 2.*pi/self.ki
        
#Reading time series
        filetype = timeSeriesFile[-3:]
        if (filetype is not "txt") or (filetype is not "csv"):
            logEvent("Timeseries must be given in txt or csv format",level=0)
            sys.exit(1)
        elif (filetype is "csv"):
            tdata = np.loadtxt(fid,skiprows=skiprows,delimiter=",")
        else:
            tdata = np.loadtxt(fid,skiprows=skiprows)
#Checks for tseries file
        ncols = len(tdata[0,:])
        if ncols is not 2:
            logEvent("Timeseries file must have two colunms for a single probe data",level=0)
            sys.exit(1)



        time_temp = tdata[:,0] 
        self.dt = (time_temp[-1]-time_temp[0])/len(time_temp)
        doInterp = False
        for i in enumerate(time_temp,start=1):
            dt_temp = time_temp[i]-time_temp[i-1]
        #check if time is at first column
            if time_temp[i]<time_temp[i-1]:
                logEvent("Found not consistent time entry between %s and %s row. Warning: time variable must be always at the first column of the file and increasing monotonically" %(i-1,i) )
                sys.exit(1)            
        #check if sampling rate is constant
            if dt is not self.dt:
                logEvent("Not constant sampling rate found, proceeding to signal interpolation to a constant sampling rate")
                doInterp = True
        if(doInterp):
            self.time = np.linspace(time_temp[0],time_temp[-1],len(time_temp))
            self.eta = np.interp(self.time,time_temp,tdata[:,1])
        else:
            self.time = time_temp
            self.eta = tdata[:,1]
        del tdata


class directionalWaves:

    """Generate approximate directiona random wave solutions

    :param Tp: peak period [T]
    :param Hs: significant wave height [L]
    :param  d: depth [L]
    :param fp: frequency [1/T]
    :param bandFactor: width factor for band  around fp [-]
    :param N: number of frequency bins [-]
    :param mwl: mean water level [L]"""
    
    def __init__(self,
                 Tp ,         #s peak period
                 Hs ,         #m significant wave height
                 d ,           #m depth
                 fp ,      #peak  frequency
                 bandFactor, #controls width of band  around fp
                 N ,          #number of frequency bins
                 M  ,         #half number of directional bins
                 mwl,        #mean water level
                 waveDir,   #main wave direction
                 normalWaveDir,# Normal to wave direction, on propagation plane
                 g,         #accelerationof gravity
                 spec_fun ,                  # spectral function
                 thetamax = pi,         #max directional band, measured from lead wave direction, defaults to pi
                 s =5 ,                              # dir function coefficient
                 dir_fun = cos2s               # directional function
                 ): #wave spectrum


        self.waveDir = waveDir/sqrt(sum(waveDir * waveDir))
        self.normalWaveDir = normalWaveDir/sqrt(sum(normalWaveDir * normalWaveDir))
        self.g = g
        self.gAbs = sqrt(sum(g * g))
        self.Tp = Tp
        self.Hs = Hs
        self.d = d
        self.fp = fp
        self.bandFactor = bandFactor
        self.N = N
        self.mwl = mwl
        self.fmax = self.bandFactor*self.fp
        self.fmin = self.fp/self.bandFactor
        self.df = (self.fmax-self.fmin)/float(self.N-1)
        self.fi=np.zeros(self.N,'d')
        for i in range(self.N):
            self.fi[i] = self.fmin+self.df*i
        self.ki = dispersion(2.0*pi*self.fi,self.d,g=self.gAbs)
        self.wi = 2.*math.pi/self.ki
        #ai = np.sqrt((Si_J[1:]+Si_J[:-1])*(fi[1:]-fi[:-1]))
        fim_tmp = (0.5*(self.fi[1:]+self.fi[:-1])).tolist()
        self.fim = np.array([fim_tmp[0]-0.5*self.df]+fim_tmp+[fim_tmp[-1]+0.5*self.df])
        self.Si_Jm = spec_fun(self.fim,f0=self.fp,Hs=self.Hs,g=self.g,gamma=3.3)
        self.ai = np.sqrt((self.Si_Jm[1:]+self.Si_Jm[:-1])*(self.fim[1:]-self.fim[:-1]))
        self.waves = MonochromaticWaves
        self.M = M
        self.thetas = np.linspace(0,thetamax,self.M+1)        
        self.dth = thetas[1]-thetas[0]
        self.spread = dir_fun(thetas,s)
        self.dirs = zeros((2*self.M + 1,3),'d')
        self.ai_d = zeros((self.N,2*M+1),'d')
        self.phi = zeros((self.N,2*M+1),'d')
        self.G_Int = normInt(self.dth,self.dir_fun,s,self.M+1)
        for ii in range(1,self.M+1):
            self.dirs[self.M+ii,:]= cos(self.thetas[ii])*waveDir + sin(self.thetas[ii])*normalWaveDir
            self.dirs[self.M-ii,:] = cos(self.thetas[ii])*waveDir - sin(self.thetas[ii])*normalWaveDir
            self.ai_d[self.M+ii,:] = self.ai*self.G_Int*spread[ii]  
            self.ai_d[self.M-ii,:] = self.ai*self.G_Int*spread[ii]  
            self.phi[self.M+ii,:] = 2.0*pi*np.random.random(self.fi.shape[0])
            self.phi[self.M-ii,:] = 2.0*pi*np.random.random(self.fi.shape[0])

        self.dirs[self.M,:] = self.waveDir
        self.phi[self.M,:] = 2.0*pi*np.random.random(self.fi.shape[0])
        self.ai_d[self.M,:] = self.ai*self.G_Int*spread[0] 
 
    def eta(self,x,y,z,t):
        """Free surface displacement
        
        :param x: floating point x coordinate
        :param t: time"""
        Eta=0.
        for jj in range(2*self.M + 1):
            for ii in range(self.N):
                Eta+=waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii,jj],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.dirs[ii,jj],wavelength=wi[ii], phi0 = self.phi[ii,jj]).eta(x,y,z,t)
        return Eta
    #        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()
    
    def u(self,x,z,t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        U=0.
        for jj in range(2*self.M + 1):
            for ii in range(self.N):
                U+=waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii,jj],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.dirs[ii,jj],wavelength=wi[ii], phi0 = self.phi[ii,jj]).u(x,y,z,t)
        return U
    
    def v(self,x,z,t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        V=0.
        for jj in range(2*self.M + 1):
            for ii in range(self.N):
                V+=waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii,jj],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.dirs[ii,jj],wavelength=wi[ii], phi0 = self.phi[ii,jj]).v(x,y,z,t)
        return V
    
    def w(self,x,z,t):
        """x-component of velocity

        :param x: floating point x coordinate
        :param z: floating point z coordinate (height above bottom)
        :param t: time
        """
        W=0.
        for jj in range(2*self.M + 1):
            for ii in range(self.N):
                W+=waves(period = 1./self.fi[ii], waveHeight = 2.*self.ai[ii,jj],mwl = self.mwl, depth = self.d,g = self.g,waveDir = self.dirs[ii,jj],wavelength=wi[ii], phi0 = self.phi[ii,jj]).w(x,y,z,t)
        return W
    
    

if __name__ == '__main__':
    from matplotlib.pyplot import *
    import os as os
    global VTIME
    VTIME = 0.

    def etaCalc(x,y,z,t,wave_fun,ttime):
        eta = np.zeros((len(x),len(y),len(z),len(t)),float)
        for n,T in enumerate(t):
            for J,xi in enumerate(x):
                for I,yi in enumerate(y):
                    for K,zi in enumerate(z):                        
                        eta[J,I,K,n] = wave_fun.eta(xi,yi,zi,T)

        return eta
    def velCalc(x,y,z,t,wave_fun,ttime,ss):
        eta = np.zeros((len(x),len(y),len(z),len(t)),float)
        for n,T in enumerate(t):
            for J,xi in enumerate(x):
                for I,yi in enumerate(y):
                    for K,zi in enumerate(z):                        
                        eta[J,I,K,n] = wave_fun.u(xi,yi,zi,T,ss)

        return eta
    def plotSeriesAlongAxis(x,t,ts,ifig,string):
       # print(x)
        fig = figure(ifig)
        for i,xi in enumerate(x):
            line1 = plot(t,ts[i,:])
            
        savefig("timeseries_%s.png" %string)


    print "Loading variables"
    
    Tp = 5.0 #s peak period
    Hs = 2.0 #m significant wave height
    mwl = 0.0 #mean water level    Hs = 2.0
    depth = 10.0 # water depth
    waveDir = np.array([0,1,0])#wave Direction
    g = np.array([0,0,-9.81]) # 
    print "Setting space and time arrays"
    li =2.*pi/dispersion(2.*pi/Tp,depth,g=9.81)
    bandFactor = 2.
    x = np.linspace(0,0,1)
    y = np.linspace(0,0,1)
    z = np.linspace(0,0,1)
    t=np.linspace(0,50.*Tp/1.1,625)
#Fenton coefficients
    
    Y = [0.04160592, #Surface elevation Fourier coefficients for non-dimensionalised solution
         0.00555874,
         0.00065892,
         0.00008144,
         0.00001078,
         0.00000151,
         0.00000023,
         0.00000007]

    B = [0.05395079,
         0.00357780,
         0.00020506,
         0.00000719,
         -0.00000016,
         -0.00000005,
         0.00000000,
         0.00000000]
    print "Calculating waves"



    waveType = "Fenton" # Change between Linear, Fenton and random
    if waveType is "Linear": 
        waves = MonochromaticWaves( Tp, Hs, mwl,  depth, g,  waveDir)
    elif waveType is "Fenton": 
        wlength = 8.12
        t=np.linspace(0,30,500)
        waves = MonochromaticWaves( 2.95, 0.109 , mwl,  0.873, g,  waveDir,wlength, waveType, Y, B) #period,waveHeight,mwl,depth,g,waveDir,wavelength=None,waveType="Linear",Ycoeff = None, Bcoeff =None, meanVelocity = 0.,phi0 = 0.):
    elif waveType is "Random":
        waves = RandomWaves(Tp = Tp,
                            Hs = Hs,
                            d = depth,
                            fp = 1./Tp,
                            bandFactor = bandFactor,
                            N = 101,
                            mwl = mwl,
                            g = g)
    print "Calculating free-surface"

    print "Start-time: %s" %tt.time()
    ttime = tt.time()
    eta = etaCalc(x,y,z,t,waves,ttime)
    v1 = velCalc(x,y,z,t,waves,ttime,"x")
    v2 = velCalc(x,y,z,t,waves,ttime,"y")
    v3 = velCalc(x,y,z,t,waves,ttime,"z")
    duration = tt.time()-ttime
    print "Duration: %s" %duration
    print "Plotting free surface elevation"

    plotSeriesAlongAxis(x,t,eta[0,0,:,:],0,"Eta")
    plotSeriesAlongAxis(z,t,v1[0,0,:,:],1,"UX")
    plotSeriesAlongAxis(y,t,v2[0,0,:,:],2,"UY")
    plotSeriesAlongAxis(z,t,v3[0,0,:,:],3,"UZ")
    print VTIME
