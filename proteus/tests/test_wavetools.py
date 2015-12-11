from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest
import random
from math import cos,sin,cosh,sinh,pi,tanh,log
import sys,os
import logging


comm = Comm.init()
Profiling.procID = comm.rank()
def getpath():
    path =str(os.getcwd())
    if "tests" in path[-6:]:
        path =""
    else:
        path = path+"/proteus/tests/"
    return path


Profiling.logEvent("Testing WaveTools")
class TestAuxFunctions(unittest.TestCase):
    def testVDir(self):
        from proteus.WaveTools import setVertDir
        self.assertTrue(np.array_equal(setVertDir(np.array([0,-9.81,0])), np.array([0,1,0])))
    def testDirVector(self):
        from proteus.WaveTools import setDirVector
        self.assertTrue(all(setDirVector(np.array([2.,2.,1.]))== np.array([2.,2.,1])/3.))
    def testDirCheck(self):
        from proteus.WaveTools import dirCheck
        dirCheck(np.array([1.,2.,3.]),np.array([7.,4.,-5.]) )# Just loading the function with two vertical vectors
        with self.assertRaises(SystemExit) as cm:
            dirCheck(np.array([9,9,9]),np.array([4,5,6]))
        self.assertEqual(cm.exception.code, 1)     
        
    def testReduceToIntervals(self):
        from proteus.WaveTools import reduceToIntervals
        fi = np.linspace(0, 100, 101)
        df = 1
        fr = reduceToIntervals(fi,df)
        fi[:-1]+=0.5
        fi_te = np.zeros(len(fi)+1,)
        fi_te[1:] = fi[:]
        self.assertTrue((fr- fi_te == 0).all())

    def testIntegrateRectangles(self): # Testing the integration fynction for y = 2*x at [0,1]. The area should be 1

        from proteus.WaveTools import reduceToIntervals,returnRectangles 
        x = np.linspace(0,1,101)
        dx = 0.01
        xim = reduceToIntervals(x,dx)
        y = 2*xim 
        A = sum(returnRectangles(y,xim))
        self.assertTrue(round(A,10) == 1.0)
    def testIntegrateRectangles3D(self): # Testing the integration fynction for y = 2*x at [0,1]. The area should be 1
        from proteus.WaveTools import reduceToIntervals,returnRectangles3D
        x = np.linspace(0,1,101)
        dx = 0.01
        z = np.linspace(0,1,201)
        dz = 0.005
        xim = reduceToIntervals(x,dx)
        zim = reduceToIntervals(z,dz)
        y1 = np.zeros((len(xim),len(zim)),)
                
        for j in range(len(zim)):
            for i in range(len(xim)):
                y1[i,j] = 2.*xim[i]*zim[j]
        
        A = sum(sum(returnRectangles3D(y1,xim,zim)))
        # Integrate function z*(2*x) over x[0,1], z[0,1] result == 0.5
        self.assertTrue(round(A,10)== 0.5)
    def testNormInt(self): # Testing the integration fynction for y = 2*x at [0,1]. The area should be 1
        from proteus.WaveTools import normIntegral, reduceToIntervals, returnRectangles
        #pickin
        
        
        x = np.linspace(0,1,101)
        dx = 0.01
        xim = reduceToIntervals(x,dx)        
        y = 5*xim
        A  = normIntegral(y,xim)
        A = sum(returnRectangles(A,xim))
        self.assertTrue(round(A,10)== 1)
       
        
        
    def testEtaMode(self):
        from proteus.WaveTools import eta_mode
        x = 10.
        y =15.
        z = 2.
        t =50.
        kDir = [0.05,0.02,0.0]
        omega = 0.5
        phi = 3.14/5.
        amplitude =0.2
        eta = amplitude*cos(kDir[0]*x+kDir[1]*y+kDir[2]*z - omega*t +phi)
        self.assertTrue((eta - eta_mode(x,y,z,t,kDir,omega,phi,amplitude)==0.))# check eta
    def testVelMode(self): # Checking particle velocities
        from proteus.WaveTools import vel_mode

        kDir = np.array([2*pi,0.0,0.0])# Wavelength == 1
        omega = 2*pi
        phi =0.
        amplitude = 1.         
        g = np.array( [0,0.,-9.81])
        depth = 2.
        mwl =3.
        x=  pi/4./kDir[0]
        y = 0.
        z= 2.
        vDir = np.array([0,0,1])
        t= 0.
        kAbs = 2*pi
        for i in range(4):
            U_x = vel_mode(x,y,z,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"x")
            U_y = vel_mode(x,y,z,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"y")
            U_z = vel_mode(x,y,z,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"z")
            x+= 0.25
            # Checking velocity signs with quadrants
            if i ==0:
                #1st Quadrant
               self.assertTrue((U_x > 0.) and (U_y == 0.) and (U_z >0.))
            if i ==1:
                #2nd Quadrant
               self.assertTrue((U_x < 0.) and (U_y == 0.) and (U_z > 0.))
            if i ==2:
                #3rd Quadrant
               self.assertTrue((U_x < 0.) and (U_y == 0.) and (U_z < 0.))
            if i ==3:
                #4th Quadrant
               self.assertTrue((U_x > 0.) and (U_y == 0.) and (U_z < 0.))
        #Checking that the code does not allow z to be outside (-d,0)
#Checking vertical coherency
# U_z = 0 at z = mwl-d
        self.assertTrue(vel_mode(x,y,1.,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"z")==0.)
    
    def testTophat(self):    
        from proteus.WaveTools import tophat
        a  = np.random.rand(100)
        filt = tophat(100,0.1)
        af = a*filt
        a[:10] = 0.
        a[-10:] =0.
        self.assertTrue( a.all() == af.all())

    def testcosTap(self):    
        from proteus.WaveTools import costap
        a  = np.random.rand(100)
        filt = costap(100,0.1)
        af = a*filt
        a[:10] = 0.5*(1.-np.cos(pi*np.linspace(0,9,10)/10.))
        a[-10:] =0.5*(1.-np.cos(pi*np.linspace(9,0,10)/10.))
        self.assertTrue( a.all() == af.all())
    def testDecomposeFFT(self):    
        from proteus.WaveTools import decompose_tseries
        dt = 0.01
        time = np.linspace(0,199,200 )
        eta = np.cos(0.2*pi*time + 0.2)
        nfft = len(time)
        # testing full decompistion
        dec = decompose_tseries(time,eta,1)
        time = np.linspace(0,199,len(dec[1]))
        eta = np.cos(0.2*pi*time + 0.2)
        rec = np.zeros(len(time),)
        
        for ii in range(len(dec[0])):
            rec[:]+=dec[1][ii]*np.cos(dec[0][ii]*time[:]+dec[2][ii])
        rec[:]+=dec[3]
        self.assertTrue( rec.all() == eta.all())


class TestWaveParameters(unittest.TestCase):
#Checking dispersion calculation for a predicted wavelenght of 5.00m
    def test_dispersion(self):
        from proteus.WaveTools import dispersion
        length = 2*pi/dispersion(2*pi/1.94,1.)
        self.assertTrue(abs(length - 5.)/5.<0.001)
        length = 2*pi/dispersion([2*pi/1.94,2*pi/1.94,],1.)
        length-=5.
        length/=5
        self.assertTrue( (all(length) <0.001) or  (all(length) > -0.001))
#Check  sigma
    def test_sigma(self):
        from proteus.WaveTools import sigma,JONSWAP
        omega0=0.01
        sigma0 = 0.07
        sigma1 = 0.09
        x = np.ones((3),'d')
        x[0] = 0.5*omega0
        x[1] = omega0
        x[2] = 2.*omega0
        sigma = sigma(x,omega0)
        self.assertTrue((sigma[0] == sigma0).all())
        self.assertTrue((sigma[1] == sigma0).all())   
        self.assertTrue((sigma[2] == sigma1).all())
    def test_Jonswap(self): #JONSWAP tests
# Test Jonswap spectrum without TMA modification
        from proteus.WaveTools import sigma, JONSWAP, dispersion
        import random
        f0 = random.random() + 1.
        f = np.linspace(f0/2.,2.*f0,100)
        sig = sigma(f,f0)
        gamma = 6.*random.random() + 1. 
        Hs = random.random()
        bj = 0.0624*(1.094 - 0.01915*log(gamma))/(0.23+0.0336*gamma-0.185/(1.9+gamma))
        r_exp = np.exp(-(f/f0 -1 )**2/(2.*sig**2))
        JON = (bj*(Hs**2)*(f0**4)/f**5)*np.exp(-1.25*(f0/f)**4)*(gamma**r_exp)
        JON2 = JONSWAP(f,f0,Hs,gamma,TMA=False, depth = None)
        
        JCOMP = JON2/JON
        self.assertTrue((np.around(JCOMP,10)==1).all())
        h = random.random()
# Checking failure mode
        with self.assertRaises(SystemExit) as cm:
            JON2 = JONSWAP(f,f0,Hs,gamma,TMA=True)
        self.assertEqual(cm.exception.code, 1)  
# Check TMA modification           
        k = dispersion(2*pi*f,h)
        TMA = np.tanh(k*h)*np.tanh(k*h)/(1.+2.*k*h/np.sinh(2*k*h))
        JON2 = JONSWAP(f,f0,Hs,gamma,TMA=True, depth=h)
        JCOMP = JON2/(TMA*JON)
        self.assertTrue((np.around(JCOMP,10)==1).all())
    def test_PM(self): #PM tests
        from proteus.WaveTools import PM_mod
        f0 = random.random() + 1.
        f = np.linspace(f0/2.,2.*f0,10)        
        Hs = random.random()
        g = 9.81
        S_PM = (5./16.) * Hs**2 * f0**4 / f**5 * np.exp(-5./4. * (f0/f)**4)
        S_PM2 =  PM_mod(f,f0,Hs)
        SCOMP = S_PM2/S_PM
        self.assertTrue((np.around(SCOMP,10)==1).all())
    def testCos2s(self):
        from proteus.WaveTools import cos2s
        f0 = random.random() + 1.
        f = np.linspace(f0/2.,2.*f0,10.)        
        thetas = np.linspace(-pi/2.,pi/2.,11)
        s = 10. + 10. * np.random.random()
        S_PM = np.zeros((len(thetas),len(f)),)
        for ii in range(len(thetas)):
            for jj in range(len(f)):
                S_PM[ii,jj]= np.cos(thetas[ii]/2.)**(2*s)
        S_PM2 =  cos2s(thetas,f,s)
        SCOMP = S_PM2/S_PM
        self.assertTrue(np.array_equal(S_PM,S_PM2))

    def testMitsuyasu(self):
        from proteus.WaveTools import mitsuyasu
        f0 = random.random() + 1.
        f = np.linspace(f0/2.,2.*f0,10.)        
        thetas = np.linspace(-pi/2.,pi/2.,11)
        s = 10 + 10. * np.random.random()
        ss = np.zeros(len(f),)
        ss = (f/f0)**5
        i = np.where(f>f0)[0][0]
        ss[i:] = (f[i:]/f0)**(-2.5)
        ss[:] *= s
        S_PM = np.zeros((len(thetas),len(f)),)
        for ii in range(len(thetas)):
            for jj in range(len(f)):
                S_PM[ii,jj]= np.cos(thetas[ii]/2.)**(2.*ss[jj])
        S_PM2 =  mitsuyasu(thetas,f,f0,s)
        self.assertTrue(np.array_equal(S_PM,S_PM2))
       
        

class CheckMonochromaticWavesFailures(unittest.TestCase):
    def testFailureModes(self):
        from proteus.WaveTools import MonochromaticWaves
#Failure 1: Give non existent waveType
        with self.assertRaises(SystemExit) as cm1:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([1,0,0]),wavelength=None,waveType="Error",Ycoeff = None, Bcoeff =None, meanVelocity = np.array([0.,0.,0.]),phi0 = 0.)   
        self.assertEqual(cm1.exception.code, 1)     
#Failure 2: Give gravity direction not vertical to wave direction
        with self.assertRaises(SystemExit) as cm2:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,0,-9.81]),wavelength=None,waveType="Linear",Ycoeff = None, Bcoeff =None, meanVelocity = np.array([0.,0.,0.]),phi0 = 0.)   
        self.assertEqual(cm2.exception.code, 1)     
# Failure 3: Give Fenton type without wavelength
        with self.assertRaises(SystemExit) as cm3:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([1,0,0]),wavelength=None,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =np.array([1.,1,1.]), meanVelocity = np.array([0.,0.,0.]),phi0 = 0.)   
        self.assertEqual(cm3.exception.code, 1) 
# Failure 4: Give Fenton type without YCoeff  
        with self.assertRaises(SystemExit) as cm4:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,1,0]),wavelength=5.,waveType="Fenton",Ycoeff = None, Bcoeff =np.array([1.,1,1.]), meanVelocity = np.array([0.,0.,0.]),phi0 = 0.)   
        self.assertEqual(cm4.exception.code, 1)   
# Failure 5: Give Fenton type without BCoeff  
        with self.assertRaises(SystemExit) as cm5:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,1.,0]),wavelength=5.,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =None, meanVelocity = np.array([0.,0.,0.]),phi0 = 0.)   
        self.assertEqual(cm5.exception.code, 1) 
  # Failure 6: Give meanVelocity a vector value but not with 3 components
        with self.assertRaises(SystemExit) as cm6:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,1,0]),wavelength=5.,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =np.array([1.,1,1.]), meanVelocity =np.array([0.,0.,0.,0.]) ,phi0 = 0.)   
        self.assertEqual(cm6.exception.code, 1) 
  # Success!: Give all parameters in correct form!
        a = MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,1,0]),wavelength=5.,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =np.array([1.,1,1.]), meanVelocity =np.array([0.,0.,0.]) ,phi0 = 0.)   
        self.assertTrue(None == None)

class VerifyMonoChromaticLinearWaves(unittest.TestCase):
    def testLinear(self):
        from proteus.WaveTools import MonochromaticWaves
        import random
# Wave direction, random in x,y plane
        period = 1 
        waveHeight = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        dir1 = 2*random.random() - 1 
        dir2 = 2*random.random() - 1 
        waveDir = np.array([dir1,dir2, 0])
        phi0 = random.random()*2.*pi        
        a = MonochromaticWaves(period,waveHeight,mwl,depth,g,waveDir,wavelength=None,waveType="Linear",Ycoeff = None, Bcoeff =None, meanVelocity = np.array([0.,0,0.]),phi0 = phi0)
        x = random.random()*200. - 100.
        y = random.random()*200. - 100.
        z = mwl - depth + random.random()*( depth)
        t =  random.random()*200. - 100.
        eta = a.eta(x,y,z,t)
        ux = a.u(x,y,z,t,"x")
        uy = a.u(x,y,z,t,"y")
        uz = a.u(x,y,z,t,"z")

        omega = 2.*pi/period
# dispersion and setDirVector are tested above
        from proteus.WaveTools import dispersion,setDirVector
        kw = dispersion(omega,depth,gAbs)
        normDir = setDirVector(waveDir)
        amp = 0.5 * waveHeight
# Flow equation from Wikipedia, Airy wave theory https://en.wikipedia.org/wiki/Airy_wave_theoryhttps://en.wikipedia.org/wiki/Airy_wave_theory
        etaRef = amp*cos(kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - omega * t +phi0)
        z0 = z - mwl
        uxRef = normDir[0]*amp*omega*cosh(kw*(z0+depth))*cos(kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - omega * t +phi0)/sinh(kw*depth)
        uyRef = normDir[1]*amp*omega*cosh(kw*(z0+depth))*cos(kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - omega * t +phi0)/sinh(kw*depth)
        uzRef = amp*omega*sinh(kw*(z0+depth))*sin(kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - omega * t +phi0)/sinh(kw*depth)
       
        self.assertTrue(round(eta,8) == round(etaRef,8) )
        self.assertTrue(round(ux,8) == round(uxRef,8) )
        self.assertTrue(round(uy,8) == round(uyRef,8) )
        self.assertTrue(round(uz,8) == round(uzRef,8) )
class VerifyMonoChromaticFentonWaves(unittest.TestCase):
#Fenton methodology equations at http://johndfenton.com/Papers/Fenton88-The-numerical-solution-of-steady-water-wave-problems.pdf
#http://johndfenton.com/Steady-waves/Fourier.html

    def testFenton(self):
        from proteus.WaveTools import MonochromaticWaves
        period = 1. 
        waveHeight = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        dir1 = 2*random.random() - 1 
        dir2 = 2*random.random() - 1 
        waveDir = np.array([dir1,dir2, 0])
        phi0 = random.random()*2.*pi 
        wl = 10.
        YC =  np.array([5.,4.,3.,2.,1.])
        BC =  np.array([1.,2.,3.,4.,5.]) 
        mv =  np.array([6.,7.,8.])
# Set-up of Y and B coeffs does not correspond to physical properties
        a = MonochromaticWaves(period,waveHeight,mwl,depth,g,waveDir,wavelength=wl,waveType="Fenton",Ycoeff = YC, Bcoeff = BC, meanVelocity = mv,phi0 = phi0) 
        x = random.random()*200. - 100.
        y = random.random()*200. - 100.
        z =  mwl - depth + random.random()*( depth)
        t =  random.random()*200. - 100.
        eta = a.eta(x,y,z,t)
        ux = a.u(x,y,z,t,"x")
        uy = a.u(x,y,z,t,"y")
        uz = a.u(x,y,z,t,"z")
        omega = 2.*pi/period
# setDirVector are tested above
        from proteus.WaveTools import setDirVector
        kw = 2*pi/wl
        z0 = z - mwl
        normDir = setDirVector(waveDir)       
        amp = 0.5 * waveHeight
        etaRef = 0.
        uxRef = 0.
        uyRef = 0.
        uzRef = 0.
        jj = 0
        uxRef= mv[0]
        uyRef= mv[1]
        uzRef= mv[2]

        for ii in range(len(YC)):
            jj+=1
            etaRef+=YC[ii]*cos(jj*kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-jj*omega*t +phi0)/kw
            uxRef += normDir[0]* np.sqrt(gAbs/kw)*jj*BC[ii]*cosh(jj*kw*(z0+depth)) *cos(jj*kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-jj*omega*t +phi0)/cosh(jj*kw*depth)
            uyRef += normDir[1]* np.sqrt(gAbs/kw)*jj*BC[ii]*cosh(jj*kw*(z0+depth)) *cos(jj*kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-jj*omega*t +phi0)/cosh(jj*kw*depth)
            uzRef +=  np.sqrt(gAbs/kw)*jj*BC[ii]*sinh(jj*kw*(z0+depth)) *sin(jj*kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-jj*omega*t +phi0)/cosh(jj*kw*depth)
            
#            uxRef+=  normDir[0]*amp*omega*cosh(kw*(z0+depth))*cos(kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - omega * t +phi0)/sinh(kw*depth)
#*jj*BC[ii]*normDir[0]*cosh(jj*kw*(z0+depth))*cos(jj*kw*(normDir[0]*x+normDir[1]*y+normDir[2]*z) - jj*omega * t +phi0)*tanh(jj*kw*depth)/sinh(jj*kw*depth)
        self.assertTrue(round(eta,8) == round(etaRef,8) )
        self.assertTrue(round(ux,8) == round(uxRef,8))
        self.assertTrue(round(uy,8) == round(uyRef,8))
        self.assertTrue(round(uz,8) == round(uzRef,8))
       
        
#========================================= RANDOM WAVES ======================================
 


class CheckRandomWavesFailures(unittest.TestCase):
    def testFailureModes(self):
        from proteus.WaveTools import RandomWaves
#Failure 1: Give a wrong spectrum name
        with self.assertRaises(SystemExit) as cm1:
            RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"blahblah", spectral_params= None )
        self.assertEqual(cm1.exception.code, 1 )     

#Failure 2: Give gravity direction not vertical to wave direction
        with self.assertRaises(SystemExit) as cm2:
            RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,0,-9.81]),100,2.,"JONSWAP", spectral_params=None )
        self.assertEqual(cm2.exception.code, 1 )
#Failure 3: Give random spectral parameters
        with self.assertRaises(SystemExit) as cm3:
            RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", spectral_params= {"random1": 3.3, "random2":True,"random3" : 10.}  )
        self.assertEqual(cm3.exception.code, 1)
#Failure 4: Give wrong type of phases
        with self.assertRaises(SystemExit) as cm4:
            RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", spectral_params= {"random1": 3.3, "random2":True,"random3" : 10.}, phi = 0.  )
        self.assertEqual(cm4.exception.code, 1)
#Failure 5: Give wrong number of phase 
        with self.assertRaises(SystemExit) as cm5:
            RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", spectral_params= {"random1": 3.3, "random2":True,"random3" : 10.}, phi = np.zeros(1,)  )
        self.assertEqual(cm5.exception.code, 1)

  # Success!: Give all parameters in correct form!
        RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,1,0]),100,2.,"JONSWAP", spectral_params=None )
        RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,1,0]),100,2.,"JONSWAP", spectral_params={"gamma": 3.3, "TMA":True,"depth": 10.} )
        RandomWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,1,0]),100,2.,"JONSWAP", spectral_params={"gamma": 3.3, "TMA":True,"depth": 10.}, phi = np.zeros(100, float) )
        self.assertTrue(None == None)
    
class VerifyRandomWaves(unittest.TestCase):
    def testRandom(self):
        from proteus.WaveTools import RandomWaves
        import random
        # Assinging a random value at a field and getting the expected output
        Tp = 1. 
        Hs = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        dir1 = 2*random.random() - 1 
        dir2 = 2*random.random() - 1 
        waveDir = np.array([dir1,dir2, 0])
        N = 100
        phi = np.random.rand(N)
        gamma = 1.2
        TMA = True
        spectName = "JONSWAP"
        bandFactor = 2.0

        a= RandomWaves(Tp,
                     Hs,
                     mwl,#m significant wave height
                     depth ,           #m depth
                     waveDir,
                     g,      #peak  frequency
                     N,
                     bandFactor,         #accelerationof gravity
                     spectName# random words will result in error and return the available spectra 
                   )
        x = random.random()*200. - 100.
        y = random.random()*200. - 100.
        z =  mwl - depth + random.random()*( depth)
        t =  random.random()*200. - 100.
        # Just loading functions
        eta = a.eta(x,y,z,t)
        ux = a.u(x,y,z,t,"x")
        uy = a.u(x,y,z,t,"y")
        uz = a.u(x,y,z,t,"z")

        # Testing with a specific phi array
        a= RandomWaves(
            Tp,
            Hs,
            mwl,#m significant wave height
            depth ,           #m depth
            waveDir,
            g,      #peak  frequency
            N,
            bandFactor,         #accelerationof gravity
            spectName, 
            spectral_params =  {"gamma": gamma, "TMA": TMA,"depth": depth}, 
            phi = phi# random words will result in error and return the available spectra 
    )
        eta = a.eta(x,y,z,t)
        ux = a.u(x,y,z,t,"x")
        uy = a.u(x,y,z,t,"y")
        uz = a.u(x,y,z,t,"z")


        # setDirVector are tested above
        from proteus.WaveTools import setDirVector, dispersion, reduceToIntervals, returnRectangles, JONSWAP
        fmin = 1./(Tp * bandFactor) 
        fmax = bandFactor/(Tp)
        fi = np.linspace(fmin,fmax,N)
        df = (fmax-fmin)/(N -1 )
        ki = dispersion(2*pi*fi,depth)
        z0 = z - mwl
        normDir = setDirVector(waveDir) 
        fim = reduceToIntervals(fi,df)
        Si_Jm = JONSWAP(fim,1./Tp,Hs,gamma,TMA, depth)
        ai = np.sqrt(2.*returnRectangles(Si_Jm,fim))
        omega = 2*pi*fi
        etaRef = 0.
        uxRef = 0.
        uyRef = 0.
        uzRef = 0.

        for ii in range(N):
            etaRef+=ai[ii]*cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[ii])
            uxRef += normDir[0]*ai[ii]*omega[ii] *cosh(ki[ii]*(z0+depth)) *cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[ii])/sinh(ki[ii]*depth)
            uyRef += normDir[1]*ai[ii]*omega[ii] *cosh(ki[ii]*(z0+depth)) * cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[ii])/sinh(ki[ii]*depth)
            uzRef +=  ai[ii]*omega[ii] *sinh(ki[ii]*(z0+depth)) * sin(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[ii])/sinh(ki[ii]*depth)
        

        self.assertTrue(round(eta,8) == round(etaRef,8) )
        self.assertTrue(round(ux,8) == round(uxRef,8))
        self.assertTrue(round(uy,8) == round(uyRef,8))
        self.assertTrue(round(uz,8) == round(uzRef,8))

class CheckMultiSpectraRandomWavesFailures(unittest.TestCase):
    def testFailureModes(self):
        

        from proteus.WaveTools import MultiSpectraRandomWaves
#Failure 1: Give parameters as float rather than list
        Tp = 1. 
        Hs = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        dir1 = 2*random.random() - 1 
        dir2 = 2*random.random() - 1 
        waveDir = np.array([dir1,dir2, 0])
        N = 100
        phi = np.random.rand(N)
        gamma = 1.2
        TMA = True
        spectName = "JONSWAP"
        bandFactor = 2.0
# Give wrong length of list
        with self.assertRaises(SystemExit) as cm1:
            MultiSpectraRandomWaves(
            2,
            Tp,
            [Hs,Hs],
            mwl,#m significant wave height
            depth ,           #m depth
            [waveDir,waveDir],
            g,      #peak  frequency
            np.array([N,N]),
            [bandFactor,bandFactor],         #accelerationof gravity
            [spectName, spectName], 
            spectral_params =[  {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth} ], 
            phi = [phi,phi]# random words will result in error and return the available spectra 
        )
        self.assertEqual(cm1.exception.code, 1 )     

        with self.assertRaises(SystemExit) as cm2:
            MultiSpectraRandomWaves(
            2,
                [Tp,Tp,Tp],
            [Hs,Hs],
            mwl,#m significant wave height
            depth ,           #m depth
            [waveDir,waveDir],
            g,      #peak  frequency
            np.array([N,N]),
            [bandFactor,bandFactor],         #accelerationof gravity
            [spectName, spectName], 
            spectral_params =[  {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth} ], 
            phi = [phi,phi]# random words will result in error and return the available spectra 
    )

        self.assertEqual(cm2.exception.code, 1 )     

        # Success!: Give all parameters in correct form!
        MultiSpectraRandomWaves(
            2,
            [Tp,Tp],
            [Hs,Hs],
            mwl,#m significant wave height
            depth ,           #m depth
            [waveDir,waveDir],
            g,      #peak  frequency
            np.array([N,N]),
            [bandFactor,bandFactor],         #accelerationof gravity
            [spectName, spectName], 
            spectral_params =[  {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth} ], 
            phi = [phi,phi]# random words will result in error and return the available spectra 
    )

        self.assertTrue(None == None)

class VerifyMultiSpectraRandomWaves(unittest.TestCase):
    def testMultiSpectraDoubleExact(self):
        from proteus.WaveTools import MultiSpectraRandomWaves, RandomWaves 
        Tp = 1. 
        Hs = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        dir1 = 2*random.random() - 1 
        dir2 = 2*random.random() - 1 
        waveDir = np.array([dir1,dir2, 0])
        N = 100
        phi = np.random.rand(N)
        gamma = 1.2
        TMA = True
        spectName = "JONSWAP"
        bandFactor = 2.0
        x = random.random()*200. - 100.
        y = random.random()*200. - 100.
        z =  mwl - depth + random.random()*( depth)
        t =  random.random()*200. - 100.
        # Testing with a specific phi array
        a= RandomWaves(
            Tp,
            Hs,
            mwl,#m significant wave height
            depth ,           #m depth
            waveDir,
            g,      #peak  frequency
            N,
            bandFactor,         #accelerationof gravity
            spectName, 
            spectral_params =  {"gamma": gamma, "TMA": TMA,"depth": depth}, 
            phi = phi
    )
        eta = a.eta(x,y,z,t)
        ux = a.u(x,y,z,t,"x")
        uy = a.u(x,y,z,t,"y")
        uz = a.u(x,y,z,t,"z")

#Doubling the spectral properties
        aa= MultiSpectraRandomWaves(
            2,
            [Tp,Tp],
            [Hs,Hs],
            mwl,#m significant wave height
            depth ,           #m depth
            [waveDir,waveDir],
            g,      #peak  frequency
            np.array([N,N]),
            [bandFactor,bandFactor],         #accelerationof gravity
            [spectName, spectName], 
            spectral_params =[  {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth} ], 
            phi = [phi,phi]# random words will result in error and return the available spectra 
    )
        eta2 = aa.eta(x,y,z,t)
        ux2 = aa.u(x,y,z,t,"x")
        uy2 = aa.u(x,y,z,t,"y")
        uz2 = aa.u(x,y,z,t,"z")
        
        self.assertTrue(round(2.*eta,8) == round(eta2,8))
        self.assertTrue(round(2.*ux,8) == round(ux2,8))
        self.assertTrue(round(2.*uy,8) == round(uy2,8))
        self.assertTrue(round(2.*uz,8) == round(uz2,8))
# Testing with 5 spectra
        aa= MultiSpectraRandomWaves(
            5,
            [Tp,Tp,Tp,Tp,Tp],
            [Hs,Hs,Hs,Hs,Hs],
            mwl,#m significant wave height
            depth ,           #m depth
            [waveDir,waveDir,waveDir,waveDir,waveDir],
            g,      #peak  frequency
            np.array([N,N,N,N,N]),
            [bandFactor,bandFactor,bandFactor,bandFactor,bandFactor],         #accelerationof gravity
            [spectName, spectName,spectName,spectName,spectName], 
            spectral_params =[  {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth}, {"gamma": gamma, "TMA": TMA,"depth": depth},  {"gamma": gamma, "TMA": TMA,"depth": depth}, {"gamma": gamma, "TMA": TMA,"depth": depth}  ], 
            phi = [phi,phi,phi,phi,phi]# random words will result in error and return the available spectra 
    )
        eta2 = aa.eta(x,y,z,t)
        ux2 = aa.u(x,y,z,t,"x")
        uy2 = aa.u(x,y,z,t,"y")
        uz2 = aa.u(x,y,z,t,"z")

        self.assertTrue(round(5.*eta,8) == round(eta2,8))
        self.assertTrue(round(5.*ux,8) == round(ux2,8))
        self.assertTrue(round(5.*uy,8) == round(uy2,8))
        self.assertTrue(round(5.*uz,8) == round(uz2,8))


class CheckDirectionalWaveFailures(unittest.TestCase):
    def testFailureModes(self):
        from proteus.WaveTools import DirectionalWaves
        DirectionalWaves(200,1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", "cos2s", spectral_params= None, spread_params = {"s":10}, phi = None, phiSymm = False  )
       # Putting a silly name
        with self.assertRaises(SystemExit) as cm1:
            DirectionalWaves(200,1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", "aaaargh", spectral_params= None, spread_params = {"s":10}, phi = None, phiSymm = False  )
        self.assertEqual(cm1.exception.code, 1 )     
       # Putting incorrect phi
        with self.assertRaises(SystemExit) as cm2:
            DirectionalWaves(200,1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", "cos2s", spectral_params= None, spread_params = {"s":10}, phi = np.zeros(15,), phiSymm = False  )
        self.assertEqual(cm2.exception.code, 1 )     
        #putting non existent parameters
        with self.assertRaises(SystemExit) as cm3:
            DirectionalWaves(200,1.,1.,0.,10.,np.array([0,0,1]),np.array([0,-9.81,0]),100,2.,"JONSWAP", "cos2s", spectral_params= None, spread_params = {"blah":10}, phi = None, phiSymm = False  )
        self.assertEqual(cm3.exception.code, 1 )     
        

            
class VerifyDirectionals(unittest.TestCase):
    def testDirectional(self):
        from proteus.WaveTools import DirectionalWaves
        import random
        # Assinging a random value at a field and getting the expected output
        Tp = 1. 
        Hs = 0.15
        mwl = 4.5
        depth = 0.9
        g = np.array([0,0,-9.81])
        gAbs = 9.81
        # Setting the angle to the first quadrant otherwise the acos function further below is confused
        theta0 = 2*pi* random.random()
        dir1 = cos(theta0)
        dir2 = sin(theta0)

        waveDir = np.array([dir1,dir2, 0])
        N = 5
        M = 10
        phi =  2.0*pi*np.random.rand(2*M+1,N)
        gamma = 1.2
        TMA = True
        spectName = "JONSWAP"
        spectral_params =  {"gamma": gamma, "TMA": TMA,"depth": depth}
        bandFactor = 2.0
        spreadName = "mitsuyasu"
        spread_params = {"f0": 1./Tp, "smax": 15.}
        phiSymm = False

        aa= DirectionalWaves(
            M,
            Tp,
            Hs,
            mwl,#m significant wave height
            depth ,           #m depth
            waveDir,
            g,      #peak  frequency
            N,
            bandFactor,         #accelerationof gravity
            spectName,
            spreadName,
            spectral_params,
            spread_params,
            phi,
            phiSymm
            )
        x = random.random()*200. - 100.
        y = random.random()*200. - 100.
        z =  mwl - depth + random.random()*( depth)
        t =  random.random()*200. - 100.
        eta = aa.eta(x,y,z,t)
        ux = aa.u(x,y,z,t,"x")
        uy = aa.u(x,y,z,t,"y")
        uz = aa.u(x,y,z,t,"z")

        # Testing with a specific phi array

        # setDirVector are tested above
        from proteus.WaveTools import setDirVector, dispersion, reduceToIntervals, returnRectangles3D, JONSWAP,mitsuyasu, normIntegral

        fmin = 1./(Tp * bandFactor) 
        fmax = bandFactor/(Tp)
        fi = np.linspace(fmin,fmax,N)
        thetas = np.linspace(theta0 - pi/2,theta0+pi/2,2*M+1)
        dth = pi/(2*M)
        df = (fmax-fmin)/(N - 1 )
        ki = dispersion(2*pi*fi,depth)
        waveDirs = np.zeros((2*M+1,3),)
        for jj in range(2*M+1):
            waveDirs[jj,:] = np.array([cos(thetas[jj]),sin(thetas[jj]),0])
        z0 = z - mwl
        fim = reduceToIntervals(fi,df)
        thetas-=theta0
        thetas_m = reduceToIntervals(thetas,dth)
        Si_Jm = JONSWAP(fim,1./Tp,Hs,gamma,TMA, depth)
        Si_dir = mitsuyasu(thetas_m,fim,1./Tp, 15.)
        for ii in range(0,N):            
            Si_dir[:,ii] = normIntegral(Si_dir[:,ii],thetas_m)
            Si_dir[:,ii]*= Si_Jm[ii] 
        ai = np.sqrt(2.*returnRectangles3D(Si_dir,thetas_m,fim))

        omega = 2*pi*fi
        etaRef = 0.
        uxRef = 0.
        uyRef = 0.
        uzRef = 0.

        for ii in range(N):
            for jj in range(2*M+1):
                normDir = waveDirs[jj,:]
                etaRef+=ai[jj,ii]*cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[jj,ii])
                uxRef += normDir[0]*ai[jj,ii]*omega[ii] *cosh(ki[ii]*(z0+depth)) *cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[jj,ii])/sinh(ki[ii]*depth)
                uyRef += normDir[1]*ai[jj,ii]*omega[ii] *cosh(ki[ii]*(z0+depth)) * cos(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[jj,ii])/sinh(ki[ii]*depth)
                uzRef +=  ai[jj,ii]*omega[ii] *sinh(ki[ii]*(z0+depth)) * sin(ki[ii]*(normDir[0]*x+normDir[1]*y+normDir[2]*z)-omega[ii]*t +phi[jj,ii])/sinh(ki[ii]*depth)
        
class CheckTimeSeriesFailureModes(unittest.TestCase):
    def testTimeSeriesFailureModes(self):
        from proteus.WaveTools import TimeSeries
        #load successfully - direct decomposition
        path = getpath()
        aa= TimeSeries(
            path+"test_timeSeries.txt",
            0,
            np.array([1.,0,0]), 
            1.  ,
            64 ,          #number of frequency bins
            1. ,        
            np.array([1.,0,0]), 
            np.array([0,0,-9.81])
            )

        with self.assertRaises(SystemExit) as cm1:
            aa= TimeSeries(
                path+"test_timeSeries.dat",
                0,
                np.array([1.,0,0]), 
                1.  ,
               64 ,          #number of frequency bins
                1. ,        
                np.array([1,0,0]), 
                np.array([0,0,-9.81])
            )
        with self.assertRaises(SystemExit) as cm2:
            aa= TimeSeries(
                path+"test_timeSeries_err1.csv",
                0,
                np.array([1.,0,0]), 
                1.  ,
                64 ,          #number of frequency bins
                1. ,        
                np.array([1,0,0]), 
                np.array([0,0,-9.81])
            )
        self.assertEqual(cm2.exception.code, 1 )     
        with self.assertRaises(SystemExit) as cm3:
            aa= TimeSeries(
                path+"test_timeSeries_err2.txt",
                0,
                np.array([1.,0,0]), 
                1.  ,
                64 ,          #number of frequency bins
                1. ,        
                np.array([1,0,0]), 
                np.array([0,0,-9.81])
            )
        self.assertEqual(cm3.exception.code, 1 )     


        with self.assertRaises(SystemExit) as cm4:
            aa= TimeSeries(
            path+"test_timeSeries.txt",
            0,
            np.array([1.,0,0]), 
            1.  ,
            64 ,          #number of frequency bins
            1. ,        
            np.array([1,0,0]), 
            np.array([1,0,0])
            )
        self.assertEqual(cm4.exception.code, 1 )     

        with self.assertRaises(SystemExit) as cm5:
            aa= TimeSeries(
            path+"test_timeSeries.txt",
            0,
            np.array([0,1.,0,0]), 
            1.  ,
            64 ,          #number of frequency bins
            1. ,        
            np.array([1,0,0]), 
            np.array([0,0,-9.81])
            )
        self.assertEqual(cm5.exception.code, 1 )     
        
class VerifyTimeSeries(unittest.TestCase):
    def testDirect(self):
        path =getpath()
        from proteus.WaveTools import TimeSeries, costap
        import random
        aa= TimeSeries(
            path+"test_timeSeries.txt",
            0,
             np.array([0.,0.,0]),            
            1.  ,
            256,          #number of frequency bins
            1. ,        
            np.array([1,0,0]), 
            np.array([0,0,-9.81])
            )
        fid = open(path+"test_timeSeries.txt","r")
        data = np.loadtxt(fid)
        timeRef = data[:,0]
        etaRef = data[:,1]
        
        timeInt = timeRef #np.linspace(timeRef[0],timeRef[-1],len(timeRef))
        etaInt = etaRef #np.interp(timeInt, timeRef, etaRef)
        etaTest = np.zeros(len(timeRef),"d")
        x = 0.
        y = 0.
        z = 0.
        ii = -1
        for tt in timeInt:
            ii+=1
            etaTest[ii] = aa.etaDirect(x,y,z,tt) 

        etaInt-=np.mean(etaInt)
        etaInt*=costap(len(data),0.025)
        norm = max(etaRef)
        err = (etaInt - etaTest)**2
        err = np.sqrt(sum(err))/len(etaInt)/max(abs(etaInt))
        self.assertTrue(err<1e-3 )     
        
        
        


if __name__ == '__main__':
    unittest.main(verbosity=2)
