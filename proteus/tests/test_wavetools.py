from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest
import random
from math import cos,sin,cosh,sinh,pi

comm = Comm.init()
Profiling.procID = comm.rank()


Profiling.logEvent("Testing WaveTools")
class TestModes(unittest.TestCase):
    def testVDir(self):
        from proteus.WaveTools import setVertDir
        self.assertTrue(np.array_equal(setVertDir(np.array([0,-9.81,0])), np.array([0,1,0])))
    def testDirVector(self):
        from proteus.WaveTools import setDirVector
        self.assertTrue(all(setDirVector(np.array([2.,2.,1.]))== np.array([2.,2.,1])/3.))
    def dirCheck(self):
        from proteus.WaveTools import dirCheck
        with self.assertRaises(SystemExit) as cm:
            dirCheck(np.array([9,9,9]),np.array([4,5,6]))
        self.assertEqual(cm.exception.code, 1)     
        self.assertTrue(dirCheck(np.array([1.,2.,3.]),np.array([7.,4.,-5.])==0))

            
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
        zB = np.array([0.5,3.5])
        for z in zB:
            with self.assertRaises(SystemExit) as cm:
                vel_mode(x,y,z,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"x")
            self.assertEqual(cm.exception.code, 1)       
#Checking vertical coherency
# U_z = 0 at z = mwl-d
        self.assertTrue(vel_mode(x,y,1.,t,kDir,kAbs,omega,phi,amplitude,mwl,depth,g,vDir,"z")==0.)
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

#Check lower sigma
    def test_lower(self):#idea was  just to test what  happens  when omega is  constant and < omega0
        from proteus.WaveTools import sigma
        omega0=0.01
        sigma0 = 0.07
        x = np.ones((10,),'d')
        x *= omega0
        sigma = sigma(x,omega0)
        self.assertTrue((sigma == sigma0).all())

class TestMonochromaticWaves(unittest.TestCase):
    def testFailureModes(self):
        from proteus.WaveTools import MonochromaticWaves
#Failure 1: Give non existent waveType
        with self.assertRaises(SystemExit) as cm1:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([1,0,0]),wavelength=None,waveType="Error",Ycoeff = None, Bcoeff =None, meanVelocity = 0.,phi0 = 0.)   
        self.assertEqual(cm1.exception.code, 1)     
#Failure 2: Give gravity direction parallel to wave direction
        with self.assertRaises(SystemExit) as cm2:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,1]),np.array([0,0,-9.81]),wavelength=None,waveType="Linear",Ycoeff = None, Bcoeff =None, meanVelocity = 0.,phi0 = 0.)   
        self.assertEqual(cm2.exception.code, 1)     
# Failure 3: Give Fenton type without wavelength
        with self.assertRaises(SystemExit) as cm3:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([1,0,0]),wavelength=None,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =np.array([1.,1,1.]), meanVelocity = 0.,phi0 = 0.)   
        self.assertEqual(cm3.exception.code, 1) 
# Failure 4: Give Fenton type without YCoeff  
        with self.assertRaises(SystemExit) as cm4:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,0,1]),wavelength=5.,waveType="Fenton",Ycoeff = None, Bcoeff =np.array([1.,1,1.]), meanVelocity = 0.,phi0 = 0.)   
        self.assertEqual(cm4.exception.code, 1)   
# Failure 5: Give Fenton type without BCoeff  
        with self.assertRaises(SystemExit) as cm5:
            MonochromaticWaves(1.,1.,0.,10.,np.array([0,0,-9.81]),np.array([0,0,1]),wavelength=5.,waveType="Fenton",Ycoeff = np.array([1.,1,1.]), Bcoeff =None, meanVelocity = 0.,phi0 = 0.)   
        self.assertEqual(cm5.exception.code, 1)   




if __name__ == '__main__':
    unittest.main()
