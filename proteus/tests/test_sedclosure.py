from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest

comm = Comm.init()
Profiling.procID = comm.rank()

Profiling.logEvent("Testing SedClosure")


class TestHsu(unittest.TestCase):
    def testGradularDrag1(self):
        from proteus.mprans.SedClosure import HsuSedStress
        # Constant params
        aDarcy = 1.
        bForch = 1.
        grain = 0.1
        packFraction = 0.2
        packMargin = 0.01
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))

        # Testing for sedF >> packFraction
        nu = 1.
        sedF = 0.5
        drag = aDarcy * nu* sedF /((1.-sedF)*grain**2) +  bForch * umag / grain
        sedSt = HsuSedStress( aDarcy, bForch, grain, packFraction, packMargin)
        # For sedF > pacFraction - > drag = a * nu* sedF /((1-sedF)*grain^2) + beta * umag / grain
        drag2 = sedSt.granularDrag(sedF, uf, us, nu)

        self.assertTrue(drag == drag2)
    def testGradularDrag2(self):
        from proteus.mprans.SedClosure import HsuSedStress
        aDarcy = 1.
        bForch = 1.
        grain = 0.1
        packFraction = 0.2
        packMargin = 0.01
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF << packFraction and Rep = (1. - sed) * umag*nu/grain = 0.9 * 5. * 0.1 / 1. = 0.45
        nu = 1.
        sedF = 0.1
        Rep = (1.- sedF)*umag*grain / nu    
        drag =  ( 24. * (1.+0.15*Rep**(0.687))/Rep) * 0.75 * umag * (1. -sedF)**(-1.65) / grain # Chen and Hsu 2014
        sedSt = HsuSedStress(aDarcy, bForch, grain, packFraction, packMargin)
        drag2 = sedSt.granularDrag(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))

    def testGradularDrag3(self):
        from proteus.mprans.SedClosure import HsuSedStress
        # Constant params
        aDarcy = 1.
        bForch = 1.
        grain = 0.1
        packFraction = 0.2
        packMargin = 0.01
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF << packFraction and Rep > 1000
        sedF = 0.1
        nu = 1e-4
        Rep = (1.- sedF) * umag * grain / nu    
        drag =  ( 0.44 * 0.75 * umag * (1. -sedF)**(-1.65) )/ grain # Chen and Hsu 2014
        sedSt = HsuSedStress( aDarcy, bForch, grain, packFraction, packMargin)
        drag2 = sedSt.granularDrag(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))
    def testGradularDrag4(self):
        from proteus.mprans.SedClosure import HsuSedStress
        # Constant params
        aDarcy = 1.
        bForch = 1.
        grain = 0.1
        packFraction = 0.2
        packMargin = 0.01
        uf = np.array([5.,4.],"d")
        us = np.array([1.,1.],"d") 
        umag_temp = (uf - us)*(uf - us)
        umag = np.sqrt(sum(umag_temp))
        # Testing for sedF =  packFraction +0.5 packmargin and Rep > 1000
        sedF = 0.205
        nu = 1e-4
        Rep = (1.- sedF) * umag * grain / nu
        draga = aDarcy * nu* sedF /((1.-sedF)*grain**2) +  bForch * umag / grain
    
        dragb =  ( 0.44 * 0.75 * umag * (1. -sedF)**(-1.65) )/ grain # Chen and Hsu 2014
        w = 0.5 + (sedF - packFraction) / (2. * packMargin)
        drag = w*draga + (1.-w) * dragb
        sedSt = HsuSedStress( aDarcy, bForch, grain, packFraction, packMargin)
        drag2 = sedSt.granularDrag(sedF, uf, us, nu)
        self.assertTrue(round(drag,10) == round(drag2,10))



if __name__ == '__main__':
    unittest.main(verbosity=2)
