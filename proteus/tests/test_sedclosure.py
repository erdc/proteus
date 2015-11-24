from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest

comm = Comm.init()
Profiling.procID = comm.rank()

Profiling.logEvent("Testing SedClosure")
class TestHsu(unittest.TestCase):
    def testVDir(self):
        from proteus.mprans.SedClosure import HsuSedStress
        closure = HsuSedStress(3.0)
        self.assertTrue(closure.M_sf_x(2.0,3.0) == 21.0)

if __name__ == '__main__':
    unittest.main(verbosity=2)
