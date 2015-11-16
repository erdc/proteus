from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import unittest

comm = Comm.init()
Profiling.procID = comm.rank()

class TestSigma(unittest.TestCase):
    def test_lower(self):#idea was  just to test what  happens  when omega is  constant and < omega0
        from proteus.WaveTools import sigma
        omega0=0.01
        sigma0 = 0.07
        x = np.ones((10,),'d')
        x *= omega0
        sigma = sigma(x,omega0)
        self.assertTrue((sigma == sigma0).all())

if __name__ == '__main__':
    unittest.main()
