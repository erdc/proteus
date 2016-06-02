#!usr/bin/env python
""" Test module for null space """

import os
import numpy
import tables
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose  = True
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace
from petsc4py import PETSc as p4pyPETSc

class TestStokesDrivenCavity():

    def setUp(self):
        """ This function tests the driven cavity with a boundary null space """
        import stokesDrivenCavity_2d_p
        import stokesDrivenCavity_2d_n
        pList = [stokesDrivenCavity_2d_p]
        nList = [stokesDrivenCavity_2d_n]
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList = [default_s]
        opts.verbose = True
        opts.profile = True
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)

    def tearDown(self):
        os.remove('drivenCavityStokesTrial.h5')
        os.remove('drivenCavityStokesTrial.xmf')

    def test_01_RunsWithNullSpace(self):
        self.ns.calculateSolution('testNullSpace')
        set_trace()
        test_path = os.path.dirname(os.path.abspath(__file__))
        expected = tables.openFile(os.path.join(test_path,
                                   'drivenCavityStokesTrial_expected.h5'),
                                   'r')
        actual = tables.openFile('drivenCavityStokesTrial.h5','r')
        assert numpy.allclose(expected.root.velocity_t1,
                              actual.root.velocity_t1)
        expected.close()
        actual.close()
            
if __name__=='__main__':
    from proteus import Comm
    comm = Comm.init()
