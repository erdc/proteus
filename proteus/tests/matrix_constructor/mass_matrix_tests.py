#!/usr/bin/env python
"""
Test module for 2D Quadrilateral Meshes
"""

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

class TestMassConstructionStokes2D():
    """ Verifty construction of 2D Mass Matrix from a Stokes Problem. """

    def __init__(self):
        pass

    def setUp(self):
        """ Initialize the test problem """
        import stokes_2d_p
        import stokes_2d_n
        pList = [stokes_2d_p]
        nList = [stokes_2d_n]
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList = [default_s]
        opts.verbose = True
        opts.profile = True
        self.ns = NumericalSolution.NS_base(so,
                                            pList,
                                            nList,
                                            so.sList,
                                            opts)

    def tearDown(self):
        """ Tear down function """
        FileList = ['Mass_matrix_test.xmf',
                    'Mass_matrix_test.h5',
                    'rdomain.ele',
                    'rdomain.neig',
                    'rdomain.node',
                    'rdomain.edge',
                    'rdomain.poly']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_1(self):
        mm = self.ns.modelList[0].levelModelList[1]
        mm.attachMassOperator()
        mass_mat = LinearAlgebraTools.superlu_sparse_2_dense(mm.MassOperator)
        #        numpy.savetxt('mass_mat_stokes_ex.txt',mass_mat,fmt='%1.5f')
        comparison_mat = numpy.loadtxt('mass_mat_stokes_ex.txt')
        assert numpy.allclose(mass_mat,comparison_mat,rtol=1e-03)

if __name__ == '__main__':
    tt = TestMassConstructionStokes2D()
    tt.setUp()
    tt.test_1()
    tt.tearDown()
