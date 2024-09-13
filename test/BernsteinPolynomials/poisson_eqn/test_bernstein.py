#!/usr/bin/env python
"""
Test module for Berstein basis FE
"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=1
Profiling.verbose=True
import os
import numpy as np
import h5py
import pytest
from . import (parameters_for_poisson, poisson_p, poisson_n)

class TestBernstein(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        reload(parameters_for_poisson)
        self.pList = [poisson_p]
        self.nList = [poisson_n]
        self.sList = [default_s]
        self.so = default_so
        self._scriptdir = os.path.dirname(__file__)
        self.sim_names = []
        self.aux_names = []

    def teardown_method(self,method):
        pass

    @pytest.mark.skipif(sys.platform == "darwin", reason="does not run on macOS")
    def test_2D_hex(self):
        # Set parameters for this test
        parameters_for_poisson.ct.nd = 2
        parameters_for_poisson.ct.useHex = True
        # Reload _p and _n modules
        reload(poisson_p)
        reload(poisson_n)
        poisson_n.conservativeFlux=None
        # Update name
        self.so.name = "2D_"+self.pList[0].name+"_hex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/'+self.so.name+'.h5'
        #expected = h5py.File(os.path.join(self._scriptdir,expected_path))
        #actual = h5py.File(self.so.name+'.h5','r')
        #assert np.allclose(expected.root.u0_t1,
        #                   actual['u0_t1'],
        #                   atol=1e-10)
        #expected.close()

        actual = h5py.File(self.so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + self.so.name + '_u0_t1.csv'
        #write comparison file
        #np.array(actual['u0_t1']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u0_t1']).flatten(),decimal=10)
        actual.close()

    def test_2D_simplex(self):
        # Set parameters for test
        parameters_for_poisson.ct.nd = 2
        parameters_for_poisson.ct.useHex = False
        # Reload _p and _n modules
        reload(poisson_p)
        reload(poisson_n)
        poisson_p.domain.MeshOptions.nnz=1
        # Update name
        self.so.name = "2D_"+self.pList[0].name+"_simplex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/'+self.so.name+'.h5'
        #expected = h5py.File(os.path.join(self._scriptdir,expected_path))
        #actual = h5py.File(self.so.name+'.h5','r')
        #assert np.allclose(expected.root.u0_t1,
        #                   actual['u0_t1'],
        #                   atol=1e-10)
        actual = h5py.File(self.so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + self.so.name + '_u0_t1.csv'
        #write comparison file
        #np.array(actual['u0_t1']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u0_t1']).flatten(),decimal=10)
        actual.close()

    def test_3D_hex(self):
        # Set parameters for test
        parameters_for_poisson.ct.nd = 3
        parameters_for_poisson.useHex = True
        # Reload _p and _n modules
        reload(poisson_p)
        reload(poisson_n)
        poisson_n.conservativeFlux=None
        # Update name
        self.so.name = "3D_"+self.pList[0].name+"_hex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/'+self.so.name+'.h5'
        #expected = h5py.File(os.path.join(self._scriptdir,expected_path))
        #actual = h5py.File(self.so.name+'.h5','r')
        #assert np.allclose(expected.root.u0_t1,
        #                   actual['u0_t1'],
        #                   atol=1e-10)
        actual = h5py.File(self.so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + self.so.name + '_u0_t1.csv'
        #write comparison file
        #np.array(actual['u0_t1']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u0_t1']).flatten(),decimal=10)
        actual.close()

    def test_3D_simplex(self):
        # Set parameters for test
        parameters_for_poisson.ct.nd = 3
        parameters_for_poisson.ct.useHex = False
        # Reload _p and _n modules
        reload(poisson_p)
        reload(poisson_n)
        # Update name
        self.so.name = "3D_"+self.pList[0].name+"_simplex_degree2"
        # NUMERICAL SOLUTION #
        ns = proteus.NumericalSolution.NS_base(self.so,
                                               self.pList,
                                               self.nList,
                                               self.sList,
                                               opts)
        self.sim_names.append(ns.modelList[0].name)
        ns.calculateSolution('poisson')
        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/'+self.so.name+'.h5'
        #expected = h5py.File(os.path.join(self._scriptdir,expected_path))
        #actual = h5py.File(self.so.name+'.h5','r')
        #assert np.allclose(expected.root.u0_t1,
        #                   actual['u0_t1'],
        #                   atol=1e-10)
        actual = h5py.File(self.so.name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + self.so.name + '_u0_t1.csv'
        #write comparison file
        #np.array(actual['u0_t1']).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual['u0_t1']).flatten(),decimal=10)
        actual.close()
