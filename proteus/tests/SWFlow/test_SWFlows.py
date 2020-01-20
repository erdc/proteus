#!/usr/bin/env python
"""
Test module for SWFlow
"""
import os
import pytest
import tables
import numpy as np

class TestSWFlow(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)

    def test_solitary(self):
        # call runSWEs
        os.system("PYTHONPATH={0} parun --SWEs solitary.py".format(self._scriptdir))

        # COMPARE VS SAVED FILES #
        #expected_path = 'comparison_files/solitary.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('solitary.h5','r')
        #assert np.allclose(expected.root.h_t3,actual.root.h_t3,atol=1e-10)
        #expected.close()

        expected_path = 'comparison_files/' + 'comparison_' + 'solitary' + '_h_t3.csv'
        #write comparison file
        #np.array(actual.root.h_t3).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.h_t3).flatten(),decimal=10)
        actual.close()

    def test_parab1D(self):
        # Call runSWEs
        os.system("PYTHONPATH={0} parun -v  --SWEs parab1D.py".format(self._scriptdir))

        # COMPARE VS SAVED FILES #
        actual = tables.open_file('parab1D.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + 'parab1D' + '_h_t11.csv'
        #write comparison file
        #np.array(actual.root.h_t11).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.h_t11).flatten(),decimal=10)


        actual.close()

    def test_dam3Bumps(self):
        # Call runSWEs
        os.system("PYTHONPATH={0} parun -v  --SWEs dam3Bumps.py".format(self._scriptdir))

        # COMPARE VS SAVED FILES #
        actual = tables.open_file('dam3Bumps.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + 'dam3Bumps' + '_h_t4.csv'
        #write comparison file
        #np.array(actual.root.h_t4).tofile(os.path.join(self._scriptdir, expected_path),sep=",")
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.h_t4).flatten(),decimal=10)
        actual.close()
