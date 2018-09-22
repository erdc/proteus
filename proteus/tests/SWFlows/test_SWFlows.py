#!/usr/bin/env python
"""
Test module for SWFlows
"""
import os
import pytest
import tables
import numpy as np

class TestSWFlows(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        
    def test_solitary(self):
        # call runSWEs
        os.system("parun --SWEs solitary.py")
        
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/solitary.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('solitary.h5','r')
        assert np.allclose(expected.root.h_t3,actual.root.h_t3,atol=1e-10)
        expected.close()
        actual.close()
        
    def test_parab1D(self):
        # Call runSWEs
        os.system("parun --SWEs parab1D.py")

        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/parab1D.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('parab1D.h5','r')
        assert np.allclose(expected.root.h_t11,actual.root.h_t11,atol=1e-10)
        expected.close()
        actual.close()
