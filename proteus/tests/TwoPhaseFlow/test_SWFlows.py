#!/usr/bin/env python
"""
Test module for TwoPhaseFlow
"""
import os
import pytest
import tables
import numpy as np

class TestTwoPhaseFlow(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)

    def test_fallingBubble(self):
        # call runSWEs
        os.system("parun --TwoPhaseFlow -f fallingBubble.py -C 'final_time=0.1 dt_output=0.1 refinement=2'")
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/fallingBubble.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('fallingBubble.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()

    def test_marin(self):
        # call runSWEs
        os.system("parun --TwoPhaseFlow -f marin.py -C 'final_time=0.1 dt_output=0.1'")
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/marin.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('marin.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()

    def test_quiescentTank(self):
        # call runSWEs
        os.system("parun --TwoPhaseFlow -f quiescentTank.py -C 'final_time=0.1 dt_output=0.1 refinement=6'")
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/quiescentTank.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('quiescentTank.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()

    def test_risingBubble(self):
        # call runSWEs
        os.system("parun --TwoPhaseFlow -f risingBubble.py -C 'final_time=0.1 dt_output=0.1 refinement=2'")
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/risingBubble.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('risingBubble.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()

    def test_TwoDimBucklingFlow(self):
        # call runSWEs
        os.system("parun --TwoPhaseFlow -f TwoDimBucklingFlow.py -C 'final_time=0.1 dt_output=0.1 refinement=4'")
        # COMPARE VS SAVED FILES #
        expected_path = 'comparison_files/TwoDimBucklingFlow.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file('TwoDimBucklingFlow.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-10)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-10)
        expected.close()
        actual.close()
