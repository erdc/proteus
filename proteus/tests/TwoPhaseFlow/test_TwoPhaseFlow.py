#!/usr/bin/env python
"""
Test module for TwoPhaseFlow
"""
import pytest
import tables
import numpy as np
import proteus.defaults
from proteus import Context
from proteus import default_so
from proteus.iproteus import *
import os
import sys
Profiling.logLevel=1
Profiling.verbose=True

class TestTwoPhaseFlow(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.path = proteus.__path__[0]+"/tests/TwoPhaseFlow/"

    def compare_vs_saved_files(self,name):
        expected_path = 'comparison_files/' + name + '.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(name+'.h5','r')
        assert np.allclose(expected.root.phi_t2,actual.root.phi_t2,atol=1e-8)
        assert np.allclose(expected.root.velocity_t2,actual.root.velocity_t2,atol=1e-8)
        expected.close()
        actual.close()

    def test_fallingBubble(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "-f fallingBubble.py -l1 -v -C 'final_time=0.1 dt_output=0.1 refinement=2'")
        self.compare_vs_saved_files("fallingBubble")

    def test_marin(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "-f marin.py -l1 -v -C 'final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("marin")

    def test_quiescentTank(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "-f quiescentTank.py -l1 -v -C 'final_time=0.1 dt_output=0.1 refinement=6'")
        self.compare_vs_saved_files("quiescentTank")

    def test_risingBubble(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "-f risingBubble.py -l1 -v -C 'final_time=0.1 dt_output=0.1 refinement=2'")
        self.compare_vs_saved_files("risingBubble")

    def test_TwoDimBucklingFlow(self):
        os.system("parun --TwoPhaseFlow --path " + self.path + " "
                  "-f TwoDimBucklingFlow.py -l1 -v -C 'final_time=0.1 dt_output=0.1 refinement=4'")
        self.compare_vs_saved_files("TwoDimBucklingFlow")
