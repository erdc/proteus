#!/usr/bin/env python
"""
Test module for SWFlow
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

class TestSWFlow(object):

    def setup_method(self,method):
        self._scriptdir = os.path.dirname(__file__)
        self.path = proteus.__path__[0]+"/tests/SWFlow/"

    def compare_vs_saved_files(self,name):
        expected_path = 'comparison_files/' + name + '.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(name+'.h5','r')
        assert np.allclose(expected.root.h_t2,actual.root.h_t2,atol=1e-8)
        expected.close()
        actual.close()

    def test_solitary_wave(self):
        # call runSWEs
        os.system("parun --SWEs -l1 -v solitary_wave.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("solitary_wave")

    def test_parab1D(self):
        # Call runSWEs
        os.system("parun --SWEs -l1 -v parab1D.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("parab1D")

    def test_dam3Bumps(self):
        # Call runSWEs
        os.system("parun --SWEs -l1 -v dam3Bumps.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("dam3Bumps")

    def test_GN_steady(self):
        os.system("parun --SWEs -l1 -v GN_steady.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("GN_steady")

    def test_solitary_reef(self):
        os.system("parun --SWEs -l1 -v solitary_reef.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("solitary_reef")

    def test_seawall(self):
        os.system("parun --SWEs -l1 -v seawall.py -C 'refinement=4 final_time=0.1'")
        self.compare_vs_saved_files("seawall")
