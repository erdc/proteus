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
        #expected_path = 'comparison_files/' + name + '.h5'
        #expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        #actual = tables.open_file(name+'.h5','r')
        #assert np.allclose(expected.root.h_t2,actual.root.h_t2,atol=1e-8)
        #expected.close()
        #actual.close()

        actual = tables.open_file(name+'.h5','r')
        expected_path = 'comparison_files/' + 'comparison_' + name + '_h_t2.csv'
        #write comparison file
        # write_path = './comparison_files/' + 'comparison_' + name + '_h_t2.csv'
        # np.array(actual.root.h_t2).tofile(os.path.join(self._scriptdir, write_path),sep=",")
        #
        np.testing.assert_almost_equal(np.fromfile(os.path.join(self._scriptdir, expected_path),sep=","),np.array(actual.root.h_t2).flatten(),decimal=7)
        actual.close()

    def test_solitary_wave(self):
        # call runSWEs
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v solitary_wave.py -C 'refinement=3 final_time=0.1'")
        self.compare_vs_saved_files("solitary_wave")

    def test_parab1D(self):
        # Call runSWEs
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v parab1D.py -C 'refinement=3 final_time=10.0 dt_output=10.0'")
        self.compare_vs_saved_files("parab1D")

    def test_dam3Bumps(self):
        # Call runSWEs
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v dam3Bumps.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("dam3Bumps")

    def test_dSWEs_steady_state(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v dSWEs_steady_state.py  -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("dSWEs_steady_state")

    def test_reef_island_runup(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v reef_island_runup.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("reef_island_runup")

    def test_seawall(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v seawall.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("seawall")

    def test_conical_island(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v  conical_island.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("conical_island")

    def test_transcritical_bump(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v transcritical_bump.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("transcritical_bump")

    def test_obstacle_flow(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v obstacle_flow.py -C 'he=4.0 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("obstacle_flow")

    def test_santos_step(self):
        os.system("parun --SWEs --path " + self.path + " "
                  "-l1 -v santos_step.py -C 'refinement=3 final_time=0.1 dt_output=0.1'")
        self.compare_vs_saved_files("santos_step")
