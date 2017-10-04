#!/usr/bin/env python
"""

Test module for BDM2 Elements

"""
import os,sys,inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0,cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() )) [0],"import_modules")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0,cmd_subfolder)

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy as np
import pytest
from post_processing.import_modules import bdm_tests_template_mesh8 as bt

@pytest.mark.PostProcessingTools
class TestBDM2Mesh8():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass


    def setup_method(self,method):
        """Initialize the test problem. """
        reload(bt)
        self.transport_obj = bt.ns.modelList[0].levelModelList[0]
        self.bdm2_obj = self.transport_obj.velocityPostProcessor.vpp_algorithms[0]
        self._setRelativePath()

    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)
        
    def teardown_method(self,method):
        filenames = ['poisson_bdm1_test.h5', 'poisson_bdm1_test.xmf','reference_triangle.ele',
                     'reference_triangle.node', 'reference_triangle.poly','proteus.log',
                     'blockDomain.poly']
        for file in filenames:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename, e.strerror ))
            else:
                pass

    def test_BDM2_reference_triangle_full_in_space(self):
        rel_path_1 = "comparison_files/bdm_bdy_func_values_mesh_8.npy"
        rel_path_2 = "comparison_files/bdm_func_values_mesh_8.npy"
        bdm_bdy_values = np.load(os.path.join(self.scriptdir,rel_path_1))
        bdm_values = np.load(os.path.join(self.scriptdir,rel_path_2))

        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values.copy()
        self.bdm2_obj.q[('velocity',0)] = bdm_values.copy()

        self.bdm2_obj.evaluateLocalVelocityRepresentation(0,True)
        assert np.allclose(self.bdm2_obj.q[('velocity',0)],bdm_values)
 

if __name__ == '__main__':
    pass
