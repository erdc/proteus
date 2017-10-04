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
from post_processing.import_modules import bdm_tests_template as bt_temp

@pytest.mark.PostProcessingTools
class TestBDM2Reference1():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        reload(bt_temp)
        self.transport_obj = bt_temp.ns.modelList[0].levelModelList[0]
        self.bdm2_obj = self.transport_obj.velocityPostProcessor.vpp_algorithms[0]
        self._setRelativePath()

    def teardown_method(self,method):
        """Tear down the test problem. """
        filenames = ['poisson_bdm1_test.h5', 'poisson_bdm1_test.xmf','reference_triangle.ele',
                     'reference_triangle.node', 'reference_triangle.poly','proteus.log']
        for file in filenames:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename, e.strerror ))
            else:
                pass

    def _setRelativePath(self):
        self.scriptdir = os.path.dirname(__file__)
    
    def test_BDM2_reference_triangle(self):
        '''
        Test the construction of a BDM2 projection matrix and rhs
        on the reference triangle
        '''

        # ******************* TEST PROJECTION MATRIX CONSTRUCTION ************
        
        # need to override factored BDM projection matrix
        self.bdm2_obj.BDMprojectionMat_element \
                         = np.zeros_like(self.bdm2_obj.BDMprojectionMat_element)

        self.bdm2_obj.buildLocalBDM2projectionMatrices \
                                   (self.bdm2_obj.degree,
                                    self.bdm2_obj.vt.ebq[('w*dS_u',0)],
                                    self.bdm2_obj.vt.ebq['n'],
                                    self.bdm2_obj.vt.ebq[('v',0)],
                                    self.bdm2_obj.q[('w',0)],     
                                    self.bdm2_obj.weightedInteriorTestGradients,  
                                    self.bdm2_obj.weightedInteriorDivFreeElement, 
                                    self.bdm2_obj.piola_trial_function,
                                    self.bdm2_obj.edgeFlags,
                                    self.bdm2_obj.BDMprojectionMat_element)        

        #    The following .savetxt command  generates the comparison output.  Be sure
        #    this is actually generating what you want before you uncomment!  The 
        #    currently stored file should be correct.

        #    np.savetxt('bdm2_ref_proj_mat.txt', bdm2_obj.BDMprojectionMat_element[0])
        rel_path = "comparison_files/bdm2_ref_proj_mat.txt"
        comparison_matrix = np.loadtxt(os.path.join(self.scriptdir,rel_path), dtype = float)
        assert np.allclose(comparison_matrix,self.bdm2_obj.BDMprojectionMat_element)

        # ******************** TEST RHS CONSTRUCTION *************************

        # construct a RHS vector from a velocity field of all 1's
        self.bdm2_obj.ebq[('velocity',0)] = np.ones_like(self.bdm2_obj.ebq[('velocity',0)])
        self.bdm2_obj.q[('velocity',0)] = np.ones_like(self.bdm2_obj.q[('velocity',0)])

        self.bdm2_obj.buildBDM2rhs(self.bdm2_obj.BDMprojectionMat_element,
                              self.bdm2_obj.BDMprojectionMatPivots_element,
                              self.bdm2_obj.vt.ebq[('w*dS_u',0)],
                              self.bdm2_obj.vt.ebq['n'],
                              self.bdm2_obj.weightedInteriorTestGradients,
                              self.bdm2_obj.weightedInteriorDivFreeElement,
                              self.bdm2_obj.ebq[('velocity',0)],
                              self.bdm2_obj.q[('velocity',0)],
                              self.bdm2_obj.q[('velocity_dofs',0)],
                              self.bdm2_obj.edgeFlags)

        test_rhs = self.bdm2_obj.q[('velocity_dofs',0)]
        comparison_rhs = np.array([ 3.33333333e-01,  3.33333333e-01,  1.33333333e+00,
                                   -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                                   -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                                   -1.00000000e+00,  5.00000000e-01,  4.33680869e-19])
        assert np.allclose(comparison_rhs,test_rhs)

    def test_BDM2_reference_triangle_full_in_space(self):
        rel_path_1 = "comparison_files/bdm_bdy_func_values.npy"
        rel_path_2 = "comparison_files/bdm_func_values.npy"
        bdm_bdy_values = np.load(os.path.join(self.scriptdir,rel_path_1))
        bdm_values = np.load(os.path.join(self.scriptdir,rel_path_2))
        
        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values.copy()
        self.bdm2_obj.q[('velocity',0)] = bdm_values.copy()
        self.bdm2_obj.evaluateLocalVelocityRepresentation(0,True)
        
        assert np.allclose(self.bdm2_obj.q[('velocity',0)],bdm_values)

    def test_BDM2_reference_triangle_full_not_in_space(self):
        rel_path_1 = "comparison_files/bdm_bdy_func_values_trig.npy"
        rel_path_2 = "comparison_files/bdm_func_values_trig.npy"
        bdm_bdy_values = np.load(os.path.join(self.scriptdir,rel_path_1))
        bdm_values = np.load(os.path.join(self.scriptdir,rel_path_2))
        
        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values
        self.bdm2_obj.q[('velocity',0)] = bdm_values

        self.bdm2_obj.evaluateLocalVelocityRepresentation(0,True)

        rel_path_3 = "comparison_files/trig_velocity_rep.npy"
        comparison_vec = np.load(os.path.join(self.scriptdir,rel_path_3))
        assert np.allclose(self.bdm2_obj.q[('velocity',0)],comparison_vec)


if __name__ == '__main__':
    pass
