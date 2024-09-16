#!/usr/bin/env python
"""

Test module for BDM2 Elements

"""
import proteus.test_utils.TestTools
import os
import sys
import inspect

proteus.test_utils.TestTools.addSubFolders( inspect.currentframe() )

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy as np
import pytest

import bdm_tests_template as bt_temp

@pytest.mark.PostProcessingTools
class TestBDM2Reference1(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        """Initialize the test problem. """
        from importlib import reload
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
                except OSError as e:
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
        import cpostprocessing
        # ******************* TEST PROJECTION MATRIX CONSTRUCTION ************

        # need to override factored BDM projection matrix
        self.bdm2_obj.BDMprojectionMat_element \
                         = np.zeros_like(self.bdm2_obj.BDMprojectionMat_element)

        cpostprocessing.buildLocalBDM2projectionMatrices \
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
        np.testing.assert_almost_equal(comparison_matrix,self.bdm2_obj.BDMprojectionMat_element[0],decimal=6)

        # ******************** TEST RHS CONSTRUCTION *************************

        # construct a RHS vector from a velocity field of all 1's
        self.bdm2_obj.ebq[('velocity',0)] = np.ones_like(self.bdm2_obj.ebq[('velocity',0)])
        self.bdm2_obj.q[('velocity',0)] = np.ones_like(self.bdm2_obj.q[('velocity',0)])

        cpostprocessing.buildBDM2rhs(self.bdm2_obj.BDMprojectionMat_element,
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
        np.testing.assert_almost_equal(comparison_rhs,test_rhs[0],decimal=6)

    def test_BDM2_reference_triangle_full_in_space(self):
        rel_path_1 = "comparison_files/bdm_bdy_func_values.npy"
        rel_path_2 = "comparison_files/bdm_func_values.npy"
        bdm_bdy_values = np.load(os.path.join(self.scriptdir,rel_path_1))
        bdm_values = np.load(os.path.join(self.scriptdir,rel_path_2))

        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values.copy()
        self.bdm2_obj.q[('velocity',0)] = bdm_values.copy()

        self.bdm2_obj.evaluateLocalVelocityRepresentation(0,True)
        #np.save(os.path.join(self.scriptdir,rel_path_2), self.bdm2_obj.q[('velocity',0)])
        np.testing.assert_almost_equal(self.bdm2_obj.q[('velocity',0)],bdm_values,decimal=6)


if __name__ == '__main__':
    pass
