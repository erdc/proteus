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
import bdm_tests_template as bt_temp

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
                                    self.bdm2_obj.w_dS[0],
                                    self.bdm2_obj.vt.ebq['n'],
                                    self.bdm2_obj.vt.ebq[('v',0)],
                                    self.bdm2_obj.q[('w',0)],     
                                    self.bdm2_obj.weightedInteriorTestGradients,  
                                    self.bdm2_obj.weightedInteriorDivFreeElement, 
                                    self.bdm2_obj.piola_trial_function,           
                                    self.bdm2_obj.BDMprojectionMat_element)        

        #    The following .savetxt command  generates the comparison output.  Be sure
        #    this is actually generating what you want before you uncomment!  The 
        #    currently stored file should be correct.

        #    np.savetxt('bdm2_ref_proj_mat.txt', bdm2_obj.BDMprojectionMat_element[0])
        comparison_matrix = np.loadtxt('./comparison_files/bdm2_ref_proj_mat.txt', dtype = float)
        assert np.allclose(comparison_matrix,self.bdm2_obj.BDMprojectionMat_element)

        # ******************** TEST RHS CONSTRUCTION *************************

        # construct a RHS vector from a velocity field of all 1's
        self.bdm2_obj.ebq[('velocity',0)] = np.ones_like(self.bdm2_obj.ebq[('velocity',0)])
        self.bdm2_obj.q[('velocity',0)] = np.ones_like(self.bdm2_obj.q[('velocity',0)])

        self.bdm2_obj.buildBDM2rhs(self.bdm2_obj.BDMprojectionMat_element,
                              self.bdm2_obj.BDMprojectionMatPivots_element,
                              self.bdm2_obj.w_dS[0],
                              self.bdm2_obj.vt.ebq['n'],
                              self.bdm2_obj.weightedInteriorTestGradients,
                              self.bdm2_obj.weightedInteriorDivFreeElement,
                              self.bdm2_obj.ebq[('velocity',0)],
                              self.bdm2_obj.q[('velocity',0)],
                              self.bdm2_obj.q[('velocity_dofs',0)])

        test_rhs = self.bdm2_obj.q[('velocity_dofs',0)]

        comparison_rhs = np.array([ 3.33333333e-01,  3.33333333e-01,  1.33333333e+00,
                                   -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                                   -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                                   -1.00000000e+00,  5.00000000e-01,  4.33680869e-19])

        assert np.allclose(comparison_rhs,test_rhs)

    def test_BDM2_reference_triangle_full_in_space(self):
        bdm_bdy_values = np.load('./comparison_files/bdm_bdy_func_values.npy')
        bdm_values = np.load('./comparison_files/bdm_func_values.npy')
        
        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values.copy()
        self.bdm2_obj.q[('velocity',0)] = bdm_values.copy()
        
        self.bdm2_obj.evaluateLocalVelocityRepresentation(0)

        assert np.allclose(self.bdm2_obj.q[('velocity',0)],bdm_values)

    def test_BDM2_reference_triangle_full_not_in_space(self):
        bdm_bdy_values = np.load('./comparison_files/bdm_bdy_func_values_trig.npy')
        bdm_values = np.load('./comparison_files/bdm_func_values_trig.npy')

        self.bdm2_obj.ebq[('velocity',0)] = bdm_bdy_values
        self.bdm2_obj.q[('velocity',0)] = bdm_values

        self.bdm2_obj.evaluateLocalVelocityRepresentation(0)

        comparison_vec = np.load('./comparison_files/trig_velocity_rep.npy')
        assert np.allclose(self.bdm2_obj.q[('velocity',0)],comparison_vec)


    def test_bdm_sshaped_region(self):
        '''
        Tests a post-processed velocity field projected
        into bdm1 space.
        '''
        # import and run a small 2D poiseulle problem
        reload(default_so)
        reload(default_p)
        reload(default_n)
        import sShaped_block_2d_p
        import sShaped_block_2d_n
        import tables
        import numpy as np
        import os

        pList = [sShaped_block_2d_p]
        nList = [sShaped_block_2d_n]
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList=[default_s]
        #opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        # create conservative velocity field
        simFlagsList = [{}]
        simFlagsList[0]['storeQuantities'] = ['u',"q:('velocity',0)"]
        # solve problem
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts,simFlagsList)
        ns.calculateSolution('test1')
        #
        test_path = os.path.dirname(os.path.abspath(__file__))
        expected = tables.openFile(os.path.join(test_path,
                                   './comparison_files/test_bdm_sshaped_region_expected.h5'),'r')
        actual = tables.openFile('poisson_bdm1_test.h5','r')

        assert np.allclose(expected.root.velocity_0_elementQuadrature_p_t1, \
                           actual.root.velocity_0_elementQuadrature_t1), \
               'post-processed velocity field is no longer producing expectout out'

        expected.close()
        actual.close()

        # delete output files
        filenames = ['poisson_bdm1_test.h5', 'poisson_bdm1_test.xmf','blockDomain.ele',
                     'blockDomain.node', 'blockDomain.poly',
                     'blockDomain.edge', 'blockDomain.neig',
                     'proteus.log',
                     'BDM2_Test_File.h5', 'BDM2_Test_File.xmf']
        for file in filenames:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename, e.strerror ))
            else:
                pass

# def test_piola_mapping():
#     '''
#     Test the construction of a BDM2 projection matrix and rhs on a 
#     triangle different from the reference triangle.
#     '''


if __name__ == '__main__':
    pass
