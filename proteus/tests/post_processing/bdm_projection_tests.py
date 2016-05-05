#!/usr/bin/env python
"""

Test module for BDM2 Elements

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

def test_BDM2_reference_triangle():
    '''
    Test the construction of a BDM2 projection matrix and rhs
    on the reference triangle
    '''
    # bdm_tests_template loads 
    import bdm_tests_template as bt
    import numpy as np

    transport_obj = bt.ns.modelList[0].levelModelList[0]
    bdm2_obj = transport_obj.velocityPostProcessor.vpp_algorithms[0]


    # ******************* TEST PROJECTION MATRIX CONSTRUCTION ************


    # need to override factored BDM projection matrix
    bdm2_obj.BDMprojectionMat_element = \
                        np.zeros_like(bdm2_obj.BDMprojectionMat_element)

    bdm2_obj.buildLocalBDM2projectionMatrices \
                               (bdm2_obj.degree,
                                bdm2_obj.w_dS[0],
                                bdm2_obj.vt.ebq['n'],
                                bdm2_obj.vt.ebq[('v',0)],
                                bdm2_obj.q[('w',0)],     
                                bdm2_obj.weightedInteriorTestGradients,  
                                bdm2_obj.weightedInteriorDivFreeElement, 
                                bdm2_obj.piola_trial_function,           
                                bdm2_obj.BDMprojectionMat_element)        

#    The following .savetxt command  generates the comparison output.  Be sure
#    this is actually generating what you want before you uncomment!  The 
#    currently stored file should be correct.

#    np.savetxt('bdm2_ref_proj_mat.txt', bdm2_obj.BDMprojectionMat_element[0])

    comparison_matrix = np.loadtxt('bdm2_ref_proj_mat.txt', dtype = float)

    assert np.allclose(comparison_matrix,bdm2_obj.BDMprojectionMat_element)


    # ******************** TEST RHS CONSTRUCTION *************************

    # construct a RHS vector from a velocity field of all 1's
    bdm2_obj.ebq[('velocity',0)] = np.ones_like(bdm2_obj.ebq[('velocity',0)])
    bdm2_obj.q[('velocity',0)] = np.ones_like(bdm2_obj.q[('velocity',0)])

    bdm2_obj.buildBDM2rhs(bdm2_obj.BDMprojectionMat_element,
                          bdm2_obj.BDMprojectionMatPivots_element,
                          bdm2_obj.w_dS[0],
                          bdm2_obj.vt.ebq['n'],
                          bdm2_obj.weightedInteriorTestGradients,
                          bdm2_obj.weightedInteriorDivFreeElement,
                          bdm2_obj.ebq[('velocity',0)],
                          bdm2_obj.q[('velocity',0)],
                          bdm2_obj.q[('velocity_dofs',0)])

    test_rhs = bdm2_obj.q[('velocity_dofs',0)]
    comparison_rhs = np.array([ 3.33333333e-01,  3.33333333e-01,  1.33333333e+00,
                               -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                               -1.66666667e-01, -1.66666667e-01, -6.66666667e-01,
                               -1.00000000e+00,  5.00000000e-01,  4.33680869e-19])

    assert np.allclose(comparison_rhs,test_rhs)

def test_bdm_sshaped_region():
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

    # start by deleting output files if they exist
    filenames = ['poisson_bdm1_test.h5', 'poisson_bdm1_test.xmf']
    for file in filenames:
        if os.path.exists(file):
            try:
                os.remove(file)
            except OSError, e:
                print ("Error: %s - %s." %(e.filename, e.strerror ))
        else:
            pass

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
    expected = tables.openFile(
        os.path.join(test_path,
                     'test_bdm_sshaped_region_expected.h5'),'r')
    actual = tables.openFile('poisson_bdm1_test.h5','r')

    assert np.allclose(expected.root.velocity_0_elementQuadrature_p_t1, \
                       actual.root.velocity_0_elementQuadrature_p_t1), \
           'post-processed velocity field is no longer producing expectout out'

    expected.close()
    actual.close()

    

def test_piola_mapping():
    '''
    Test the construction of a BDM2 projection matrix and rhs on a 
    triangle different from the reference triangle.
    '''


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_BDM2_reference_triangle()
    test_bdm_sshaped_region()
   # import nose
   # nose.main()

