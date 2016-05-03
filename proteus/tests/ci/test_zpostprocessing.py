#!/usr/bin/env python
"""

Test module for Post-Processing routines.

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq


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

    assert np.allclose(expected.root.velocity_0_elementQuadrature_t1, \
                       actual.root.velocity_0_elementQuadrature_t1), \
           'post-processed velocity field is no longer producing expectout out'

    expected.close()
    actual.close()
