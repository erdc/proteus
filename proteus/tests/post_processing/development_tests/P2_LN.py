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

def test_BDM2_P2():
    '''
    Test the construction of a BDM2 projection matrix and rhs
    on the reference triangle
    '''
    # bdm_tests_template loads 
    import example1 as ex
    import numpy as np

    transport_obj = ex.ns.modelList[0].levelModelList[0]
    bdm2_obj = transport_obj.velocityPostProcessor.vpp_algorithms[0]
    
    import pdb
    pdb.set_trace()
    # Test correct node star is being created


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_BDM2_P2()
   # import nose
   # nose.main()

