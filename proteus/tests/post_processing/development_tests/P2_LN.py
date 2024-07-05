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

def test_BDM2_P2():
    '''
    Test the construction of a BDM2 projection matrix and rhs
    on the reference triangle
    '''
    # bdm_tests_template loads 
    import example1 as ex
    import numpy as np
    import expected_output as eo

    transport_obj = ex.ns.modelList[0].levelModelList[0]
    bdm2_obj = transport_obj.velocityPostProcessor.vpp_algorithms[0]

    assert eo.globalDOF2Element_test == bdm2_obj.globalDOF2globalElementList
    assert eo.globalDOFGlobalElement2StarElement_test == bdm2_obj.globalDOFGlobalElement2StarElement
    assert np.allclose(eo.dofStarElementsArray_cmp,bdm2_obj.dofStarElementsArray)

    import pdb
    pdb.set_trace()
    # Test correct node star is being created

def test_BDM_P1():
    '''
    Tests some BDM stuff
    '''
    # bdm_tests_template loads 
    import example1 as ex
    import numpy as np
    import expected_output as eo

    import pdb
    pdb.set_trace()


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_BDM_P1()
#    test_BDM2_P2()
