from pytest import approx

def test_wdot():
    """
    test_wdot

    Verifies that the proteus.LinearAlgebraTools.wdot computes dot product correctly.
    """
    import numpy as np
    import numpy.testing as npt
    from proteus.LinearAlgebraTools import wDot
    from proteus import Comm
    comm = Comm.init()
    x = np.array([-1, 2, 3])
    y = np.array([5, 6, 10])
    h = np.array([0.5, 1.2, 6.0])
    t1 = 191.9        
    t2 = wDot(x, y, h)
    test = npt.assert_almost_equal
    test.description = 'test_wdot'
    assert t1 == approx(t2)

if __name__ == '__main__':
    for test, v1, v2 in test_wdot():
        print(test.description)
        test(v1,v2)
