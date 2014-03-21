import numpy as np
import numpy.linalg
import numpy.random
from nose.tools import eq_ as eq
from nose.tools import ok_ as ok

def test_l2Norm():
    """ Test to see if l2Norm function observes rules of norms
    """
    from proteus.LinearAlgebraTools import l2Norm
    eq(l2Norm(np.zeros(10,))<1.0e-30, True)
    x= numpy.random.rand(100)
    y= numpy.random.rand(100)
    eq(l2Norm(x+y)<= (l2Norm(x) + l2Norm(y)), True)
    eq(abs(l2Norm(100.0*x)- 100.0*l2Norm(x))< 1.0e-30,True)

test_l2Norm()

