from proteus import (Comm, Profiling, Quadrature)
from proteus.Profiling import logEvent

comm = Comm.get()
Profiling.logLevel = 2
Profiling.verbose = True
import numpy as np


def ex2(N, hk, x0, y0):
    r"""
    Compute the integral of :math:`1_{x0*y+y0*x<= x0*y0}`
    :math:`0\leqx0\leq1` and :math:`0\leqy0\leq1`

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    x0 : integrand parameter
    y0 : integrand parameter
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(1 / hk))

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:, 0] * y0 + comp_quad.points[:, 1] * x0, x0 * y0)])
    logEvent("hk=%f\t true-int=%f\t comp-int=%f error=%f" % (1.0 / N, x0 * y0 * 0.5, ii, np.abs(ii - x0 * y0 * 0.5)))
    return 1.0 / N, ii - x0 * y0


def ex3(N, hk, r0):
    """
    Compute the integral of 1_{|[x,y]|<= r0}
    r0<=0.5*sqrt(2)
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(1 / hk))

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:, 0] ** 2 + comp_quad.points[:, 1] ** 2, r0 * r0)])
    print "hk=%f\t true-int=%f\t comp-int=%f error=%f" % (1.0 / N, r0**2 * 0.25 * np.pi, ii, np.abs(ii - r0**2 * 0.25 * np.pi))

M = 10
def test_line():
    #=========================================================================
    # example 1: straight line
    #=========================================================================
    for i in range(M):
        ex2(1, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(2, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(3, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"

def test_circle():
    #=========================================================================
    # example 2: circle
    #=========================================================================

    for i in range(M):
        ex3(1, 1.0 / 2**i, 0.5)
    print "\n \n"

    for i in range(M):
        ex3(2, 1.0 / 2**i, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(3, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"
    

if __name__ == '__main__':

    M = 10
    #=========================================================================
    # example 1: straight line
    #=========================================================================
    for i in range(M):
        ex2(1, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(2, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(3, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"
    #=========================================================================
    # example 2: circle
    #=========================================================================

    for i in range(M):
        ex3(1, 1.0 / 2**i, 0.5)
    print "\n \n"

    for i in range(M):
        ex3(2, 1.0 / 2**i, 0.3)
    print "\n \n"

    for i in range(M):
        ex2(3, 1.0 / 2**i, 0.2, 0.3)
    print "\n \n"
