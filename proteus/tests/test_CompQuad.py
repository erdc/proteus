from __future__ import division
from builtins import range
from past.utils import old_div
from proteus import (Comm, Profiling, Quadrature)
from proteus.Profiling import logEvent

import os
import sys
comm = Comm.get()
Profiling.procID = comm.rank()
Profiling.logFile = sys.stdout
Profiling.logLevel = 2
Profiling.verbose = False

import numpy as np
import numpy.testing as npt
import unittest


logEvent("Testing Composite quadrature rule")


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

    N = int(np.ceil(old_div(1, hk)))

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:, 0] * y0 + comp_quad.points[:, 1] * x0, x0 * y0)])

    ee = np.abs(ii - x0 * y0 * 0.5)

    logEvent("hk=%f\t true-int=%f\t comp-int=%f error=%f" %
             (old_div(1.0, N), x0 * y0 * 0.5, ii, ee))
    return old_div(1.0, N), ee


def ex3(N, hk, r0):
    r"""
    Compute the integral of  :math:`1_{|[x,y]|<= r0}`
    :math: `r0<=0.5*sqrt(2)`

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    r0 : integrand parameter
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(old_div(1, hk)))

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:, 0] ** 2 + comp_quad.points[:, 1] ** 2, r0 * r0)])
    ee = np.abs(ii - r0**2 * 0.25 * np.pi)
    logEvent("hk=%f\t true-int=%f\t comp-int=%f error=%f" %
             (old_div(1.0, N), r0**2 * 0.25 * np.pi, ii, ee))
    return old_div(1.0, N), ee


def ex4(N, hk):
    r"""
    Compute the integral of :math:`x` over the reference triangle

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(old_div(1, hk)))

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(quad_weights * quad_points[:, 0])
    ii_comp_quad = np.sum(comp_quad.weights * comp_quad.points[:, 0])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (old_div(1.0, N), ii_quad, ii_comp_quad, ee))
    return old_div(1.0, N), ee


def ex5(N, hk):
    r"""
    Compute the integral of :math:`x*y` over the reference triangle

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(old_div(1, hk)))

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(quad_weights * quad_points[:, 0] * quad_points[:, 1])
    ii_comp_quad = np.sum(comp_quad.weights *
                          comp_quad.points[:, 0] * comp_quad.points[:, 1])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (old_div(1.0, N), ii_quad, ii_comp_quad, ee))
    return old_div(1.0, N), ee


def ex6(N, hk):
    r"""
    Compute the integral of :math:`x*y*y` over the reference triangle

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    """

    quad = Quadrature.GaussTriangle(N)
    comp_quad = Quadrature.CompositeTriangle(quad, hk)

    N = int(np.ceil(old_div(1, hk)))

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(
        quad_weights * quad_points[:, 0] * quad_points[:, 1] * quad_points[:, 1])
    ii_comp_quad = np.sum(comp_quad.weights *
                          comp_quad.points[:, 0] * comp_quad.points[:, 1] * comp_quad.points[:, 1])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (old_div(1.0, N), ii_quad, ii_comp_quad, ee))
    return old_div(1.0, N), ee


def get_convergence_rate(hh, ee, cc):
    r"""
    Compute the convergence cc from mesh size hh and error ee

    Parameters
    ----------
    hh : numpy double array, with shape (M,) 
        mesh size
    e : numpy double array, with shape (M,) 
        error
    cc : numpy double array, with shape (M,) 
        convergence rate: :math:`cc[i]=log(ee[i]/ee[i-1])/log(hh[i]/hh[i-1]), i>0`
    """
    for i in range(1, hh.shape[0]):
        cc[i] = old_div(np.log(old_div(ee[i], (ee[i - 1] + 1e-15))), np.log(old_div(hh[i], hh[i - 1])))


class TestCompQuad(unittest.TestCase):

    def test_constant_exact(self):
        #======================================================================
        # example 0: 1
        #======================================================================
        M = 10
        for i in range(M):
            hk, error = ex2(1, old_div(1.0, 2**i), 1.0, 1.0)
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_1st_poly_exact(self):
        #======================================================================
        # example 1: x
        #======================================================================
        M = 10
        for i in range(M):
            hk, error = ex4(1, old_div(1.0, 2**i))
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_2nd_poly_exact(self):
        #======================================================================
        # example 2: x*y
        #======================================================================
        M = 10
        for i in range(M):
            hk, error = ex5(2, old_div(1.0, 2**i))
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_3rd_poly_exact(self):
        #======================================================================
        # example 3: x*y^2
        #======================================================================
        M = 10
        for i in range(M):
            hk, error = ex6(3, old_div(1.0, 2**i))
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_line(self):
        M = 10
        #======================================================================
        # example 5: straight line
        #======================================================================
        cell_size = np.zeros((M,), 'd')
        error = np.zeros((M,), 'd')
        convergence_rate = np.zeros((M,), 'd')

        for i in range(M):
            cell_size[i], error[i] = ex2(1, old_div(1.0, 2**i), 0.6, 1.0)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(2, old_div(1.0, 2**i), 0.6, 1.0)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(3, old_div(1.0, 2**i), 0.6, 1.0)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(4, old_div(1.0, 2**i), 0.6, 1.0)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

    def test_circle(self):
        M = 10
        #======================================================================
        # example 6: circle
        #======================================================================
        cell_size = np.zeros((M,), 'd')
        error = np.zeros((M,), 'd')
        convergence_rate = np.zeros((M,), 'd')

        for i in range(M):
            cell_size[i], error[i] = ex3(1, old_div(1.0, 2**i), 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(2, old_div(1.0, 2**i), 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(3, old_div(1.0, 2**i), 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(4, old_div(1.0, 2**i), 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")


if __name__ == '__main__':

    unittest.main(verbosity=2)
