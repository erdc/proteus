from proteus import (Comm, Profiling, Quadrature)
from proteus.Profiling import logEvent

import os
import sys
comm = Comm.get()
Profiling.procID = comm.rank()
Profiling.logFile = sys.stdout
Profiling.logLevel = 2
Profiling.verbose = True

import numpy as np
import numpy.testing as npt
import unittest


logEvent("Testing Composite quadrature rule")


def ex2(N, hk, x0=1, y0=1, z0=1):
    r"""
    Compute the integral of :math:`1_{x/x0+y/y0+z/z0<= 1}`
    :math:`0< x0\leq 1` and :math:`0< y0\leq 1` and :math:`0< z0\leq 1`

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    x0 : integrand parameter
    y0 : integrand parameter
    z0 : integrand parameter
    """

    quad = Quadrature.GaussTetrahedron(N)
    comp_quad = Quadrature.CompositeTetrahedron(quad, hk)

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:,0]/x0 + comp_quad.points[:,1]/y0 + +comp_quad.points[:,2]/z0, 1.0)])

    ee = np.abs(ii - x0 * y0 * z0 / 6.0)

    logEvent("hk=%f\t true-int=%f\t comp-int=%f error=%f" %
             (comp_quad.h, x0 * y0 * z0 / 6.0, ii, ee))
    return comp_quad.h, ee


def ex3(N, hk, r0):
    r"""
    Compute the integral of  :math:`1_{|[x,y,z]|<= r0}`
    :math: `r0<=0.5*sqrt(2)`

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    r0 : integrand parameter
    """

    quad = Quadrature.GaussTetrahedron(N)
    comp_quad = Quadrature.CompositeTetrahedron(quad, hk)

    ii = np.sum(comp_quad.weights[np.less_equal(
        comp_quad.points[:, 0] ** 2 + comp_quad.points[:, 1] ** 2 + comp_quad.points[:, 2] ** 2, r0 * r0)])
    ee = np.abs(ii - r0**3 * np.pi / 6.0)
    logEvent("hk=%f\t true-int=%f\t comp-int=%f error=%f" %
             (comp_quad.h, r0**3 * np.pi / 6.0, ii, ee))
    return comp_quad.h, ee


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

    quad = Quadrature.GaussTetrahedron(N)
    comp_quad = Quadrature.CompositeTetrahedron(quad, hk)

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(quad_weights * quad_points[:, 0])
    ii_comp_quad = np.sum(comp_quad.weights * comp_quad.points[:, 0])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (comp_quad.h, ii_quad, ii_comp_quad, ee))
    return comp_quad.h, ee


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

    quad = Quadrature.GaussTetrahedron(N)
    comp_quad = Quadrature.CompositeTetrahedron(quad, hk)

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(quad_weights * quad_points[:, 0] * quad_points[:, 1])
    ii_comp_quad = np.sum(comp_quad.weights *
                          comp_quad.points[:, 0] * comp_quad.points[:, 1])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (comp_quad.h, ii_quad, ii_comp_quad, ee))
    return comp_quad.h, ee


def ex6(N, hk):
    r"""
    Compute the integral of :math:`x*y*z` over the reference triangle

    Parameters
    ----------
    N : int
        Order of base quadrature rule
    hk : double
         submesh element diameter
    """

    quad = Quadrature.GaussTetrahedron(N)
    comp_quad = Quadrature.CompositeTetrahedron(quad, hk)

    quad_points = np.asarray(quad.points, 'd')
    quad_weights = np.asarray(quad.weights, 'd')

    ii_quad = np.sum(
        quad_weights * quad_points[:, 0] * quad_points[:, 1] * quad_points[:, 2])
    ii_comp_quad = np.sum(comp_quad.weights *
                          comp_quad.points[:, 0] * comp_quad.points[:, 1] * comp_quad.points[:, 2])

    ee = np.abs(ii_quad - ii_comp_quad)
    logEvent("hk=%f\t quad-int=%f\t comp-quad-int=%f error=%f" %
             (comp_quad.h, ii_quad, ii_comp_quad, ee))
    return comp_quad.h, ee


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
        cc[i] = np.log(ee[i]/(ee[i - 1] + 1e-15))/np.log(hh[i]/hh[i - 1])


class TestCompQuad(unittest.TestCase):

    def test_constant_exact(self):
        #======================================================================
        # example 0: 1
        #======================================================================
        M = 5
        for i in range(M):
            hk, error = ex2(1, 1.0/2**i, 1.0, 1.0, 1.0)
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_1st_poly_exact(self):
        #======================================================================
        # example 1: x
        #======================================================================
        M = 5
        for i in range(M):
            hk, error = ex4(1, 1.0/2**i)
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_2nd_poly_exact(self):
        #======================================================================
        # example 2: x*y
        #======================================================================
        M = 5
        for i in range(M):
            hk, error = ex5(2, 1.0/2**i)
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_3rd_poly_exact(self):
        #======================================================================
        # example 3: x*y*z
        #======================================================================
        M = 5
        for i in range(M):
            hk, error = ex6(3, 1.0/2**i)
            assert np.allclose(error, 0.0, atol=1e-10)

    def test_plane(self):
        M = 5
        #======================================================================
        # example 5: plane
        #======================================================================
        cell_size = np.zeros((M,), 'd')
        error = np.zeros((M,), 'd')
        convergence_rate = np.zeros((M,), 'd')

        for i in range(M):
            cell_size[i], error[i] = ex2(1, 1.0/2**i, 0.5, 0.6, 0.7)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(2, 1.0/2**i, 0.5, 0.6, 0.7)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(3, 1.0/2**i, 0.5, 0.6, 0.7)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex2(4, 1.0/2**i, 0.5, 0.6, 0.7)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

    def test_sphere(self):
        M = 5
        #======================================================================
        # example 6: circle
        #======================================================================
        cell_size = np.zeros((M,), 'd')
        error = np.zeros((M,), 'd')
        convergence_rate = np.zeros((M,), 'd')

        for i in range(M):
            cell_size[i], error[i] = ex3(1, 1.0/2**i, 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(2, 1.0/2**i, 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(3, 1.0/2**i, 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")

        for i in range(M):
            cell_size[i], error[i] = ex3(4, 1.0/2**i, 0.5)

        get_convergence_rate(cell_size, error, convergence_rate)
        logEvent("average convergence rate is %f" %
                 np.average(convergence_rate[1:]))
        self.assertGreater(np.average(
            convergence_rate[1:]), 1.0, "convergence should be > 1")


if __name__ == '__main__':

    unittest.main(verbosity=2)
