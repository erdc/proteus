#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus import default_n,default_s,default_so
Profiling.logLevel=7
Profiling.verbose=True
def test_c0p1():
    import poisson_3d_p
    import poisson_3d_c0p1_n
    pList = [poisson_3d_p]
    nList = [poisson_3d_c0p1_n]
    so = default_so
    so.name = pList[0].name = "poisson_3d_c0p1"
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_3d_c0p1')
    assert(True)

# def test_c0p2():
#     import poisson_3d_p
#     import poisson_3d_c0p2_n
#     pList = [poisson_3d_p]
#     nList = [poisson_3d_c0p2_n]
#     so = default_so
#     so.name = pList[0].name = "poisson_3d_c0p2"
#     so.sList=[default_s]
#     nList[0].linearSolver=default_n.LU
#     nList[0].multilevelLinearSolver=default_n.LU
#     ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
#     ns.calculateSolution('poisson_3d_c0p2')
#     assert(True)

# def test_c0q1():
#     import poisson_3d_p
#     import poisson_3d_c0q1_n
#     pList = [poisson_3d_p]
#     nList = [poisson_3d_c0q1_n]
#     so = default_so
#     so.name = pList[0].name = "poisson_3d_c0q1"
#     so.sList=[default_s]
#     nList[0].linearSolver=default_n.LU
#     nList[0].multilevelLinearSolver=default_n.LU
#     ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
#     ns.calculateSolution('poisson_3d_c0q1')
#     assert(True)

# def test_c0q2():
#     import poisson_3d_p
#     import poisson_3d_c0q2_n
#     pList = [poisson_3d_p]
#     nList = [poisson_3d_c0q2_n]
#     so = default_so
#     so.name = pList[0].name = "poisson_3d_c0q2"
#     so.sList=[default_s]
#     nList[0].linearSolver=default_n.LU
#     nList[0].multilevelLinearSolver=default_n.LU
#     ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
#     ns.calculateSolution('poisson_3d_c0q2')
#     assert(True)

if __name__ == '__main__':
    test_c0p1()
#    test_c0p2()
#    test_c0q1()
#    test_c0q2()
