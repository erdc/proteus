#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import poisson_3d_tetgen_p
import poisson_3d_tetgen_c0p1_n
import poisson_3d_tetgen_c0p2_n
reload(poisson_3d_tetgen_p)
reload(poisson_3d_tetgen_c0p1_n)
reload(poisson_3d_tetgen_c0p2_n)
import pytest

@pytest.mark.modelTest
@pytest.mark.poissonTest
class TestPoissonTetgen():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        reload(poisson_3d_tetgen_p)

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ['meshNoVessel.1.neigh',
                    'meshNoVessel.edge',
                    'meshNoVessel.ele',
                    'meshNoVessel.face',
                    'meshNoVessel.node',
                    'meshNoVessel.poly',
                    'poisson_3d_tetgen_c0p1pe1.xmf',
                    'poisson_3d_tetgen_c0p1pe1.h5',
                    'poisson_3d_tetgen_c0p2pe1.xmf',
                    'poisson_3d_tetgen_c0p2pe1.h5' ]
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_c0p1(genMesh=True):
        reload(poisson_3d_tetgen_c0p1_n)
        pList = [poisson_3d_tetgen_p]
        nList = [poisson_3d_tetgen_c0p1_n]
        so = default_so
        so.name = pList[0].name = "poisson_3d_tetgen_c0p1"+"pe"+`comm.size()`
        so.sList=[default_s]
        Profiling.logLevel=7
        Profiling.verbose=False
        opts.generatePartitionedMeshFromFiles = True
        opts.gatherArchive=True
        pList[0].genMesh=genMesh
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        #nList[0].linearSolver=default_n.LU
        #nList[0].multilevelLinearSolver=default_n.LU
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0p1')
        assert(True)

    @pytest.mark.slowTest
    def test_c0p2(genMesh=True):
        reload(poisson_3d_tetgen_c0p2_n)
        pList = [poisson_3d_tetgen_p]
        nList = [poisson_3d_tetgen_c0p2_n]
        so = default_so
        so.name = pList[0].name = "poisson_3d_tetgen_c0p2"+"pe"+`comm.size()`
        so.sList=[default_s]
        Profiling.logLevel=7
        Profiling.verbose=False
        opts.generatePartitionedMeshFromFiles = True
        opts.gatherArchive=True
        pList[0].genMesh=genMesh
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0p2')
        assert(True)

# def test_c0q1():
#     import poisson_3d_p
#     import poisson_3d_c0q1_n
#     pList = [poisson_3d_p]
#     nList = [poisson_3d_c0q1_n]
#     so = default_so
#     so.name = pList[0].name = "poisson_3d_c0q1"+"pe"+`comm.size()`
#     so.sList=[default_s]
#     opts.logLevel=7
#     opts.verbose=True
#     nList[0].linearSolver=default_n.KSP_petsc4py
#     nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
#     ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
#     ns.calculateSolution('poisson_3d_c0q1')
#     assert(True)

# def test_c0q2():
#     import poisson_3d_p
#     import poisson_3d_c0q2_n
#     pList = [poisson_3d_p]
#     nList = [poisson_3d_c0q2_n]
#     so = default_so
#     so.name = pList[0].name = "poisson_3d_c0q2"+"pe"+`comm.size()`
#     so.sList=[default_s]
#     opts.logLevel=7
#     opts.verbose=True
#     nList[0].linearSolver=default_n.KSP_petsc4py
#     nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
#     ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
#     ns.calculateSolution('poisson_3d_c0q2')
#     assert(True)

if __name__ == '__main__':
    pass
