#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
import os, pytest
from proteus.iproteus import *
from proteus import Comm, defaults
comm = Comm.get()
modulepath = os.path.dirname(os.path.abspath(__file__))

@pytest.mark.modelTest
@pytest.mark.poissonTest
class TestPoissonTetgen(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        pass

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ['meshNoVessel.neigh',
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

    @pytest.mark.slowTest
    def test_c0p1(genMesh=True):
        pList = [defaults.load_physics('poisson_3d_tetgen_p',
                                       modulepath)]
        nList = [defaults.load_numerics('poisson_3d_tetgen_c0p1_n',
                                        modulepath)]
        so = defaults.System_base()
        so.name = pList[0].name = "poisson_3d_tetgen_c0p1"+"pe"+repr(comm.size())
        so.sList=[default_s]
        Profiling.logLevel=7
        Profiling.verbose=False
        opts.generatePartitionedMeshFromFiles = True
        opts.gatherArchive=True
        pList[0].genMesh=genMesh
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0p1')
        assert(True)

    @pytest.mark.slowTest
    def test_c0p2(genMesh=True):
        pList = [defaults.load_physics('poisson_3d_tetgen_p',
                                       modulepath)]
        nList = [defaults.load_numerics('poisson_3d_tetgen_c0p2_n',
                                        modulepath)]
        so = defaults.System_base()
        so.name = pList[0].name = "poisson_3d_tetgen_c0p2"+"pe"+repr(comm.size())
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

if __name__ == '__main__':
    pass
