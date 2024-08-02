#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus.defaults import load_physics, load_numerics, System_base
from petsc4py import PETSc
import os
import pytest
modulepath = os.path.dirname(os.path.abspath(__file__))

@pytest.mark.modelTest
@pytest.mark.poissonTest
class TestPoisson(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        OptDB = PETSc.Options()
        OptDB.setValue("ksp_type", "cg")
        OptDB.setValue("pc_type", "gamg")

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ['tetgen.ele',
                    'tetgen.face',
                    'tetgen.node',
                    'hexmesh_3x3.mesh',
                    'poisson_3d_c0p1pe1.xmf',
                    'poisson_3d_c0p1pe1.h5',
                    'poisson_3d_c0p2pe1.xmf',
                    'poisson_3d_c0p2pe1.h5',
                    'poisson_3d_c0q1_hexMesh_pe1.xmf',
                    'poisson_3d_c0q1_hexMesh_pe1.h5',
                    'poisson_3d_c0q1_proteusMesh_pe1.xmf',
                    'poisson_3d_c0q1_proteusMesh_pe1.h5',
                    'poisson_3d_tetgen_c0p1pe1.xmf',
                    'poisson_3d_tetgen_c0p1pe1.h5',
                    'poisson_3d_c0q2pe1.xmf',
                    'poisson_3d_c0q2pe1.h5',
                    'reference_triangle_2d.node',
                    'reference_triangle_2d.ele',
                    'reference_triangle_2d.poly',
                    'reference_triangle_3d.ele',
                    'reference_triangle_3d.node',
                    'reference_triangle_3d.poly',
                    'reference_triangle_3d.face']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_c0p1(self):
        pList = [load_physics('poisson_3d_p',
                              modulepath)]
        nList = [load_numerics('poisson_3d_c0p1_n',
                               modulepath)]
        so = System_base()
        so.name = pList[0].name = "poisson_3d_c0p1"+"pe"+repr(comm.size())
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        #nList[0].linearSolver=default_n.LU
        #nList[0].multilevelLinearSolver=default_n.LU
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0p1')
        del ns
        assert(True)

    @pytest.mark.slowTest
    def test_c0p2(self):
        pList = [load_physics('poisson_3d_p',
                              modulepath)]
        nList = [load_numerics('poisson_3d_c0p2_n',
                               modulepath)]
        so = System_base()
        so.name = pList[0].name = "poisson_3d_c0p2"+"pe"+repr(comm.size())
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.gatherArchive=True
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0p2')
        del ns
        assert(True)

    def check_c0q1(self,test_hexMesh_3x3=False,use_petsc=True,name="_hexMesh_"):
        poisson_3d_p = load_physics('poisson_3d_p',
                                    modulepath)
        poisson_3d_c0q1_n = load_numerics('poisson_3d_c0q1_n',
                                          modulepath)
        poisson_3d_c0q1_n.hex=True
        if test_hexMesh_3x3 == True:
            poisson_3d_p.meshfile='hexMesh_3x3'
            poisson_3d_p.domain = Domain.MeshHexDomain(poisson_3d_p.meshfile)
            poisson_3d_p.x0 = (-3.,-3.,-3.)
            poisson_3d_p.L  = ( 6., 6., 6.)
        pList = [poisson_3d_p]
        nList = [poisson_3d_c0q1_n]
        so = System_base()
        so.name = pList[0].name = "poisson_3d_c0q1"+name+"pe"+repr(comm.size())
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True

        if use_petsc == True:
            nList[0].linearSolver=default_n.KSP_petsc4py
            nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0q1')
        del ns
        assert(True)

    def test_c0q1(self):
        """
        Test Hexes with direct solver, regular domain
        """
        self.check_c0q1(test_hexMesh_3x3=False,use_petsc=True, name="_proteusMesh_")

    def test_c0q1_hex_mesh(self):
        from proteus import MeshTools
        xmf_archive_base=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..','hex_cube_3x3')
        heavy_file_base = xmf_archive_base
        mesh_info = MeshTools.readMeshXdmf(xmf_archive_base,heavy_file_base)
        hex_meshfile_base = 'hexmesh_3x3'
        MeshTools.writeHexMesh(mesh_info,hex_meshfile_base,index_base=1)
        self.check_c0q1(test_hexMesh_3x3=False,use_petsc=True, name="_hexMesh_")

    @pytest.mark.slowTest
    def test_c0q2(self):
        pList = [load_physics('poisson_3d_p',
                              modulepath)]
        nList = [load_numerics('poisson_3d_c0q2_n',
                               modulepath)]
        so = System_base()
        so.name = pList[0].name = "poisson_3d_c0q2"+"pe"+repr(comm.size())
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('poisson_3d_c0q2')
        del ns
        assert(True)

if __name__ == '__main__':
    pass
