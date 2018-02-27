#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from petsc4py import PETSc
import os
import poisson_3d_p
import poisson_3d_c0p1_n
import poisson_3d_c0p2_n
import poisson_3d_c0q1_n
import poisson_3d_c0q2_n
import pytest

@pytest.mark.modelTest
@pytest.mark.poissonTest
class TestPoisson():

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
        reload(poisson_3d_p)

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
        reload(poisson_3d_p)
        reload(poisson_3d_c0p1_n)
        pList = [poisson_3d_p]
        nList = [poisson_3d_c0p1_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_3d_c0p1"+"pe"+`comm.size()`
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
        reload(poisson_3d_p)
        reload(poisson_3d_c0p2_n)
        pList = [poisson_3d_p]
        nList = [poisson_3d_c0p2_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_3d_c0p2"+"pe"+`comm.size()`
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

    def check_c0q1(self,test_hexMesh_3x3=False,use_petsc=False,name="_hexMesh_"):
        reload(poisson_3d_p)
        reload(poisson_3d_c0q1_n)
        poisson_3d_c0q1_n.hex=True
        if test_hexMesh_3x3 == True:
            poisson_3d_p.meshfile='hexMesh_3x3'
            poisson_3d_p.domain = Domain.MeshHexDomain(poisson_3d_p.meshfile)
            poisson_3d_p.x0 = (-3.,-3.,-3.)
            poisson_3d_p.L  = ( 6., 6., 6.)
        pList = [poisson_3d_p]
        nList = [poisson_3d_c0q1_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_3d_c0q1"+name+"pe"+`comm.size()`
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
        self.check_c0q1(test_hexMesh_3x3=True,use_petsc=True, name="_hexMesh_")

    @pytest.mark.slowTest
    def test_c0q2(self):
        reload(poisson_3d_p)
        reload(poisson_3d_c0q2_n)
        pList = [poisson_3d_p]
        nList = [poisson_3d_c0q2_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "poisson_3d_c0q2"+"pe"+`comm.size()`
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
