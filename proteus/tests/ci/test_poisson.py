#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
def test_c0p1():
    import poisson_3d_p
    import poisson_3d_c0p1_n
    pList = [poisson_3d_p]
    nList = [poisson_3d_c0p1_n]
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
    assert(True)

def test_c0p2():
    import poisson_3d_p
    import poisson_3d_c0p2_n
    pList = [poisson_3d_p]
    nList = [poisson_3d_c0p2_n]
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
    assert(True)

def check_c0q1(test_hexMesh_3x3=False,use_petsc=False):
    import poisson_3d_p
    import poisson_3d_c0q1_n
    poisson_3d_c0q1_n.hex=True
    if test_hexMesh_3x3 == True:
        poisson_3d_p.meshfile='hexMesh_3x3'
        poisson_3d_p.domain = Domain.MeshHexDomain(poisson_3d_p.meshfile)
        poisson_3d_p.x0 = (-3.,-3.,-3.)
        poisson_3d_p.L  = ( 6., 6., 6.)
    pList = [poisson_3d_p]
    nList = [poisson_3d_c0q1_n]
    so = default_so
    so.name = pList[0].name = "poisson_3d_c0q1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True

    if use_petsc == True:
        nList[0].linearSolver=default_n.KSP_petsc4py
        nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_3d_c0q1')
    assert(True)

def test_c0q1():
    """
    Test Hexes with direct solver, regular domain
    """
    check_c0q1(test_hexMesh_3x3=False,use_petsc=False)

def test_c0q1_hex_mesh():
    from proteus import MeshTools
    xmf_archive_base=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..','hex_cube_3x3')
    heavy_file_base = xmf_archive_base
    mesh_info = MeshTools.readMeshXdmf(xmf_archive_base,heavy_file_base)
    hex_meshfile_base = 'hexmesh_3x3'
    MeshTools.writeHexMesh(mesh_info,hex_meshfile_base,index_base=1)
    check_c0q1(test_hexMesh_3x3=True,use_petsc=False)

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
    test_c0p1()
    test_c0p2()
    test_c0q1()
    test_c0q1_hex_mesh()
    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass
