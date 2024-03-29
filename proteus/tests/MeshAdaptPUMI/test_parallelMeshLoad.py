from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus import Domain
from proteus.MeshAdaptPUMI import MeshAdapt
from proteus import Comm
import os
from petsc4py import PETSc
import pytest

#Run these tests with `mpirun -n 2 py.test --boxed test_parallelMeshLoad.py`

@pytest.mark.skip(reason="need to run in parallel,also files are too large")
def test_3DparallelLoadPUMI(verbose=0):
    """Test to load 3D parallel PUMI model and mesh"""
    comm = Comm.init()
    eq(comm.size(),2)
    testDir=os.path.dirname(os.path.abspath(__file__))
    domain = Domain.PUMIDomain()
    Model=testDir+ '/Prism.dmg'
    Mesh=testDir + '/Prism.smb'
    domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI()
    domain.PUMIMesh.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    mesh.convertFromPUMI(domain.PUMIMesh, domain.faceList, domain.regList,parallel = comm.size() > 1, dim = domain.nd)
    eq(mesh.nElements_global,8148)
    eq(mesh.nNodes_global,1880)
    eq(mesh.nEdges_global,11001)
    eq(mesh.nElementBoundaries_global,17270)
    #Ideally, we can assert whether each rank gets the proper number of mesh entities/
    #With the present setup, that information doesn't seem accessible
    #if(comm.rank()==0):
    #    eq(mesh.nElements_owned,4074)
    #    eq(mesh.nNodes_owned,994)
    #    eq(mesh.nEdges_owned,8729)
    #    eq(mesh.nElementBoundaries_owned,5648)

@pytest.mark.skip(reason="need to run in parallel")
def test_2DparallelLoadPUMI(verbose=0):
    """Test to load 2D parallel PUMI model and mesh"""
    comm = Comm.init()
    eq(comm.size(),2)
    testDir=os.path.dirname(os.path.abspath(__file__))
    domain = Domain.PUMIDomain(dim=2)
    Model=testDir+ '/Rectangle.dmg'
    Mesh=testDir + '/Rectangle.smb'
    domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI()
    domain.PUMIMesh.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    mesh = MeshTools.TriangularMesh()
    mesh.cmesh = cmeshTools.CMesh()
    mesh.convertFromPUMI(domain,domain.PUMIMesh, domain.faceList, domain.regList,parallel = comm.size() > 1, dim = domain.nd)
    assert mesh.nElements_global == 8
    assert mesh.nNodes_global == 10
    assert mesh.nEdges_global == 17
    assert mesh.nElementBoundaries_global == 17

