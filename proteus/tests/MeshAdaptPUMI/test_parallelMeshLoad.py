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

def test_3DparallelLoadPUMI(verbose=0):
    """Test to load 3D parallel PUMI model and mesh"""
    comm = Comm.init()
    if comm.size()!= 2:
        print("Skipping 3DprallelLoadPUMI,must be run on 2 MPI tasks")
        return
    testDir=os.path.dirname(os.path.abspath(__file__))
    domain = Domain.PUMIDomain(manager=MeshAdapt.AdaptManager())
    Model=testDir+ '/cubepar.dmg'
    Mesh=testDir + '/cubepar.smb'
    #domain.PUMIMesh=MeshAdapt.MeshAdaptPUMI()
    domain.AdaptManager.PUMIAdapter.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    mesh.convertFromPUMI(domain, domain.AdaptManager.PUMIAdapter, domain.faceList, domain.regList,parallel = comm.size() > 1, dim = domain.nd)
    assert mesh.nElements_global == 670
    assert mesh.nNodes_global == 190
    assert mesh.nEdges_global == 977
    assert mesh.nElementBoundaries_global == 1458
    #Ideally, we can assert whether each rank gets the proper number of mesh entities/
    #With the present setup, that information doesn't seem accessible
    #if(comm.rank()==0):
    #    eq(mesh.nElements_owned,4074)
    #    eq(mesh.nNodes_owned,994)
    #    eq(mesh.nEdges_owned,8729)
    #    eq(mesh.nElementBoundaries_owned,5648)

def test_2DparallelLoadPUMI(verbose=0):
    """Test to load 2D parallel PUMI model and mesh"""
    comm = Comm.init()
    if comm.size()!= 2:
        print("Skipping 3DprallelLoadPUMI,must be run on 2 MPI tasks")
        return
    testDir=os.path.dirname(os.path.abspath(__file__))
    domain = Domain.PUMIDomain(dim=2,manager=MeshAdapt.AdaptManager())
    Model=testDir+ '/squarepar.dmg'
    Mesh=testDir + '/squarepar.smb'
    #domain.PUMIMesh=MeshAdapt.MeshAdaptPUMI()
    domain.AdaptManager.PUMIAdapter.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    mesh = MeshTools.TriangularMesh()
    mesh.cmesh = cmeshTools.CMesh()
    mesh.convertFromPUMI(domain, domain.AdaptManager.PUMIAdapter, domain.faceList, domain.regList,parallel = comm.size() > 1, dim = domain.nd)
    assert mesh.nElements_global == 40
    assert mesh.nNodes_global == 29
    assert mesh.nEdges_global == 68
    assert mesh.nElementBoundaries_global == 68
