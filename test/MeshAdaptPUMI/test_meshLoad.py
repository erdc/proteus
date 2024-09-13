#!/usr/bin/env python
from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdapt
from proteus import Domain
import os
import pytest

#@pytest.mark.skip(reason=".smb file was removed during cleaning")
def test_meshLoadPUMI(verbose=0):
    """Test to load serial PUMI model and mesh"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    cubeMdl=testDir + '/cube.dmg'
    cube670p1=testDir + '/cube.smb'

    manager=MeshAdapt.AdaptManager()
    PUMIAdapter=manager.PUMIAdapter
    PUMIAdapter.loadModelAndMesh(cubeMdl.encode('ascii'), cube670p1.encode('ascii'))

    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    PUMIAdapter.constructFromSerialPUMIMesh(mesh.cmesh)
    cmeshTools.allocateGeometricInfo_tetrahedron(mesh.cmesh)
    cmeshTools.computeGeometricInfo_tetrahedron(mesh.cmesh)
    mesh.buildFromC(mesh.cmesh)
    assert mesh.nElements_global == 670
    assert mesh.nNodes_global == 190
    assert mesh.nEdges_global == 977
    assert mesh.nElementBoundaries_global == 1458
