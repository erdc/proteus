#!/usr/bin/env python
from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdaptPUMI
from nose.tools import eq_ as eq
from nose.tools import ok_ as ok
import os

def test_meshLoadPUMI(verbose=0):
    """Test to load serial PUMI model and mesh"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    cubeMdl=testDir + '/cube.dmg'
    cube670p1=testDir + '/cube.smb'
    meshAdaptInstance = MeshAdaptPUMI.MeshAdaptPUMI()
    meshAdaptInstance.loadModelAndMesh(cubeMdl, cube670p1)
    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    meshAdaptInstance.constructFromSerialPUMIMesh(mesh.cmesh)
    cmeshTools.allocateGeometricInfo_tetrahedron(mesh.cmesh)
    cmeshTools.computeGeometricInfo_tetrahedron(mesh.cmesh)
    mesh.buildFromC(mesh.cmesh)
    eq(mesh.nElements_global,670)
    eq(mesh.nNodes_global,190)
    eq(mesh.nEdges_global,977)
    eq(mesh.nElementBoundaries_global,1458)

if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_meshLoad:test_meshLoadPUMI')

