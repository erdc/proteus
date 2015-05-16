#!/usr/bin/env python
from proteus.MeshTools import readMeshXdmf,writeHexMesh

import os
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

def test_3x3_cube(verbose=0):
    """
    Read sample openfoam mesh from aggelos and check that the basic information is correct
    """
    xmf_base=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'hex_cube_3x3')
    h5_base = xmf_base
    mesh_info = readMeshXdmf(xmf_base,h5_base,verbose=0)

    eq(mesh_info.nElements_global,27)
    eq(mesh_info.nNodes_global,64)
    
    eq(mesh_info.nElements_global,mesh_info.nElements_owned)
    eq(mesh_info.nNodes_global,mesh_info.nNodes_owned)

    eq(mesh_info.nodeArray.shape,(mesh_info.nNodes_owned,3))
    eq(mesh_info.elementNodesArray.shape,(mesh_info.nElements_owned,8))
    
    eq(mesh_info.elementTopologyName,'Hexahedron')

    eq(len(mesh_info.nodeMaterialTypes),mesh_info.nNodes_owned)
    eq(len(mesh_info.elementMaterialTypes),mesh_info.nElements_owned)

def test_write_3x3_cube(verbose=0):
    """
    Read sample openfoam mesh from aggelos and try to write in Ido's hex format
    """
    xmf_base=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'hex_cube_3x3')
    h5_base = xmf_base
    mesh_info = readMeshXdmf(xmf_base,h5_base,verbose=0)
    
    writeHexMesh(mesh_info,'hexMesh_3x3',index_base=0)

if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_read_mesh_xml:test_3x3_cube,test_read_mesh_xml:test_write_3x3_cube')
