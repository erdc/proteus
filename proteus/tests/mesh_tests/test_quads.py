#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

"""
from __future__ import division
from builtins import range
from builtins import object
from past.utils import old_div
from proteus.iproteus import *
from proteus.test_utils import TestTools
TestTools.addSubFolders( inspect.currentframe() )

import stokes_2d_p
import stokes_2d_n
import pytest
import inspect
import numpy as np

@pytest.fixture(scope="module")
def simple_mesh():
    nnx = 4
    nny = 4
    x = (-1.0,-1.0)
    L = (2.0,2.0)
    refinementLevels = 2
    nLayersOfOverlap = 1
    parallelPartitioningType = 0
    skipInit = False
    mlMesh = MeshTools.MultilevelQuadrilateralMesh(nnx,nny,1,
                                                   x[0],x[1],0.0,
                                                   L[0],L[1],1,
                                                   refinementLevels,
                                                   skipInit,
                                                   nLayersOfOverlap,
                                                   parallelPartitioningType,
                                                   useC=False)
    yield mlMesh, nnx, nny

@pytest.fixture(scope="module")
def simple_mesh_with_c():
    nnx = 4
    nny = 4
    x = (-1.0,-1.0)
    L = (2.0,2.0)
    refinementLevels = 1
    nLayersOfOverlap = 1
    parallelPartitioningType = 0
    skipInit = False
    mlMesh = MeshTools.MultilevelQuadrilateralMesh(nnx,nny,1,
                                                   x[0],x[1],0.0,
                                                   L[0],L[1],1,
                                                   refinementLevels,
                                                   skipInit,
                                                   nLayersOfOverlap,
                                                   parallelPartitioningType,
                                                   useC=True)
    yield mlMesh, nnx, nny

@pytest.mark.MeshTools
def test_mesh_build(simple_mesh):
    """  Test mesh generation and refinment """
    mlMesh,nnx,nny = simple_mesh
    assert mlMesh.meshList[0].nElements_global == (nnx-1)*(nny-1), 'Mesh generator has built incorrect number of quads'
    assert mlMesh.meshList[1].nElements_global == 4*(nnx-1)*(nny-1), 'Mesh generator has built incorrect number of quads'

@pytest.mark.MeshTools
def test_mesh_build_1(simple_mesh_with_c):
    """ Test mesh generation with c """
    mlMesh,nnx,nny = simple_mesh_with_c
    assert mlMesh.meshList[0].nElements_global == (nnx-1)*(nny-1), 'Mesh generator has built incorrect number of quads'
    # TODO (ARB) - add an additional test for refinement when its ready

@pytest.mark.MeshTools
def test_calc_quad_area(simple_mesh):
    mlMesh, nnx, nny = simple_mesh
    for i in range(9):
        assert mlMesh.meshList[0]._calc_quad_area(i) == 4. / 9.

@pytest.mark.MeshTools
def test_calc_quad_area(simple_mesh_with_c):
    mlMesh, nnx, nny = simple_mesh_with_c
    for i in range(9):
        assert mlMesh.meshList[0]._calc_quad_area(i) == 4. / 9.
        
@pytest.mark.MeshTools
def test_calc_hmax(simple_mesh):
    mlMesh, nnx, nny = simple_mesh
    quad_mesh = mlMesh.meshList[0]
    for i in range(9):
        hmax_i = quad_mesh._calc_hmax(i)
        h = quad_mesh._calc_pt_distance((-1.0,old_div(-1.,3.)),(old_div(-1.,3.),old_div(1.,3.)))
        assert abs(h-hmax_i) < 1e-12

@pytest.mark.MeshTools
def test_buildNodeDiameterArray_1(simple_mesh):
    mlMesh, nnx, nny = simple_mesh
    quad_mesh = mlMesh.meshList[0]
    quad_mesh.buildNodeDiameterArray()
    expected = np.full((16),0.94280904)
    assert np.allclose(quad_mesh.nodeDiametersArray,expected)
    assert abs(quad_mesh.volume-4.0) < 1e-12

@pytest.mark.modelTest
@pytest.mark.moderateTest
@pytest.mark.MeshTools
class Test2DStokesOnQuads(object):
    """ Runs a 2D Poiseulle Stokes problem on Quads with TH elements """

    @classmethod    
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        reload(stokes_2d_p)
        reload(stokes_2d_n)
        pList = [stokes_2d_p]
        nList = [stokes_2d_n]    
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList=[default_s]
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        self._scriptdir = os.path.dirname(__file__)
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)

    def teardown_method(self,method):
        Filelist = ["rdomain.edge",
                    "rdomain.ele",
                    "rdomain.neig",
                    "rdomain.node",
                    "rdomain.poly",
                    "reference_triangle.poly",
                    "reference_simplex.poly",
                    "proteus.log",
#                    "poiseulleFlow.xmf",
                    "poiseulleFlow.h5",
                    "poiseulleFlow0.h5"]
        TestTools.removeFiles(Filelist)
       
    def test_01_FullRun(self):
        import filecmp
        self.ns.calculateSolution('test1')
        if self.ns.ar[0].global_sync:
            relpath = "comparison_files/poiseulle_global_xmf.output"
        else:
            relpath = "comparison_files/poiseulle_xmf.output"
        xmf_file = filecmp.cmp('poiseulleFlow.xmf',os.path.join(self._scriptdir,relpath))
        assert xmf_file == True, '******** xmf_file compare failed **********'

if __name__ == '__main__':
    pass
