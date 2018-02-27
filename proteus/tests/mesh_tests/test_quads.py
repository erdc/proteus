#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

"""
from proteus.iproteus import *
from proteus.test_utils import TestTools
TestTools.addSubFolders( inspect.currentframe() )

import stokes_2d_p
import stokes_2d_n
import pytest
import inspect

@pytest.mark.MeshTools
def test_mesh_build():
    """  Test mesh generation and refinment """
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
                                                   parallelPartitioningType)

    assert mlMesh.meshList[0].nElements_global == (nnx-1)*(nny-1), 'Mesh generator has built incorrect number of quads'
    assert mlMesh.meshList[1].nElements_global == 4*(nnx-1)*(nny-1), 'Mesh generator has built incorrect number of quads'

@pytest.mark.modelTest
@pytest.mark.moderateTest
@pytest.mark.MeshTools
class Test2DStokesOnQuads():
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
                    "poiseulleFlow.xmf",
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
