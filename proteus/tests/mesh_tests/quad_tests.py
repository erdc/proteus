#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#Profiling.logLevel=7
#Profiling.verbose=True
import numpy.testing as npt
import nose
from nose.plugins.attrib import attr
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

@attr('fast','mesh')
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

@attr('fast','stokes','mesh')
class Test2DStokesOnQuads():
    """ Runs a simple 2D Poiseulle Stokes problem on Quads with TH elements """

    def setUp(self):
        """ Initialize the test problem """
        import stokes_2d_p
        import stokes_2d_n
        pList = [stokes_2d_p]
        nList = [stokes_2d_n]    
        so = default_so
        so.tnList = [0.,1.]
        so.name = pList[0].name
        so.sList=[default_s]
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        self.ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)

    def tearDown(self):
        if os.path.exists('proteus.log'):
            os.remove('proteus.log')
        else:
            pass
        os.remove('poiseulleFlow.xmf')
        os.remove('poiseulleFlow.h5')

    def test_01_FullRun(self):
        import filecmp
        self.ns.calculateSolution('test1')
        xmf_file = filecmp.cmp('poiseulleFlow.xmf','proteus/tests/mesh_tests/poiseulle_xmf.output')
        assert xmf_file == True, '******** xmf_file compare failed **********'

if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    nose.main()

