#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

"""
from proteus.iproteus import *
from proteus.test_utils import TestTools

TestTools.addSubFolders( inspect.currentframe() )

import h5py
import numpy
import pytest
import inspect
import tables

import poiseulle_stokes_2d_p
import poiseulle_stokes_2d_n
import stokesDrivenCavity_2d_p
import stokesDrivenCavity_2d_n

pytestmark = pytest.mark.meshtest

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

class Test2DPoiseulleStokesOnQuads(proteus.test_utils.TestTools.SimulationTest):
    """ Runs a 2D Poiseulle Stokes problem on Quads elements """

    def setup_method(self,method):
        reload(poiseulle_stokes_2d_p)
        reload(poisuelle_stokes_2d_n)
        pList = [poisuelle_stokes_2d_p]
        nList = [poisuelle_stokes_2d_n]    
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
        self.remove_files(Filelist)

    @pytest.mark.skip(reason="WIP unknown regression")
    def test_01_FullRun(self):
        self.ns.calculateSolution('test1')
        if self.ns.ar[0].global_sync:
            relpath = "comparison_files/poiseulleFlow_expected.h5"
        else:
            relpath = "comparison_files/poiseulleFlow_expected.h5"
#        xmf_file = failecmp.cmp('poiseulleFlow.xmf',os.path.join(self._scriptdir,relpath))
        import pdb ;  pdb.set_trace()
        f_expected = h5py.File(os.path.join(self._scriptdir,relpath),'r')
        f_actual = h5py.File('poiseulleFlow.h5','r')
        pressure_expected = f_expected[u'p_analytical1'].value
        pressure_actual = f_actual[u'p_analytical1'].value
        assert xmf_file == True, '******** xmf_file compare failed **********'

class Test2DDrivenCavityStokesOnQuads(proteus.test_utils.TestTools.SimulationTest):

    def setup_method(self,method):
        reload(stokesDrivenCavity_2d_n)
        reload(stokesDrivenCavity_2d_p)
        pList = [stokesDrivenCavity_2d_p]
        nList = [stokesDrivenCavity_2d_n]    
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
        extens = ('edge','ele','neig','node','poly','prof0','log','xmf','h5')
        for currentFile in os.listdir('.'):
            if any(currentFile.endswith(ext) for ext in extens):
                os.remove(currentFile)

    def test_01_C0Q1(self):
        self.ns.calculateSolution('test1')
        f_actual = tables.openFile(os.path.join(self._scriptdir,
                                                'drivenCavityStokesTrial.h5'),
                                   'r')
        f_expected = tables.openFile(os.path.join(self._scriptdir,
                                                  "comparison_files/drivenCavityStokesC0Q1C0Q1.h5"),
                                     'r')

        if (numpy.allclose(f_actual.root.v_t1,
                           f_expected.root.v_t1,atol=1e-4)) == False:
            raise Exception, 'Unexpected C0Q1C0Q1 Driven Cavity Output!'

        if (numpy.allclose(f_actual.root.u_t1,
                           f_expected.root.u_t1,atol=1e-4)) == False:
            raise Exception, 'Unexpected C0Q1C0Q1 Driven Cavity Output!'

        f_actual.close()
        f_expected.close()

        
if __name__ == '__main__':
    pass
