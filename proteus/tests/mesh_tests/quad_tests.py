#!/usr/bin/env python
"""

Test module for 2D Quadrilateral Meshes

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from nose.tools import set_trace

def test_mesh_build():
    '''
    Test mesh generation and refinment
    '''
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

def test_2D_stokes():
    import stokes_2d_p
    import stokes_2d_n
    pList = [stokes_2d_p]
    nList = [stokes_2d_n]    
    so = default_so
    so.tnList = [0.,1.]
    so.name = pList[0].name
    so.sList=[default_s]
    #opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    # nList[0].femSpaces[0]  = default_n.C0_AffineLinearOnSimplexWithNodalBasis
    # nList[0].linearSolver=default_n.KSP_petsc4py
    # nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    # nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #nList[0].linearSolver=default_n.LU
    #nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('test1')
  #  pass


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    import nose
    nose.main()

