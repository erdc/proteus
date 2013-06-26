#!/usr/bin/env python
"""
Test module for linear boundary value problems in 2d

This module solves equations of the form

.. _math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=True
def test_c0p1_serial():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    nList[0].femSpaces = {0:default_n.C0_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = False
    nList[0].nn = 21
    nList[0].nLevels = 2
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p1_serial')
    assert(True)
def test_c0p1():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p1"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    #
    nList[0].femSpaces = {0:default_n.C0_AffineLinearOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nn = 21
    nList[0].nLevels = 1
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.node
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    nList[0].cfluxtag  = 'pwl'
    #
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('cg')
    linearSolver.ksp.max_it = 5000
    ns.calculateSolution('poisson_2d_c0p1')
    assert(True)
def test_c0p2_serial():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p2_serial"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    nList[0].femSpaces = {0:default_n.C0_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = False
    nList[0].nn = 21
    nList[0].nLevels = 2
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.element
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('poisson_2d_c0p1_serial')
    assert(True)

def test_c0p2():
    import poisson_het_2d_p
    import poisson_het_2d_c0pk_n
    pList = [poisson_het_2d_p]
    nList = [poisson_het_2d_c0pk_n]
    so = default_so
    so.name = pList[0].name = "poisson_2d_c0p2"+"pe"+`comm.size()`
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    #
    nList[0].femSpaces = {0:default_n.C0_AffineQuadraticOnSimplexWithNodalBasis}
    nList[0].parallel = True
    nList[0].nn = 21
    nList[0].nLevels = 1
    nList[0].linearSolver=default_n.KSP_petsc4py
    nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
    nList[0].parallelPartitioningType = default_n.MeshParallelPartitioningTypes.node
    nList[0].numericalFluxType = default_n.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    nList[0].cfluxtag  = 'pwl-bdm'
    #
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    #configure linear solver
    linearSolver = ns.nlsList[-1].solverList[-1].linearSolver
    assert 'ksp' in dir(linearSolver), "couldn't find ksp in {} with dir {}".format(linearSolver,dir(linearSolver))
    linearSolver.ksp.setType('fgmres')
    linearSolver.ksp.max_it = 10000
    linearSolver.ksp.getPC().setType('asm')
    linearSolver.ksp.getPC().setASMOverlap(2)
    linearSolver.ksp.getPC().setASMType(1) #3--basic, 1 -- restrict, 2 -- interpolate
    ns.calculateSolution('poisson_2d_c0p2')
    assert(True)


if __name__ == '__main__':
    import optparse

    usage = "usage: %prog [options] "
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--parallel",
                      help="testing in parallel? ",
                      action="store_true",
                      dest="parallel",
                      default=False)
    
    (run_opts,args) = parser.parse_args()
    
    if run_opts.parallel:
        test_c0p1()
        test_c0p2()
    else:
        test_c0p1_serial()
        test_c0p2_serial()

    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass
