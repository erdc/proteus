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
from petsc4py import PETSc as p4pyPETSc

def test_Qp_mat():
    '''
    First, verify that Qp returns the correct pressure mass matrix
    Second, verfity that QpShell performs matrix multiplies correctly
    '''
    # import and run a small 2D poiseulle problem
    import stokes_2d_p
    import stokes_2d_n
    import numpy as np
    from scipy.sparse import csr_matrix
    pList = [stokes_2d_p]
    nList = [stokes_2d_n]    
    so = default_so
    so.tnList = [0.,1.]
    so.name = pList[0].name
    so.sList=[default_s]
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    # *** 1 *** : test whether the pressure mass matrix is being created properly.
    matrix = LinearSolvers.NavierStokes3D(L=ns.modelList[0].par_jacobianList[1],schurPC='Qp')
    matrix.setUp()
    # Maybe this matrix should be stored in a seperate file that is loaded.
    pressure_mass_matrix = np.array([[0.08333333, 0., 0., 0., 0.04166667, 0., 0.04166667, 0., 0.],
                                     [0.,0.16666667, 0., 0., 0., 0.08333333, 0.04166667, 0.04166667, 0.],
                                     [0., 0., 0.08333333, 0., 0., 0., 0., 0.04166667,  0.04166667],
                                     [0., 0., 0., 0.16666667, 0.04166667, 0.08333333, 0., 0., 0.04166667],
                                     [0.04166667, 0., 0., 0.04166667, 0.25, 0.08333333, 0.08333333, 0., 0.],
                                     [0., 0.08333333, 0., 0.08333333, 0.08333333, 0.5, 0.08333333, 0.08333333, 0.08333333],
                                     [0.04166667, 0.04166667, 0., 0., 0.08333333, 0.08333333, 0.25, 0., 0.],
                                     [0., 0.04166667, 0.04166667, 0., 0., 0.08333333, 0., 0.25, 0.08333333],
                                     [0., 0., 0.04166667, 0.04166667, 0., 0.08333333, 0., 0.08333333, 0.25]])

    Qp_raw = matrix.Qp.getValuesCSR()
    Qp_dense = LinearAlgebraTools.petsc4py_sparse_2_dense(Qp_raw)
    assert np.allclose(pressure_mass_matrix,Qp_dense)
    # *** 2 *** verify QpShell performs a matrix vector product correctly.
    # create a vector of 1's to multiply the vector by and an empty vector the
    # solution in in.
    vector_x = p4pyPETSc.Vec().create()
    vector_y = p4pyPETSc.Vec().create()
    vector_x.setType('standard')
    vector_y.setType('standard')
    vector_x.setSizes(len(pressure_mass_matrix[0]),len(pressure_mass_matrix[0]))
    vector_y.setSizes(len(pressure_mass_matrix[0]),len(pressure_mass_matrix[0]))
    vector_x.set(1)
    matrix.Qp_shell.mult(vector_x,vector_y)
    # TODO - verify true_y
    true_y = np.array([[0.16666667,0.333333333,0.166666667,0.33333333,0.5,1.,0.5,0.5,0.5]])
    assert np.allclose(vector_y,true_y)


if __name__ == '__main__':
    from proteus import Comm
    comm = Comm.init()
    test_Qp_mat()

    # test = LinearSolvers.KSP_petsc4py(ns.modelList[0].jacobianList[0],
    #                                   ns.modelList[0].par_jacobianList[0],
    #                                   Preconditioner=SimpleNavierStokes3D)
