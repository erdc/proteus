import math
import pytest
import numpy as np
from proteus import superluWrappers

# def test_SparseFactor_initialization():
#     test_mat = superluWrappers.SparseFactor(1)
#     assert type(test_mat) == superluWrappers.SparseFactor

@pytest.fixture(scope='module')
def simple_sparse_mat():
    nr = 4 ; nc = 4; nnz = 4
    nzval = np.array([5., 8., 3., 6.])
    rowptr = np.array([0, 0, 2, 3, 4])
    colind = np.array([0, 1, 2, 1])
    A = superluWrappers.SparseMatrix(nr,
                                     nc,
                                     nnz,
                                     nzval,
                                     colind,
                                     rowptr)
    yield A

def test_SparseMatrix_1(simple_sparse_mat):
    A = simple_sparse_mat
    nzval = np.array([5., 8., 3., 6.])
    rowptr = np.array([0, 0, 2, 3, 4])
    colind = np.array([0, 1, 2, 1])    
    # test full csr_rep
    A_csr_rep = A.getCSRrepresentation()
    assert np.array_equal(A_csr_rep[0], rowptr)
    assert np.array_equal(A_csr_rep[1], colind)
    assert np.array_equal(A_csr_rep[2], nzval)
    # test sub csr rep using all rows
    A_sub_csr_rep = A.getSubMatCSRrepresentation(0,4)
    assert np.array_equal(A_sub_csr_rep[0], rowptr)
    assert np.array_equal(A_sub_csr_rep[1], colind)
    assert np.array_equal(A_sub_csr_rep[2], nzval)
    # test sub csr for last two rows
    A_sub_csr_rep = A.getSubMatCSRrepresentation(3,4)
    assert np.array_equal(A_sub_csr_rep[0], np.array([3,4]))
    assert np.array_equal(A_sub_csr_rep[1], np.array([1]))
    assert np.array_equal(A_sub_csr_rep[2], np.array([6]))
    # test sub csr for middle two rows
    A_sub_csr_rep = A.getSubMatCSRrepresentation(1,3)
    assert np.array_equal(A_sub_csr_rep[0], np.array([0,2,3]))
    assert np.array_equal(A_sub_csr_rep[1], np.array([0,1,2]))
    assert np.array_equal(A_sub_csr_rep[2], np.array([5.,8.,3.]))
    # TODO - ARB : there are different forms in which one
    # may expect the output of the submatrix to appear
    # I believe this is correct, but this should be verified
    # before merging

def test_SparseMatrix_matvecmult_1(simple_sparse_mat):
    A = simple_sparse_mat
    xp = np.array([1.,1.,1.,1.])
    yp = np.zeros(4)
    A.matvec(xp,yp)
    assert np.array_equal(yp, np.array([0., 13., 3., 6.]))
    xp = np.array([0., 2., 1., 3.])
    A.matvec(xp,yp)
    assert np.array_equal(yp, np.array([0., 16., 3., 12.]))

