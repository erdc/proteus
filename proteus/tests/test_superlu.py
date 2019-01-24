import os
import math
import pytest
import numpy as np
from numpy import linalg
import scipy.io as sio

from proteus import superluWrappers

@pytest.fixture(scope='module')
def simple_singular_sparse_mat():
    nr = 3 ; nc = 3; nnz = 8
    nzval = np.array([1., 2., 3., 6., -1., 1., 2., 1.], dtype=np.float64)
    rowptr = np.array([0, 2, 5, 8], dtype=np.int32)
    colind = np.array([0, 1, 0, 1, 2, 0, 1, 2], dtype=np.int32)
    A = superluWrappers.SparseMatrix(nr,
                                     nc,
                                     nnz,
                                     nzval,
                                     colind,
                                     rowptr)
    yield A

@pytest.fixture(scope='module')
def simple_nonsingular_sparse_mat():
    nr = 3 ; nc = 3; nnz = 7
    nzval = np.array([1., 2., -1., 2., 4., 1., -1.], dtype=np.float64)
    rowptr = np.array([0, 3, 5, 7], dtype=np.int32)
    colind = np.array([0, 1, 2, 0, 1, 1, 2], dtype=np.int32)
    A = superluWrappers.SparseMatrix(nr,
                                     nc,
                                     nnz,
                                     nzval,
                                     colind,
                                     rowptr)
    yield A
    
@pytest.fixture(scope='module')
def simple_sparse_mat():
    nr = 4 ; nc = 4; nnz = 4
    nzval = np.array([5., 8., 3., 6.], dtype=np.float64)
    rowptr = np.array([0, 0, 2, 3, 4], dtype=np.int32)
    colind = np.array([0, 1, 2, 1], dtype=np.int32)
    A = superluWrappers.SparseMatrix(nr,
                                     nc,
                                     nnz,
                                     nzval,
                                     colind,
                                     rowptr)
    yield A

@pytest.fixture(scope='module')
def small_sparse_mat():
    script_dir = os.path.dirname(__file__)
    mat = sio.mmread(os.path.join(script_dir,'sparse_mat_ex.mtx'))
    mat_csr = mat.tocsr()
    nr = mat_csr.shape[0] ; nc = mat_csr.shape[1] ; nnz = mat_csr.nnz
    nzvals = mat_csr.data ; colind = mat_csr.indices
    rowptr = mat_csr.indptr
    A = superluWrappers.SparseMatrix(nr,
                                     nc,
                                     nnz,
                                     nzvals,
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

def test_SparseMatrix_matvecmult_2(small_sparse_mat):
    A = small_sparse_mat
    xp = np.ones(A.nc) ; yp = np.zeros(A.nr)
    A.matvec(xp,yp)
    assert abs(linalg.norm(yp)-1460.0312) < 0.00001

def test_fwrite_w(simple_sparse_mat):
    A = simple_sparse_mat
    A.fwrite('test.t', base=0)
    s = '''4 4 4 
1 0 5.00000000e+00
1 1 8.00000000e+00
2 2 3.00000000e+00
3 1 6.00000000e+00'''

    with open('test.t', 'r') as file:
        content_as_string = file.read()

    os.remove('test.t')
    assert s.strip()==content_as_string.strip()

def test_sparseFactorPrepare_1(simple_nonsingular_sparse_mat):
    A = simple_nonsingular_sparse_mat
    sparseFactor = superluWrappers.SparseFactor(A.nr)
    superluWrappers.sparseFactorPrepare(A, sparseFactor)
    x = np.ones(A.nr)
    superluWrappers.sparseFactorSolve(sparseFactor, x)
    assert np.array_equal(x, np.array([-.5, .5, -.5]))
