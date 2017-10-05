from proteus import Comm, Profiling
import numpy as np
import numpy.testing as npt
import pytest
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

comm = Comm.init()
Profiling.procID = comm.rank()

class MockMat():
    """ petsc4py-based mock SuperLU Matrix for testing.  Filled like this:

    [ 1 1     ]
    [   1 2   ]
    [     1 3 ]
    [       1 ]
    """


    def __init__(self, n):
        from petsc4py import PETSc
        A = PETSc.Mat().create()
        A.setSizes(n)
        A.setType('aij')
        A.setUp()
        A.assemblyBegin()
        for i in range(n):
            A[i, i] = 1
            if i + 1 < n:
                A[i, i + 1] = i
        A.assemblyEnd()
        self.A = A

    def getCSRrepresentation(self):
        return self.A.getValuesCSR()

    def getSubMatCSRrepresentation(self, start, end):
        from petsc4py import PETSc
        ids = range(start, end)
        isg = PETSc.IS().createGeneral(ids)
        B = self.A.getSubMatrix(isg)
        return B.getValuesCSR()

@pytest.mark.LinearAlgebraTools
def test_vec_create():
    """test_vec_create

    Verifies that the proteus.LinearAlgebraTools.Vec constructor
    correctly creates one-dimensional arrays of the given length of
    type double precision and with entries set to zero for several
    trials.
    """
    from proteus.LinearAlgebraTools import Vec
    for n in [1, 10, 100, 1000]:
        x = Vec(n)
        # Vector of length n
        eq(x.size, n)
        # One-dimensional
        eq(x.shape, (n,))
        # Of type double-precision
        eq(x.dtype, np.double)
        # All entries are zero
        eq(np.count_nonzero(x), 0)
        # Verify assignment works
        x[:] = range(1, n+1)
        eq(np.count_nonzero(x), n)

@pytest.mark.LinearAlgebraTools
def test_mat_create():
    """test_mat_create

    Verifies that the proteus.LinearAlgebraTools.Mat constructor
    correctly creates double-precision matrices of the given
    dimensions and with entries set to zero for several trials.
    """

    from proteus.LinearAlgebraTools import Mat
    for m, n in [(1, 1), (1, 10), (100, 1), (500, 500)]:
        x = Mat(m, n)
        # Matrix containing m*n entries
        eq(x.size, m*n)
        # Two-dimensional
        eq(x.shape, (m, n))
        # Of type double-precision
        eq(x.dtype, np.double)
        # All entries are zero
        eq(np.count_nonzero(x), 0)
        # Assign a row
        x[0, :] = range(1, n+1)
        # Assign a column
        x[:, 0] = range(1, m+1)
        eq(np.count_nonzero(x), m+n-1)

@pytest.mark.LinearAlgebraTools
def test_vec_scalar_math():
    """test_vec_scalar_math

    Verifies that the Vec object created by
    proteus.LinearAlgebraTools.Vec behaves as expected when computing
    basic linear algebra for several trials.
    """

    from proteus.LinearAlgebraTools import Vec

    v1 = Vec(2)
    v1[:] = [1, 2]

    npt.assert_almost_equal(v1.sum(), 3)
    npt.assert_almost_equal(v1.prod(), 2)

    v2 = Vec(2)
    v2[:] = [2, 0.5]

    # v1 + v2
    v_sum = np.asarray([3, 2.5])
    npt.assert_almost_equal(v_sum, v1+v2)

    # v1 + 2*v2
    v_scaled_sum = np.asarray([5, 3])
    npt.assert_almost_equal(v_scaled_sum, v1+2*v2)

    # 0.5*v1 + 1
    v_scalar_sum = np.asarray([1.5, 2])
    npt.assert_almost_equal(v_scalar_sum, 0.5*v1+1)

@pytest.mark.LinearAlgebraTools
def test_mat_vec_math():
    """test_mat_vec_math

    Verifies that the Mat and Vec objects from
    proteus.LinearAlgebraTools behave as expected together when
    computing basic linear algebra for one trial.
    """

    from proteus.LinearAlgebraTools import Vec
    from proteus.LinearAlgebraTools import Mat

    v1 = Vec(2)
    v1[:] = [1, 1]

    m1 = Mat(3, 2)
    m1[0, :] = [1, 2]
    m1[1, :] = [3, 2]
    m1[2, :] = [4, 5]

    # m1*v1
    dot_product = np.asarray([3, 5, 9])
    npt.assert_almost_equal(dot_product, m1.dot(v1))

@pytest.mark.LinearAlgebraTools
def test_superlu_mat():
    """test_superlu_mat

    Simple verification of functionality of SuperLU SparseMat operators.
    """


def compute_norms(h, A, vecs):
    from proteus.LinearAlgebraTools import l2Norm, l1Norm, lInfNorm, rmsNorm
    from proteus.LinearAlgebraTools import wl2Norm, wl1Norm, wlInfNorm
    from proteus.LinearAlgebraTools import energyNorm
    for norm in [l2Norm, l1Norm, lInfNorm, rmsNorm]:
        yield norm.__name__, [norm(x) for x in vecs]
    for norm in [wl2Norm, wl1Norm, wlInfNorm]:
        yield norm.__name__, [norm(x, h) for x in vecs]
    for norm in [energyNorm]:
        yield norm.__name__, [energyNorm(x, A) for x in vecs]

@pytest.mark.LinearAlgebraTools
def test_norm_correctness():
    from math import sqrt
    from proteus.LinearAlgebraTools import Mat
    from proteus.LinearAlgebraTools import l2Norm, l1Norm, lInfNorm, rmsNorm
    from proteus.LinearAlgebraTools import wl2Norm, wl1Norm, wlInfNorm
    from proteus.LinearAlgebraTools import energyNorm

    x = np.ones(2)
    h = np.ones(2)
    A = Mat(2, 2)
    A[0, 0] = 1
    A[1, 1] = 1

    assert l1Norm(x)==2
    assert l2Norm(x)==sqrt(2)
    assert lInfNorm(x)==1
    assert rmsNorm(x)==1
    assert wl1Norm(x,h)==2
    assert wl2Norm(x,h)==sqrt(2)
    assert wlInfNorm(x,h)==1
    assert energyNorm(x,A)==sqrt(2)

@pytest.mark.LinearAlgebraTools
def test_norm_zero():
    """test_norm_zero

    Norm of a zero vector better be close to zero.
    """

    from proteus.LinearAlgebraTools import Vec
    from proteus.LinearAlgebraTools import Mat

    h = Vec(2)
    h[:] = [0.5, 0.5]

    A = Mat(2, 2)
    #needs to be SPD
    A[0, :] = [5, 2]
    A[1, :] = [2, 6]

    v = Vec(2)
    v[:] = [0, 0]

    for name, norms in compute_norms(h, A, [v]):
        n = norms[0]
        assert np.allclose(n,0)

@pytest.mark.LinearAlgebraTools
def test_norm_homogeneity():
    """test_norm_homogeneity

    Test if defined norms in proteus.LinearAlgebraTools obey
    absolute homoegeneity for several trials.
    """

    from proteus.LinearAlgebraTools import Vec
    from proteus.LinearAlgebraTools import Mat

    h = Vec(2)
    h[:] = [0.5,0.5]

    A = Mat(2,2)
    #needs to be SPD
    A[0, :] = [5, 2]
    A[1, :] = [2, 6]

    v1 = Vec(2)
    v1[:] = [1, 1]

    for a in [0.5, -2]:
        for name, norms in compute_norms(h, A, [v1, a*v1]):
            t1, t2 = norms
            assert np.allclose(abs(a)*t1,t2)

@pytest.mark.LinearAlgebraTools
def test_norm_triangle_inequality():
    """test_norm_triangle_inequality

    Test if defined norms in proteus.LinearAlgebraTools obey
    the triangle inequality for several trials.
    """

    from proteus.LinearAlgebraTools import Vec
    from proteus.LinearAlgebraTools import Mat

    h = Vec(2)
    h[:] = [1.0, 1.0]

    #need an SPD matrix for this to be a norm
    A = Mat(2, 2)
    A[0, :] = [4, 2]
    A[1, :] = [2, 5]

    v1 = Vec(2)
    v1[:] = [1, 1]

    v2 = Vec(2)
    v2[:] = [-1, 4]

    for name, norms in compute_norms(h, A, [v1+v2, v1, v2]):
        t1, t2, t3 = norms
        assert t1 <= t2 + t3

@pytest.mark.LinearAlgebraTools
def test_petsc_binary_mat_io():
    """test_petsc_binary_mat_io

    Verifies that binary save/loads on PETSc Mats works correctly.
    """

    from proteus.LinearAlgebraTools import ParMat_petsc4py

    n = 10
    par_bs = 1
    par_n = n
    par_N = n
    par_nghost = 0
    subdomain2global = range(n)
    ghosted_csr_mat = MockMat(n)

    M = ParMat_petsc4py(ghosted_csr_mat, par_bs, par_n, par_N, par_nghost, subdomain2global)

    M.save('M.petsc_mat')

    from petsc4py import PETSc

    viewer = PETSc.Viewer().createBinary('M.petsc_mat', 'r')
    A = PETSc.Mat().load(viewer)

    mi, mj, mv = M.csr_rep
    ai, aj, av = A.getValuesCSR()
    npt.assert_equal(mi, ai)
    npt.assert_equal(mj, aj)
    npt.assert_equal(mv, av)

if __name__ == '__main__':
    import nose
    nose.main()
