import numpy as np
import cython
cimport numpy as np
from proteus cimport superluWrappers

ctypedef int PROTEUS_LAPACK_INTEGER
# ARB - the compiler macro does not seem to be picking this up...
ctypedef double [:] DDATA
ctypedef int [:] IDATA
ctypedef superluWrappers._SuperMatrix SuperMatrix

cdef extern from "smoothers.h":
    void cjacobi_NR_prepare "jacobi_NR_prepare"(superluWrappers._SuperMatrix *A, double w, double tol, double* M)
    void cjacobi_NR_solve "jacobi_NR_solve"(superluWrappers._SuperMatrix *A, double* M, double* R, int* node_order, double* dX)
    void cnl_jacobi_NR_solve "nl_jacobi_NR_solve"(superluWrappers._SuperMatrix *A, double* R, int* node_order, double w, double tol, double* dX)
    void cgauss_seidel_NR_prepare "gauss_seidel_NR_prepare"(superluWrappers._SuperMatrix *A, double w, double tol, double* M)
    void cgauss_seidel_NR_solve "gauss_seidel_NR_solve"(superluWrappers._SuperMatrix *A, double *M, double *R, int *node_order, double *dX)
    void cnl_gauss_seidel_NR_solve "nl_gauss_seidel_NR_solve"(superluWrappers._SuperMatrix *A, double *R, int *node_order, double w, double tol, double *dX)
    int casm_NR_init "asm_NR_init"(superluWrappers._SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p)
    void casm_NR_free "asm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots)
    void casm_NR_prepare "asm_NR_prepare"(superluWrappers._SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    void casm_NR_solve "asm_NR_solve"(superluWrappers._SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    int cbasm_NR_init "basm_NR_init"(int rowBlocks, superluWrappers._SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p, PROTEUS_LAPACK_INTEGER*** subdomain_col_pivots_p)
    void cbasm_NR_free "basm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots, PROTEUS_LAPACK_INTEGER** subdomain_col_pivots)
    void cbasm_NR_prepare "basm_NR_prepare"(int rowBlocks, int N, superluWrappers._SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)
    void cbasm_NR_solve "basm_NR_solve"(int rowBlocks, int N, superluWrappers._SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)

class ASMFactor(object):

    def __init__(self, L):
        self.L = L
        self._cASMFactor = cASMFactor(self.L._cSparseMatrix)

cdef class cASMFactor(object):
    cdef int N
    cdef int *subdomain_dim
    cdef int **l2g_L
    cdef double **subdomain_L
    cdef double **subdomain_R
    cdef double **subdomain_dX
    cdef PROTEUS_LAPACK_INTEGER **subdomain_pivots

    def __cinit__(self,
                  superluWrappers.cSparseMatrix L):
        cdef SuperMatrix AS
        AS.Stype = superluWrappers._SLU_NR
        AS.Dtype = superluWrappers._SLU_D
        AS.Mtype = superluWrappers._SLU_GE
        AS.nrow = L.nr
        AS.ncol = L.nc
        AS.Store = &L.A
        if casm_NR_init(&AS,
                        &self.subdomain_dim,
                        &self.l2g_L,
                        &self.subdomain_L,
                        &self.subdomain_R,
                        &self.subdomain_dX,
                        &self.subdomain_pivots):
            assert 1==0 

class BASMFactor(object):

    def __init__(self, L, bs):
        self.L = L
        self.bs = bs
        self._cBASMFactor = cBASMFactor(self.L._cSparseMatrix,
                                        self.bs)

cdef class cBASMFactor(object):
    cdef int N
    cdef int bs
    cdef int *subdomain_dim
    cdef int **l2g_L
    cdef double **subdomain_L
    cdef double **subdomain_R
    cdef double **subdomain_dX
    cdef PROTEUS_LAPACK_INTEGER **subdomain_pivots
    cdef PROTEUS_LAPACK_INTEGER **subdomain_col_pivots

    def __cinit__(self,
                  superluWrappers.cSparseMatrix L,
                  int bs):
        cdef SuperMatrix AS
        AS.Stype = superluWrappers._SLU_NR
        AS.Dtype = superluWrappers._SLU_D
        AS.Mtype = superluWrappers._SLU_GE
        AS.nrow = L.nr
        AS.ncol = L.nc
        AS.Store = &L.A
        if cbasm_NR_init(bs,
                         &AS,
                         &self.subdomain_dim,
                         &self.l2g_L,
                         &self.subdomain_L,
                         &self.subdomain_R,
                         &self.subdomain_dX,
                         &self.subdomain_pivots,
                         &self.subdomain_col_pivots):
            assert 1==0
    
def jacobi_NR_prepare(A, w, tol, M):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    w : double
    tol : double
    M : np.array
    """
    smootherWrappersjacobi_NR_prepare(A._cSparseMatrix, w, tol, M)

cdef void smootherWrappersjacobi_NR_prepare(superluWrappers.cSparseMatrix sm,
                                            double w,
                                            double tol,
                                            DDATA M):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cjacobi_NR_prepare(&AS, w, tol, &M[0])

def jacobi_NR_solve(A, M, R, node_order, dX):
    """
    
    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    M : np.array double
    R : np.array double
    node_order : np.array int
    dX : np.array double
    """
    smootherWrappersjacobi_NR_solve(A._cSparseMatrix, M, R, node_order, dX)

cdef void smootherWrappersjacobi_NR_solve(superluWrappers.cSparseMatrix sm,
                                          DDATA M,
                                          DDATA R,
                                          IDATA node_order,
                                          DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cjacobi_NR_solve(&AS, &M[0], &R[0], &node_order[0], &dX[0])

def nl_jacobi_NR_solve(A, R, node_order, w, tol, dX):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    R : np.array double
    node_order : np.array int
    w : np.float
    tol : np.float
    dX : np.array double
    """
    smootherWrappersnl_jacobi_NR_solve(A._cSparseMatrix, R, node_order, w, tol, dX)

cdef void smootherWrappersnl_jacobi_NR_solve(superluWrappers.cSparseMatrix sm,
                                             DDATA R,
                                             IDATA node_order,
                                             double w,
                                             double tol,
                                             DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cnl_jacobi_NR_solve(&AS, &R[0], &node_order[0], w, tol, &dX[0])

def gauss_seidel_NR_preare(A, w, tol, M):
    """

    Arguments
    ---------
    A :
    w :
    tol :
    M :
    """
    smootherWrappersgauss_seidel_NR_prepare(A._cSparseMatrix, w, tol, M)

cdef void smootherWrappersgauss_seidel_NR_prepare(superluWrappers.cSparseMatrix sm,
                                                  double w,
                                                  double tol,
                                                  DDATA M):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cgauss_seidel_NR_prepare(&AS, w, tol, &M[0])

def gauss_seidel_NR_solve(A, M, R, node_order, dX):
    """

    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    M : np.array double
    R : np.array double
    node_order : np.array int
    dX : np.array double
    """
    smootherWrappersgauss_seidel_NR_solve(A._cSparseMatrix, M, R, node_order, dX)

cdef void smootherWrappersgauss_seidel_NR_solve(superluWrappers.cSparseMatrix sm,
                                                DDATA M,
                                                DDATA R,
                                                IDATA node_order,
                                                DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cgauss_seidel_NR_solve(&AS, &M[0], &R[0], &node_order[0], &dX[0])

def nl_gauss_seidel_NR_solve(A, R, node_order, w, tol, dX):
    """
    
    Arguments
    ---------
    A : superluWrappers.SparseMatrix
    R : np.array double
    node_order : np.array int
    w : double
    tol : double
    dX : np.array double
    """
    smootherWrappers_nl_gauss_seidel_NR_solve(A._cSparseMatrix,
                                              R,
                                              node_order,
                                              w,
                                              tol,
                                              dX)

cdef smootherWrappers_nl_gauss_seidel_NR_solve(superluWrappers.cSparseMatrix sm,
                                               DDATA R,
                                               IDATA node_order,
                                               double w,
                                               double tol,
                                               DDATA dX):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cnl_gauss_seidel_NR_solve(&AS, &R[0], &node_order[0], w, tol, &dX[0])

def asm_NR_prepare(A, asmFactor):
    """

    Arguments
    ---------
    A :
    asmFactor :
    """
    smootherWrappers_asm_NR_prepare(A._cSparseMatrix,
                                    asmFactor._cASMFactor)
    
cdef void smootherWrappers_asm_NR_prepare(superluWrappers.cSparseMatrix sm,
                                          cASMFactor asmFactor):
    cdef SuperMatrix AS
    AS.Stype = superluWrappers._SLU_NR
    AS.Dtype = superluWrappers._SLU_D
    AS.Mtype = superluWrappers._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    casm_NR_prepare(&AS,
                    asmFactor.subdomain_dim,
                    asmFactor.l2g_L,
                    asmFactor.subdomain_L,
                    asmFactor.subdomain_pivots)

def asm_NR_solve(A, w, asmFactor, node_order, R, dX):
    """
    
    Arguments
    ---------
    A :
    w :
    asmFactor :
    node_order :
    R :
    dX :
    """
    pass

def basm_NR_prepare(A, baseFactor):
    """

    Arguments
    ---------
    A :
    basmFactor :
    """
    pass

def basm_NR_solve(A, w, basmFactor, node_order, R, dX):
    """
    
    Arguments
    ---------
    A :
    w :
    basmFactor :
    node_order :
    R :
    dX :
    """
    pass
