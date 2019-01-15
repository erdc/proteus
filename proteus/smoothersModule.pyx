import numpy as np
import cython
cimport numpy as np
cimport superluWrappers as slw

ctypedef int PROTEUS_LAPACK_INTEGER
# ARB - the compiler macro does not seem to be picking this up...

cdef extern from "smoothers.h":
    void cjacobi_NR_prepare "jacobi_NR_prepare"(slw._SuperMatrix *A, double w, double tol, double* M)
    void cjacobi_NR_solve "jacobi_NR_solve"(slw._SuperMatrix *A, double* M, double* R, int* node_order, double* dX)
    void cnl_jacobi_NR_solve "nl_jacobi_NR_solve"(slw._SuperMatrix *A, double* R, int* node_order, double w, double tol, double* dX)
    void cgauss_seidel_NR_prepare "gauss_seidel_NR_prepare"(slw._SuperMatrix *A, double w, double tol, double* M)
    void cgauss_seidel_NR_solve "gauss_seidel_NR_solve"(slw._SuperMatrix *A, double *M, double *R, int *node_order, double *dX)
    void cnl_gauss_seidel_NR_solve "nl_gauss_seidel_NR_solve"(slw._SuperMatrix *A, double *R, int *node_order, double w, double tol, double *dX)
    void casm_NR_init "asm_NR_init"(slw._SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p)
    void casm_NR_free "asm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots)
    void casm_NR_prepare "asm_NR_prepare"(slw._SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    void casm_NR_solve "asm_NR_solve"(slw._SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    int cbasm_NR_init "basm_NR_init"(int rowBlocks, slw._SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p, PROTEUS_LAPACK_INTEGER*** subdomain_col_pivots_p)
    void cbasm_NR_free "basm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots, PROTEUS_LAPACK_INTEGER** subdomain_col_pivots)
    void cbasm_NR_prepare "basm_NR_prepare"(int rowBlocks, int N, slw._SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)
    void cbasm_NR_solve "basm_NR_solve"(int rowBlocks, int N, slw._SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)

cdef struct ASMFactor:
    np.int_t N
    np.int_t *subdomain_dim
    np.int_t **l2g_L
    np.float64_t **subdomain_L, **subdomain_R, **subdomain_dX
    PROTEUS_LAPACK_INTEGER **subdomain_pivots

cdef struct BASMFactor:
    np.int_t N
    np.int_t bs
    np.int_t *subdomain_dim
    np.int_t **l2g_L
    np.float64_t **subdomain_L, **subdomain_R, **subdomain_dX
    PROTEUS_LAPACK_INTEGER **subdomain_pivots
    PROTEUS_LAPACK_INTEGER **subdomain_col_pivots

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

cdef void smootherWrappersjacobi_NR_prepare(slw.cSparseMatrix sm,
                                            double w,
                                            double tol,
                                            np.float64_t [:] M):
    cdef slw._SuperMatrix* AS
    AS.Stype = slw._SLU_NR
    AS.Dtype = slw._SLU_D
    AS.Mtype = slw._SLU_GE
    AS.nrow = sm.nr
    AS.ncol = sm.nc
    AS.Store = &sm.A
    cjacobi_NR_prepare(AS, w, tol, &M[0])

def jacobi_NR_solve(A, M, R, node_order, dX):
    """
    
    Arguments
    ---------
    A :
    M :
    R :
    node_order :
    dX :
    """
    pass

def nl_jacobi_NR_solve(A, R, node_order, w, tol, dX):
    """

    Arguments
    ---------
    A :
    R :
    node_order :
    w :
    tol :
    dX :
    """
    pass

def gauss_seidel_NR_preare(A, w, tol, M):
    """

    Arguments
    ---------
    A :
    w :
    tol :
    M :
    """
    pass

def gauss_seidel_NR_solve(A, M, R, node_order, dX):
    """

    Arguments
    ---------
    A :
    M :
    R :
    node_order :
    dX :
    """
    pass

def nl_gauss_seidel_NR_solve(A, R, node_order, w, tol, dX):
    """
    
    Arguments
    ---------
    A :
    R :
    node_order :
    w :
    tol :
    dX :
    """
    pass

def asm_NR_prepare(A, asmFactor):
    """

    Arguments
    ---------
    A :
    asmFactor :
    """
    pass

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
