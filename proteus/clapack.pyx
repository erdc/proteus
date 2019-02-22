import numpy as np
cimport numpy as np

cdef extern from "proteus_lapack.h":
    extern int cdgetrf_ "dgetrf_"(int *m, int *n, double *a, int *lda, int *ipiv, int *info)
    extern int cdgetrs_ "dgetrs_"(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info)
    extern int cdgetc2_ "dgetc2_"(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info)
    extern int cdgesc2_ "dgesc2_"(int *n, double *a, int *lda, double* rhs, int *ipiv, int *jpiv, double* scale)
    extern int cdgeev_ "dgeev_"(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr,
                                double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork,int* info)
    extern int cdgetri_ "dgetri_"(int* N,double* A,int* LDA,int* IPIV,double* WORK,int* LWORK,int* INFO )

cdef extern from "proteus_blas.h":
    extern void cdcopy_ "dcopy_"(const int* N, const double *X, const int* incX, double *Y, const int* incY)

class DenseFactor:
    def __init__(self, n):
        self.n = n
        self.lu = np.zeros(shape=(n,n), dtype=float)
        self.pivots = np.zeros(shape=n, dtype='int32')

cdef lapackWrappersDCopy(int n,
                         np.ndarray mat,
                         np.ndarray lu):
    cdef int dim=n, incr1=1, incr2=1, dim2=n**2

    cdcopy_(&dim2,
            <double *> mat.data,
            &incr1,
            <double *> lu.data,
            &incr2)

cdef lapackWrappersDenseFactorPrepare(int dim,
                                      np.ndarray lu,
                                      np.ndarray pivot):
    cdef int info=0, N=dim
    cdgetrf_(&N,
             &N,
             <double*> lu.data,
             &N,
             <int *> pivot.data,
             &info)
#     # ARB - this int should be a PROTEUS_LAPACK_INTEGER, but the complier does not seem to be picking it up

cdef lapackWrappersDenseFactorSolve(int dim,
                                    np.ndarray lu,
                                    np.ndarray pivots,
                                    np.ndarray b):
    cdef char trans ='T'
    # ARB Note - relative to how I would think about this, it seems this needs a transpose?!
    cdef int N=dim, info=0, nrhs=1
    cdgetrs_(&trans,
             &N,
             &nrhs,
             <double *> lu.data,
             &N,
             <int *> pivots.data,
             <double *> b.data,
             &N,
             &info)

cdef lapackWrappersDenseCalculateEigenvalues(char* jobvl,
                                             char* jobvr,
                                             int n,
                                             np.ndarray A,
                                             int lda,
                                             np.ndarray wr,
                                             np.ndarray wi,
                                             np.ndarray vl,
                                             int ldvl,
                                             np.ndarray vr,
                                             int ldvr,
                                             np.ndarray work,
                                             int nwork):
    cdef int info = 0
    cdgeev_(jobvl,
            jobvr,
            &n,
            <double *> A.data,
            &lda,
            <double *> wr.data,
            <double *> wi.data,
            <double *> vl.data,
            &ldvl,
            <double *> vr.data,
            &ldvr,
            <double *> work.data,
            &nwork,
            &info)             

def blasCopy(n,
             mat,
             denseFactor):
    DF_lu = denseFactor.lu

    lapackWrappersDCopy(n,
                        mat,
                        DF_lu)

def denseFactorPrepare(n,
                       mat,
                       denseFactor):
    DF_pivots = denseFactor.pivots
    DF_lu = denseFactor.lu
    
    lapackWrappersDCopy(n,
                        mat,
                        DF_lu)
    
    lapackWrappersDenseFactorPrepare(n,
                                     DF_lu,
                                     DF_pivots)

def denseFactorSolve(n,
                     mat,
                     denseFactor,
                     b):
    DF_pivots = denseFactor.pivots
    DF_lu = denseFactor.lu

    lapackWrappersDenseFactorSolve(n,
                                   DF_lu,
                                   DF_pivots,
                                   b)

def denseCalculateEigenvalues(jobvl,
                              jobvr,
                              n,
                              leig,
                              lda,
                              eigenvalues_r,
                              eigenvalues_i,
                              leftEigenvectors,
                              ldvl,
                              rightEigenvectors,
                              ldvr,
                              work,
                              lwork):
    """
    Parameters
    ----------
    jobvl (input) : str
        'N' indicates left eigenvectors are not computed
        'V' indicates left eigenvectors are computed
    jobvr (input) : str
        'N' indicates right eigenvectors are not computed
        'V' indicates right eigenvectors are computed
    n (input) : int
        order of matrix A
    leig (input / output) : np array
        on entry matrix A, on exit overwritten
    lda (input) : int
        leading dimension of the array A
    eigenvalues_r (output) : array
        real part of computed eigenvalues
    eigenvalues_i (output) : array
        imaginary part of computed eigenvalues
    vl (output) : array
        storage container for left eigenvectors
    ldvl (input) : int
        leading dimension of the array vl
    vr (output) : array
        storage container for right eigenvectos
    ldvr (input) : int
        leading dimension of array vr
    work (workspace / output) : array
    lwork (input) : int
    """

    lapackWrappersDenseCalculateEigenvalues(jobvl,
                                            jobvr,
                                            n,
                                            leig,
                                            lda,
                                            eigenvalues_r,
                                            eigenvalues_i,
                                            leftEigenvectors,
                                            ldvl,
                                            rightEigenvectors,
                                            ldvr,
                                            work,
                                            lwork)
    
