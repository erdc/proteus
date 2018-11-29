import numpy as np
cimport numpy as np
cdef extern from "proteus_lapack.h":
    extern int cdgetrf_ "dgetrf_"(int *m, int *n, double *a, int *lda, int *ipiv, int *info)
    extern int cdgetrs_ "dgetrs_"(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info)
    extern int cdgetc2_ "dgetc2_"(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info)
    extern int cdgesc2_ "dgesc2_"(int *n, double *a, int *lda, double* rhs, int *ipiv, int *jpiv, double* scale)
    extern int cdgeev_ "dgeev_"(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork,int* info)
    extern int cdgetri_ "dgetri_"(int* N,double* A,int* LDA,int* IPIV,double* WORK,int* LWORK,int* INFO )