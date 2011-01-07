#ifndef PYADH_LAPACK_PROTO_H
#define PYADH_LAPACK_PROTO_H

#ifdef __cplusplus
extern "C"
{
#endif

extern int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern int dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern int dgetc2_(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info);
extern int dgesc2_(int *n, double *a, int *lda, double* rhs, int *ipiv, int *jpiv, double* scale);
extern int dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork,int* info);
extern int dgetri_(int* N,double* A,int* LDA,int* IPIV,double* WORK,int* LWORK,int* INFO );

#ifdef __cplusplus
}
#endif

#define __CLPK_integer int
#endif
