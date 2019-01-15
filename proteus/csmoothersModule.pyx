import numpy as np
import cython
cimport numpy as np

cdef extern from "smoothers.h":
    void cjacobi_NR_prepare "jacobi_NR_prepare"(_SuperMatrix *A, double w, double tol, double* M)
    void cjacobi_NR_solve "jacobi_NR_solve"(_SuperMatrix *A, double* M, double* R, int* node_order, double* dX)
    void cnl_jacobi_NR_solve "nl_jacobi_NR_solve"(_SuperMatrix *A, double* R, int* node_order, double w, double tol, double* dX)
    void cgauss_seidel_NR_prepare "gauss_seidel_NR_prepare"(_SuperMatrix *A, double w, double tol, double* M)
    void cgauss_seidel_NR_solve "gauss_seidel_NR_solve"(_SuperMatrix *A, double *M, double *R, int *node_order, double *dX)
    void cnl_gauss_seidel_NR_solve "nl_gauss_seidel_NR_solve"(_SuperMatrix *A, double *R, int *node_order, double w, double tol, double *dX)
    void casm_NR_init "asm_NR_init"(_SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p)
    void casm_NR_free "asm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots)
    void casm_NR_prepare "asm_NR_prepare"(_SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    void casm_NR_solve "asm_NR_solve"(_SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots)
    int cbasm_NR_init "basm_NR_init"(int rowBlocks, _SuperMatrix *A, int** subdomain_dim_p, int*** l2g_L_p, double*** subdomain_L_p, double*** subdomain_R_p, double*** subdomain_dX_p, PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p, PROTEUS_LAPACK_INTEGER*** subdomain_col_pivots_p)
    void cbasm_NR_free "basm_NR_free"(int N, int* subdomain_dim, int** l2g_L, double** subdomain_L, double** subdomain_R, double** subdomain_dX, PROTEUS_LAPACK_INTEGER** subdomain_pivots, PROTEUS_LAPACK_INTEGER** subdomain_col_pivots)
    void cbasm_NR_prepare "basm_NR_prepare"(int rowBlocks, int N, _SuperMatrix *A, int* subdomain_dim, int** l2g_L, double** subdomainL, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)
    void cbasm_NR_solve "basm_NR_solve"(int rowBlocks, int N, _SuperMatrix *A, double w, double** subdomainL, int* subdomain_dim, int** l2g_L, double* R, double** subdomainR, int* node_order, double** subdomain_dX, double* dX, PROTEUS_LAPACK_INTEGER** subdomainPivots, PROTEUS_LAPACK_INTEGER** subdomainColPivots)
