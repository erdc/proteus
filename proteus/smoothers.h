#ifndef SMOOTHERS_H
#define SMOOTHERS_H
/*!
 \file smoothers.h
 \brief C implementations of multilevel smoother algorithms
*/
/*!
 \defgroup smoothers smoothers
 \brief C implementations of multilevel smoother algorithms
 @{
*/

#include "slu_ddefs.h"
#include PROTEUS_LAPACK_H

void jacobi_NR_prepare(SuperMatrix *A, 
                       double w, 
                       double tol, 
                       double* M);
void jacobi_NR_solve(SuperMatrix *A, 
                     double* M, 
                     double* R, 
                     int *node_order, 
                     double* dX);
void nl_jacobi_NR_solve(SuperMatrix *A,  
                        double* R, 
                        int *node_order, 
                        double w, 
                        double tol, 
                        double* dX);
void gauss_seidel_NR_prepare(SuperMatrix *A, 
                             double w, 
                             double tol, 
                             double* M);
void gauss_seidel_NR_solve(SuperMatrix *A, 
                           double* M, 
                           double* R, 
                           int *node_order, 
                           double* dX);
void nl_gauss_seidel_NR_solve(SuperMatrix *A,  
                              double* R, 
                              int *node_order, 
                              double w, 
                              double tol, 
                              double* dX);
int asm_NR_init(SuperMatrix *A, 
                 int** subdomain_dim_p, 
                 int*** l2g_L_p,
                 double*** subdomain_L_p, 
                 double*** subdomain_R_p, 
                 double*** subdomain_dX_p,
                PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p);
void asm_NR_free(int N, 
                 int* subdomain_dim, 
                 int** l2g_L,
                 double** subdomain_L, 
                 double** subdomain_R, 
                 double** subdomain_dX,
                 PROTEUS_LAPACK_INTEGER** subdomain_pivots);
void asm_NR_prepare(SuperMatrix *A, 
                    int* subdomain_dim,
                    int** l2g_L,
                    double** subdomainL, 
                    PROTEUS_LAPACK_INTEGER** subdomainPivots);
void asm_NR_solve(SuperMatrix *A, 
                  double w,
                  double** subdomainL, 
                  int* subdomain_dim, 
                  int** l2g_L,  
                  double* R, 
                  double** subdomainR,
                  int *node_order, 
                  double** subdomain_dX,
                  double* dX, 
                  PROTEUS_LAPACK_INTEGER** subdomainPivots); 
/*repeat for blocks and use full pivoting*/
int basm_NR_init(int rowBlocks,
		 SuperMatrix *A, 
                 int** subdomain_dim_p, 
                 int*** l2g_L_p,
                 double*** subdomain_L_p, 
                 double*** subdomain_R_p, 
                 double*** subdomain_dX_p,
                 PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p,
		 PROTEUS_LAPACK_INTEGER*** subdomain_col_pivots_p);

void basm_NR_free(int N,
		  int* subdomain_dim, 
		  int** l2g_L,
		  double** subdomain_L, 
		  double** subdomain_R, 
		  double** subdomain_dX,
		  PROTEUS_LAPACK_INTEGER** subdomain_pivots,
		  PROTEUS_LAPACK_INTEGER** subdomain_col_pivots);


void basm_NR_prepare(int rowBlocks,
		     int N,
		     SuperMatrix *A, 
		     int* subdomain_dim,
		     int** l2g_L,
		     double** subdomainL, 
		     PROTEUS_LAPACK_INTEGER** subdomainPivots,
		     PROTEUS_LAPACK_INTEGER** subdomainColPivots);

void basm_NR_solve(int rowBlocks,
		   int N,
		   SuperMatrix *A, 
		   double w,
		   double** subdomainL, 
		   int* subdomain_dim, 
		   int** l2g_L,  
		   double* R, 
		   double** subdomainR,
		   int *node_order, 
		   double** subdomain_dX,
		   double* dX, 
		   PROTEUS_LAPACK_INTEGER** subdomainPivots,
		   PROTEUS_LAPACK_INTEGER** subdomainColPivots);
/** @} */
#endif
