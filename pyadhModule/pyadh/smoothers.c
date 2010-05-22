#include <math.h>
#include <stdio.h>
#include <string.h>
#include "smoothers.h"
#include PYADH_LAPACK_H
/** \file smoothers.c
    \ingroup  smoothers
    @{
*/
void jacobi_NR_prepare(SuperMatrix *A, double w, double tol, double* M)
{
    gauss_seidel_NR_prepare(A, w, tol, M);
}

void jacobi_NR_solve(SuperMatrix *A, double* M, double* R, int *node_order, double* dX)
{
  /* Purpose:  jacobi method with the additional assumption that A is of 
   * Stype NRformat
   */
  int i;
    int N = A->nrow;
  for (i=0; i<N; i++) 
    dX[i] = M[i]*R[i];   
}

void nl_jacobi_NR_solve(SuperMatrix *A,  double* R, int *node_order, double w, double tol, double* dX)
{
/* Purpose:  jacobi method with the additional assumption that A is of Stype NRformat
 */

    int i;
    int nnz;
    double *nzval;
    int *colind;
    int *rowptr;
    NRformat *AStore = (NRformat*) A->Store;
    int N = A->nrow;

    nnz = AStore->nnz;
    nzval = (double*)AStore->nzval;
    colind = AStore->colind;
    rowptr = AStore->rowptr;
    memset(dX,0,sizeof(double)*N);
    for (i=0; i<N; i++) 
    {
      int cnode = (int)node_order[i];
      int j;
      double diag = 1.0;
      int diag_found = 0;
      double val = R[cnode];
      for (j=rowptr[cnode]; j<rowptr[cnode+1]; j++)
        {
          /*check for diagonal element*/
          if (colind[j] == cnode)
            {
              diag = nzval[j];
              if (fabs(nzval[j]) >= tol)
                diag_found = 1;
            }
        }
      if (diag_found) 
        dX[cnode] = w*val/diag;
      else 
        {
          printf("D[%d] = %12.5e \n",cnode,diag);
          ABORT("Diagonal element is 0 within tol in Jacobi or is not in sparse matrix.");
        }
    }
}

void gauss_seidel_NR_prepare(SuperMatrix *A, double w, double tol, double* M)
{
  /* Purpose:  invert diagonal for gauss_seidel method with the additional assumption that A is of Stype NRformat
   */

  int i;
  int nnz;
  double *nzval;
  int *colind;
  int *rowptr;
  NRformat *AStore = (NRformat*) A->Store;
  int N = A->nrow;

  nnz = AStore->nnz;
  nzval = (double*)AStore->nzval;
  colind = AStore->colind;
  rowptr = AStore->rowptr;

  for (i=0; i<N; i++) 
    {
      int j;
      double diag=1.0;
      int diag_found = 0;
      for (j=rowptr[i]; j<rowptr[i+1]; j++)
        {
          /*check for diagonal element*/
          if (colind[j] == i)
            {
              diag = nzval[j];
              if (fabs(nzval[j]) >= tol)
                diag_found = 1;
            }
        }
      if (diag_found)
        M[i] = w/diag;
      else
        {
          printf("D[%d] = %12.5e \n",i,diag);
          ABORT("Diagonal element is 0 within tol in Gauss-Seidel or is not in sparse matrix.");
        }
    }
}
  
void gauss_seidel_NR_solve(SuperMatrix *A, double* M, double* R, int *node_order, double* dX)
{
/* Purpose:  gauss_seidel method with the additional assumption that A is of Stype NRformat
 */
    int i;
    int nnz;
    double *nzval;
    int *colind;
    int *rowptr;
    NRformat *AStore = (NRformat*) A->Store;
    int N = A->nrow;

    nnz = AStore->nnz;
    nzval = (double*)AStore->nzval;
    colind = AStore->colind;
    rowptr = AStore->rowptr;
    memset(dX,0,sizeof(double)*N);
    for (i=0; i<N; i++) 
    {
        int cnode = (int)node_order[i];
        int j;
        double val = R[cnode];
        for (j=rowptr[cnode]; j<rowptr[cnode+1]; j++)
          val -= nzval[j]*dX[colind[j]];
        dX[cnode] = M[cnode]*val;
    }
}

void nl_gauss_seidel_NR_solve(SuperMatrix *A,  double* R, int *node_order, double w, double tol, double* dX)
{
/* Purpose:  gauss_seidel method with the additional assumption that A is of Stype NRformat
 */

    int i;
    int nnz;
    double *nzval;
    int *colind;
    int *rowptr;
    NRformat *AStore = (NRformat*) A->Store;
    int N = A->nrow;

    nnz = AStore->nnz;
    nzval = (double*)AStore->nzval;
    colind = AStore->colind;
    rowptr = AStore->rowptr;
    memset(dX,0,sizeof(double)*N);
    for (i=0; i<N; i++) 
    {
        int cnode = (int)node_order[i];
        int j;
        double diag=1.0;
        int diag_found = 0;
        double val = R[cnode];

        for (j=rowptr[cnode]; j<rowptr[cnode+1]; j++)
        {
          /*check for diagonal element*/
          if (colind[j] == cnode)
            {
              diag = nzval[j];
              if (fabs(nzval[j]) >= tol)
                diag_found = 1;
            }
          val -= nzval[j]*dX[colind[j]];
        }
        if (diag_found) 
          dX[cnode] = w*val/diag;
        else 
          {
            printf("D[%d] = %12.5e \n",cnode,diag);
            ABORT("Diagonal element is 0 within tol in Gauss-Seidel or is not in sparse matrix.");
          }
    }
}

int asm_NR_init(SuperMatrix *A, 
                 int** subdomain_dim_p, 
                 int*** l2g_L_p,
                 double*** subdomain_L_p, 
                 double*** subdomain_R_p, 
                 double*** subdomain_dX_p,
                 PYADH_LAPACK_INTEGER*** subdomain_pivots_p)
{ 
  int i,N; 
  int* subdomain_dim; 
  int** l2g_L;
  double** subdomain_L; 
  double** subdomain_R; 
  double** subdomain_dX;
  PYADH_LAPACK_INTEGER** subdomain_pivots;
  NRformat* ANR = (NRformat*)A->Store;
  N = A->nrow;
  *subdomain_dim_p = (int*)malloc(N*sizeof(int));
  *subdomain_pivots_p = (PYADH_LAPACK_INTEGER**)malloc(N*sizeof(PYADH_LAPACK_INTEGER*));
  *l2g_L_p = (int**)malloc(N*sizeof(int*));
  *subdomain_R_p = (double**)malloc(N*sizeof(double*));
  *subdomain_dX_p = (double**)malloc(N*sizeof(double*));
  *subdomain_L_p = (double**)malloc(N*sizeof(double*));
  if ( (*subdomain_dim_p == NULL) ||
       (*l2g_L_p == NULL) ||
       (*subdomain_R_p == NULL) ||
       (*subdomain_dX_p == NULL) ||
       (*subdomain_L_p == NULL) ||
       (*subdomain_pivots_p == NULL) )
    {
      return 1;
    }
  subdomain_dim = *subdomain_dim_p;
  subdomain_pivots = *subdomain_pivots_p;
  l2g_L = *l2g_L_p;
  subdomain_R = *subdomain_R_p;
  subdomain_dX =*subdomain_dX_p;
  subdomain_L =*subdomain_L_p;
   /* extract subdomain sizes and allocate storage */
  for (i=0; i<N; i++)  
    { 
      subdomain_dim[i] = ANR->rowptr[i+1]-ANR->rowptr[i];
      subdomain_pivots[i] = (PYADH_LAPACK_INTEGER*)malloc(subdomain_dim[i]*sizeof(PYADH_LAPACK_INTEGER));
      l2g_L[i] = (int*)malloc(subdomain_dim[i]*subdomain_dim[i]*sizeof(int)); 
      subdomain_R[i] = (double*)malloc(subdomain_dim[i]*sizeof(double)); 
      subdomain_dX[i] = (double*)malloc(subdomain_dim[i]*sizeof(double)); 
      subdomain_L[i] = (double*)malloc(subdomain_dim[i]*subdomain_dim[i]*sizeof(double)); 
      if ((l2g_L[i] == NULL) ||
          (subdomain_R[i] == NULL) ||
          (subdomain_dX[i] == NULL) ||
          (subdomain_L[i] == NULL))
        {
          return 1;
        }
    }
   /* extract the local to global mappings */ 
  for (i=0;i<N;i++) 
    { 
      int j, jj, k, kk;
      for (kk=0; kk < subdomain_dim[i]; kk++) 
        { 
          k = ANR->colind[ANR->rowptr[i] + kk]; 
          for (jj=0;jj<subdomain_dim[i]; jj++) 
            { 
              l2g_L[i][kk*subdomain_dim[i]+jj] = -1; 
              for (j=ANR->rowptr[k];j<ANR->rowptr[k+1];j++) 
                { 
                  if (ANR->colind[j] == ANR->colind[ANR->rowptr[i] + jj]) 
                  { 
                    l2g_L[i][kk*subdomain_dim[i]+jj] = j;
                    break; 
                  } 
                } 
            } 
        } 
    } 
  return 0;
} 
  
void asm_NR_free(int N, 
                 int* subdomain_dim, 
                 int** l2g_L,
                 double** subdomain_L, 
                 double** subdomain_R, 
                 double** subdomain_dX,
                 PYADH_LAPACK_INTEGER** subdomain_pivots)
{
  int i;
  free(subdomain_dim);
  for (i=0;i<N;i++)
    {
      free(subdomain_pivots[i]);
      free(l2g_L[i]);
      free(subdomain_R[i]);
      free(subdomain_dX[i]);
      free(subdomain_L[i]);
    }
  free(subdomain_pivots);
  free(l2g_L);
  free(subdomain_R);
  free(subdomain_dX);
  free(subdomain_L);
}

void asm_NR_prepare(SuperMatrix *A, 
                    int* subdomain_dim,
                    int** l2g_L,
                    double** subdomainL, 
                    PYADH_LAPACK_INTEGER** subdomainPivots) 
{ 
  /* Purpose: extract subdomain matrices and factor  in place 
   */ 
  int i; 
  NRformat *ANR = (NRformat*) A->Store;
  double *nzval = (double*) ANR->nzval;
  for (i=0; i<A->nrow; i++)  
    { 
      int ii,jj; 
      for (ii=0;ii<subdomain_dim[i];ii++) 
        for (jj=0;jj<subdomain_dim[i];jj++)  
          { 
            int index = ii*subdomain_dim[i]  + jj; 
            if (l2g_L[i][index] != -1) 
              { 
                subdomainL[i][ii + jj*subdomain_dim[i]] =  
                  nzval[l2g_L[i][index]]; 
              } 
            else 
              subdomainL[i][ii+jj*subdomain_dim[i]] = 0.0;
          } 
      PYADH_LAPACK_INTEGER La_N=((PYADH_LAPACK_INTEGER)subdomain_dim[i]),INFO=0; 
      dgetrf_(&La_N,&La_N,subdomainL[i],&La_N,subdomainPivots[i],&INFO);  
    } 
} 

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
                  PYADH_LAPACK_INTEGER** subdomainPivots) 
{ 
/* Purpose:  asm method with the additional assumption that A is of Stype NRformat 
 */ 
    int i; 
    int nnz; 
    double *nzval; 
    int *colind; 
    int *rowptr; 
    NRformat *AStore = (NRformat*) A->Store; 
    int N = A->nrow; 
    
    PYADH_LAPACK_INTEGER La_N,INFO=0,NRHS=1; 
    char TRANS='N'; 

    nnz = AStore->nnz; 
    nzval = (double*)AStore->nzval; 
    colind = AStore->colind; 
    rowptr = AStore->rowptr; 
    memset(dX,0,sizeof(double)*N); 
    for (i=0; i<N; i++)  
    { 
        int cnode = (int)node_order[i]; 
        int j, k,jj, ii; 
        /* extract and update the subdomain residual */ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          { 
            subdomainR[cnode][jj] = R[colind[rowptr[cnode] + jj]]; 
          } 
        for (ii=0;ii<subdomain_dim[cnode];ii++)
          {
            k = colind[rowptr[cnode]+ii];
            for (j=rowptr[k];j<rowptr[k+1];j++)
              subdomainR[cnode][ii] -= nzval[j]*dX[colind[j]];
          }
        /* copy R into dX because lapack wants it that way*/ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          { 
            subdomain_dX[cnode][jj] = subdomainR[cnode][jj]; 
          }         
        /* solve  subdomain problem*/
        La_N = (PYADH_LAPACK_INTEGER) subdomain_dim[cnode];
        dgetrs_(&TRANS, 
                &La_N, 
                &NRHS, 
                subdomainL[cnode], 
                &La_N, 
                subdomainPivots[cnode], 
                subdomain_dX[cnode], 
                &La_N, 
                &INFO); 
        /* set the global correction from the subdomain correction */ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          dX[colind[rowptr[cnode] + jj]] += w*subdomain_dX[cnode][jj];
    } 
} 
/** @} */
