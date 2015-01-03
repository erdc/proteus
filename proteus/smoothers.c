#include <math.h>
#include <stdio.h>
#include <string.h>
#include "smoothers.h"
#include PROTEUS_LAPACK_H
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
      int cnode = node_order[i];
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
      int cnode = node_order[i];
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
        int cnode = node_order[i];
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
                 PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p)
{ 
  int i,N; 
  int* subdomain_dim; 
  int** l2g_L;
  double** subdomain_L; 
  double** subdomain_R; 
  double** subdomain_dX;
  PROTEUS_LAPACK_INTEGER** subdomain_pivots;
  NRformat* ANR = (NRformat*)A->Store;
  N = A->nrow;
  *subdomain_dim_p = (int*)malloc(N*sizeof(int));
  *subdomain_pivots_p = (PROTEUS_LAPACK_INTEGER**)malloc(N*sizeof(PROTEUS_LAPACK_INTEGER*));
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
      subdomain_pivots[i] = (PROTEUS_LAPACK_INTEGER*)malloc(subdomain_dim[i]*sizeof(PROTEUS_LAPACK_INTEGER));
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
                 PROTEUS_LAPACK_INTEGER** subdomain_pivots)
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
                    PROTEUS_LAPACK_INTEGER** subdomainPivots) 
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
      PROTEUS_LAPACK_INTEGER La_N=((PROTEUS_LAPACK_INTEGER)subdomain_dim[i]),INFO=0; 
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
                  PROTEUS_LAPACK_INTEGER** subdomainPivots) 
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
    
    PROTEUS_LAPACK_INTEGER La_N,INFO=0,NRHS=1; 
    char TRANS='N'; 

    nnz = AStore->nnz; 
    nzval = (double*)AStore->nzval; 
    colind = AStore->colind; 
    rowptr = AStore->rowptr; 
    memset(dX,0,sizeof(double)*N); 
    for (i=0; i<N; i++)  
    { 
        int cnode = node_order[i]; 
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
        La_N = (PROTEUS_LAPACK_INTEGER) subdomain_dim[cnode];
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

int basm_NR_init(int rowBlocks,
		 SuperMatrix *A, 
                 int** subdomain_dim_p, 
                 int*** l2g_L_p,
                 double*** subdomain_L_p, 
                 double*** subdomain_R_p, 
                 double*** subdomain_dX_p,
                 PROTEUS_LAPACK_INTEGER*** subdomain_pivots_p,
		 PROTEUS_LAPACK_INTEGER*** subdomain_col_pivots_p)
{ 
  int i,N; 
  int* subdomain_dim; 
  int** l2g_L;
  double** subdomain_L; 
  double** subdomain_R; 
  double** subdomain_dX;
  PROTEUS_LAPACK_INTEGER** subdomain_pivots;
  PROTEUS_LAPACK_INTEGER** subdomain_col_pivots;
  NRformat* ANR = (NRformat*)A->Store;
  /*local variables for computing dimensions etc*/
  int max_local_dim=0;  /*size of largest possible subdomain system (assuming all the column entries are unique)*/
  int local_nonzero_entries = 0;    /*total number of nonzero entries in local row system*/
  int n_unique_local_dofs=0;/*number of unique unknowns in a local
			      system*/
  int *unique_local_dofs; /*use to keep track of which unknowns
			    have already been accounted for */
  assert((A->nrow % rowBlocks) == 0);
  N = A->nrow/rowBlocks;

  for (i=0; i < N; i++)
    {
      local_nonzero_entries = ANR->rowptr[(i+1)*rowBlocks]-ANR->rowptr[i*rowBlocks];
      if (local_nonzero_entries > max_local_dim)
	max_local_dim = local_nonzero_entries;
    }
  unique_local_dofs = (int*) malloc(max_local_dim*sizeof(int));
  for (i=0; i < max_local_dim; i++)
    unique_local_dofs[i]=-1;

  *subdomain_dim_p = (int*)malloc(N*sizeof(int));
  *subdomain_pivots_p = (PROTEUS_LAPACK_INTEGER**)malloc(N*sizeof(PROTEUS_LAPACK_INTEGER*));
  *subdomain_col_pivots_p = (PROTEUS_LAPACK_INTEGER**)malloc(N*sizeof(PROTEUS_LAPACK_INTEGER*));
  *l2g_L_p = (int**)malloc(N*sizeof(int*));
  *subdomain_R_p = (double**)malloc(N*sizeof(double*));
  *subdomain_dX_p = (double**)malloc(N*sizeof(double*));
  *subdomain_L_p = (double**)malloc(N*sizeof(double*));
  if ( (*subdomain_dim_p == NULL) ||
       (*l2g_L_p == NULL) ||
       (*subdomain_R_p == NULL) ||
       (*subdomain_dX_p == NULL) ||
       (*subdomain_L_p == NULL) ||
       (*subdomain_pivots_p == NULL) ||
       (*subdomain_col_pivots_p == NULL))
    {
      /*clean up*/
      free(unique_local_dofs);
      return 1;
    }
  subdomain_dim = *subdomain_dim_p;
  subdomain_pivots = *subdomain_pivots_p;
  subdomain_col_pivots = *subdomain_col_pivots_p;
  l2g_L = *l2g_L_p;
  subdomain_R = *subdomain_R_p;
  subdomain_dX =*subdomain_dX_p;
  subdomain_L =*subdomain_L_p;
   /* extract subdomain sizes and allocate storage */
  for (i=0; i<N; i++)  
    { 
      /*count number of unique dofs in the local system*/
      int jj,kk,j,found_j;
      local_nonzero_entries = ANR->rowptr[(i+1)*rowBlocks]-ANR->rowptr[i*rowBlocks];
      for (jj=0; jj < local_nonzero_entries; jj++)
	unique_local_dofs[jj] = -12345;
      n_unique_local_dofs = 0;
      for (jj=0; jj < local_nonzero_entries; jj++)
	{
	  /*global column id*/
	  j = ANR->colind[ANR->rowptr[i*rowBlocks]+jj];
	  /*is this entry unique*/
	  found_j = -1;
	  for (kk=0; kk < n_unique_local_dofs; kk++)
	    {
	      if (unique_local_dofs[kk] == j)
		{
		  found_j = 1;
		  break;
		}
	    }
	  if (found_j == -1)
	    {
	      unique_local_dofs[n_unique_local_dofs] = j;
	      n_unique_local_dofs++;
	    }
	  assert(n_unique_local_dofs <= local_nonzero_entries);
	}
      
      subdomain_dim[i] = n_unique_local_dofs;
      subdomain_pivots[i] = (PROTEUS_LAPACK_INTEGER*)malloc(subdomain_dim[i]*sizeof(PROTEUS_LAPACK_INTEGER));
      subdomain_col_pivots[i] = (PROTEUS_LAPACK_INTEGER*)malloc(subdomain_dim[i]*sizeof(PROTEUS_LAPACK_INTEGER));
      l2g_L[i] = (int*)malloc(subdomain_dim[i]*subdomain_dim[i]*sizeof(int)); 
      subdomain_R[i] = (double*)malloc(subdomain_dim[i]*sizeof(double)); 
      subdomain_dX[i] = (double*)malloc(subdomain_dim[i]*sizeof(double)); 
      subdomain_L[i] = (double*)malloc(subdomain_dim[i]*subdomain_dim[i]*sizeof(double)); 
      if ((l2g_L[i] == NULL) ||
          (subdomain_R[i] == NULL) ||
          (subdomain_dX[i] == NULL) ||
          (subdomain_L[i] == NULL))
        {
	  /*clean up*/
	  free(unique_local_dofs);
          return 1;
        }
    }
   /* extract the local to global mappings */ 
  for (i=0;i<N;i++) 
    { 
      int j, jj, k, kk, found_j;
      /*again determine unique column indeces*/
      local_nonzero_entries = ANR->rowptr[(i+1)*rowBlocks]-ANR->rowptr[i*rowBlocks];
      for (jj=0; jj < local_nonzero_entries; jj++)
	unique_local_dofs[jj] = -12345;
      n_unique_local_dofs = 0;
      for (jj=0; jj < local_nonzero_entries; jj++)
	{
	  /*global column id*/
	  j = ANR->colind[ANR->rowptr[i*rowBlocks]+jj];
	  /*is this entry unique*/
	  found_j = -1;
	  for (kk=0; kk < n_unique_local_dofs; kk++)
	    {
	      if (unique_local_dofs[kk] == j)
		{
		  found_j = 1;
		  break;
		}
	    }
	  if (found_j == -1)
	    {
	      unique_local_dofs[n_unique_local_dofs] = j;
	      n_unique_local_dofs++;
	    }
	  assert(n_unique_local_dofs <= local_nonzero_entries);
	}
      assert(subdomain_dim[i] == n_unique_local_dofs);
      for (kk=0; kk < subdomain_dim[i]; kk++) 
        { 
          k = unique_local_dofs[kk];
          for (jj=0;jj<subdomain_dim[i]; jj++) 
            { 
              l2g_L[i][kk*subdomain_dim[i]+jj] = -1; 
              for (j=ANR->rowptr[k];j<ANR->rowptr[k+1];j++)/*just look through k's row*/ 
                { 
                  if (ANR->colind[j] == unique_local_dofs[jj])
                  { 
                    l2g_L[i][kk*subdomain_dim[i]+jj] = j;
                    break; 
                  } 
                } 
            } 
        } 
    } 
  /*mwf debug*/
  printf("BASM_init N= %d rowBlocks=%d \n",N,rowBlocks);
  for (i=0; i < N; i++)
    {
      int kk,jj;
      printf("l2g_L[%d] dim= %d\n",i,subdomain_dim[i]);
      for (kk=0; kk < subdomain_dim[i]; kk++)
	{
	  printf("\t");
	  for (jj=0; jj < subdomain_dim[i]; jj++)
	    {
	      printf(" %d ",l2g_L[i][kk*subdomain_dim[i] + jj]);
	    }
	  printf("\n");
	}
    }
  /*clean up*/
  free(unique_local_dofs);

  return 0;
} 

void basm_NR_free(int N,
		  int* subdomain_dim, 
		  int** l2g_L,
		  double** subdomain_L, 
		  double** subdomain_R, 
		  double** subdomain_dX,
		  PROTEUS_LAPACK_INTEGER** subdomain_pivots,
		  PROTEUS_LAPACK_INTEGER** subdomain_col_pivots)
{
  int i;
  free(subdomain_dim);
  for (i=0;i<N;i++)
    {
      free(subdomain_pivots[i]);
      free(subdomain_col_pivots[i]);
      free(l2g_L[i]);
      free(subdomain_R[i]);
      free(subdomain_dX[i]);
      free(subdomain_L[i]);
    }
  free(subdomain_pivots);
  free(subdomain_col_pivots);
  free(l2g_L);
  free(subdomain_R);
  free(subdomain_dX);
  free(subdomain_L);
}

void basm_NR_prepare(int rowBlocks,
		     int N,
		     SuperMatrix *A, 
		     int* subdomain_dim,
		     int** l2g_L,
		     double** subdomainL, 
		     PROTEUS_LAPACK_INTEGER** subdomainPivots,
		     PROTEUS_LAPACK_INTEGER** subdomainColPivots) 
{ 
  /* Purpose: extract subdomain matrices and factor  in place 
   */ 
  int i; 
  NRformat *ANR = (NRformat*) A->Store;
  double *nzval = (double*) ANR->nzval;
  assert (N*rowBlocks == A->nrow);
  /*zero things for safety*/
  for (i=0; i< N; i++)
    {
      int ii,jj;
      for (ii=0; ii<subdomain_dim[i]; ii++)
	{
	  subdomainPivots[i][ii] = 0; subdomainColPivots[i][ii] = 0;
	  for (jj=0; jj<subdomain_dim[i]; jj++)
	    subdomainL[i][jj*subdomain_dim[i]+ii] = 0.0;
	}
    }
  for (i=0; i<N; i++)  
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
      PROTEUS_LAPACK_INTEGER La_N=((PROTEUS_LAPACK_INTEGER)subdomain_dim[i]),INFO=0; 
      dgetc2_(&La_N,subdomainL[i],&La_N,subdomainPivots[i],subdomainColPivots[i],&INFO);  
      /*doesn't seem to be handling trivial case of dim=1 and L_00 = 0 well*/
      /*mwf debug*/
      if (INFO > 0)
	{
	  printf("basm prepare jac dgetc2 INFO=%d nN=%d \n",(int)(INFO),subdomain_dim[i]);
	  for (ii=0;ii<subdomain_dim[i];ii++)
	    {
	      for(jj=0;jj<subdomain_dim[i];jj++)
		{
		  
		  printf("%12.5e \t",subdomainL[i][ii*subdomain_dim[i] + jj]);
		}
	      printf("\n");
	    }
	  for (ii=0; ii<subdomain_dim[i];ii++)
	    {
	      printf("colPivot[%d]= %dl \t",ii,subdomainColPivots[i][ii]);
	      printf("pivot[%d]= %dl \t",ii,subdomainPivots[i][ii]);
	    }
	}
    } 
} 

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
		   PROTEUS_LAPACK_INTEGER** subdomainColPivots) 
{ 
/* Purpose:  asm method with the additional assumption that A is of Stype NRformat 
 */ 
    int i; 
    int nnz; 
    double *nzval; 
    int *colind; 
    int *rowptr; 
    NRformat *AStore = (NRformat*) A->Store; 
    assert(N*rowBlocks == A->nrow);
    double scale = 1.0;
    PROTEUS_LAPACK_INTEGER La_N; 
    /*have to keep track of unique dofs in local system now*/
    int max_local_dim=0;  /*size of largest possible subdomain system (assuming all the column entries are unique)*/
    int local_nonzero_entries = 0;    /*total number of nonzero entries in local row system*/
    int n_unique_local_dofs=0;/*number of unique unknowns in a local
				system*/
    int *unique_local_dofs; /*use to keep track of which unknowns
			      have already been accounted for */
    

    for (i=0; i < N; i++)
      {
	local_nonzero_entries = AStore->rowptr[(i+1)*rowBlocks]-AStore->rowptr[i*rowBlocks];
	if (local_nonzero_entries > max_local_dim)
	  max_local_dim = local_nonzero_entries;
      }
    unique_local_dofs = (int*) malloc(max_local_dim*sizeof(int));
    for (i=0; i < max_local_dim; i++)
      unique_local_dofs[i]=-1;
    
    nnz = AStore->nnz; 
    nzval = (double*)AStore->nzval; 
    colind = AStore->colind; 
    rowptr = AStore->rowptr; 
    memset(dX,0,sizeof(double)*N*rowBlocks); 
    for (i=0; i<N; i++)  
    { 
        int cnode = node_order[i]; 
        int j, k,kk, jj, ii, found_j=-1; 
	/*again determine unique column indeces maybe go ahead and store these*/
	local_nonzero_entries = AStore->rowptr[(cnode+1)*rowBlocks]-AStore->rowptr[cnode*rowBlocks];
	for (jj=0; jj < local_nonzero_entries; jj++)
	  unique_local_dofs[jj] = -12345;
	n_unique_local_dofs = 0;
	for (jj=0; jj < local_nonzero_entries; jj++)
	  {
	    /*global column id*/
	    j = AStore->colind[AStore->rowptr[cnode*rowBlocks]+jj];
	    /*is this entry unique*/
	    found_j = -1;
	    for (kk=0; kk < n_unique_local_dofs; kk++)
	      {
		if (unique_local_dofs[kk] == j)
		  {
		    found_j = 1;
		    break;
		  }
	      }
	    if (found_j == -1)
	      {
		unique_local_dofs[n_unique_local_dofs] = j;
		n_unique_local_dofs++;
	      }
	    assert(n_unique_local_dofs <= local_nonzero_entries);
	  }
	assert(subdomain_dim[cnode] == n_unique_local_dofs);
        /* extract and update the subdomain residual */ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          { 
            subdomainR[cnode][jj] = R[unique_local_dofs[jj]]; /*colind[rowptr[cnode*rowBlocks] + jj]];*/ 
          } 
        for (ii=0;ii<subdomain_dim[cnode];ii++)
          {
            k = unique_local_dofs[ii]; /*colind[rowptr[cnode*rowBlocks]+ii];*/
            for (j=rowptr[k];j<rowptr[k+1];j++)
              subdomainR[cnode][ii] -= nzval[j]*dX[colind[j]];/*check this if needs to be only entries in nodestar*/
          }
        /* copy R into dX because lapack wants it that way*/ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          { 
            subdomain_dX[cnode][jj] = subdomainR[cnode][jj]; 
          }         
        /* solve  subdomain problem*/
        La_N = (PROTEUS_LAPACK_INTEGER) subdomain_dim[cnode];
	/*doesn't seem to be handling trivial case of dim=1 and L_00 = 0 well*/
	if (La_N == 1 && fabs(subdomainL[cnode][0]) <= 1.0e-64)
	  subdomain_dX[cnode][0] = 0.0;
	else
	  {
	    dgesc2_(&La_N, 
		    subdomainL[cnode], 
		    &La_N, 
		    subdomain_dX[cnode], 
		    subdomainPivots[cnode], 
		    subdomainColPivots[cnode],
		    &scale);
	  }
        /* set the global correction from the subdomain correction */ 
        for (jj = 0;jj<subdomain_dim[cnode];jj++) 
          {
	    /*dX[colind[rowptr[cnode*rowBlocks] + jj]] += w*subdomain_dX[cnode][jj];*/
	    dX[unique_local_dofs[jj]] += w*subdomain_dX[cnode][jj];
	  }
    } 
} 

/** @} */
