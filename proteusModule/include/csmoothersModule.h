#ifndef SMOOTHERSMODULE_H
#define SMOOTHERSMODULE_H

#include "Python.h"
#include PROTEUS_LAPACK_H

/**
 \file csmoothersModule.h
 \brief The python interface to smoothers
*/

/**
  \defgroup csmoothers csmoothers
  \brief The python interface to the smoothers library
   @{
*/  
typedef struct 
{
  PyObject_HEAD
  int N;
  int *subdomain_dim;
  int **l2g_L;
  double **subdomain_L,
    **subdomain_R,
    **subdomain_dX;
  PROTEUS_LAPACK_INTEGER** subdomain_pivots;
} ASMFactor;
typedef struct 
{
  PyObject_HEAD
  int N; int bs;
  int *subdomain_dim;
  int **l2g_L;
  double **subdomain_L,
    **subdomain_R,
    **subdomain_dX;
  PROTEUS_LAPACK_INTEGER** subdomain_pivots;
  PROTEUS_LAPACK_INTEGER** subdomain_col_pivots;
} BASMFactor;
/** @} */
#endif
