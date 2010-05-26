#ifndef SUPERLUWRAPPERSMODULE_H
#define SUPERLUWRAPPERSMODULE_H 
#include PROTEUS_SUPERLU_H

/*!
 \file superluWrappersModule.h
 \brief Python interface to superlu
*/
/*!
 \defgroup superluWrappers superluWrappers
 \brief Python interface to superlu
 @{
*/

typedef struct
{
  PyObject_HEAD
  int dim[2];
  NRformat A;
} SparseMatrix;

/** @} */
#endif
