#ifndef FLCBDFWRAPPERSMODULE_H
#define FLCBDFWRAPPERSMODULE_H 
#include <cstddef>
#include <complex>

extern "C"
{
#include "Python.h"
#include "numpy/arrayobject.h"
#include "superluWrappersModule.h"
}

#include "cmeshToolsModule.h"
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "petscsys.h"
#include "petsc.h"
#include "petscmat.h"
#include "petscao.h"
#include "petscbt.h"
#include "petscksp.h"
#include "petscconf.h"
#include <stdio.h>
#include <valarray>
#include "partitioning.h"

/**
   \defgroup flcbdfWrappersOld flcbdfWrappersOld
   \brief Python interface to adaptive BDF code (daetk)
   @{ 
*/

extern "C"
{
#ifdef FLCBDF_WRAPPERS_MODULE
  static MPI_Comm PROTEUS_COMM_WORLD;
#else
  static void **PyFLCBDFWrappersOld_API;
  #define PROTEUS_COMM_WORLD \
    *(MPI_Comm*)(PyFLCBDFWrappersOld_API[0])
  static int
  import_flcbdfWrappersOld(void)
  {
    PyObject* module = PyImport_ImportModule("proteus.flcbdfWrappersOld");
    
    if (module != NULL)
      {
	PyObject *c_api_object = PyObject_GetAttrString(module,"_C_API");
	if (c_api_object == NULL)
	  return -1;
	if (PyCObject_Check(c_api_object))
	  PyFLCBDFWrappersOld_API = (void **)PyCObject_AsVoidPtr(c_api_object);
	Py_DECREF(c_api_object);
      }
    return 0;
  }
#endif
}
/** @} */
#endif
