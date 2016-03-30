#ifndef FLCBDFWRAPPERSMODULE_H
#define FLCBDFWRAPPERSMODULE_H 
#include <cstddef>

extern "C"
{
#include "Python.h"
#include "numpy/arrayobject.h"
#include "superluWrappersModule.h"
}

#include "cmeshToolsModule.h"
#include "mesh.h"

namespace Daetk 
{
  namespace Petsc
  {
    namespace cc
    {
      extern "C"
      {
/* hack on diamond */
#define MYCPLUSPLUS __cplusplus
#ifdef __cplusplus
#undef __cplusplus
#define PROTEUSREDEFINECPP
#endif
#include "mpi.h"
#ifdef PETSC_INCLUDE_AS_C
#define PETSC_RESTRICT __restrict__
#define PETSC_C_RESTRICT __restrict__
#define PETSC_CXX_RESTRICT __restrict__
#include "petscsys.h"
#define PETSC_RESTRICT __restrict__
#define PETSC_C_RESTRICT __restrict__
#define PETSC_CXX_RESTRICT __restrict__
#include "petsc.h"
#include "petscmat.h"
#include "petscao.h"
#include "petscbt.h"
#include "petscksp.h"
#include "petscconf.h"
#endif
#ifdef PROTEUSREDEFINECPP
#define __cplusplus MYCPLUSPLUS
#endif
#ifndef PETSC_INCLUDE_AS_C
#include "petsc.h"
#include "petscmat.h"
#include "petscao.h"
#include "petscbt.h"
#include "petscksp.h"
#include "petscconf.h"
#endif
	/*cek try adding forward declarations since I can't find the header for these intel functions */
#ifdef PETSC_HAVE__INTEL_FAST_MEMSET
          #include <string.h>
          extern void* _intel_fast_memset(void *b, int c, size_t len);
#endif
#ifdef PETSC_HAVE__INTEL_FAST_MEMCPY
          #include <string.h>
          extern void* _intel_fast_memcpy(void* s1, const void* s2, size_t n);
#endif
/*mwf debug*/
/*#include "parmetis.h" cek hack on lonestar*/
      }
    }
  }
}
#include "Definitions.h"
#include "FullDataFile.h"
#include "WeightedRMSNorm.h"
#include "FLCBDF_lite.h"
#include <stdio.h>
#include <valarray>
/*!
 \file flcbdfWrappersModule.h
 \brief Python interface to adpative BDF code (daetk)
*/

/**
   \defgroup flcbdfWrappers flcbdfWrappers
   \brief Python interface to adaptive BDF code (daetk)
   @{ 
*/

extern "C"
{
#ifdef FLCBDF_WRAPPERS_MODULE
  static Daetk::Petsc::cc::MPI_Comm PROTEUS_COMM_WORLD;
#else
  static void **PyFLCBDFWrappers_API;
  #define PROTEUS_COMM_WORLD \
    *(Daetk::Petsc::cc::MPI_Comm*)(PyFLCBDFWrappers_API[0])
  static int
  import_flcbdfWrappers(void)
  {
    PyObject* module = PyImport_ImportModule("proteus.flcbdfWrappers");
    
    if (module != NULL)
      {
	PyObject *c_api_object = PyObject_GetAttrString(module,"_C_API");
	if (c_api_object == NULL)
	  return -1;
	if (PyCObject_Check(c_api_object))
	  PyFLCBDFWrappers_API = (void **)PyCObject_AsVoidPtr(c_api_object);
	Py_DECREF(c_api_object);
      }
    return 0;
  }
#endif

typedef struct
{
  PyObject_HEAD
  Daetk::Petsc::Vec* sizeVec,*yVec,*DyVec,*yprimeVec,*DyprimeVec;
  Daetk::Petsc::Sys* petscSys;
  Daetk::WeightedL2Norm* wNorm;
  Daetk::FullDataFile* data;
  //mwf could make this base class to allow dummy data collection on some processors
  //Daetk::DataCollector * data;
  Daetk::FLCBDF_lite* flcbdf;
} FLCBDF_integrator;

typedef struct
{
  PyObject_HEAD
  Daetk::Petsc::Sys* petscSys;
} DaetkPetscSys;

}
/** @} */
#endif
