#include "Python.h"
#include "numpy/arrayobject.h"
#include PROTEUS_LAPACK_H
#include PROTEUS_BLAS_H
#include "proteus_blas.h"
/** \file lapackWrappersModule.c
    \defgroup lapackWrappers lapackWrappers
    \brief Python interface to lapack
    @{
*/
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))

typedef struct 
{
  PyObject_HEAD
  double* lu;
  PROTEUS_LAPACK_INTEGER* pivots;
} DenseFactor;

#define DFP(p) ((DenseFactor*)p)

static PyObject*
DenseFactor_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  DenseFactor *self;
  self = (DenseFactor *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
DenseFactor_init(DenseFactor *self, PyObject *args, PyObject *kwds)
{
  int dim;
  if(!PyArg_ParseTuple(args,
                       "i",
                       &dim))
    return -1;
  self->pivots=(PROTEUS_LAPACK_INTEGER*)malloc(dim*sizeof(PROTEUS_LAPACK_INTEGER));
  self->lu=(double*)malloc(dim*dim*sizeof(double));
  if ( (self->pivots == NULL) || (self->lu == NULL))
    {
      PyErr_NoMemory();
      return -1;
    }
  return 0;
}

static  void
DenseFactor_dealloc(DenseFactor* self)
{
  free(self->lu);
  free(self->pivots);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject DenseFactorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "lapackWrappers.DenseFactor",             /*tp_name*/
  sizeof(DenseFactor), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)DenseFactor_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "dense factor objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)DenseFactor_init,      /* tp_init */
  0,                         /* tp_alloc */
  DenseFactor_new,                 /* tp_new */
};


static PyObject* lapackWrappersDenseFactorPrepare(PyObject* self,
                                                  PyObject* args)
{
  int dim,dim2;
  int incr1,incr2;
  PROTEUS_LAPACK_INTEGER N,INFO=0;
  PyObject *mat,*denseFactor;
  if(!PyArg_ParseTuple(args,"iOO",
                       &dim,
                       &mat,
                       &denseFactor))
    return NULL;
  N = (PROTEUS_LAPACK_INTEGER)(dim);
  dim2=dim*dim; incr1=1; incr2=1;
  dcopy_(&dim2,DDATA(mat),&incr1,DFP(denseFactor)->lu,&incr2);
  dgetrf_(&N,
          &N,
          DFP(denseFactor)->lu,
          &N,
          DFP(denseFactor)->pivots,
          &INFO);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* lapackWrappersDenseFactorSolve(PyObject* self,
                                                PyObject* args)
{
  int dim,use_transpose=0;
  PROTEUS_LAPACK_INTEGER N,INFO=0,NRHS=1;
  char TRANS='N';
  PyObject *mat,*b,*denseFactor;
  if(!PyArg_ParseTuple(args,"iOOO|i",
                       &dim,
                       &mat,
                       &denseFactor,
                       &b,
		       &use_transpose))
    return NULL;
  if (use_transpose == 1)
    TRANS='T';
  else
    TRANS='N';

  N = (PROTEUS_LAPACK_INTEGER)(dim);
  dgetrs_(&TRANS,
          &N,
          &NRHS,
          DFP(denseFactor)->lu,
          &N,
          DFP(denseFactor)->pivots,
          DDATA(b),
          &N,
          &INFO);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* lapackWrappersDenseCalculateEigenvalues(PyObject* self,
                                                         PyObject* args)
{
  char JOBVL,JOBVR;
  int N,LDA,LDVL,LDVR,LWORK,INFO=0;
  PyObject *A,*VL,*VR,*WI,*WORK,*WR;
  if (!PyArg_ParseTuple(args,
                        "cciOiOOOiOiOi",
                        &JOBVL,
                        &JOBVR,
                        &N,
                        &A,
                        &LDA,
                        &WR,
                        &WI,
                        &VL,
                        &LDVL,
                        &VR,
                        &LDVR,
                        &WORK,
                        &LWORK))
    return NULL;
  PROTEUS_LAPACK_INTEGER laN=((PROTEUS_LAPACK_INTEGER)N),
    laLDA=((PROTEUS_LAPACK_INTEGER)LDA),
    laLDVL=((PROTEUS_LAPACK_INTEGER)LDVL),
    laLDVR=((PROTEUS_LAPACK_INTEGER)LDVR),
    laLWORK=((PROTEUS_LAPACK_INTEGER)LWORK),
    laINFO=((PROTEUS_LAPACK_INTEGER)INFO);
  dgeev_(&JOBVL,
         &JOBVR,
         &laN,
         DDATA(A),
         &laLDA,
         DDATA(WR),
         DDATA(WI),
         DDATA(VL),
         &laLDVL,
         DDATA(VR),
         &laLDVR,
         DDATA(WORK),
         &laLWORK,
         &laINFO);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef lapackWrappersMethods[] = {
  { "denseFactorPrepare",
    lapackWrappersDenseFactorPrepare,
    METH_VARARGS,
    "factor the dense  matrix"},
  {"denseFactorSolve",
   lapackWrappersDenseFactorSolve,
   METH_VARARGS,
   "solve the dense linear system that was previously factorized"},
  {"denseCalculateEigenvalues",
   lapackWrappersDenseCalculateEigenvalues,
   METH_VARARGS,
   "calculate eigenvalues and  eigenvectors of a dense  matrix"},
  {NULL}
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initlapackWrappers(void) 
{
    PyObject *m,*d;

    if (PyType_Ready(&DenseFactorType) < 0)
      return;
    
    m = Py_InitModule3("lapackWrappers", 
                       lapackWrappersMethods,
                       "lapack wrappers module");
    d = PyModule_GetDict(m);
    import_array();
    Py_INCREF(&DenseFactorType);
    PyModule_AddObject(m, "DenseFactor", (PyObject *)&DenseFactorType);
}
/** @} */
