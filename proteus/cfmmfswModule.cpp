#include "Python.h"
#include "numpy/arrayobject.h"
#include "cfmmfsw.h"
#include "mesh.h"
#include "cmeshToolsModule.h"
#include "FMMandFSW.h"
#include <iostream>

#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

extern "C"
{
static PyObject*
FMMEikonalSolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FMMEikonalSolver_wrapper *self;
  self = (FMMEikonalSolver_wrapper *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
FMMEikonalSolver_init(FMMEikonalSolver_wrapper* self, PyObject *args, PyObject *kwds)
{
  PyObject *Cmesh;
  int nSpace;
  if(!PyArg_ParseTuple(args,
                       "iO",
                       &nSpace,
		       &Cmesh))

    return -1;
  
  if (FMMEikonalSolver_wrapper_init(&(MESH(Cmesh)),
				    nSpace,
				    self))
    {
      PyErr_NoMemory();
      return -1;
    }
  return 0;
}
static  void
FMMEikonalSolver_dealloc(FMMEikonalSolver_wrapper* self)
{
  FMMEikonalSolver_wrapper_del(self);

  self->ob_type->tp_free((PyObject*)self);
}

static PyObject*
FMMEikonalSolver_copy(FMMEikonalSolver_wrapper* self, PyObject* args)
{

  PyObject *other;
  if (!PyArg_ParseTuple(args,
			"O",
			&other))
    return NULL;
  assert(self);
  assert(other);

  FMMEikonalSolver_wrapper_copy((FMMEikonalSolver_wrapper *)other,
				self);
  
  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject*
FMMEikonalSolver_solve(FMMEikonalSolver_wrapper* self, PyObject* args,
		       PyObject* keywds)
{

  PyObject *phi0, *nodalSpeeds, *T;
  int failed;
  double zeroTol=1.0e-4,trialTol=1.0e-1;
  int initFlag=-1,//ignored if -1, 1 --> frontIntersection, 0 --> magnitude only
    verbose=0; 
  /*requires keyword entry in list for all arguments, but still seems to enforce positional
    ones using | */
  static char *kwlist[] = {"phi0","nodalSpeeds","T","zeroTol","trialTol","initFlag",
			   "verbose",NULL};

  if (!PyArg_ParseTupleAndKeywords(args,keywds,
				   "OOO|ddii",kwlist,
				   &phi0,
				   &nodalSpeeds,
				   &T,
				   &zeroTol,
				   &trialTol,
				   &initFlag,
				   &verbose))
    return NULL;
  assert(self);
  //mwf hack
  //printf("\n cfmmfswModule about to call FMMEikonalSOlver_wrapper verbose= %d \n",verbose);

  failed = FMMEikonalSolver_wrapper_solve(self,
					  DDATA(phi0),
					  DDATA(nodalSpeeds),
					  DDATA(T),
					  zeroTol,
					  trialTol,
					  initFlag,
					  verbose);

  
  return Py_BuildValue("i",failed);

}

static PyMethodDef FMMEikonalSolverMethods[] = {
  {"copy",
   (PyCFunction)FMMEikonalSolver_copy,
   METH_VARARGS,
   "copy other solver"},
  {"solve",
   (PyCFunction)FMMEikonalSolver_solve,
   METH_VARARGS | METH_KEYWORDS,
   "solve eikonal equation with variable speed"},
  {NULL} /*end*/
};

static PyTypeObject FMMEikonalSolverType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cfmmfsw.FMMEikonalSolver",             /*tp_name*/
  sizeof(FMMEikonalSolver_wrapper), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)FMMEikonalSolver_dealloc,                         /*tp_dealloc*/
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
  "FMMSolver wrapper ",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  FMMEikonalSolverMethods,     /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)FMMEikonalSolver_init,      /* tp_init */
  0,                         /* tp_alloc */
  FMMEikonalSolver_new,                 /* tp_new */
};

/**********************************************************************
 FSW wrapper 
***********************************************************************/

static PyObject*
FSWEikonalSolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FSWEikonalSolver_wrapper *self;
  self = (FSWEikonalSolver_wrapper *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
FSWEikonalSolver_init(FSWEikonalSolver_wrapper* self, PyObject *args, PyObject *keywds)
{
  PyObject *Cmesh,*refPoints;
  int nSpace,maxIts,initFlag,nRefPoints;
  double atol,rtol;

  initFlag  = 0;/*MAGNITUDE*/
  atol = 1.0e-8; rtol = 1.0e-8; maxIts=1000;
  refPoints = NULL;
  nRefPoints= 0;
  /*requires keyword entry in list for all arguments, but still seems to enforce positional
    ones using | */
  static char *kwlist[] = {"nSpace","mesh","atol","rtol","maxIts","initFlag",
			   "nRefPoints","refPoints",NULL};
  
  if(!PyArg_ParseTupleAndKeywords(args,keywds,
				  "iO|ddiiiO",kwlist,
				  &nSpace,
				  &Cmesh,
				  &atol,
				  &rtol,
				  &maxIts,
				  &initFlag,
				  &nRefPoints,
				  &refPoints))

    return -1;
  if (nRefPoints > 0 && refPoints == NULL)
    return -1;
  if (nRefPoints == 0 && refPoints != NULL)
    return -1;
  if (nRefPoints > 0)
    {
      if (FSWEikonalSolver_wrapper_init(&(MESH(Cmesh)),
					nSpace,
					atol,
					rtol,
					maxIts,
					initFlag,
					nRefPoints,
					DDATA(refPoints),
					self))
	{
	  PyErr_NoMemory();
	  return -1;
	}
    }
  else 
    {
      if (FSWEikonalSolver_wrapper_init(&(MESH(Cmesh)),
					nSpace,
					atol,
					rtol,
					maxIts,
					initFlag,
					nRefPoints,
					NULL,
					self))
	{
	  PyErr_NoMemory();
	  return -1;
	}
    }

  return 0;
}
static  void
FSWEikonalSolver_dealloc(FSWEikonalSolver_wrapper* self)
{
  FSWEikonalSolver_wrapper_del(self);

  self->ob_type->tp_free((PyObject*)self);
}

static PyObject*
FSWEikonalSolver_copy(FSWEikonalSolver_wrapper* self, PyObject* args)
{

  PyObject *other;
  if (!PyArg_ParseTuple(args,
			"O",
			&other))
    return NULL;
  assert(self);
  assert(other);

  FSWEikonalSolver_wrapper_copy((FSWEikonalSolver_wrapper *)other,
				self);
  
  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject*
FSWEikonalSolver_solve(FSWEikonalSolver_wrapper* self, PyObject* args,
		       PyObject* keywds)
{

  PyObject *phi0, *nodalSpeeds, *T;
  int failed;
  double zeroTol=1.0e-4,trialTol=1.0e-1;
  int initFlag=-1,//ignored if -1, 1 --> frontIntersection, 0 --> magnitude only
    verbose=0; 
  /*requires keyword entry in list for all arguments, but still seems to enforce positional
    ones using | */
  static char *kwlist[] = {"phi0","nodalSpeeds","T","zeroTol","trialTol","initFlag",
			   "verbose",NULL};

  if (!PyArg_ParseTupleAndKeywords(args,keywds,
				   "OOO|ddii",kwlist,
				   &phi0,
				   &nodalSpeeds,
				   &T,
				   &zeroTol,
				   &trialTol,
				   &initFlag,
				   &verbose))
    return NULL;
  assert(self);

  failed = FSWEikonalSolver_wrapper_solve(self,
					  DDATA(phi0),
					  DDATA(nodalSpeeds),
					  DDATA(T),
					  zeroTol,
					  trialTol,
					  initFlag,
					  verbose);

  
  return Py_BuildValue("i",failed);

}

static PyMethodDef FSWEikonalSolverMethods[] = {
  {"copy",
   (PyCFunction)FSWEikonalSolver_copy,
   METH_VARARGS,
   "copy other solver"},
  {"solve",
   (PyCFunction)FSWEikonalSolver_solve,
   METH_VARARGS | METH_KEYWORDS,
   "solve eikonal equation with variable speed"},
  {NULL} /*end*/
};

static PyTypeObject FSWEikonalSolverType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cfmmfsw.FSWEikonalSolver",             /*tp_name*/
  sizeof(FSWEikonalSolver_wrapper), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)FSWEikonalSolver_dealloc,                         /*tp_dealloc*/
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
  "FMMSolver wrapper ",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  FSWEikonalSolverMethods,     /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)FSWEikonalSolver_init,      /* tp_init */
  0,                         /* tp_alloc */
  FSWEikonalSolver_new,                 /* tp_new */
};

/***********************************************************************
 module functions
 **********************************************************************/
static PyObject*
cfmmfswLocalPWLreconstruction(PyObject* self, PyObject* args)
{
  PyObject *phi0,*phi0R,*Cmesh;
  double zeroTol;
  int nSpace,verbose;
  if (!PyArg_ParseTuple(args,"iOOOdi",
			&nSpace,
			&Cmesh,
			&phi0,
			&phi0R,
			&zeroTol,
			&verbose))
    return NULL;

  localPWL_wrapper(nSpace,
		   &(MESH(Cmesh)),
		   DDATA(phi0),
		   DDATA(phi0R),
		   zeroTol,
		   verbose);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfmmfswCopyOverOutsideBand(PyObject* self, PyObject* args)
{
  PyObject *out,*in;
  double tolerance;
  if (!PyArg_ParseTuple(args,"dOO",
			&tolerance,
			&in,
			&out))
    return NULL;

  bool failed = NearFrontInitialization::copyOverOutsideBand(SHAPE(in)[0],
							     tolerance,
							     DDATA(in),
							     DDATA(out));
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef cfmmfswMethods[] = {
  { "localPWLreconstruction", 
    cfmmfswLocalPWLreconstruction,
    METH_VARARGS, 
    "simple pwl reconstruction algorithm" },
  { "copyOverOutsideBand", 
    cfmmfswCopyOverOutsideBand,
    METH_VARARGS, 
    "copy from in to out where |out| > tolerance" },
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcfmmfsw(void)
{
  PyObject *m,*d;
  if (PyType_Ready(&FMMEikonalSolverType) < 0)
    return;
  if (PyType_Ready(&FSWEikonalSolverType) < 0)
    return;
  m = Py_InitModule3("cfmmfsw",cfmmfswMethods ,"unstructured FMM and FSW module");

  Py_INCREF(&FMMEikonalSolverType);
  PyModule_AddObject(m, "FMMEikonalSolver", (PyObject *)&FMMEikonalSolverType);
  Py_INCREF(&FSWEikonalSolverType);
  PyModule_AddObject(m, "FSWEikonalSolver", (PyObject *)&FSWEikonalSolverType);
  d = PyModule_GetDict(m);
  import_array();
}
}
