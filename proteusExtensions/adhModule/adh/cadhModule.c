#include "Python.h"
#include "numpy/arrayobject.h"
#include "share_declare.h"
#define ADH_SHARE_INCLUDED
#include "cadhimpl.h"

#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

/** Run adh "natively" the body of this function is a very slightly modified adh.c*/
static PyObject* cadhModule_cadhRun(PyObject* self,
				    PyObject* args)
{
  PyObject *filename,*runname;
  if (!PyArg_ParseTuple(args,
                        "OO",
                        &filename,
			&runname))
    return NULL;
  char *argv[3];
  int argc=3;
  int failure = 0;
  argv[0] = "adh";
  argv[1] = PyString_AsString(filename);/* allow a list of runnames later */
  argv[2] = PyString_AsString(runname);/* allow a list of runnames later */
  
  printf("preadh\n");
  /*preadh*/
  /*could do this in python in an ADH.py module but I'm trying to keep all the cadh stuff in this module*/
  /*note, you can't link the pre_adh code with the adh code because the use incompatible libraries, so I just do a system call*/
  const char* fs = PyString_AsString(filename);
  const char* preadh = "${PROTEUS_PACKAGES}/adh/bin/pre_adh ";
  size_t clen = strlen(preadh)+strlen(fs);
  char command[clen];
  strcpy(command,preadh);
  strcat(command,fs);
  printf("Running %s \n",command);
  system(command);

  failure =  cadhMain(argc,argv);
  return Py_BuildValue("i",failure);
}

/** The C representation of of ADH_NumericalSolution, a class for
    representing the entire ADH numerical solution, which may consist
    of multiple equation sets. This is basically the local variables from adh.c */
typedef struct
{
  PyObject_HEAD
  ADH_NumericalSolution adh;
} cADH_NumericalSolution;

/**Do some stuff before taking a step*/
static PyObject*
cADH_NumericalSolution_step(cADH_NumericalSolution *self, PyObject *args, PyObject *kwds)
{
  int failure = 0;
  failure = ADH_NumericalSolution_step(&self->adh);
  return Py_BuildValue("i",failure);
}

/**Do some things after taking a step*/
static PyObject*
cADH_NumericalSolution_stepTaken(cADH_NumericalSolution *self, PyObject *args, PyObject *kwds)
{
  int notdone = 0;
  notdone = ADH_NumericalSolution_stepTaken(&self->adh);
  return Py_BuildValue("i",notdone);
}

/**Instatiate a new cADH_NumericalSolution object*/
static PyObject*
cADH_NumericalSolution_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  cADH_NumericalSolution *self;
  self = (cADH_NumericalSolution *)type->tp_alloc(type,0);
  return (PyObject*)self;
}


/**Initialize a NumericalSolution object */
static int
cADH_NumericalSolution_init(cADH_NumericalSolution* self,PyObject* args,PyObject *kwds)
{
  PyObject *filename,*runname;
  if (!PyArg_ParseTuple(args,
                        "OO",
                        &filename,
			&runname))
    return -1;
  char *argv[3];
  int argc=3;
  int failure = 0;
  argv[0] = "adh";
  argv[1] = PyString_AsString(filename);/* allow a list of runnames later */
  argv[2] = PyString_AsString(runname);/* allow a list of runnames later */
  
  printf("preadh\n");
  /*preadh*/
  /*could do this in python in an ADH.py module but I'm trying to keep all the cadh stuff in this module*/
  /*note, you can't link the pre_adh code with the adh code because the use incompatible libraries, so I just do a system call*/
  const char* fs = PyString_AsString(filename);
  const char* preadh = "${PROTEUS_PACKAGES}/adh/bin/pre_adh ";
  size_t clen = strlen(preadh)+strlen(fs);
  char command[clen];
  strcpy(command,preadh);
  strcat(command,fs);
  printf("Running %s \n",command);
  system(command);

  failure = ADH_NumericalSolution_init(&self->adh,argc,argv);
  return failure;
}

/** Deallocate the memory for the NumericalSolution object,including file scope variables */
static int 
cADH_NumericalSolution_dealloc(cADH_NumericalSolution *self, PyObject *args, PyObject *kwds)
{
  int failure = ADH_NumericalSolution_dealloc(&self->adh);
  self->ob_type->tp_free((PyObject*)self);
  return failure;
}

/** Loop over the time steps from the input file and calculation solution. This is the central part of adh.c */
static cADH_NumericalSolution* 
cADH_NumericalSolution_calculateSolution(cADH_NumericalSolution* self, 
					 PyObject *args)
{
  int failure = ADH_NumericalSolution_calculateSolution(&self->adh);
  return self;
}

/** I was getting ready to try to evaluate individual residuals and jacobians */
/* static PyObject*  */
/* ADH_gw_resid(ADH*self,  */
/* 	     PyObject *args) */
/* { */
/*   int iwhich_flag; */
/*   PyObject *bc_mask,*residual; */
/*   if (!PyArg_ParseTuple(args, */
/*                         "iOO", */
/* 			iwhich_flag, */
/*                         &bc_mask, */
/* 			&residual)) */
/*     return -1; */
/*   fe_gw_resid(iwhich_flag,   /\* dummy flag - should be 0 *\/ */
/* 	      IDATA(bc_mask),  /\* boundary condition mask *\/ */
/*               DDATA(residual));   /\* the residual *\/ */
/*   return self */
/* } */

/* static PyObject*  */
/* ADH_gw_load(ADH*self,  */
/* 	    PyObject *args) */
/* { */
/*   int iwhich_flag; */
/*   PyObject *bc_mask,*jacobian,*residual; */
/*   if (!PyArg_ParseTuple(args, */
/*                         "iOO", */
/* 			iwhich_flag, */
/*                         &bc_mask, */
/* 			&residual)) */
/*     return -1; */
/*   fe_gw_load(iwhich_flag,		/\* dummy flag - should be 0 *\/ */
/* 	     bc_mask,			/\* boundary condition mask *\/ */
/* 	     self.diagonal,		/\* the diagonal of the matrix *\/ */
/* 	     self.matrix,		/\* the matrix *\/ */
/* 	     DDATA(residual));   /\* the residual *\/ */
/*   return self */
/* } */

/** The methods of the cADH_NumericalSolution objects */
static PyMethodDef cADH_NumericalSolution_methods[] = {
  {"calculateSolution", 
   (PyCFunction)cADH_NumericalSolution_calculateSolution, 
   METH_VARARGS, 
   "calculates the solutions for the time steps in the ADH input files"},
  {"step", 
   (PyCFunction)cADH_NumericalSolution_step, 
   METH_VARARGS, 
   "prepare for next step or retry"},
  {"stepTaken", 
   (PyCFunction)cADH_NumericalSolution_stepTaken, 
   METH_VARARGS, 
   "clean up after solver step"},
  {NULL,NULL}
};

/** The type definition of the cADH_NumericalSolution objects */
static PyTypeObject cADH_NumericalSolutionType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cadh.cADH_NumericalSolution",             /*tp_name*/
  sizeof(cADH_NumericalSolution), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)cADH_NumericalSolution_dealloc,                         /*tp_dealloc*/
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
  "sparse factor objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  cADH_NumericalSolution_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)cADH_NumericalSolution_init,      /* tp_init */
  0,                         /* tp_alloc */
  cADH_NumericalSolution_new,                 /* tp_new */
};

/**The C representation of an ADH_OneLevelTransport object. These are
   basically the variables in fe_main, except that I added a pointer
   to a cADH_NumericalSolution object to make sure none of the global
   variables get deleted before where done with the transport
   object. */ 
typedef struct
{
  PyObject_HEAD
  ADH_OneLevelTransport transport;
  /*cek additions */
  cADH_NumericalSolution* numericalSolution;
} cADH_OneLevelTransport;

/**Instantiate a cADH_OneLevelTransport object */
static PyObject*
cADH_OneLevelTransport_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  cADH_OneLevelTransport *self;
  self = (cADH_OneLevelTransport *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

/**Initialize a cADH_OneLevelTransport object. This is the first part of fe_main*/
static int
cADH_OneLevelTransport_init(cADH_OneLevelTransport* self,PyObject* args,PyObject *kwds)
{
  PyObject *nsPyObject;
  int failure = 0;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &nsPyObject))
    return -1;
  self->numericalSolution = (cADH_NumericalSolution*)(nsPyObject);
  Py_INCREF(nsPyObject);
  failure = ADH_OneLevelTransport_init(&self->transport);
  return failure;
}

/**Allocate rellocations after adaption */
static int
cADH_OneLevelTransport_realloc(cADH_OneLevelTransport* self,PyObject* args,PyObject *kwds)
{
  /* PyObject *nsPyObject; */
  /* if (!PyArg_ParseTuple(args, */
  /*                       "O", */
  /*                       &nsPyObject)) */
  /*   return -1; */
  /* cADH_NumericalSolution* numericalSolution = (cADH_NumericalSolution*)(nsPyObject); */
  realloc_fe_main(&self->transport);
  return 0;
}

static int 
cADH_OneLevelTransport_dealloc(cADH_OneLevelTransport *self, PyObject *args, PyObject *kwds)
{
  int failure = ADH_OneLevelTransport_dealloc(&self->transport);
  return failure;
}

/** Call the ADH Newton solves for all the models. This will kick out
    with a failure for any model so the caller needs handle nonlinear
    solver failures */
static PyObject*
cADH_OneLevelTransport_solve(cADH_OneLevelTransport* self, PyObject* args)
{
  int failure = 0;
  failure = ADH_OneLevelTransport_solve(&self->transport);
  return Py_BuildValue("i",failure);
}

static PyMethodDef cADH_OneLevelTransport_methods[] = {
  {"solve", 
   (PyCFunction)cADH_OneLevelTransport_solve, 
   METH_VARARGS, 
   "solves the nonlinear systems at a time step"},
  {"realloc", 
   (PyCFunction)cADH_OneLevelTransport_realloc, 
   METH_VARARGS, 
   "reallocs memory"},
  {NULL,NULL}
};


static PyTypeObject cADH_OneLevelTransportType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cadh.cADH_OneLevelTransport",             /*tp_name*/
  sizeof(cADH_OneLevelTransport), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)cADH_OneLevelTransport_dealloc,                         /*tp_dealloc*/
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
  "sparse factor objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  cADH_OneLevelTransport_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)cADH_OneLevelTransport_init,      /* tp_init */
  0,                         /* tp_alloc */
  cADH_OneLevelTransport_new,                 /* tp_new */
};

static PyMethodDef cadhMethods[] = {
 { "cadhRun",
    cadhModule_cadhRun,
    METH_VARARGS, 
    "Call cadh like you would on the command line"},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcadh(void)
{
  PyObject *m,*d;
  if (PyType_Ready(&cADH_NumericalSolutionType) < 0)
    return;
  if (PyType_Ready(&cADH_OneLevelTransportType) < 0)
    return;
  m = Py_InitModule("cadh", cadhMethods);
  d = PyModule_GetDict(m);
  import_array();
  Py_INCREF(&cADH_NumericalSolutionType);
  PyModule_AddObject(m, "cADH_NumericalSolution", (PyObject *)&cADH_NumericalSolutionType);
  Py_INCREF(&cADH_OneLevelTransportType);
  PyModule_AddObject(m, "cADH_OneLevelTransport", (PyObject *)&cADH_OneLevelTransportType);
}
