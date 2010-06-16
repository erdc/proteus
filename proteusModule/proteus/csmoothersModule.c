#include "Python.h"
#include "smoothers.h"
#include "numpy/arrayobject.h"
#include "supermatrix.h"
#include <stdio.h>
#include "superluWrappersModule.h"
#include "csmoothersModule.h"
/** \file csmoothersModule.c
    \ingroup csmoothers csmoothers
    @{
*/
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

static PyObject*
csmoothers_jacobi_NR_prepare(PyObject* self, 
                            PyObject* args)
{
  PyObject *A,*M;
  double w,tol;
  SuperMatrix AS;
  if(!PyArg_ParseTuple(args,"OddO",
                       &A,
                       &w,
                       &tol,
                       &M))
    return NULL;
  AS.Stype = SLU_NC;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  jacobi_NR_prepare(&AS,
                    w,
                    tol,
                    DDATA(M));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_jacobi_NR_solve(PyObject* self, 
                          PyObject* args)
{
  PyObject *A,*M,*dX,*R,*node_order;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &A,
                       &M,
                       &R,
                       &node_order,
                       &dX))
    return NULL;
  SuperMatrix AS;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  jacobi_NR_solve(&AS,
                  DDATA(M),
                  DDATA(R),
                  IDATA(node_order),
                  DDATA(dX));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_nl_jacobi_NR_solve(PyObject* self, 
                             PyObject* args)
{
  PyObject *A,*dX,*R,*node_order;
  double w,tol;
  if(!PyArg_ParseTuple(args,"OOOddO",
                       &A,
                       &R,
                       &node_order,
                       &w,
                       &tol,
                       &dX))
    return NULL;
  SuperMatrix AS;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  nl_jacobi_NR_solve(&AS,
                     DDATA(R),
                     IDATA(node_order),
                     w,
                     tol,
                     DDATA(dX));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_gauss_seidel_NR_prepare(PyObject* self, 
                                  PyObject* args)
{
  PyObject *A,*M;
  double w,tol;
  if(!PyArg_ParseTuple(args,"OddO",
                       &A,
                       &w,
                       &tol,
                       &M))
    return NULL;
  SuperMatrix AS;
  AS.Stype = SLU_NC;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  gauss_seidel_NR_prepare(&AS,
                          w,
                          tol,
                          DDATA(M));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_gauss_seidel_NR_solve(PyObject* self, 
                                PyObject* args)
{
  PyObject *A,*M,*dX,*R,*node_order;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &A,
                       &M,
                       &R,
                       &node_order,
                       &dX))
    return NULL;
  SuperMatrix AS;
  AS.Stype = SLU_NC;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  gauss_seidel_NR_solve(&AS,
                        DDATA(M),
                        DDATA(R),
                        IDATA(node_order),
                        DDATA(dX));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_nl_gauss_seidel_NR_solve(PyObject* self, 
                                   PyObject* args)
{
  PyObject *A,*dX,*R,*node_order;
  double w,tol;
  if(!PyArg_ParseTuple(args,"OOOddO",
                       &A,
                       &R,
                       &node_order,
                       &w,
                       &tol,
                       &dX))
    return NULL;
  SuperMatrix AS;
  AS.Stype = SLU_NC;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  nl_gauss_seidel_NR_solve(&AS,
                           DDATA(R),
                           IDATA(node_order),
                           w,
                           tol,
                           DDATA(dX));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_asm_NR_prepare(PyObject* self, 
                         PyObject* args)
{
  PyObject *A,*asmFactor;
  SuperMatrix AS;
  if(!PyArg_ParseTuple(args,"OO",
                       &A,
                       &asmFactor))
    return NULL;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)A)->dim[0];
  AS.ncol = ((SparseMatrix*)A)->dim[1];
  AS.Store = &(((SparseMatrix*)A)->A);
  asm_NR_prepare(&AS, 
                 ((ASMFactor*)asmFactor)->subdomain_dim, 
                 ((ASMFactor*)asmFactor)->l2g_L, 
                 ((ASMFactor*)asmFactor)->subdomain_L, 
                 ((ASMFactor*)asmFactor)->subdomain_pivots);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_asm_NR_solve(PyObject* self, 
                       PyObject* args)
{
  PyObject *A,*asmFactor,*node_order,*R,*dX;
  double w;
  SuperMatrix AS;
  if(!PyArg_ParseTuple(args,"OdOOOO",
                       &A,
                       &w,
                       &asmFactor,
                       &node_order,
                       &R,
                       &dX))
    return NULL;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  asm_NR_solve(&AS, 
               w,
               ((ASMFactor*)asmFactor)->subdomain_L, 
               ((ASMFactor*)asmFactor)->subdomain_dim, 
               ((ASMFactor*)asmFactor)->l2g_L, 
               DDATA(R), 
               ((ASMFactor*)asmFactor)->subdomain_R, 
               IDATA(node_order),
               ((ASMFactor*)asmFactor)->subdomain_dX,
               DDATA(dX), 
               ((ASMFactor*)asmFactor)->subdomain_pivots);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ASMFactor_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ASMFactor *self;
  self = (ASMFactor *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
ASMFactor_init(ASMFactor *self, PyObject *args, PyObject *kwds)
{
  SuperMatrix AS;
  PyObject *A;
  if(!PyArg_ParseTuple(args,
                       "O",
                       &A))
    return -1;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)A)->dim[0];
  AS.ncol = ((SparseMatrix*)A)->dim[1];
  AS.Store = &(((SparseMatrix*)A)->A);
  if (asm_NR_init(&AS,
                  &self->subdomain_dim,
                  &self->l2g_L,
                  &self->subdomain_L,
                  &self->subdomain_R,
                  &self->subdomain_dX,
                  &self->subdomain_pivots))
    {
      PyErr_NoMemory();
      return -1;
    }
  return 0;
}

static  void
ASMFactor_dealloc(ASMFactor* self)
{
  asm_NR_free(self->N, 
              self->subdomain_dim,
              self->l2g_L,
              self->subdomain_L,
              self->subdomain_R,
              self->subdomain_dX,
              self->subdomain_pivots);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject ASMFactorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "asmFactor.ASMFactor",             /*tp_name*/
  sizeof(ASMFactor), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)ASMFactor_dealloc,                         /*tp_dealloc*/
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
  "asm factor objects",           /* tp_doc */
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
  (initproc)ASMFactor_init,      /* tp_init */
  0,                         /* tp_alloc */
  ASMFactor_new,                 /* tp_new */
};

static PyObject*
csmoothers_basm_NR_prepare(PyObject* self, 
			   PyObject* args)
{
  PyObject *A,*basmFactor;
  SuperMatrix AS;
  int N,bs;
  if(!PyArg_ParseTuple(args,"OO",
                       &A,
                       &basmFactor))
    return NULL;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)A)->dim[0];
  AS.ncol = ((SparseMatrix*)A)->dim[1];
  AS.Store = &(((SparseMatrix*)A)->A);
  N = ((BASMFactor*)basmFactor)->N;
  bs= ((BASMFactor*)basmFactor)->bs;
  basm_NR_prepare(bs,
		  N,
		  &AS, 
		  ((BASMFactor*)basmFactor)->subdomain_dim, 
		  ((BASMFactor*)basmFactor)->l2g_L, 
		  ((BASMFactor*)basmFactor)->subdomain_L, 
		  ((BASMFactor*)basmFactor)->subdomain_pivots,
		  ((BASMFactor*)basmFactor)->subdomain_col_pivots);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
csmoothers_basm_NR_solve(PyObject* self, 
			 PyObject* args)
{
  PyObject *A,*basmFactor,*node_order,*R,*dX;
  double w;
  SuperMatrix AS;
  int N,bs;
  if(!PyArg_ParseTuple(args,"OdOOOO",
                       &A,
                       &w,
                       &basmFactor,
                       &node_order,
                       &R,
                       &dX))
    return NULL;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)(A))->dim[0];
  AS.ncol = ((SparseMatrix*)(A))->dim[1];
  AS.Store = &((SparseMatrix*)(A))->A;
  N = ((BASMFactor*)basmFactor)->N;
  bs= ((BASMFactor*)basmFactor)->bs;
  basm_NR_solve(bs,
		N,
		&AS, 
		w,
		((BASMFactor*)basmFactor)->subdomain_L, 
		((BASMFactor*)basmFactor)->subdomain_dim, 
		((BASMFactor*)basmFactor)->l2g_L, 
		DDATA(R), 
		((BASMFactor*)basmFactor)->subdomain_R, 
		IDATA(node_order),
		((BASMFactor*)basmFactor)->subdomain_dX,
		DDATA(dX), 
		((BASMFactor*)basmFactor)->subdomain_pivots,
		((BASMFactor*)basmFactor)->subdomain_col_pivots);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
BASMFactor_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  BASMFactor *self;
  self = (BASMFactor *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
BASMFactor_init(BASMFactor *self, PyObject *args, PyObject *kwds)
{
  SuperMatrix AS;
  PyObject *A;
  int bs=1;
  if(!PyArg_ParseTuple(args,
                       "O|i",
                       &A,
		       &bs))
    return -1;
  AS.Stype = SLU_NR;
  AS.Dtype = SLU_D;
  AS.Mtype = SLU_GE;
  AS.nrow = ((SparseMatrix*)A)->dim[0];
  AS.ncol = ((SparseMatrix*)A)->dim[1];
  AS.Store = &(((SparseMatrix*)A)->A);
  if (basm_NR_init(bs,
		   &AS,
		   &self->subdomain_dim,
		   &self->l2g_L,
		   &self->subdomain_L,
		   &self->subdomain_R,
		   &self->subdomain_dX,
		   &self->subdomain_pivots,
		   &self->subdomain_col_pivots))
    {
      PyErr_NoMemory();
      return -1;
    }
  self->bs = bs;
  self->N  = AS.nrow/bs;
  assert (self->N*self->bs == AS.nrow);
  return 0;
}

static  void
BASMFactor_dealloc(BASMFactor* self)
{
  basm_NR_free(self->N, 
	       self->subdomain_dim,
	       self->l2g_L,
	       self->subdomain_L,
	       self->subdomain_R,
	       self->subdomain_dX,
	       self->subdomain_pivots,
	       self->subdomain_col_pivots);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject BASMFactorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "basmFactor.BASMFactor",             /*tp_name*/
  sizeof(BASMFactor), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)BASMFactor_dealloc,                         /*tp_dealloc*/
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
  "block asm factor objects",           /* tp_doc */
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
  (initproc)BASMFactor_init,      /* tp_init */
  0,                         /* tp_alloc */
  BASMFactor_new,                 /* tp_new */
};

static PyMethodDef csmoothersMethods[] = {
  { "jacobi_NR_prepare",
    csmoothers_jacobi_NR_prepare, 
    METH_VARARGS, 
    "add doc"}, 
  { "jacobi_NR_solve",
    csmoothers_jacobi_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { "nl_jacobi_NR_solve",
    csmoothers_nl_jacobi_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { "gauss_seidel_NR_prepare",
    csmoothers_gauss_seidel_NR_prepare, 
    METH_VARARGS, 
    "add doc"}, 
  { "gauss_seidel_NR_solve",
    csmoothers_gauss_seidel_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { "nl_gauss_seidel_NR_solve",
    csmoothers_nl_gauss_seidel_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { "asm_NR_prepare",
    csmoothers_asm_NR_prepare, 
    METH_VARARGS, 
    "add doc"}, 
  { "asm_NR_solve",
    csmoothers_asm_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { "basm_NR_prepare",
    csmoothers_basm_NR_prepare, 
    METH_VARARGS, 
    "add doc"}, 
  { "basm_NR_solve",
    csmoothers_basm_NR_solve, 
    METH_VARARGS, 
    "Add doc"}, 
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcsmoothers(void)
{
  PyObject *m,*d;
  if (PyType_Ready(&ASMFactorType) < 0)
    return;
  if (PyType_Ready(&BASMFactorType) < 0)
    return;
  m = Py_InitModule3("csmoothers", 
                     csmoothersMethods,
                     "csmoothers  module");
  Py_INCREF(&ASMFactorType);
  PyModule_AddObject(m, "ASMFactor", (PyObject *)&ASMFactorType);
  Py_INCREF(&BASMFactorType);
  PyModule_AddObject(m, "BASMFactor", (PyObject *)&BASMFactorType);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
