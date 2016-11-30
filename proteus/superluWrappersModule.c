#include "Python.h"
#include <stdio.h>
#include "numpy/arrayobject.h"
#include "superluWrappersModule.h"
/** \file superluWrappersModule.c
    \ingroup superluWrappers
    @{
*/
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject* 
SparseMatrix_matvec(SparseMatrix *self, 
                    PyObject *args)
{
  register int i,k;
  register double s;
  int_t *rowptr,*colind;
  double *a,*x,*y;
  PyObject *xp, *yp;  
  if(!PyArg_ParseTuple(args,
                       "OO",
                       &xp,
                       &yp))
    return NULL;
  a = (double*)self->A.nzval;
  rowptr = self->A.rowptr;
  colind = self->A.colind;
  x = DDATA(xp);
  y = DDATA(yp);
  for (i = 0; i < self->dim[0]; i ++) 
    {
      s = 0.0;
      for (k=rowptr[i]; k<rowptr[i+1]; k++)
        s += a[k]* x[colind[k]];
      y[i] = s;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
SparseMatrix_fwrite(SparseMatrix *self, 
                    PyObject *args)
{
  register int i,k;
  int base=0;
  int_t *rowptr,*colind;
  double *a;
  char *filename;
  FILE *file;
  if(!PyArg_ParseTuple(args,
                       "s|i",
                       &filename,
                       &base))
    return NULL;
  a = (double*)self->A.nzval;
  rowptr = self->A.rowptr;
  colind = self->A.colind;
  file = fopen(filename,"w");
  fprintf(file,"%i %i %i \n",self->dim[0],self->dim[1],rowptr[self->dim[0]]);
  /* printf("%i %i %i \n",self->dim[0],self->dim[1],rowptr[self->dim[0]]); */
  for (i = 0; i < self->dim[0]; i ++) 
    {
      for (k=rowptr[i]; k<rowptr[i+1]; k++)
	fprintf(file,"%d %d %13.8e\n",i+base,colind[k]+base,a[k]);
    }
  fclose(file);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
SparseMatrix_getCSRrepresentation(SparseMatrix *self, 
				  PyObject *args)
{
  PyObject *rowout,*colout,*aout;
  int nnz,nr;
  npy_intp dims[1];
  int_t *rowptr,*colind;
  double *a;
  a = (double*)self->A.nzval;
  rowptr = self->A.rowptr;
  colind = self->A.colind;
  nr  = self->dim[0];
  nnz = self->A.nnz;

  dims[0] = nr+1;
  rowout = PyArray_SimpleNewFromData(1,dims,NPY_INT,
				   (void*)rowptr);
  dims[0] = nnz;
  colout = PyArray_SimpleNewFromData(1,dims,NPY_INT,
				     (void*)colind);

  dims[0] = nnz;
  aout = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,
				   (void*)a);

  return Py_BuildValue("OOO",
		       rowout,
		       colout,
		       aout);
}

static PyObject* 
SparseMatrix_getSubMatCSRrepresentation(SparseMatrix *self, 
					PyObject *args)
{
  PyObject *rowout,*colout,*aout;
  int nnz,nr;
  npy_intp dims[1];
  int_t *rowptr,*colind;
  double *a;
  /*rows for submatrix */
  int range_start,range_end;
  if(!PyArg_ParseTuple(args,
                       "ii",
                       &range_start,
		       &range_end))
    return NULL;
  assert(range_start >= 0);
  assert(range_end <= self->dim[0]);
  assert(range_end > range_start);
  nr = range_end-range_start;
  assert(nr <= self->dim[0]);
  
  rowptr = self->A.rowptr + range_start;
  colind = self->A.colind + rowptr[0];
  a = (double*)(self->A.nzval) + rowptr[0];
  nnz = self->A.rowptr[range_end]-self->A.rowptr[range_start];

  dims[0] = nr+1;
  rowout = PyArray_SimpleNewFromData(1,dims,NPY_INT,
				     (void*)rowptr);
  dims[0] = nnz;
  colout = PyArray_SimpleNewFromData(1,dims,NPY_INT,
				     (void*)colind);

  dims[0] = nnz;
  aout = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,
				   (void*)a);

  return Py_BuildValue("OOO",
		       rowout,
		       colout,
		       aout);
}

static PyMethodDef SparseMatrix_methods[] = {
  {"matvec", 
   (PyCFunction)SparseMatrix_matvec, 
   METH_VARARGS, 
   "a.matvec(x, y)\n \
\n \
compute the sparse matrix-vector product y := a * x. \n \
a is a d1 by d2 sparse matrix.\n \
x and y are two 1-dimensional Numeric arrays of appropriate size.\n"},
  {"fwrite",
   (PyCFunction)SparseMatrix_fwrite,
   METH_VARARGS,
   "a.fwrite('filename') \
\n \
write the sparse matrix to a file. \n"},
  {"getCSRrepresentation",
   (PyCFunction)SparseMatrix_getCSRrepresentation,
   METH_NOARGS,
   "return CSR representation (rowptr,colind,nnzval) as Numpy arrays"},
  {"getSubMatCSRrepresentation",
   (PyCFunction)SparseMatrix_getSubMatCSRrepresentation,
   METH_VARARGS,
   "return CSR representation (rowptr,colind,nnzval) as Numpy arrays for sub matrix corresponding to local rows [istart,iend)"},
  {NULL,NULL}
};

#define SMP(p) ((SparseMatrix*))p)

static PyObject*
SparseMatrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  SparseMatrix *self;
  self = (SparseMatrix *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
SparseMatrix_init(SparseMatrix *self, PyObject *args, PyObject *kwds)
{
  int i,nr,nc,nnz;
  PyObject *nzval,*colind,*rowptr;
  if(!PyArg_ParseTuple(args,
                       "iiiOOO",
                       &nr,
                       &nc,
                       &nnz,
                       &nzval,
                       &colind,
                       &rowptr))
    return -1;
  self->dim[0]=nr;
  self->dim[1]=nc;
  self->A.nnz = nnz;
  self->A.nzval = DDATA(nzval);
  self->A.colind = (int_t*)malloc(nnz*sizeof(int_t));
  for (i=0;i<nnz;i++)
    self->A.colind[i] = IDATA(colind)[i];
  self->A.rowptr = (int_t*)malloc((nr+1)*sizeof(int_t));
  for (i=0;i<nr+1;i++)
    self->A.rowptr[i] = IDATA(rowptr)[i];
  SuperMatrix Aprint;
  Aprint.Stype = SLU_NC;
  Aprint.Dtype = SLU_D;
  Aprint.Mtype = SLU_GE;
  Aprint.nrow = nr;
  Aprint.ncol = nc;
  Aprint.Store = &self->A;
  return 0;
}

static  void
SparseMatrix_dealloc(SparseMatrix* self)
{
  free(self->A.colind);
  free(self->A.rowptr);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
SparseMatrix_getattr(SparseMatrix *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->dim[0], self->dim[1]);
  if (strcmp(name, "nnz") == 0)
    return PyInt_FromLong(self->A.nnz);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
	PyList_SetItem(list, i, PyString_FromString(members[i]));
      if (PyErr_Occurred()) {
	Py_DECREF(list);
	list = NULL;
      }
    }
    return list;
  }
  return Py_FindMethod(SparseMatrix_methods, (PyObject *)self, name);
}

static PyTypeObject SparseMatrixType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "superluWrappers.SparseMatrix",             /*tp_name*/
  sizeof(SparseMatrix), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)SparseMatrix_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  (getattrfunc)SparseMatrix_getattr,                         /*tp_getattr*/
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
  SparseMatrix_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)SparseMatrix_init,      /* tp_init */
  0,                         /* tp_alloc */
  SparseMatrix_new,                 /* tp_new */
};

typedef struct 
{
  PyObject_HEAD
  superlu_options_t options;
  SuperMatrix A,AC,L,U,X;
  GlobalLU_t Glu;
  NCformat storeA;
  DNformat storeX;
  SuperLUStat_t stat;
  int use_same_perm_c,use_same_sparsity;
  int *perm_c,*perm_r,*etree;
} SparseFactor;

#define SFP(p) ((SparseFactor*)p)

static PyObject*
SparseFactor_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  SparseFactor *self;
  self = (SparseFactor *)type->tp_alloc(type,0);
  self->A.Stype = SLU_NC;
  self->A.Dtype = SLU_D;
  self->A.Mtype = SLU_GE;
  self->A.Store = &self->storeA;

  self->AC.Stype = SLU_NCP;
  self->AC.Dtype = SLU_D;
  self->AC.Mtype = SLU_GE;
  self->AC.Store = NULL;

  self->L.Stype = SLU_SC;
  self->L.Dtype = SLU_D;
  self->L.Mtype = SLU_TRLU;
  self->L.Store = NULL;

  self->U.Stype = SLU_NC;
  self->U.Dtype = SLU_D;
  self->U.Mtype = SLU_TRU;
  self->U.Store = NULL;

  self->X.Stype = SLU_DN;
  self->X.Dtype = SLU_D;
  self->X.Mtype = SLU_GE;
  self->X.Store = &self->storeX;
  return (PyObject*)self;
}

static int
SparseFactor_init(SparseFactor *self, PyObject *args, PyObject *kwds)
{
  int dim;
  if(!PyArg_ParseTuple(args,
                       "i",
                       &dim))
    return -1;
  self->A.nrow = dim;
  self->A.ncol = dim;

  self->L.nrow = dim;
  self->L.ncol = dim;

  self->U.nrow = dim;
  self->U.ncol = dim;

  self->storeX.lda = dim;
  self->X.nrow = dim;
  self->X.ncol = 1;

  self->perm_c = (int*)malloc(dim*sizeof(int));
  self->perm_r = (int*)malloc(dim*sizeof(int));
  self->etree  = (int*)malloc(dim*sizeof(int));
  StatInit(&(self->stat));
  set_default_options(&(self->options));

  self->use_same_perm_c = 0;
  self->use_same_sparsity = 0;
  if ( (self->perm_c == NULL) || (self->perm_r == NULL) || (self->etree == NULL))
    {
      PyErr_NoMemory();
      return -1;
    }
  return 0;
}

static  void
SparseFactor_dealloc(SparseFactor* self)
{
  free(self->perm_c);
  free(self->perm_r);
  free(self->etree);
  if (self->AC.Store != NULL)
    Destroy_CompCol_Permuted(&self->AC);
  if (self->L.Store != NULL)
    Destroy_SuperNode_Matrix(&self->L);
  if (self->U.Store != NULL)
    Destroy_CompCol_Matrix(&self->U);
  StatFree(&self->stat);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject SparseFactorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "superluWrappers.SparseFactor",             /*tp_name*/
  sizeof(SparseFactor), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)SparseFactor_dealloc,                         /*tp_dealloc*/
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
  0,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)SparseFactor_init,      /* tp_init */
  0,                         /* tp_alloc */
  SparseFactor_new,                 /* tp_new */
};


static PyObject* superluWrappersSparseFactorPrepare(PyObject* self,
                                                    PyObject* args)
{
  int i,n,relax=1,panel_size=10,lwork=0,info=0,permc_spec=3;
  double drop_tol=-1.0;/* not used by superlu */
  void *work=NULL;
  PyObject *mat,*sparseFactor;
  if(!PyArg_ParseTuple(args,"OO",
                       &mat,
                       &sparseFactor))
    return NULL;
  SFP(sparseFactor)->storeA.nnz = ((SparseMatrix*)mat)->A.nnz;
  SFP(sparseFactor)->storeA.nzval = ((SparseMatrix*)mat)->A.nzval;
  SFP(sparseFactor)->storeA.colptr = ((SparseMatrix*)mat)->A.rowptr;
  SFP(sparseFactor)->storeA.rowind = ((SparseMatrix*)mat)->A.colind;
  /* calc column permutation */
  if ( SFP(sparseFactor)->use_same_perm_c == 0)
    {
      get_perm_c(permc_spec, 
                 &SFP(sparseFactor)->A, 
                 SFP(sparseFactor)->perm_c);
      SFP(sparseFactor)->use_same_perm_c = 1;
    }
  if ( SFP(sparseFactor)->use_same_sparsity == 0)
    {
      if (SFP(sparseFactor)->AC.Store  != NULL)
        {
          Destroy_CompCol_Permuted(&SFP(sparseFactor)->AC);
          Destroy_SuperNode_Matrix(&SFP(sparseFactor)->L);
          Destroy_CompCol_Matrix(&SFP(sparseFactor)->U);
        }
      /* apply column permutation and build AC and etree*/
      sp_preorder(&SFP(sparseFactor)->options, 
                  &SFP(sparseFactor)->A, 
                  SFP(sparseFactor)->perm_c, 
                  SFP(sparseFactor)->etree, 
                  &SFP(sparseFactor)->AC);
      SFP(sparseFactor)->use_same_sparsity = 1;
    }
  else
    {
      /* apply column permutation */
      SFP(sparseFactor)->options.Fact = SamePattern_SameRowPerm;
      n = SFP(sparseFactor)->A.ncol;
      for (i = 0; i < n; i++) 
        {
          ((NCPformat*)SFP(sparseFactor)->AC.Store)->colbeg[SFP(sparseFactor)->perm_c[i]] = ((NCformat*)SFP(sparseFactor)->A.Store)->colptr[i]; 
          ((NCPformat*)SFP(sparseFactor)->AC.Store)->colend[SFP(sparseFactor)->perm_c[i]] = ((NCformat*)SFP(sparseFactor)->A.Store)->colptr[i+1];
        }
    }
  dgstrf(&SFP(sparseFactor)->options,
         &SFP(sparseFactor)->AC,
         relax,
         panel_size,
         SFP(sparseFactor)->etree,
         work,
         lwork,
         SFP(sparseFactor)->perm_c,
         SFP(sparseFactor)->perm_r,
         &SFP(sparseFactor)->L,
         &SFP(sparseFactor)->U,
	 &SFP(sparseFactor)->Glu,
         &SFP(sparseFactor)->stat,
         &info);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* superluWrappersSparseFactorSolve(PyObject* self,
                                                PyObject* args)
{
  trans_t trans=TRANS;
  int info=0;
  PyObject *sparseFactor,*x;
  if(!PyArg_ParseTuple(args,"OO",
                       &sparseFactor,
                       &x))
    return NULL;
  SFP(sparseFactor)->storeX.nzval = DDATA(x);
  dgstrs(trans,
         &SFP(sparseFactor)->L,
         &SFP(sparseFactor)->U,
         SFP(sparseFactor)->perm_c,
         SFP(sparseFactor)->perm_r,
         &SFP(sparseFactor)->X,
         &SFP(sparseFactor)->stat,
         &info);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef superluWrappersMethods[] = {
  { "sparseFactorPrepare",
    superluWrappersSparseFactorPrepare,
    METH_VARARGS,
    "factor the sparse matrix"},
  {"sparseFactorSolve",
   superluWrappersSparseFactorSolve,
   METH_VARARGS,
   "solve the sparse linear system that was previously factorized"},
  {NULL}
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initsuperluWrappers(void) 
{
    PyObject* m,*d;
    if (PyType_Ready(&SparseMatrixType) < 0)
      return;
    if (PyType_Ready(&SparseFactorType) < 0)
      return;
    m = Py_InitModule3("superluWrappers", 
                       superluWrappersMethods,
                       "superlu wrappers module");
    d = PyModule_GetDict(m);
    import_array();
    Py_INCREF(&SparseMatrixType);
    PyModule_AddObject(m, "SparseMatrix", (PyObject *)&SparseMatrixType);
    Py_INCREF(&SparseFactorType);
    PyModule_AddObject(m, "SparseFactor", (PyObject *)&SparseFactorType);
}
/* @} */
