#define FLCBDF_WRAPPERS_MODULE
#include "flcbdfWrappersModule.h"
#include <algorithm>
#include "meshio.h"

//extern "C"
//{
//#include "metis.h"
//}
using namespace Daetk::Petsc::cc;
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define ND(p) ((PyArrayObject *)p)->nd
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define PARVEC_PETSCVEC(p) ((ParVec*)p)->v
#define PARVEC_ARRAY(p) ((ParVec*)p)->array
#define SMP(p) ((SparseMatrix*)p)
#define MESH(p) ((CMesh*)p)->mesh

//--memory profiling
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
  return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
  /* AIX and Solaris ------------------------------------------ */
  struct psinfo psinfo;
  int fd = -1;
  if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
    return (size_t)0L;/* Can't open? */
  if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
      close( fd );
      return (size_t)0L;/* Can't read? */
    }
  close( fd );
  return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
  return (size_t)rusage.ru_maxrss;
#else
  return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
  /* Unknown OS ----------------------------------------------- */
  return (size_t)0L;/* Unsupported. */
#endif
}





/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
  return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
  /* OSX ------------------------------------------------------ */
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		  (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
    return (size_t)0L;/* Can't access? */
  return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE* fp = NULL;
  if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
    return (size_t)0L;/* Can't open? */
  if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
      fclose( fp );
      return (size_t)0L;/* Can't read? */
    }
  fclose( fp );
  return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
  /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
  return (size_t)0L;/* Unsupported. */
#endif
}
//--memory profiling
typedef struct
{
  PyObject_HEAD
  double* array;
  Vec v;
  PyArrayObject* numpy_array;
} ParVec;


extern "C"
{
  
// DAETK global Sys object, needs to be instantiated to instantiate DAETK
static Daetk::Petsc::Sys* DAETK_SYS;

static int
ensure_comm()
{
  if (PROTEUS_COMM_WORLD == MPI_COMM_NULL) {
    PyErr_SetString(PyExc_RuntimeError, "flcbdfWrappersModule is not initialized!");
    return 0;
  }
  return 1;
}

static PyObject*
ParVec_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ParVec *self;
  self = (ParVec *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
ParVec_init(ParVec *self, PyObject *args, PyObject *kwds)
{
  int bs,n,N,nghost,useBlockVec;
  PyObject *subdomain2global,*array;

  if (!ensure_comm()) {
    return -1;
  }

  if(!PyArg_ParseTuple(args,
                       "iiiiOOi",
                       &bs,
                       &n,
                       &N,
                       &nghost,
                       &subdomain2global,
                       &array,
                       &useBlockVec))
    
    return -1;
  self->numpy_array = ((PyArrayObject *)array);
  Py_INCREF(self->numpy_array);
  self->array = DDATA(array);
  if (nghost >= 0)
    {
      std::valarray<int> ghosts(bs*nghost);
      for (int i=0;i<nghost;i++)
        {
          for (int j=0;j<bs;j++)
            {
              ghosts[i*bs+j] = IDATA(subdomain2global)[n+i]*bs+j;
            }
        }
      if (bs==1)
        {
          VecCreateGhostWithArray(PROTEUS_COMM_WORLD,n,N,nghost,&ghosts[0],self->array,&self->v);
        }
      else
        {
          if (useBlockVec)
            VecCreateGhostBlockWithArray(PROTEUS_COMM_WORLD,bs,bs*n,bs*N,nghost,IDATA(subdomain2global)+n,self->array,&self->v);
          else
            VecCreateGhostWithArray(PROTEUS_COMM_WORLD,bs*n,bs*N,bs*nghost,&ghosts[0],self->array,&self->v);
        }
    }
  else
    {
      if (bs==1)
        {
          VecCreateMPIWithArray(PROTEUS_COMM_WORLD,bs,n,N,self->array,&self->v);
        }
      else
        {
          if (useBlockVec)
            VecCreateMPIWithArray(PROTEUS_COMM_WORLD,bs,n,N,self->array,&self->v);
          else
            VecCreateMPIWithArray(PROTEUS_COMM_WORLD,1,bs*n,bs*N,self->array,&self->v);
        }
      if (useBlockVec)
        VecSetBlockSize(self->v,bs);
    }
  if (bs > 1 && useBlockVec)
    {
      int bstmp;
      VecGetBlockSize(self->v,&bstmp);
      assert(bstmp==bs);
    }
  VecGetArray(self->v,&self->array);
  assert(self->array==DDATA(array));
  return 0;
}

static  void
ParVec_dealloc(ParVec* self)
{
  if (self->numpy_array != NULL)//shouldn't be NULL, here in any case
    {
      VecRestoreArray(self->v,&self->array);
      VecDestroy(&self->v);
      Py_DECREF(self->numpy_array);
    }
}

static PyObject*
ParVec_scatter_forward_insert(ParVec *self, PyObject* args)
{
  VecRestoreArray(self->v,&self->array);
  VecGhostUpdateBegin(self->v,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(self->v,INSERT_VALUES,SCATTER_FORWARD);
  VecGetArray(self->v,&self->array);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ParVec_scatter_reverse_add(ParVec *self, PyObject* args)
{
  VecRestoreArray(self->v,&self->array);
  VecGhostUpdateBegin(self->v,ADD_VALUES,SCATTER_REVERSE);
  VecGhostUpdateEnd(self->v,ADD_VALUES,SCATTER_REVERSE);
  VecGetArray(self->v,&self->array);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ParVec_scatter_forward_insertAll(ParVec *self, PyObject* args)
{
  PyObject *subdomain2global,*arrayAll;
  if(!PyArg_ParseTuple(args,
                       "OO",
                       &subdomain2global,
                       &arrayAll))    
    return NULL;
  VecRestoreArray(self->v,&self->array);
  VecGetValues(self->v,SHAPE(subdomain2global)[0],IDATA(subdomain2global),DDATA(arrayAll));
  VecGetArray(self->v,&self->array);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ParVec_scatter_reverse_addAll(ParVec *self, PyObject* args)
{
  PyObject *subdomain2global,*arrayAll;
  if(!PyArg_ParseTuple(args,
                       "OO",
                       &subdomain2global,
                       &arrayAll))    
    return NULL;
  VecRestoreArray(self->v,&self->array);
  VecSetValues(self->v,SHAPE(subdomain2global)[0],IDATA(subdomain2global),DDATA(arrayAll),ADD_VALUES);
  VecAssemblyBegin(self->v);
  VecAssemblyEnd(self->v);
  VecGetArray(self->v,&self->array);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ParVec_printType(ParVec *self, PyObject* args)
{
  VecType type;
  VecGetType(self->v,const_cast<VecType *> (&type));
  //VecGetType(self->v,(const VecType *) (&type));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef ParVec_methods[] = {
  {"scatter_forward_insert", 
   (PyCFunction)ParVec_scatter_forward_insert, 
   METH_NOARGS, 
   "update the ghost values"},
  {"scatter_reverse_add", 
   (PyCFunction)ParVec_scatter_reverse_add, 
   METH_NOARGS, 
   "sum the local and ghost values"},
  {"scatter_forward_insertAll", 
   (PyCFunction)ParVec_scatter_forward_insertAll, 
   METH_VARARGS, 
   "update the ghost values"},
  {"scatter_reverse_addAll", 
   (PyCFunction)ParVec_scatter_reverse_addAll, 
   METH_VARARGS, 
   "sum the local and ghost values"},
  {"printType", 
   (PyCFunction)ParVec_printType, 
   METH_NOARGS, 
   "print the PETSc Vec type"},
  {NULL}
};

static PyTypeObject ParVecType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "flcbdfWrappers.ParVec",             /*tp_name*/
  sizeof(ParVec), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)ParVec_dealloc,                         /*tp_dealloc*/
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
  "ParVec objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  ParVec_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)ParVec_init,      /* tp_init */
  0,                         /* tp_alloc */
  ParVec_new,                 /* tp_new */
};


}
typedef struct
{
  PyObject_HEAD
  Mat m;//,m2;
  ISLocalToGlobalMapping subdomain2globalIS;
} ParMat;


#define PETSCMAT(p) ((ParMat*)p)->m
#define PETSCMAT2(p) ((ParMat*)p)->m2

extern "C"
{
static PyObject*
ParMat_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ParMat *self;
  self = (ParMat *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
ParMat_init(ParMat *self, PyObject *args, PyObject *kwds)
{
  int bs,n,N,nghost,max_nNeighbors;
  PyObject *subdomain2global,*L;

  if (!ensure_comm()) {
    return -1;
  }

  if(!PyArg_ParseTuple(args,
                       "iiiiiOO",
                       &bs,
                       &n,
                       &N,
                       &nghost,
                       &max_nNeighbors,
                       &subdomain2global,
                       &L))
    
    return -1;
  //set up a local 2 global mapping
  std::vector<int> indices(bs*SHAPE(subdomain2global)[0]);
  for (int i=0;i<SHAPE(subdomain2global)[0];i++)
    for (int bI=0;bI<bs;bI++)
      {
        indices[i*bs+bI] = IDATA(subdomain2global)[i]*bs+bI;
      }
  if (bs==1)
    {
      //        MatCreateMPIAIJ(PROTEUS_COMM_WORLD,n,n,N,N,1,PETSC_NULL,max_nNeighbors,PETSC_NULL,&self->m2);

      MatCreate(PROTEUS_COMM_WORLD,&self->m);
      MatSetSizes(self->m,n,n,N,N);
      MatSetFromOptions(self->m);
      //try putting indeces in global numbering
      std::vector<int> j(SMP(L)->A.rowptr[n]);
      for (int i=0;i<n;i++)
        for (int k=SMP(L)->A.rowptr[i];k<SMP(L)->A.rowptr[i+1];k++)
          {
            j[k] = indices[SMP(L)->A.colind[k]];
          }
      MatSeqAIJSetPreallocationCSR(self->m,SMP(L)->A.rowptr,&j[0],(double*)(SMP(L)->A.nzval));
      MatMPIAIJSetPreallocationCSR(self->m,SMP(L)->A.rowptr,&j[0],(double*)(SMP(L)->A.nzval));
    }
  else
    {
      //      MatCreateMPIAIJ(PROTEUS_COMM_WORLD,bs*n,bs*n,bs*N,bs*N,1,PETSC_NULL,bs*max_nNeighbors,PETSC_NULL,&self->m2);

      MatCreate(PROTEUS_COMM_WORLD,&self->m);
      MatSetSizes(self->m,bs*n,bs*n,bs*N,bs*N);
      MatSetFromOptions(self->m);
      std::vector<int> j(SMP(L)->A.rowptr[bs*n]);
      for (int i=0;i<bs*n;i++)
        for (int k=SMP(L)->A.rowptr[i];k<SMP(L)->A.rowptr[i+1];k++)
          {
            j[k] = indices[SMP(L)->A.colind[k]];
          }
      MatSeqAIJSetPreallocationCSR(self->m,SMP(L)->A.rowptr,&j[0],(double*)(SMP(L)->A.nzval));
      MatMPIAIJSetPreallocationCSR(self->m,SMP(L)->A.rowptr,&j[0],(double*)(SMP(L)->A.nzval));
    }
  //cek hack
  //PetscOptionsPrint(stdout);
  ISLocalToGlobalMappingCreate(PROTEUS_COMM_WORLD,1, bs*SHAPE(subdomain2global)[0],&indices[0],PETSC_COPY_VALUES,&self->subdomain2globalIS);
  MatSetLocalToGlobalMapping(self->m,self->subdomain2globalIS,self->subdomain2globalIS);
  return 0;
}

static  void
ParMat_dealloc(ParMat* self)
{
  MatDestroy(&self->m);
  //  MatDestroy(self->m2);
  ISLocalToGlobalMappingDestroy(&self->subdomain2globalIS);
}

static PyTypeObject ParMatType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "flcbdfWrappers.ParMat",             /*tp_name*/
  sizeof(ParMat), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)ParMat_dealloc,                         /*tp_dealloc*/
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
  "ParMat objects",           /* tp_doc */
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
  (initproc)ParMat_init,      /* tp_init */
  0,                         /* tp_alloc */
  ParMat_new,                 /* tp_new */
};
}//extern c

//--- make petsc linear solvers converge based on true residual---
//storage for work vectors
typedef struct
{
  int max_it;
  double atol;
  double rtol;
  double divtol;
  double ttol;
  double rnorm0;
  Vec t;
  Vec v;
  Vec residual;
} TrueResidualTestCtx;
#define TRUERESCTX(p) ((TrueResidualTestCtx*)p)

static 
int destroy_TrueResidualTestCtx(void* ctx)
{
  int ierr = 0;
  if (TRUERESCTX(ctx)->t)
    ierr = VecDestroy(&TRUERESCTX(ctx)->t);
  if (ierr) return ierr;
  if (TRUERESCTX(ctx)->v)
    ierr = VecDestroy(&TRUERESCTX(ctx)->v);
  if (ierr) return ierr;
  if (TRUERESCTX(ctx)->residual)
    ierr = VecDestroy(&TRUERESCTX(ctx)->residual);
  return ierr;
}

static 
int converged_true_residual(KSP ksp, int it, double rnorm, KSPConvergedReason *reason, void* ctx)
{
  int max_it = TRUERESCTX(ctx)->max_it;
  double truenorm=0.0,atol = 1.0e-10,rtol = 0.0,divtol = 1.0e-10;

  atol = TRUERESCTX(ctx)->atol;  rtol = TRUERESCTX(ctx)->rtol;  divtol = TRUERESCTX(ctx)->divtol;
  //default not converged
  *reason = KSP_CONVERGED_ITERATING;

  KSPBuildResidual(ksp,TRUERESCTX(ctx)->t,TRUERESCTX(ctx)->v,&TRUERESCTX(ctx)->residual);
  VecNorm(TRUERESCTX(ctx)->residual,NORM_2,&truenorm);
  if (it == 0)
    {
      TRUERESCTX(ctx)->rnorm0=truenorm;
      TRUERESCTX(ctx)->ttol = fmax(rtol*truenorm,atol);
    }
  if (truenorm <= TRUERESCTX(ctx)->ttol)
    {
      if (truenorm <= atol)
	*reason = KSP_CONVERGED_ATOL;
      else
	*reason = KSP_CONVERGED_RTOL;
    }
  else if (truenorm >= divtol*TRUERESCTX(ctx)->rnorm0)
    {
      *reason = KSP_DIVERGED_DTOL;
    }
  else if (it > max_it)
    {
      *reason = KSP_DIVERGED_ITS;
    }
  return 0;
}

typedef struct
{
  PyObject_HEAD
  KSP ksp; 
  TrueResidualTestCtx tres_ctx;
} CKSP;

#define KSP(p) ((CKSP*)p)->ksp

extern "C"
{
static PyObject*
CKSP_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  CKSP *self;
  self = (CKSP *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
CKSP_init(CKSP *self, PyObject *args, PyObject *kwds)
{
  PyObject *par_L;
  char* prefix(0);

  if (!ensure_comm()) {
    return -1;
  }

  if(!PyArg_ParseTuple(args,
                       "O|s",
                       &par_L,
                       &prefix))
    
    return -1;

  KSPCreate(PROTEUS_COMM_WORLD,&self->ksp);
  if (prefix)
    KSPSetOptionsPrefix(self->ksp,prefix);
  KSPSetFromOptions(self->ksp);
  PC pc;
  KSPGetPC(self->ksp,&pc);
  PCSetFromOptions(pc);
  /*mwf no array defined in ksp so comment out?
    assert(self->array==DDATA(array));
  */
  /*now include manual true residual test*/
  self->tres_ctx.t = PETSC_NULL;
  self->tres_ctx.v = PETSC_NULL;
  self->tres_ctx.residual = PETSC_NULL;
  self->tres_ctx.atol = 1.0e-10;
  self->tres_ctx.rtol = 0.0;
  self->tres_ctx.divtol= 1.0e50;
  self->tres_ctx.max_it= 10000;


  PetscReal  *res;
  int its;
  
  its = 250; // Ido HACK ...
  PetscMalloc(its*sizeof(PetscReal),&res); 

  KSPSetResidualHistory(self->ksp, res, its, PETSC_TRUE);

  return 0;
}

static PyObject*
CKSP_useTrueResidualConvergence(CKSP *self, PyObject *args, PyObject *kwds)
{
  PyObject *par_v;
  double atol,rtol,divtol;
  int max_it;
  PetscBool found;
  if(!PyArg_ParseTuple(args,
                       "O",
                       &par_v))
    
    return NULL;
  //now set residual monitor? try to get sizing info
  if (!self->tres_ctx.v)
    VecDuplicate(PARVEC_PETSCVEC(par_v),&self->tres_ctx.v);
  if (!self->tres_ctx.t)
    VecDuplicate(PARVEC_PETSCVEC(par_v),&self->tres_ctx.t);
  if (!self->tres_ctx.residual)
    VecDuplicate(PARVEC_PETSCVEC(par_v),&self->tres_ctx.residual);

  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-ksp_atol",&atol,&found);
  if (found) self->tres_ctx.atol  = atol;
  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-ksp_rtol",&rtol,&found);
  if (found) self->tres_ctx.rtol  = rtol;
  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-ksp_divtol",&divtol,&found);
  if (found) self->tres_ctx.divtol  = divtol;
  PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-ksp_max_it",&max_it,&found);
  if (found) self->tres_ctx.max_it= max_it;

  KSPSetConvergenceTest(self->ksp,converged_true_residual,(void *)&self->tres_ctx,PETSC_NULL);
  Py_INCREF(Py_None); 
  return Py_None;

}


static  void
CKSP_dealloc(CKSP* self)
{
  KSPDestroy(&self->ksp);
  destroy_TrueResidualTestCtx(&(self->tres_ctx));
}

static PyObject*
CKSP_prepare(CKSP *self, PyObject* args)
{
  PyObject *L,*par_L;
  int overlap(1);
  if(!PyArg_ParseTuple(args,
                       "OO|i",
                       &L,
                       &par_L,
                       &overlap))
    return NULL;
  //to do, get rid of copy and load into petsc storage directly
  //copy L into par_L
  //assume compressed row for now
  int *rowptr,*colind;
  double *a;
  a = (double*) (SMP(L)->A.nzval);
  rowptr = SMP(L)->A.rowptr;
  colind = SMP(L)->A.colind;
  MatZeroEntries(PETSCMAT(par_L));
  int irow[1];
  int offset_rank,offset_rankP1,ldim;
  MatGetOwnershipRange(PETSCMAT(par_L),&offset_rank,&offset_rankP1);
  ldim = offset_rankP1-offset_rank;
  if (overlap <= 0)
    {
      for (int i=0;i<SMP(L)->dim[0];i++)
        {
          irow[0] = i;
          MatSetValuesLocal(PETSCMAT(par_L),1,irow,rowptr[i+1]-rowptr[i],&colind[rowptr[i]],&a[rowptr[i]],ADD_VALUES);
        }
    }
  else
    {
      for (int i=0;i<ldim;i++)
        {
          irow[0] = i;
          MatSetValuesLocal(PETSCMAT(par_L),1,irow,rowptr[i+1]-rowptr[i],&colind[rowptr[i]],&a[rowptr[i]],INSERT_VALUES);
        }
    }
  MatAssemblyBegin(PETSCMAT(par_L),MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(PETSCMAT(par_L),MAT_FINAL_ASSEMBLY);
  KSPSetOperators(self->ksp,PETSCMAT(par_L),PETSCMAT(par_L));
  KSPSetUp(self->ksp);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
CKSP_solve(CKSP *self, PyObject* args)
{
  PyObject *par_u,*par_b;
  if(!PyArg_ParseTuple(args,
                       "OO",
                       &par_u,
                       &par_b))
    return NULL;
  VecRestoreArray(PARVEC_PETSCVEC(par_b),&PARVEC_ARRAY(par_b));
  VecRestoreArray(PARVEC_PETSCVEC(par_u),&PARVEC_ARRAY(par_u));
  //VecView(PARVEC_PETSCVEC(par_b),PETSC_VIEWER_STDOUT_WORLD);
  //VecView(PARVEC_PETSCVEC(par_u),PETSC_VIEWER_STDOUT_WORLD);
  KSPSolve(self->ksp,PARVEC_PETSCVEC(par_b),PARVEC_PETSCVEC(par_u));
  VecGetArray(PARVEC_PETSCVEC(par_b),&PARVEC_ARRAY(par_b));
  VecGetArray(PARVEC_PETSCVEC(par_u),&PARVEC_ARRAY(par_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
CKSP_info(CKSP *self, PyObject* args)
{
  PetscInt     its;
  PetscReal*   res;
  double       first,last,rel_res; 

  if (!ensure_comm()) {
    return NULL;
  }

  KSPGetResidualHistory(self->ksp, &res, &its);

  first = *res;
  for (int i=0; i<its; ++i)
    {
      last = *res;
      res++;
    } 

  rel_res = 100.0*last/first;
  PetscPrintf(PROTEUS_COMM_WORLD,"\n       Iterations: %D    Error reduction: %g (%%)\n\n", its,rel_res);
 
  
  Py_INCREF(Py_None); 
  return Py_None;
}



 
static PyMethodDef CKSP_methods[] = {
  {"prepare", 
   (PyCFunction)CKSP_prepare, 
   METH_VARARGS, 
   "prepare the linear solver and preconditioner"},
  {"solve", 
   (PyCFunction)CKSP_solve, 
   METH_VARARGS, 
   "solve the problem"},
  {"info", 
   (PyCFunction)CKSP_info, 
   METH_VARARGS, 
   "print solver performance info"},
  {"useTrueResidualConvergence", 
   (PyCFunction)CKSP_useTrueResidualConvergence, 
   METH_VARARGS, 
   "use true residual convergence test"},
  {NULL,NULL}
};

static PyTypeObject CKSPType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "flcbdfWrappers.CKSP",             /*tp_name*/
  sizeof(CKSP), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)CKSP_dealloc,                         /*tp_dealloc*/
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
  "CKSP objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  CKSP_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)CKSP_init,      /* tp_init */
  0,                         /* tp_alloc */
  CKSP_new,                 /* tp_new */
};

}//extern C

/** \file flcbdfWrappersModule.cpp
    \ingroup flcbdfWrappers
    @{
*/
extern "C"
{


static PyObject* 
FLCBDF_integrator_choose_dt(FLCBDF_integrator *self, 
                            PyObject *args)
{
  double t,tout,DT;
  if(!PyArg_ParseTuple(args,
                       "dd",
                       &t,
                       &tout))
    return NULL;
  DT = self->flcbdf->chooseDT(t,tout);
  return Py_BuildValue("d",DT);
}

static PyObject* 
FLCBDF_integrator_set_dt(FLCBDF_integrator *self, 
                         PyObject *args)
{
  double DT;
  if(!PyArg_ParseTuple(args,
                       "d",
                       &DT))
    return NULL;
  self->flcbdf->setDT(DT);
  return Py_BuildValue("d",DT);
}
static PyObject* 
FLCBDF_integrator_set_order(FLCBDF_integrator *self, 
                            PyObject *args)
{
  int k;
  if(!PyArg_ParseTuple(args,
                       "i",
                       &k))
    return NULL;
  self->flcbdf->useFixedOrder(k);
  return Py_BuildValue("i",k);
}

static PyObject* 
FLCBDF_integrator_initialize_dt(FLCBDF_integrator *self, 
                                PyObject *args)
{
  double t0,tOut,DT;
  PyObject *y,*yPrime;
  if(!PyArg_ParseTuple(args,
                       "ddOO",
                       &t0,
                       &tOut,
                       &y,
                       &yPrime))
    return NULL;
  Daetk::Petsc::Vec yVec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_),
    yPrimeVec(Daetk::Vec::REF,DDATA(yPrime),self->flcbdf->yn.ldim_);
  DT = self->flcbdf->chooseInitialStepSize(t0,tOut,yVec,yPrimeVec);
  return Py_BuildValue("d",DT);
}

static PyObject* 
FLCBDF_integrator_setInitialGuess(FLCBDF_integrator *self, 
                                  PyObject *args)
{
  PyObject *y;
  if(!PyArg_ParseTuple(args,
                       "O",
                       &y))
    return NULL;
  Daetk::Petsc::Vec yVec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_);
  yVec = self->flcbdf->yn;
  return Py_None;
}

static PyObject* 
FLCBDF_integrator_lastStepErrorOk(FLCBDF_integrator *self, 
                                  PyObject *args)
{
  bool Ok;
  PyObject *y;
  if(!PyArg_ParseTuple(args,
                       "O",
                       &y))
    return NULL;
  Daetk::Petsc::Vec yVec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_);
  Ok = !self->flcbdf->errorForStepTooLarge(yVec);
  return Py_BuildValue("i",Ok);
}

static PyObject* 
FLCBDF_integrator_calculate_yprime(FLCBDF_integrator *self, 
                                   PyObject *args)
{
  PyObject *y,*Dy,*yprime,*Dyprime;
  if(!PyArg_ParseTuple(args,
                       "OOOO",
                       &y,
                       &Dy,
                       &yprime,
                       &Dyprime))
    return NULL;
  if (self->yVec==0)
    {
      self->yVec = new Daetk::Petsc::Vec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_);
      self->DyVec = new Daetk::Petsc::Vec(Daetk::Vec::REF,DDATA(Dy),self->flcbdf->yn.ldim_);
      self->yprimeVec = new Daetk::Petsc::Vec(Daetk::Vec::REF,DDATA(yprime),self->flcbdf->yn.ldim_);
      self->DyprimeVec = new Daetk::Petsc::Vec(Daetk::Vec::REF,DDATA(Dyprime),self->flcbdf->yn.ldim_);
    }
  self->flcbdf->calculate_yprime(*self->yVec,*self->DyVec,*self->yprimeVec,*self->DyprimeVec);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
FLCBDF_integrator_stepTaken(FLCBDF_integrator *self, 
                            PyObject *args)
{
  PyObject *y;
  if(!PyArg_ParseTuple(args,
                       "O",
                       &y))
    return NULL;
  Daetk::Petsc::Vec yVec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_);
  self->flcbdf->estimateError(yVec);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
FLCBDF_integrator_retryStep_errorFailure(FLCBDF_integrator *self, 
                                         PyObject *args)
{
  double h;
  h = self->flcbdf->retryStep_errorFailure();
  Py_INCREF(Py_None);
  return Py_BuildValue("d",h);
}

static PyObject* 
FLCBDF_integrator_retryStep_solverFailure(FLCBDF_integrator *self, 
                                          PyObject *args)
{
  double h;
  h = self->flcbdf->retryStep_solverFailure();
  Py_INCREF(Py_None);
  return Py_BuildValue("d",h);
}

static PyObject* 
FLCBDF_integrator_initializeTimeHistory(FLCBDF_integrator *self, 
                                        PyObject *args)
{
  PyObject *y,*yPrime;
  if(!PyArg_ParseTuple(args,
                       "OO",
                       &y,
                       &yPrime))
    return NULL;
  Daetk::Petsc::Vec yVec(Daetk::Vec::REF,DDATA(y),self->flcbdf->yn.ldim_),
    yPrimeVec(Daetk::Vec::REF,DDATA(yPrime),self->flcbdf->yn.ldim_);
  self->flcbdf->initializeTimeHistory(yVec,yPrimeVec);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
FLCBDF_integrator_setTolerances(FLCBDF_integrator *self, 
                                PyObject *args)
{
  double atol,rtol;
  PyObject* dV;
  if(!PyArg_ParseTuple(args,
                       "ddO",
                       &atol,&rtol,&dV))
    return NULL;
  Daetk::Petsc::Vec dV_Vec(Daetk::Vec::REF,DDATA(dV),self->flcbdf->yn.ldim_);
  self->wNorm->setTolerances(atol,rtol,dV_Vec);
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject* 
FLCBDF_integrator_getCurrentAlpha(FLCBDF_integrator *self, 
                                  PyObject *args)
{
  double alpha;
  assert(self->flcbdf);
  alpha = self->flcbdf->getCurrentAlpha();
  //do I need Py_INCREF(Py_None);
  return Py_BuildValue("d",alpha);

}

static PyMethodDef FLCBDF_integrator_methods[] = {
  {"lastStepErrorOk",
   (PyCFunction)FLCBDF_integrator_lastStepErrorOk,
   METH_VARARGS,
   "check the temporal error estimate"},
  {"choose_dt",
   (PyCFunction)FLCBDF_integrator_choose_dt,
   METH_VARARGS,
   "choose the next time step"},
  {"set_dt",
   (PyCFunction)FLCBDF_integrator_set_dt,
   METH_VARARGS,
   "set the next time step"},
  {"set_order",
   (PyCFunction)FLCBDF_integrator_set_order,
   METH_VARARGS,
   "set the order of flcbdf"},
  {"initialize_dt",
   (PyCFunction)FLCBDF_integrator_initialize_dt,
   METH_VARARGS,
   "choose the initial time step"},
  {"setInitialGuess",
   (PyCFunction)FLCBDF_integrator_setInitialGuess,
   METH_VARARGS,
   "choose the initial guess for nonlinear solve"},
  {"calculate_yprime",
   (PyCFunction)FLCBDF_integrator_calculate_yprime,
   METH_VARARGS,
   "calculate the discrete time derivative"},
  {"stepTaken", 
   (PyCFunction)FLCBDF_integrator_stepTaken, 
   METH_VARARGS, 
   "accept step and calculate error estimates"},
  {"retryStep_errorFailure", 
   (PyCFunction)FLCBDF_integrator_retryStep_errorFailure, 
   METH_VARARGS, 
   "retry step when temporal error is too large"},
  {"retryStep_solverFailure", 
   (PyCFunction)FLCBDF_integrator_retryStep_solverFailure, 
   METH_VARARGS, 
   "retry step when nonlinear solver fails"},
  {"initializeTimeHistory", 
   (PyCFunction)FLCBDF_integrator_initializeTimeHistory, 
   METH_VARARGS, 
   "accept step and calculate error estimates"},
  {"setTolerances",
   (PyCFunction)FLCBDF_integrator_setTolerances,
   METH_VARARGS,
   "set (scalar) absolute and relative integration tolerances"},
  {"getCurrentAlpha", 
   (PyCFunction)FLCBDF_integrator_getCurrentAlpha, 
   METH_VARARGS, 
   "get alpha from FLCBDF approximation for jacobians (off-diagonal)"},
  {NULL,NULL}
};

static PyObject*
FLCBDF_integrator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FLCBDF_integrator *self;
  self = (FLCBDF_integrator *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
FLCBDF_integrator_init(FLCBDF_integrator *self, PyObject *args, PyObject *kwds)
{
  int i,ldim=1,dim;
  const char* yName;
  char dataFilename[80];
  PyObject* y;
  if(!PyArg_ParseTuple(args,
                       "Os",
                       &y,
                       &yName))
    return -1;
  for (i=0;i<ND(y);i++)
    ldim*=SHAPE(y)[i];
  self->petscSys = new Daetk::Petsc::Sys();
  self->sizeVec = new Daetk::Vec(Daetk::Vec::REF,DDATA(y),ldim);//this is so we have a parallel vector in the registry
  dim=self->sizeVec->dim();
  if (self->petscSys->getSize() > 1)
    self->sizeVec->setExample();
  self->wNorm = new Daetk::WeightedL2Norm(dim);
  sprintf(dataFilename,"data_%i_%s.txt",dim,yName);
  self->data = new Daetk::FullDataFile(0.0,dataFilename);
  //mwf if wanted to only write from proc 0 
  //if (self->petscSys->getRank() == 0)
  //self->data = new Daetk::FullDataFile(0.0,dataFilename);
  //else
  // self->data = new Daetk::DataCollector();
  self->flcbdf = new Daetk::FLCBDF_lite(dim,*self->wNorm,*self->data);
  self->yVec=0;
  self->DyVec=0;
  self->yprimeVec=0;
  self->DyprimeVec=0;
  return 0;
}

static  void
FLCBDF_integrator_dealloc(FLCBDF_integrator* self)
{
  delete self->flcbdf;
  delete self->wNorm;
  delete self->data;
  delete self->petscSys;
  delete self->sizeVec;
  if (self->yVec!=0)
    {
      delete self->yVec;
      delete self->DyVec;
      delete self->yprimeVec;
      delete self->DyprimeVec;
    }
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject FLCBDF_integratorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "flcbdfWrappers.FLCBDF_integrator",             /*tp_name*/
  sizeof(FLCBDF_integrator), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)FLCBDF_integrator_dealloc,                         /*tp_dealloc*/
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
  "FLCBDF integtrator objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  FLCBDF_integrator_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)FLCBDF_integrator_init,      /* tp_init */
  0,                         /* tp_alloc */
  FLCBDF_integrator_new,                 /* tp_new */
};

static PyObject* 
DaetkPetscSys_barrier(DaetkPetscSys *self, 
                      PyObject *args)
{
  self->petscSys->barrier();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
DaetkPetscSys_isMaster(DaetkPetscSys *self, 
                       PyObject *args)
{
  return Py_BuildValue("b",self->petscSys->master());
}

static PyObject* 
DaetkPetscSys_isInitialized(DaetkPetscSys *self, 
                            PyObject *args)
{
  return Py_BuildValue("b",self->petscSys->isInitialized());
}

static PyObject* 
DaetkPetscSys_beginSequential(DaetkPetscSys *self, 
                              PyObject *args)
{
  self->petscSys->beginSequential();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
DaetkPetscSys_endSequential(DaetkPetscSys *self, 
                            PyObject *args)
{
  self->petscSys->endSequential();
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
DaetkPetscSys_catchError(DaetkPetscSys *self, 
                         PyObject *args)
{
  bool error;
  if(!PyArg_ParseTuple(args,
                       "b",
                       &error))
    return NULL;
  return Py_BuildValue("b",self->petscSys->catchError(error));
}

static PyObject* 
DaetkPetscSys_rank(DaetkPetscSys *self, 
                   PyObject *args)
{
  return Py_BuildValue("i",self->petscSys->getRank());
}

static PyObject* 
DaetkPetscSys_size(DaetkPetscSys *self, 
                   PyObject *args)
{
  return Py_BuildValue("i",self->petscSys->getSize());
}

//todo add overlap for element based partitions
int partitionElementsOriginal(Mesh& mesh, int nElements_overlap)
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return NULL;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  //Contents 
  //
  //1. Partition the elements in the "default" partition (contiguous
  //chunks in given ordering) 
  //
  //2. Partition the elementNeighbors based on this partition
  //
  //3. Pass to Parmetis to build a better partition of the elements
  //
  //4. Tag a subset of the nodes on the subdomain elements as owned
  //using a mark and pass approach.** 
  //
  //5. Extract the nodes in the
  //overlapping elements.** 
  //
  //6. Build the subdomain mesh from the
  //subdomain elements 
  //
  //**To be more general we could get all the support (i.e. faces
  //and edges) and partitiong them, but the main reason for
  //partitioning is to keep track of a global numbering for degrees
  //of freedom that live on each type of geometric entity. We only
  //have node and element based DOFs so I just rebuild the other
  //information once we have elements and nodes partitioned.
  //
  // \todo check that I restore all data that PETSc expects to have
  // back, add PETSc error checking macros
  //
  //1. Build default partitioning
  //
  //get offsets so we can calculate the processor to global mapping
  //for elements in the old (default) partitioning
  valarray<int> elementOffsets_old(size+1);
  elementOffsets_old[0] = 0;
  for(int sdN=0;sdN<size;sdN++)
    elementOffsets_old[sdN+1] = elementOffsets_old[sdN] + 
      int(mesh.nElements_global)/size + 
      (int(mesh.nElements_global)%size > sdN);
    
  //2. Extract subdomain element adjacency information could read
  //only the required portion from a file
  int nElements_subdomain = (elementOffsets_old[rank+1] - elementOffsets_old[rank]);
  PetscInt *elementNeighborsOffsets_subdomain,*elementNeighbors_subdomain,*weights_subdomain;
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain+1),&elementNeighborsOffsets_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain*mesh.nElementBoundaries_element),&elementNeighbors_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain*mesh.nElementBoundaries_element),&weights_subdomain);
  //this wastes a little space
  elementNeighborsOffsets_subdomain[0] = 0;
  for (int eN=0,offset=0; eN < nElements_subdomain; eN++)
    {
      int eN_global = elementOffsets_old[rank] + eN;
      int offsetStart=offset;
      for (int ebN=0; ebN< mesh.nElementBoundaries_element;ebN++)
        {
          int eN_neighbor_global = mesh.elementNeighborsArray[eN_global*mesh.nElementBoundaries_element + ebN];
          if (eN_neighbor_global >= 0 )
            elementNeighbors_subdomain[offset++] = eN_neighbor_global;
        }
      elementNeighborsOffsets_subdomain[eN+1]=offset;
      sort(&elementNeighbors_subdomain[offsetStart],&elementNeighbors_subdomain[offset]);
      int weight = (elementNeighborsOffsets_subdomain[eN+1] - elementNeighborsOffsets_subdomain[eN]);
      for (int k=elementNeighborsOffsets_subdomain[eN];k < elementNeighborsOffsets_subdomain[eN+1];k++)
        weights_subdomain[k] = weight;
      // for (int o = offsetStart; o < offset; o++)
      //   {
      //     std::cout<<elementNeighbors_subdomain[o]<<'\t';
      //   }
      // std::cout<<std::endl;
    }
  //3. Generate the  new partitiong using PETSc, this is done in parallel using parmetis
  Mat petscAdjacency;
  //     MatCreateMPIAdj(PROTEUS_COMM_WORLD,
  //                     nElements_subdomain, mesh.nElements_global, 
  //                     &elementNeighborsOffsets_subdomain[0], &elementNeighbors_subdomain[0], 
  //                     &weights_subdomain[0], 
  //                     &petscAdjacency);
  ierr = MatCreateMPIAdj(PROTEUS_COMM_WORLD,
                         nElements_subdomain, 
                         mesh.nElements_global,
                         elementNeighborsOffsets_subdomain, 
                         elementNeighbors_subdomain,
                         PETSC_NULL,//weights_subdomain,
                         &petscAdjacency);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  MatPartitioning petscPartition;
  MatPartitioningCreate(PROTEUS_COMM_WORLD,&petscPartition);
  MatPartitioningSetAdjacency(petscPartition,petscAdjacency);
  MatPartitioningSetFromOptions(petscPartition);
  
  //get a petsc index set that has the new submdomain number for each element
  IS elementPartitioningIS_new;
  MatPartitioningApply(petscPartition,&elementPartitioningIS_new); 
  MatPartitioningDestroy(&petscPartition);
  //MatDestroy(petscAdjacency);
  //ISView(elementPartitioningIS_new,PETSC_VIEWER_STDOUT_WORLD);
    
  //experiment with metis
  //mwf set some defaults and not call if size == 1 since metis crashes
  //cek commenting out for now
  //int etype=1,edgecut=0,base=0;
  //epart assign everything to processor zero by default
  //valarray<int> epart(0,mesh.nElements_global),npart(mesh.nNodes_global);
  //if (size > 1)
  //  METIS_PartMeshNodal(&mesh.nElements_global,&mesh.nNodes_global,mesh.elementNodesArray,&etype,&base,&size,&edgecut,&epart[0],&npart[0]);
  //ISCreateGeneralWithArray(PETSC_COMM_SELF,mesh.nElements_global,&epart[0],&elementPartitioningIS_new);
  //write mesh to view with showme
  //     std::ofstream nodeout("mesh.node"),eleout("mesh.ele"),partout("mesh.part");
  //     eleout<<mesh.nElements_global<<" 3 0"<<std::endl;
  //     partout<<mesh.nElements_global<<"\t"<<size<<std::endl;
  //     for (int eN=0;eN<mesh.nElements_global;eN++)
  //       {
  // 	partout<<(eN+1)<<"\t"<<(1+epart[eN])<<std::endl;
  // 	eleout<<(eN+1)<<"\t"<<(1+mesh.elementNodesArray[eN*3+0])
  // 	      <<"\t"<<(1+mesh.elementNodesArray[eN*3+1])
  // 	      <<"\t"<<(1+mesh.elementNodesArray[eN*3+2])
  // 	      <<std::endl;
  //       }
  //     nodeout<<mesh.nNodes_global<<" 2 0 0"<<std::endl;
  //     for (int nN=0;nN<mesh.nNodes_global;nN++)
  //       {
  // 	nodeout<<(nN+1)<<"\t"<<mesh.nodeArray[nN*3+0]<<"\t"<<mesh.nodeArray[nN*3+1]<<std::endl;
  //       }
  //     eleout.close();
  //     partout.close();
  //count the new number of elements on each subdomain
  valarray<int> nElements_subdomain_new(size);
  ISPartitioningCount(elementPartitioningIS_new,size,&nElements_subdomain_new[0]);

  //get the new offsets for the subdomain to global numbering
  valarray<int> elementOffsets_new(size+1);
  elementOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    elementOffsets_new[sdN+1] = elementOffsets_new[sdN] + nElements_subdomain_new[sdN];

  //get the new element numbers for the elements on this  subdomain
  IS elementNumberingIS_subdomain_old2new;
  ISPartitioningToNumbering(elementPartitioningIS_new,&elementNumberingIS_subdomain_old2new);

  //now get the new element numbers for the whole mesh so that we
  //can just read this processors elements, reorder, and renumber**
  //
  //**We could do this in parallel by scattering all the element
  //information
  IS elementNumberingIS_global_old2new;
  ISAllGather(elementNumberingIS_subdomain_old2new,&elementNumberingIS_global_old2new);
  //ISView(elementNumberingIS_global_old2new,PETSC_VIEWER_STDOUT_SELF);
  const PetscInt *elementNumbering_global_old2new;
  ISGetIndices(elementNumberingIS_global_old2new,&elementNumbering_global_old2new);
  valarray<int> elementNumbering_global_new2old(mesh.nElements_global);
  for(int eN=0;eN<mesh.nElements_global;eN++)
    elementNumbering_global_new2old[elementNumbering_global_old2new[eN]] = eN;

  //Sort element based arrays, maybe I don't need to do this, maybe
  //I just need to start writing into the subdomain mesh here and
  //preserve subdomain2old and subdomain2global mappings
  valarray<int> elementNodesArray_new(mesh.nElements_global*mesh.nNodes_element),
    elementNeighborsArray_new(mesh.nElements_global*mesh.nElementBoundaries_element),
    elementMaterialTypes_new(mesh.nElements_global),
    elementBoundaryElementsArray_new(mesh.nElementBoundaries_global*2),
    elementBoundaryNodesArray_new(mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary),
    edgeNodesArray_new(mesh.nEdges_global*2);
  valarray<int> elementBoundaryMaterialTypes_new(mesh.nElementBoundaries_global);
  for (int eN=0;eN<mesh.nElements_global;eN++)
    {
      for (int nN=0;nN<mesh.nNodes_element;nN++)
        elementNodesArray_new[eN*mesh.nNodes_element + nN] = 
          mesh.elementNodesArray[elementNumbering_global_new2old[eN]*mesh.nNodes_element+nN];
      for (int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
        elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN] = 
          mesh.elementNeighborsArray[elementNumbering_global_new2old[eN]*mesh.nElementBoundaries_element+ebN];
      elementMaterialTypes_new[eN] = mesh.elementMaterialTypes[elementNumbering_global_new2old[eN]];
    }
  //renumber references to element numbers
  for (int eN=0;eN<mesh.nElements_global;eN++)
    {
      for (int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
        {
          int eN_ebN = elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN];
          if (eN_ebN >= 0)
            elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN] = 
              elementNumbering_global_old2new[eN_ebN];
        }
    }
  for(int ebN=0;ebN<mesh.nElementBoundaries_global;ebN++)
    {     
      int eN_L_old = mesh.elementBoundaryElementsArray[ebN*2+0],
        eN_R_old = mesh.elementBoundaryElementsArray[ebN*2+1];
      elementBoundaryElementsArray_new[ebN*2+0] = elementNumbering_global_old2new[eN_L_old];
      if(eN_R_old >= 0)
        elementBoundaryElementsArray_new[ebN*2+1] = elementNumbering_global_old2new[eN_R_old];
      //mwf assume same numbering scheme for element boundaries for now?
      elementBoundaryMaterialTypes_new[ebN] = mesh.elementBoundaryMaterialTypes[ebN];
    }
  //4. now we need to build a new node ordering with better data locality for C0 finite elements
  //otherwise we could just grab the nodes on the subdomain and not worry about ownership
  //in the long run it wouldn't be bad to do a global repartition of faces and edges for mixed hybrid
  //and non-conforming finite elements
  MPI_Status status;
  PetscBT nodeMask; 
  PetscBTCreate(mesh.nNodes_global,&nodeMask);
  if (rank > 0) 
    {
      MPI_Recv(nodeMask,PetscBTLength(mesh.nNodes_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status);
    }
  //mark the unmarked nodes on this subdomain and store the node numbers
  set<int> nodes_subdomain_owned;
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    for(int nN=0;nN<mesh.nNodes_element;nN++)
      {
        int nN_global = elementNodesArray_new[eN*mesh.nNodes_element+nN];
        if (!PetscBTLookupSet(nodeMask,nN_global))
          nodes_subdomain_owned.insert(nN_global);
      }
  //ship off the mask
  if (rank < size-1)
    MPI_Send(nodeMask,PetscBTLength(mesh.nNodes_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&nodeMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy"<<endl;
  //get the number of nodes on each processor
  valarray<int> nNodes_subdomain_new(size),
    nodeOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nNodes_subdomain_new[sdN] = nodes_subdomain_owned.size();
    else
      nNodes_subdomain_new[sdN] = 0;
  valarray<int> nNodes_subdomain_new_send=nNodes_subdomain_new;
  MPI_Allreduce(&nNodes_subdomain_new_send[0],&nNodes_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  nodeOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    nodeOffsets_new[sdN+1] = nodeOffsets_new[sdN]+nNodes_subdomain_new[sdN];
  //Now as with elements build a global node numbering, sort node
  //based information, and renumber references to node numbers
  valarray<int> nodeNumbering_new2old(nodes_subdomain_owned.size());
  set<int>::iterator nN_ownedp=nodes_subdomain_owned.begin();
  for (int nN=0;nN<int(nodes_subdomain_owned.size());nN++)
    {
      nodeNumbering_new2old[nN] = *nN_ownedp++;
    }
  IS nodeNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nodes_subdomain_owned.size(),&nodeNumbering_new2old[0],PETSC_COPY_VALUES,&nodeNumberingIS_new2old);
  IS nodeNumberingIS_global_new2old;
  ISAllGather(nodeNumberingIS_new2old,&nodeNumberingIS_global_new2old);
  const PetscInt *nodeNumbering_global_new2old;
  valarray<int> nodeNumbering_old2new_global(mesh.nNodes_global);
  ISGetIndices(nodeNumberingIS_global_new2old,&nodeNumbering_global_new2old);
  for (int nN=0;nN<mesh.nNodes_global;nN++)
    {
      nodeNumbering_old2new_global[nodeNumbering_global_new2old[nN]] = nN;
    }
  for (int eN=0;eN < mesh.nElements_global; eN++)
    {
      int nN_old;
      for (int nN=0;nN < mesh.nNodes_element; nN++)
        {
          nN_old = elementNodesArray_new[eN*mesh.nNodes_element+nN];
          elementNodesArray_new[eN*mesh.nNodes_element+nN] = nodeNumbering_old2new_global[nN_old];
        }
    }
  for (int i=0;i<mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary;i++)
    {
      int nN_old = mesh.elementBoundaryNodesArray[i];
      elementBoundaryNodesArray_new[i] = nodeNumbering_old2new_global[nN_old];
    }
  for (int i=0;i<mesh.nEdges_global*2;i++)
    {
      int nN_old = mesh.edgeNodesArray[i];
      edgeNodesArray_new[i] = nodeNumbering_old2new_global[nN_old];
    }
  valarray<int> nodeStarArray_new(mesh.nodeStarOffsets[mesh.nNodes_global]);
  for (int i=0;i<mesh.nodeStarOffsets[mesh.nNodes_global];i++)
    {
      int nN_old = mesh.nodeStarArray[i];
      nodeStarArray_new[i] = nodeNumbering_old2new_global[nN_old];
    }
  valarray<double> nodeArray_new(mesh.nNodes_global*3);
  valarray<int>  nodeMaterialTypes_new(mesh.nNodes_global);
  for (int nN=0;nN<mesh.nNodes_global;nN++)
    {
      int nN_new = nodeNumbering_old2new_global[nN];
      nodeArray_new[nN_new*3+0] = mesh.nodeArray[nN*3+0];
      nodeArray_new[nN_new*3+1] = mesh.nodeArray[nN*3+1];
      nodeArray_new[nN_new*3+2] = mesh.nodeArray[nN*3+2];
      nodeMaterialTypes_new[nN_new] = mesh.nodeMaterialTypes[nN];
    }
  //write partitioned mesh to view with "showme"
  //     std::ofstream nodeout("mesh.node"),eleout("mesh.ele"),partout("mesh.part");
  //     eleout<<mesh.nElements_global<<" 3 0"<<std::endl;
  //     partout<<mesh.nElements_global<<"\t"<<size<<std::endl;
  //     for (int eN=0;eN<mesh.nElements_global;eN++)
  //       {
  // 	partout<<(eN+1)<<"\t"<<(1+epart[elementNumbering_global_new2old[eN]])<<std::endl;
  // 	//partout<<(eN+1)<<"\t"<<(1+epart[eN])<<std::endl;
  // 	eleout<<(eN+1)<<"\t"<<(1+elementNodesArray_new[eN*3+0])
  // 	      <<"\t"<<(1+elementNodesArray_new[eN*3+1])
  // 	      <<"\t"<<(1+elementNodesArray_new[eN*3+2])
  // 	      <<std::endl;
  //       }
  //     nodeout<<mesh.nNodes_global<<" 2 0 0"<<std::endl;
  //     for (int nN=0;nN<mesh.nNodes_global;nN++)
  //       {
  // 	nodeout<<(nN+1)<<"\t"<<nodeArray_new[nN*3+0]<<"\t"<<nodeArray_new[nN*3+1]<<std::endl;
  //       }
  //     eleout.close();
  //     partout.close();
  //5. At this point we have new, renumbered and sorted the global element and node based information, and
  //we have it all on each processor so we can add overlap
  set<int> elements_overlap,nodes_overlap;
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    for (int nN=0;nN<mesh.nNodes_element;nN++)
      {
        int nN_global = elementNodesArray_new[eN*mesh.nNodes_element+nN];
        if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
          nodes_overlap.insert(nN_global);
      }
    

  if (nElements_overlap > 0)
    {
      //get all elements in the node stars and their nodes
      // for (int nN_new=nodeOffsets_new[rank]; nN_new < nodeOffsets_new[rank+1]; nN_new++)
      //   {
      //     int nN = nodeNumbering_global_new2old[nN_new];
      //     for (int offset =mesh.nodeElementOffsets[nN];offset<mesh.nodeElementOffsets[nN+1];offset++)
      //       {
      //         int eN = mesh.nodeElementsArray[offset];
      //         int eN_new = elementNumbering_global_old2new[eN];
      //         if (eN_new < elementOffsets_new[rank] or eN_new >= elementOffsets_new[rank+1])
      //           {
      //             elements_overlap.insert(eN_new);
      //             for (int nN_element=0;nN_element<mesh.nNodes_element;nN_element++)
      //               {
      //                 int nN_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN_element];
      //                 if (nN_global < nodeOffsets_new[rank] or nN_global >= nodeOffsets_new[rank+1])
      //                   nodes_overlap.insert(nN_global);
      //               }
      //           }
      //       }
      //   }
      //get all elements in the node stars of owned nodes and those elements' nodes
      for (set<int>::iterator nN=nodes_subdomain_owned.begin();nN != nodes_subdomain_owned.end();nN++)
        {
          for (int offset =mesh.nodeElementOffsets[*nN];offset<mesh.nodeElementOffsets[(*nN)+1];offset++)
            {
              int eN = mesh.nodeElementsArray[offset];
              int eN_new = elementNumbering_global_old2new[eN];
              if (eN_new < elementOffsets_new[rank] or eN_new >= elementOffsets_new[rank+1])
                {
                  elements_overlap.insert(eN_new);
                  for (int nN_element=0;nN_element<mesh.nNodes_element;nN_element++)
                    {
                      int nN_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN_element];
                      if (nN_global < nodeOffsets_new[rank] or nN_global >= nodeOffsets_new[rank+1])
                        nodes_overlap.insert(nN_global);
                    }
                }
            }
        }
      //get all the element neighbors
      for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
        for(int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
          {
            int eN_ebN = elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN];
            if (eN_ebN >= 0 && 
                (eN_ebN < elementOffsets_new[rank] || eN_ebN >= elementOffsets_new[rank+1]))
              {
                elements_overlap.insert(eN_ebN);
                for (int nN=0;nN<mesh.nNodes_element;nN++)
                  {
                    int nN_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN];
                    if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
                      nodes_overlap.insert(nN_global);
                  }
              }
          }
    }
  for (int layer=1;layer<nElements_overlap;layer++)
    {
      for (set<int>::iterator eN_p=elements_overlap.begin();eN_p != elements_overlap.end();eN_p++)
        {
          int eN_global = *eN_p;
          for(int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
            {
              int eN_ebN = elementNeighborsArray_new[eN_global*mesh.nElementBoundaries_element+ebN];
              if (eN_ebN >= 0 &&
                  (eN_ebN < elementOffsets_new[rank] || eN_ebN >= elementOffsets_new[rank+1]))
                {
                  elements_overlap.insert(eN_ebN);
                  for (int nN=0;nN<mesh.nNodes_element;nN++)
                    {
                      int nN_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN];
                      if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
                        nodes_overlap.insert(nN_global);
                    }
                }
            }
        }            
    }
  //6. Now build subdomain mesh
  //
  //set what we know
  if (mesh.subdomainp == NULL)
    mesh.subdomainp = new Mesh();
  mesh.subdomainp->nElements_global = nElements_subdomain_new[rank] + elements_overlap.size();
  mesh.subdomainp->nNodes_global = nNodes_subdomain_new[rank] + nodes_overlap.size();
  mesh.subdomainp->nNodes_element = mesh.nNodes_element;
  mesh.subdomainp->nNodes_elementBoundary = mesh.nNodes_elementBoundary;
  mesh.subdomainp->nElementBoundaries_element = mesh.nElementBoundaries_element;
  //load the elements and nodes
  valarray<int> nodeNumbering_subdomain2global(mesh.subdomainp->nNodes_global);
  map<int,int> nodeNumbering_global2subdomain;
  mesh.subdomainp->nodeArray = new double[mesh.subdomainp->nNodes_global*3];
  mesh.subdomainp->nodeMaterialTypes = new int[mesh.subdomainp->nNodes_global];
  for(int nN=0;nN<nNodes_subdomain_new[rank];nN++)
    {
      int nN_global = nN + nodeOffsets_new[rank];
      nodeNumbering_subdomain2global[nN] = nN_global;
      nodeNumbering_global2subdomain[nN_global] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = nodeArray_new[nN_global*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = nodeArray_new[nN_global*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = nodeArray_new[nN_global*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN]= nodeMaterialTypes_new[nN_global];
    }
  //note: sets in C++ are sorted so the overlap is laid out in
  //contiguous chunks corresponding to the partitions
  set<int>::iterator nN_p=nodes_overlap.begin();
  for(int nN=nNodes_subdomain_new[rank];nN < nNodes_subdomain_new[rank] + int(nodes_overlap.size()); nN++)
    {
      int nN_global = *nN_p++;
      nodeNumbering_subdomain2global[nN] = nN_global;
      nodeNumbering_global2subdomain[nN_global] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = nodeArray_new[nN_global*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = nodeArray_new[nN_global*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = nodeArray_new[nN_global*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN]= nodeMaterialTypes_new[nN_global];
    }
  mesh.subdomainp->elementNodesArray = new int[mesh.subdomainp->nElements_global*mesh.subdomainp->nNodes_element];
  mesh.subdomainp->elementMaterialTypes = new int[mesh.subdomainp->nElements_global];
  valarray<int> elementNumbering_subdomain2global(mesh.subdomainp->nElements_global);
  for (int eN=0;eN<nElements_subdomain_new[rank];eN++)
    {
      int eN_global = eN+elementOffsets_new[rank];
      elementNumbering_subdomain2global[eN] = eN_global;
      mesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypes_new[eN_global];
      for (int nN=0;nN<mesh.subdomainp->nNodes_element;nN++)
        mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+nN] = 
          nodeNumbering_global2subdomain[elementNodesArray_new[eN_global*mesh.nNodes_element + nN]];
    }
  set<int>::iterator eN_p=elements_overlap.begin();
  for(int eN=nElements_subdomain_new[rank];eN < nElements_subdomain_new[rank]+int(elements_overlap.size());eN++)
    {
      int eN_global = *eN_p++;
      mesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypes_new[eN_global];
      elementNumbering_subdomain2global[eN] = eN_global;
      for (int nN=0;nN<mesh.subdomainp->nNodes_element;nN++)
        mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+nN] = 
          nodeNumbering_global2subdomain[elementNodesArray_new[eN_global*mesh.nNodes_element + nN]];
    }

  if (mesh.subdomainp->nNodes_element == 2)
    {
      constructElementBoundaryElementsArray_edge(*mesh.subdomainp);
      allocateGeometricInfo_edge(*mesh.subdomainp);
      computeGeometricInfo_edge(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 3)
    {
      constructElementBoundaryElementsArray_triangle(*mesh.subdomainp);
      allocateGeometricInfo_triangle(*mesh.subdomainp);
      computeGeometricInfo_triangle(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 4)
    {
      constructElementBoundaryElementsArray_tetrahedron(*mesh.subdomainp);
      allocateGeometricInfo_tetrahedron(*mesh.subdomainp);
      computeGeometricInfo_tetrahedron(*mesh.subdomainp);
    }
  if (mesh.elementBoundaryMaterialTypes != NULL)
    {
      assert(mesh.elementBoundariesArray != NULL);
      //mwftodo need to copy over elementBoundaryMaterialTypes now
      for (int eN=0;eN<mesh.subdomainp->nElements_global;eN++)
        {
          int eN_global_new = elementNumbering_subdomain2global[eN];
          int eN_global_old = elementNumbering_global_new2old[eN_global_new];
          for (int ebN_element = 0; ebN_element < mesh.nElementBoundaries_element; ebN_element++)
            {
              int ebN_global_old = mesh.elementBoundariesArray[eN_global_old*mesh.nElementBoundaries_element+ebN_element];
              int ebN_subdomain = mesh.subdomainp->elementBoundariesArray[eN*mesh.nElementBoundaries_element+ebN_element];
              mesh.subdomainp->elementBoundaryMaterialTypes[ebN_subdomain] = mesh.elementBoundaryMaterialTypes[ebN_global_old];
            }
        }
    }
  //now we've got the old mesh in the old ordering and the subdomain mesh in the new ordering
  //the first chunk of nodes and elements are the owned elements so we need to know how many of those there are
  //and the offset of the first one so we can compute subdomain2global
  mesh.nodeOffsets_subdomain_owned = new int[size+1];
  mesh.elementOffsets_subdomain_owned = new int[size+1];
  for (int sdN=0;sdN<size+1;sdN++)
    {
      mesh.nodeOffsets_subdomain_owned[sdN] = nodeOffsets_new[sdN];
      mesh.elementOffsets_subdomain_owned[sdN] = elementOffsets_new[sdN];
    }
  //we also need the subdomain 2 new global mappings
  mesh.nodeNumbering_subdomain2global = new int[mesh.subdomainp->nNodes_global];
  for (int nN=0;nN<mesh.subdomainp->nNodes_global;nN++)
    mesh.nodeNumbering_subdomain2global[nN] = nodeNumbering_subdomain2global[nN];
  mesh.elementNumbering_subdomain2global = new int[mesh.subdomainp->nElements_global];
  for (int eN=0;eN<mesh.subdomainp->nElements_global;eN++)
    mesh.elementNumbering_subdomain2global[eN] = elementNumbering_subdomain2global[eN];

  ISRestoreIndices(elementNumberingIS_global_old2new,&elementNumbering_global_old2new);

  ISDestroy(&elementPartitioningIS_new);
  ISDestroy(&elementNumberingIS_subdomain_old2new);
  ISDestroy(&elementNumberingIS_global_old2new);

  ISRestoreIndices(nodeNumberingIS_global_new2old,&nodeNumbering_global_new2old);

  ISDestroy(&nodeNumberingIS_new2old);
  ISDestroy(&nodeNumberingIS_global_new2old);

  return 0;
}

int partitionNodes(Mesh& mesh, int nNodes_overlap)
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
  /***********************************************************************
    partition domain based on nodes rather than elements, basically repeats
    partitionElements with this one modification

    right now generates equivalent of 1 layer of overlap regardless of input

    1. Partition nodes in the default partition of contiguous chunks with
       input ordering

    2. Determine nodal connectivity on local processor

    3. Generate new nodal partition using PETSc interface

    4. Collect elements containing locally owned nodes and assign ownership

    5. Generate global element numbering for new subdomain ownership

    6. Create overlap (ghost) information for nodes and elements 

    7. March through additional layers of overlap if requested,

    8. Build subdomain meshes in new numbering
  ***********************************************************************/
  //
  //1. Build default nodal partition
  //
  //compute offsets to build processor (local) to global ordering for nodes
  //in default partitioning 
  valarray<int> nodeOffsets_old(size+1);
  nodeOffsets_old[0] = 0;
  for (int sdN=0; sdN < size; sdN++)
    {
      nodeOffsets_old[sdN+1] = nodeOffsets_old[sdN] + 
        int(mesh.nNodes_global)/size + (int(mesh.nNodes_global)%size > sdN);
    }
  //
  //2. Determine nodal connectivity on local processor, (local node star array)
  //
  int nNodes_subdomain = (nodeOffsets_old[rank+1] - nodeOffsets_old[rank]);
  PetscInt *nodeNeighborsOffsets_subdomain,*nodeNeighbors_subdomain,*weights_subdomain;
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain+1),&nodeNeighborsOffsets_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain*mesh.max_nNodeNeighbors_node),&nodeNeighbors_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain*mesh.max_nNodeNeighbors_node),&weights_subdomain);
  nodeNeighborsOffsets_subdomain[0] = 0;
  for (int nN = 0,offset=0; nN < nNodes_subdomain; nN++)
    {
      int nN_global = nodeOffsets_old[rank] + nN;
      for (int offset_global = mesh.nodeStarOffsets[nN_global]; 
           offset_global < mesh.nodeStarOffsets[nN_global+1]; offset_global++)
        {
          nodeNeighbors_subdomain[offset++] = mesh.nodeStarArray[offset_global];
        }
      nodeNeighborsOffsets_subdomain[nN+1]=offset;
      sort(&nodeNeighbors_subdomain[nodeNeighborsOffsets_subdomain[nN]],&nodeNeighbors_subdomain[nodeNeighborsOffsets_subdomain[nN+1]]);
      int weight= (nodeNeighborsOffsets_subdomain[nN+1] - nodeNeighborsOffsets_subdomain[nN]);
      for (int k=nodeNeighborsOffsets_subdomain[nN];k<nodeNeighborsOffsets_subdomain[nN+1];k++)
        weights_subdomain[k] = weight;
    }
  //
  //3. Generate new nodal partition using PETSc interface
  //
  Mat petscAdjacency;
  //   MatCreateMPIAdj(PROTEUS_COMM_WORLD,
  // 		  nNodes_subdomain, mesh.nNodes_global,
  // 		  &nodeNeighborsOffsets_subdomain[0], &nodeNeighbors_subdomain[0],
  // 		  &weights_subdomain[0],//PETSC_NULL,//ignore weighting for now
  // 		  &petscAdjacency);
  ierr = MatCreateMPIAdj(PROTEUS_COMM_WORLD,
			 nNodes_subdomain, 
			 mesh.nNodes_global,
			 nodeNeighborsOffsets_subdomain, 
			 nodeNeighbors_subdomain,
			 PETSC_NULL,//weights_subdomain,
			 &petscAdjacency);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  MatPartitioning petscPartition;
  MatPartitioningCreate(PROTEUS_COMM_WORLD,&petscPartition);
  MatPartitioningSetAdjacency(petscPartition,petscAdjacency);
  MatPartitioningSetFromOptions(petscPartition);
  
  //get petsc index set that has the new subdomain number for each node
  IS nodePartitioningIS_new;
  MatPartitioningApply(petscPartition,&nodePartitioningIS_new);
  MatPartitioningDestroy(&petscPartition); //gets petscAdjacency too I believe

  //determine the number of nodes per subdomain in new partitioning
  valarray<int> nNodes_subdomain_new(size);
  ISPartitioningCount(nodePartitioningIS_new,size,&nNodes_subdomain_new[0]);

  //need new offsets for subdomain to global numbering
  valarray<int> nodeOffsets_new(size+1);
  nodeOffsets_new[0] = 0;
  for (int sdN = 0; sdN < size; sdN++)
    nodeOffsets_new[sdN+1] = nodeOffsets_new[sdN] + nNodes_subdomain_new[sdN];

  //get the new node numbers for nodes on this subdomain
  IS nodeNumberingIS_subdomain_old2new;
  ISPartitioningToNumbering(nodePartitioningIS_new,&nodeNumberingIS_subdomain_old2new);

  //collect new node numbers for whole mesh so that subdomain reordering and renumbering
  //can be done easily

  IS nodeNumberingIS_global_old2new;
  ISAllGather(nodeNumberingIS_subdomain_old2new,&nodeNumberingIS_global_old2new);
  //mwf original and correct I believe 
  const PetscInt * nodeNumbering_global_old2new;//needs restore call
  ISGetIndices(nodeNumberingIS_global_old2new,&nodeNumbering_global_old2new);

  //reverse mapping for node numbers as well
  valarray<int> nodeNumbering_global_new2old(mesh.nNodes_global);
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
    nodeNumbering_global_new2old[nodeNumbering_global_old2new[nN]] = nN;

  //   //mwf debug
  //   if (rank == 0)
  //     {
  //       for (int sdN = 0; sdN < size+1; sdN++)
  // 	std::cout<<"partitionNodes rank= "<<rank<<" nodeOffset["<<sdN<<"]= "<<nodeOffsets_new[sdN]<<std::endl;
  //       for (int nN = 0; nN < mesh.nNodes_global; nN++)
  // 	{
  // 	  std::cout<<"partitionNodes rank= "<<rank<<" nN= "<<nN<<" old2new= "<<nodeNumbering_global_old2new[nN]<<" new2old= "<<nodeNumbering_global_new2old[nN]<<std::endl;
  // 	}
  //       for (int sdN = 0; sdN < size; sdN++)
  // 	{
  // 	  std::cout<<"============"<<std::endl;
  // 	  std::cout<<"partitionNodes rank= "<<sdN<<" nNodes_owned= "<<nodeOffsets_new[sdN+1]-nodeOffsets_new[sdN]<<" = "<<std::endl;
  // 	  for (int nN = nodeOffsets_new[sdN]; nN < nodeOffsets_new[sdN+1]; nN++)
  // 	    {
  // 	      std::cout<<"new number= "<<nN<<" <--> old number "<<nodeNumbering_global_new2old[nN]<<" x,y,z= "
  // 		       <<mesh.nodeArray[nodeNumbering_global_new2old[nN]*3+0]<<" , "
  // 		       <<mesh.nodeArray[nodeNumbering_global_new2old[nN]*3+1]<<" , "
  // 		       <<mesh.nodeArray[nodeNumbering_global_new2old[nN]*3+2]<<std::endl;
  // 	    }
  // 	}
      
  //     }
  //
  //4. To build subdomain meshes, go through and collect elements containing
  //   the locally owned nodes. Assign processor ownership of elements 
  //   
  MPI_Status status;
  PetscBT elementMask;
  PetscBTCreate(mesh.nElements_global,&elementMask);
  //get the owned element information 
  if (rank > 0)
    {
      MPI_Recv(elementMask,PetscBTLength(mesh.nElements_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status);
    }
  //mark the unmarked elements on this subdomain and store element numbers (in old numbering)
  set<int> elements_subdomain_owned;
  for (int nN = nodeOffsets_new[rank]; nN < nodeOffsets_new[rank+1]; nN++)
    {
      int nN_global_old = nodeNumbering_global_new2old[nN];
      for (int eN_star_offset = mesh.nodeElementOffsets[nN_global_old]; 
           eN_star_offset < mesh.nodeElementOffsets[nN_global_old+1]; eN_star_offset++)
        {
          int eN_star_old = mesh.nodeElementsArray[eN_star_offset];
          if (!PetscBTLookupSet(elementMask,eN_star_old))
            {
              elements_subdomain_owned.insert(eN_star_old);
            }
        }
    }
  //mwf debug
  //   for (int nN =  nodeOffsets_new[rank]; nN < nodeOffsets_new[rank+1]; nN++)
  //     {
  //       int nN_global_old = nodeNumbering_global_new2old[nN];
  //       int nElements_owned_nN =0;
  //       for (int eN_star_offset = mesh.nodeElementOffsets[nN_global_old]; 
  // 	   eN_star_offset < mesh.nodeElementOffsets[nN_global_old+1]; eN_star_offset++)
  // 	{
  // 	  int eN_star_old = mesh.nodeElementsArray[eN_star_offset];
  // 	  if (elements_subdomain_owned.find(eN_star_old) != elements_subdomain_owned.end())
  // 	    nElements_owned_nN++;
  // 	}
    
  //       if (nElements_owned_nN <= 0)
  // 	{
  // 	  std::cout<<"Problem? proc "<<rank<<" nN_new = "<<nN<<" nN_old= "<<nN_global_old<<" nElements_owned_for_nN = "<<nElements_owned_nN<<std::endl;
  // 	  //find out processor owners for node neighbors
  // 	  for (int offset = mesh.nodeStarOffsets[nN_global_old]; offset < mesh.nodeStarOffsets[nN_global_old+1]; offset++)
  // 	    {
  // 	      int nN_neig_old = mesh.nodeStarArray[offset];
  // 	      std::cout<<"\t neig node old "<<nN_neig_old<<" neig node new "<<nodeNumbering_global_old2new[nN_neig_old]<<" this proc offsets= ["
  // 		       <<nodeOffsets_new[rank]<<","<<nodeOffsets_new[rank+1]<<"];"<<std::endl;
  // 	    }
  // 	}
  //     }
  //pass off newly marked info
  if (rank < size-1)
    MPI_Send(elementMask,PetscBTLength(mesh.nElements_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&elementMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy"<<endl;

  //
  //5. Generate global element numbering corresponding to new subdomain ownership
  //
  valarray<int> nElements_subdomain_new(size),
    elementOffsets_new(size+1);
  for (int sdN = 0; sdN < size; sdN++)
    {
      if (sdN == rank)
        nElements_subdomain_new[sdN] = int(elements_subdomain_owned.size());
      else
        nElements_subdomain_new[sdN] = 0;
    }
  valarray<int> nElements_subdomain_new_send = nElements_subdomain_new;
  MPI_Allreduce(&nElements_subdomain_new_send[0],&nElements_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  //new size info
  elementOffsets_new[0] = 0;
  for (int sdN = 0; sdN < size; sdN++)
    elementOffsets_new[sdN+1] = elementOffsets_new[sdN] + nElements_subdomain_new[sdN];
  
  //map to old element numbering
  valarray<int> elementNumbering_subdomain_new2old(elements_subdomain_owned.size());
  set<int>::iterator eN_ownedp = elements_subdomain_owned.begin();
  for (int eN = 0; eN < int(elements_subdomain_owned.size()); eN++)
    {
      elementNumbering_subdomain_new2old[eN] = *eN_ownedp++;
    }
  //use Petsc IS to get global new2old numbering
  IS elementNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,elements_subdomain_owned.size(),&elementNumbering_subdomain_new2old[0],PETSC_COPY_VALUES,
                  &elementNumberingIS_subdomain_new2old);
  IS elementNumberingIS_global_new2old;
  ISAllGather(elementNumberingIS_subdomain_new2old,&elementNumberingIS_global_new2old);
  
  const PetscInt *elementNumbering_global_new2old;//needs to be restored
  ISGetIndices(elementNumberingIS_global_new2old,&elementNumbering_global_new2old);
  //reverse mapping
  valarray<int> elementNumbering_global_old2new(mesh.nElements_global);
  for (int eN = 0; eN < mesh.nElements_global; eN++)
    {
      elementNumbering_global_old2new[elementNumbering_global_new2old[eN]] = eN;
    }



  //4b,5b. repeat process to build global face numbering
  //first get element --> element boundaries array for new element but old element boundary numbering
  valarray<int> elementBoundariesArray_new(mesh.nElements_global*mesh.nElementBoundaries_element);
  for (int eN=0; eN < mesh.nElements_global; eN++)
    for (int ebN=0; ebN < mesh.nElementBoundaries_element; ebN++)
      {
        elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN] = 
          mesh.elementBoundariesArray[elementNumbering_global_new2old[eN]*mesh.nElementBoundaries_element+ebN];
      }
  MPI_Status status_elementBoundaries;
  PetscBT elementBoundaryMask; 
  PetscBTCreate(mesh.nElementBoundaries_global,&elementBoundaryMask);
  if (rank > 0) 
    {
      MPI_Recv(elementBoundaryMask,PetscBTLength(mesh.nElementBoundaries_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status_elementBoundaries);
    }
  //mark the unmarked faces on this subdomain and store the global face numbers
  //going through owned elements can pick up owned elementBoundaries on "outside" of owned nodes nodeStars 
  set<int> elementBoundaries_subdomain_owned;
  if (mesh.nNodes_element == 8)
    {

      int lface[6][4] = {{0,1,2,3},
                         {0,1,5,4},
                         {1,2,6,5},
                         {2,3,7,6},
                         {3,0,4,7},
                         {4,5,6,7}};
		       

      for (int nN = nodeOffsets_new[rank]; nN < nodeOffsets_new[rank+1]; nN++)
        {
          int nN_global_old = nodeNumbering_global_new2old[nN];
          //now get elements in node star
          for (int offset_old = mesh.nodeElementOffsets[nN_global_old]; 
               offset_old < mesh.nodeElementOffsets[nN_global_old+1]; offset_old++)
            {
              int eN_star_old = mesh.nodeElementsArray[offset_old];
              int eN_star_new = elementNumbering_global_old2new[eN_star_old];

              //loop through element boundaries on each element, want but want to skip
              //the element boundary across from the owned node
              for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
                {	      
                  bool foundNode = false;
                  for (int nNl=0; nNl<mesh.nNodes_elementBoundary ;nNl++)
                    {		  	 		 		    
                      int nN_global_old_across = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element + lface[ebN][nNl]];
                      if (nN_global_old_across == nN_global_old) foundNode = true;		  
                    } 
                  if (foundNode)
                    { 
                      int ebN_global=elementBoundariesArray_new[eN_star_new*mesh.nElementBoundaries_element+ebN];
                      if (!PetscBTLookupSet(elementBoundaryMask,ebN_global))
                        elementBoundaries_subdomain_owned.insert(ebN_global);
                    }
                }
	      	
            } 
        }
      
    }  
  else 
    {

      for (int nN = nodeOffsets_new[rank]; nN < nodeOffsets_new[rank+1]; nN++)
        {
          int nN_global_old = nodeNumbering_global_new2old[nN];
          //now get elements in node star
          for (int offset_old = mesh.nodeElementOffsets[nN_global_old]; 
               offset_old < mesh.nodeElementOffsets[nN_global_old+1]; offset_old++)
            {
              int eN_star_old = mesh.nodeElementsArray[offset_old];
              int eN_star_new = elementNumbering_global_old2new[eN_star_old];

              //loop through element boundaries on each element, want but want to skip
              //the element boundary across from the owned node
              for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
                {
                  int nN_global_old_across = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element+ebN];
                  if (nN_global_old_across != nN_global_old)
                    {
                      int ebN_global=elementBoundariesArray_new[eN_star_new*mesh.nElementBoundaries_element+ebN];
                      if (!PetscBTLookupSet(elementBoundaryMask,ebN_global))
                        elementBoundaries_subdomain_owned.insert(ebN_global);
                    }
                }
            }
        }
    }

  //ship off the mask
  if (rank < size-1)
    MPI_Send(elementBoundaryMask,PetscBTLength(mesh.nElementBoundaries_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&elementBoundaryMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy for elementBoundaries"<<endl;
  //get the number of elementBoundaries on each processor
  valarray<int> nElementBoundaries_subdomain_new(size),
    elementBoundaryOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nElementBoundaries_subdomain_new[sdN] = elementBoundaries_subdomain_owned.size();
    else
      nElementBoundaries_subdomain_new[sdN] = 0;
  valarray<int> nElementBoundaries_subdomain_new_send=nElementBoundaries_subdomain_new;
  MPI_Allreduce(&nElementBoundaries_subdomain_new_send[0],&nElementBoundaries_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  elementBoundaryOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    elementBoundaryOffsets_new[sdN+1] = elementBoundaryOffsets_new[sdN]+nElementBoundaries_subdomain_new[sdN];
  //Now as with elements and nodes build a global face numbering
  //resetting the face based information is a little different since much of this is currently built below based
  //on the element and node information
  //
  valarray<int> elementBoundaryNumbering_new2old(elementBoundaries_subdomain_owned.size());
  set<int>::iterator ebN_ownedp=elementBoundaries_subdomain_owned.begin();
  for (int ebN=0;ebN<int(elementBoundaries_subdomain_owned.size());ebN++)
    {
      elementBoundaryNumbering_new2old[ebN] = *ebN_ownedp++;
    }
  IS elementBoundaryNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,elementBoundaries_subdomain_owned.size(),&elementBoundaryNumbering_new2old[0],PETSC_COPY_VALUES,&elementBoundaryNumberingIS_subdomain_new2old);
  IS elementBoundaryNumberingIS_global_new2old;
  ISAllGather(elementBoundaryNumberingIS_subdomain_new2old,&elementBoundaryNumberingIS_global_new2old);
  const PetscInt *elementBoundaryNumbering_global_new2old;
  valarray<int> elementBoundaryNumbering_old2new_global(mesh.nElementBoundaries_global);
  ISGetIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);
  for (int ebN=0;ebN<mesh.nElementBoundaries_global;ebN++)
    {
      elementBoundaryNumbering_old2new_global[elementBoundaryNumbering_global_new2old[ebN]] = ebN;
    }
  for (int eN=0;eN < mesh.nElements_global; eN++)
    {
      int ebN_old;
      for (int ebN=0;ebN < mesh.nElementBoundaries_element; ebN++)
        {
          ebN_old = elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN];
          elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN] = elementBoundaryNumbering_old2new_global[ebN_old];
        }
    }
  
  //4c,5c. Build global edge numbering as well ownership is determined by who owns the left (0) node 
  //    of the edge
  map<NodeTuple<2>, int> nodesEdgeMap_global; //new global node numbers --> original edge numbering
  set<int> edges_subdomain_owned;
  for (int ig = 0; ig < mesh.nEdges_global; ig++)
    {
      const int nN0_global_old = mesh.edgeNodesArray[2*ig];
      const int nN1_global_old = mesh.edgeNodesArray[2*ig+1];
      const int nN0_global     = nodeNumbering_global_old2new[nN0_global_old];
      const int nN1_global     = nodeNumbering_global_old2new[nN1_global_old];
      int nodes[2];
      nodes[0] = nN0_global;
      nodes[1] = nN1_global;
      NodeTuple<2> et(nodes);
      nodesEdgeMap_global[et] = ig;
      if (nodeOffsets_new[rank] <= et.nodes[0] && et.nodes[0] < nodeOffsets_new[rank+1])
        edges_subdomain_owned.insert(ig);
    }

  valarray<int> nEdges_subdomain_new(size),
    edgeOffsets_new(size+1);
  
  for (int sdN=0; sdN < size; sdN++)
    if (sdN == rank)
      nEdges_subdomain_new[sdN] = edges_subdomain_owned.size();
    else
      nEdges_subdomain_new[sdN] = 0;
  //collect ownership info
  valarray<int> nEdges_subdomain_new_send=nEdges_subdomain_new;
  MPI_Allreduce(&nEdges_subdomain_new_send[0],&nEdges_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  edgeOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    edgeOffsets_new[sdN+1] = edgeOffsets_new[sdN]+nEdges_subdomain_new[sdN];
  
  //build new petsc numbering and global maps from old2new and new2old 
  valarray<int> edgeNumbering_new2old(edges_subdomain_owned.size());
  set<int>::iterator edges_ownedp = edges_subdomain_owned.begin();
  for (int i=0; i < int(edges_subdomain_owned.size());i++)
    edgeNumbering_new2old[i] = *edges_ownedp++;
    
  IS edgeNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,edges_subdomain_owned.size(),&edgeNumbering_new2old[0],PETSC_COPY_VALUES,&edgeNumberingIS_subdomain_new2old);
  IS edgeNumberingIS_global_new2old;
  ISAllGather(edgeNumberingIS_subdomain_new2old,&edgeNumberingIS_global_new2old);
  const PetscInt *edgeNumbering_global_new2old;
  valarray<int> edgeNumbering_old2new_global(mesh.nEdges_global);
  ISGetIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);
  for (int ig=0;ig<mesh.nEdges_global;ig++)
    {
      edgeNumbering_old2new_global[edgeNumbering_global_new2old[ig]] = ig;
    }
  
  
  //create  array with (new edge) --> (new node 0, new node 1)
  //and map from (new node 0, new node 1) --> (new global edge)
  valarray<int> edgeNodesArray_newNodesAndEdges(2*mesh.nEdges_global);
  map<NodeTuple<2>, int > nodesEdgeMap_global_new;
  for (int ig = 0; ig < mesh.nEdges_global; ig++)
    {
      const int nN0_global_old = mesh.edgeNodesArray[2*ig];
      const int nN1_global_old = mesh.edgeNodesArray[2*ig+1];
      const int nN0_global     = nodeNumbering_global_old2new[nN0_global_old];
      const int nN1_global     = nodeNumbering_global_old2new[nN1_global_old];
      
      const int edge_new = edgeNumbering_old2new_global[ig];
      edgeNodesArray_newNodesAndEdges[edge_new*2+0] = nN0_global;
      edgeNodesArray_newNodesAndEdges[edge_new*2+1] = nN1_global;
      int nodes[2];
      nodes[0] = nN0_global;
      nodes[1] = nN1_global;
      NodeTuple<2> et(nodes);
      nodesEdgeMap_global_new[et] = edge_new;
    }
  





  //
  //6. Figure out which elements are in node stars but are not locally owned, create ghost information
  //   for these, do the same for elements
  
  set<int> elements_overlap,nodes_overlap,elementBoundaries_overlap,edges_overlap;
  for (int nN = nodeOffsets_new[rank]; nN < nodeOffsets_new[rank+1]; nN++)
    {
      int nN_global_old = nodeNumbering_global_new2old[nN];
      //nodes in node star
      for (int offset_old = mesh.nodeStarOffsets[nN_global_old]; 
           offset_old < mesh.nodeStarOffsets[nN_global_old+1]; offset_old++)
        {
          int nN_neig_old = mesh.nodeStarArray[offset_old];
          int nN_neig_new = nodeNumbering_global_old2new[nN_neig_old];
	  
          bool offproc = nN_neig_new >=  nodeOffsets_new[rank+1] || nN_neig_new < nodeOffsets_new[rank];
          if (offproc)
            nodes_overlap.insert(nN_neig_new);
        }
      //now get elements, elementBoundaries, and edges in node star
      for (int offset_old = mesh.nodeElementOffsets[nN_global_old]; 
           offset_old < mesh.nodeElementOffsets[nN_global_old+1]; offset_old++)
        {
          int eN_star_old = mesh.nodeElementsArray[offset_old];
          int eN_star_new = elementNumbering_global_old2new[eN_star_old];
          bool offproc = eN_star_new >= elementOffsets_new[rank+1] || eN_star_new < elementOffsets_new[rank];
          if (offproc)
            elements_overlap.insert(eN_star_new);
          for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
            {
              int ebN_global=elementBoundariesArray_new[eN_star_new*mesh.nElementBoundaries_element+ebN];
              // 	      //mwf debug
              // 	      std::cout<<"partitionNode default overlap rank= "<<rank<<" nN_new= "<<nN<<" eN_star_new= "<<eN_star_new<<" ebN= "<<ebN
              // 		       <<" ebN_global= "<<ebN_global<<" ghost= "<<(ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
              // 		       <<" offsets= ["<<elementBoundaryOffsets_new[rank]<<","<<elementBoundaryOffsets_new[rank+1]<<"]"<<std::endl;
              if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                {
                  elementBoundaries_overlap.insert(ebN_global);
                }
            }
          for (int nN0=0;nN0<mesh.nNodes_element;nN0++)
            for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
              {
                const int nN0_global_old = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element+nN0];
                const int nN1_global_old = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element+nN1];
                const int nN0_global     = nodeNumbering_global_old2new[nN0_global_old];
                const int nN1_global     = nodeNumbering_global_old2new[nN1_global_old];
                bool foundEdge = false;
                int nodes[2];
                nodes[0] = nN0_global;
                nodes[1] = nN1_global;
                NodeTuple<2> et(nodes);
                const int edge_global = nodesEdgeMap_global_new[et];
                if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                  edges_overlap.insert(edge_global);
              }//edges
        }//elements in node star
    }//nodes on this processor

  //
  //7. If we want more layers of overlap, do we have to build connectivity info in new numbering then march out
  //   or can we just march through nodes in node_overlap and grab all of their elements that aren't owned?
  //

  int overlap_remaining = nNodes_overlap -1; //default gives 1 layer overlap
  //last set of overlap nodes added
  set<int> last_nodes_added2overlap = nodes_overlap;
  while (overlap_remaining > 0)
    {
      set<int> new_nodes_overlap,new_elements_overlap,new_elementBoundaries_overlap,new_edges_overlap;
      set<int>::iterator nN_p = last_nodes_added2overlap.begin();
      while (nN_p != last_nodes_added2overlap.end())
        {
          int nN_global_new = *nN_p;
          int nN_global_old = nodeNumbering_global_new2old[nN_global_new];//need old numbering for connectivity
          for (int offset_old = mesh.nodeStarOffsets[nN_global_old]; 
               offset_old < mesh.nodeStarOffsets[nN_global_old+1]; offset_old++)
            {
              int nN_neig_old = mesh.nodeStarArray[offset_old];
              int nN_neig_new = nodeNumbering_global_old2new[nN_neig_old];
              //just need to check if neighbor is offprocessor, may already be in overlap
              //but set merge will take care of that
              bool offproc = nN_neig_new >=  nodeOffsets_new[rank+1] || nN_neig_new < nodeOffsets_new[rank];
              if (offproc)
                new_nodes_overlap.insert(nN_neig_new);
            }//node neighbor loop
          nN_p++;
        }//loop adding new nodes

      //loop through added nodes, grab elements only check if not on processor
      set<int>::iterator nN_newp = last_nodes_added2overlap.begin();
      while (nN_newp != last_nodes_added2overlap.end())
        {
          int nN_global_new = *nN_newp;
          int nN_global_old = nodeNumbering_global_new2old[nN_global_new];//need old numbering for connectivity
          for (int offset_old = mesh.nodeElementOffsets[nN_global_old];
               offset_old < mesh.nodeElementOffsets[nN_global_old+1]; offset_old++)
            {
              int eN_star_old = mesh.nodeElementsArray[offset_old];
              int eN_star_new = elementNumbering_global_old2new[eN_star_old];
              bool offproc = eN_star_new >= elementOffsets_new[rank+1] || eN_star_new < elementOffsets_new[rank];
              if (offproc)
                new_elements_overlap.insert(eN_star_new);
              //element boundaries too
              for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
                {
                  int ebN_global=elementBoundariesArray_new[eN_star_new*mesh.nElementBoundaries_element+ebN];
                  // 		  //mwf debug
                  // 		  std::cout<<"partitionNode overlap_remaining= "<<overlap_remaining<<" rank= "<<rank<<" nN_global_new= "<<nN_global_new<<" eN_star_new= "<<eN_star_new<<" ebN= "<<ebN
                  // 			   <<" ebN_global= "<<ebN_global<<" ghost= "<<(ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                  // 			   <<" offsets= ["<<elementBoundaryOffsets_new[rank]<<","<<elementBoundaryOffsets_new[rank+1]<<std::endl;
                  if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                    new_elementBoundaries_overlap.insert(ebN_global);
                }//element boundaries
              //edges
              for (int nN0=0;nN0<mesh.nNodes_element;nN0++)
                for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
                  {
                    const int nN0_global_old = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element+nN0];
                    const int nN1_global_old = mesh.elementNodesArray[eN_star_old*mesh.nNodes_element+nN1];
                    const int nN0_global     = nodeNumbering_global_old2new[nN0_global_old];
                    const int nN1_global     = nodeNumbering_global_old2new[nN1_global_old];
                    bool foundEdge = false;
                    int nodes[2];
                    nodes[0] = nN0_global;
                    nodes[1] = nN1_global;
                    NodeTuple<2> et(nodes);
                    const int edge_global = nodesEdgeMap_global_new[et];
                    if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                      new_edges_overlap.insert(edge_global);
                  }//edges
            }//elements in node star
          nN_newp++;
        }//new nodes
      last_nodes_added2overlap.clear();
      set_difference(new_nodes_overlap.begin(),new_nodes_overlap.end(),
                     nodes_overlap.begin(),nodes_overlap.end(),
                     insert_iterator<set<int> >(last_nodes_added2overlap,
                                                last_nodes_added2overlap.begin()));

      //could do a set_merge 
      for (set<int>::iterator nN_addedp = new_nodes_overlap.begin();
           nN_addedp != new_nodes_overlap.end();
           nN_addedp++)
        {
          nodes_overlap.insert(*nN_addedp);
        }
      for (set<int>::iterator eN_addedp = new_elements_overlap.begin();
           eN_addedp != new_elements_overlap.end();
           eN_addedp++)
        {
          elements_overlap.insert(*eN_addedp);
        }
      for (set<int>::iterator ebN_addedp = new_elementBoundaries_overlap.begin();
           ebN_addedp != new_elementBoundaries_overlap.end();
           ebN_addedp++)
        {
          elementBoundaries_overlap.insert(*ebN_addedp);
        }
      for (set<int>::iterator edge_addedp = new_edges_overlap.begin();
           edge_addedp != new_edges_overlap.end();
           edge_addedp++)
        {
          edges_overlap.insert(*edge_addedp);
        }
      //example calls

      overlap_remaining--;
    }//ovelap loop

  //
  //8. Build subdomain meshes in new numbering, assumes memory not allocated in subdomain mesh
  //   
  if (mesh.subdomainp == NULL)
    mesh.subdomainp = new Mesh();
  mesh.subdomainp->nElements_global = nElements_subdomain_new[rank] + elements_overlap.size();
  mesh.subdomainp->nNodes_global    = nNodes_subdomain_new[rank] + nodes_overlap.size();
  mesh.subdomainp->nElementBoundaries_global = nElementBoundaries_subdomain_new[rank]+elementBoundaries_overlap.size();
  mesh.subdomainp->nEdges_global   = nEdges_subdomain_new[rank]+edges_overlap.size();
  mesh.subdomainp->nNodes_element   = mesh.nNodes_element;
  mesh.subdomainp->nNodes_elementBoundary = mesh.nNodes_elementBoundary;
  mesh.subdomainp->nElementBoundaries_element = mesh.nElementBoundaries_element;
  //subdomain 2 global mappings (including ghost info)
  valarray<int> nodeNumbering_subdomain2global(mesh.subdomainp->nNodes_global);
  valarray<int> elementNumbering_subdomain2global(mesh.subdomainp->nElements_global);
  valarray<int> elementBoundaryNumbering_subdomain2global(mesh.subdomainp->nElementBoundaries_global);
  valarray<int> edgeNumbering_subdomain2global(mesh.subdomainp->nEdges_global);
  map<int,int> nodeNumbering_global2subdomain;
  map<int,int> elementBoundaryNumbering_global2subdomain;
  mesh.subdomainp->nodeArray = new double[mesh.subdomainp->nNodes_global*3];
  mesh.subdomainp->nodeMaterialTypes = new int[mesh.subdomainp->nNodes_global];
  //locally owned
  for (int nN = 0; nN < nNodes_subdomain_new[rank]; nN++)
    {
      int nN_global_new = nN + nodeOffsets_new[rank];
      int nN_global_old = nodeNumbering_global_new2old[nN_global_new];
      nodeNumbering_subdomain2global[nN] = nN_global_new;
      nodeNumbering_global2subdomain[nN_global_new] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = mesh.nodeArray[nN_global_old*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = mesh.nodeArray[nN_global_old*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = mesh.nodeArray[nN_global_old*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN] = mesh.nodeMaterialTypes[nN_global_old];
    }
  //ghost
  //note: sets in C++ are sorted so the overlap is laid out in
  //contiguous chunks corresponding to the partitions
  set<int>::iterator nN_p = nodes_overlap.begin();
  for (int nN = nNodes_subdomain_new[rank]; nN < nNodes_subdomain_new[rank] + int(nodes_overlap.size()); nN++)
    {
      int nN_global_new = *nN_p++;
      int nN_global_old = nodeNumbering_global_new2old[nN_global_new];
      nodeNumbering_subdomain2global[nN] = nN_global_new;
      nodeNumbering_global2subdomain[nN_global_new] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = mesh.nodeArray[nN_global_old*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = mesh.nodeArray[nN_global_old*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = mesh.nodeArray[nN_global_old*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN] = mesh.nodeMaterialTypes[nN_global_old];
    }
  mesh.subdomainp->elementNodesArray = new int[mesh.subdomainp->nElements_global*mesh.subdomainp->nNodes_element];
  mesh.subdomainp->elementMaterialTypes = new int[mesh.subdomainp->nElements_global];
  //locally owned
  for (int eN = 0; eN < nElements_subdomain_new[rank]; eN++)
    {
      int eN_global_new = elementOffsets_new[rank] + eN;
      int eN_global_old = elementNumbering_global_new2old[eN_global_new];
      elementNumbering_subdomain2global[eN] = eN_global_new;
      mesh.subdomainp->elementMaterialTypes[eN] = mesh.elementMaterialTypes[eN_global_old];
      for (int nN =  0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          int nN_global_old = mesh.elementNodesArray[eN_global_old*mesh.nNodes_element + nN];
          int nN_global_new = nodeNumbering_global_old2new[nN_global_old];
          int nN_subdomain  = nodeNumbering_global2subdomain[nN_global_new];
          mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN]= nN_subdomain;
        }

    }
  //ghost
  set<int>::iterator eN_p = elements_overlap.begin();
  for (int eN = nElements_subdomain_new[rank]; eN < nElements_subdomain_new[rank] + int(elements_overlap.size()); eN++)
    {
      int eN_global_new = *eN_p++;
      int eN_global_old = elementNumbering_global_new2old[eN_global_new];
      elementNumbering_subdomain2global[eN] = eN_global_new;
      mesh.subdomainp->elementMaterialTypes[eN] = mesh.elementMaterialTypes[eN_global_old];
      for (int nN =  0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          int nN_global_old = mesh.elementNodesArray[eN_global_old*mesh.nNodes_element + nN];
          int nN_global_new = nodeNumbering_global_old2new[nN_global_old];
          int nN_subdomain  = nodeNumbering_global2subdomain[nN_global_new];
          mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN]= nN_subdomain;
        }
    }
  //element boundaries
  //locally owned
  for (int ebN=0; ebN < nElementBoundaries_subdomain_new[rank]; ebN++)
    {
      int ebN_global = ebN + elementBoundaryOffsets_new[rank];
      elementBoundaryNumbering_subdomain2global[ebN]=ebN_global;
      elementBoundaryNumbering_global2subdomain[ebN_global] = ebN;
    }
  //ghost
  set<int>::iterator ebN_p = elementBoundaries_overlap.begin();
  for(int ebN=nElementBoundaries_subdomain_new[rank];ebN < nElementBoundaries_subdomain_new[rank] + int(elementBoundaries_overlap.size()); ebN++)
    {
      int ebN_global = *ebN_p++;
      elementBoundaryNumbering_subdomain2global[ebN] = ebN_global;
      elementBoundaryNumbering_global2subdomain[ebN_global] = ebN;
    }
  //need elementBoundariesArray to assign consistent numbering on subdomain
  mesh.subdomainp->elementBoundariesArray = 
    new int[mesh.subdomainp->nElements_global*mesh.subdomainp->nElementBoundaries_element];
  for (int eN=0;eN<nElements_subdomain_new[rank];eN++)
    {
      int eN_global = eN+elementOffsets_new[rank];
      for (int ebN=0;ebN<mesh.subdomainp->nElementBoundaries_element;ebN++)
        mesh.subdomainp->elementBoundariesArray[eN*mesh.subdomainp->nElementBoundaries_element+ebN] =
          elementBoundaryNumbering_global2subdomain[elementBoundariesArray_new[eN_global*mesh.nElementBoundaries_element + ebN]];
    }
  //ghost elements
  set<int>::iterator eN_p2 = elements_overlap.begin();
  for (int eN = nElements_subdomain_new[rank]; eN < nElements_subdomain_new[rank] + int(elements_overlap.size()); eN++)
    {
      int eN_global_new = *eN_p2++;
      for (int ebN=0;ebN<mesh.subdomainp->nElementBoundaries_element;ebN++)
        mesh.subdomainp->elementBoundariesArray[eN*mesh.subdomainp->nElementBoundaries_element+ebN] =
          elementBoundaryNumbering_global2subdomain[elementBoundariesArray_new[eN_global_new*mesh.nElementBoundaries_element + ebN]];

    }      

  //edges
  mesh.subdomainp->edgeNodesArray = new int[mesh.subdomainp->nEdges_global*2];
  //locally owned
  for (int i=0; i < nEdges_subdomain_new[rank]; i++)
    {
      const int ig = i+edgeOffsets_new[rank];
      const int nN0_global = edgeNodesArray_newNodesAndEdges[ig*2+0];
      const int nN1_global = edgeNodesArray_newNodesAndEdges[ig*2+1];
      //mwf todo double check can always count on having nodes on this processor
      const int nN0_subdomain = nodeNumbering_global2subdomain[nN0_global];
      const int nN1_subdomain = nodeNumbering_global2subdomain[nN1_global];
      mesh.subdomainp->edgeNodesArray[2*i+0]=nN0_subdomain;
      mesh.subdomainp->edgeNodesArray[2*i+1]=nN1_subdomain;
      edgeNumbering_subdomain2global[i] = ig;
    }
  //ghost
  set<int>::iterator edge_p = edges_overlap.begin();
  for (int i=nEdges_subdomain_new[rank]; i < nEdges_subdomain_new[rank] + int(edges_overlap.size()); i++)
    {
      const int ig =*edge_p++;
      const int nN0_global = edgeNodesArray_newNodesAndEdges[ig*2+0];
      const int nN1_global = edgeNodesArray_newNodesAndEdges[ig*2+1];
      //mwf todo make sure always have nodes for the edge on this processor 
      const int nN0_subdomain = nodeNumbering_global2subdomain[nN0_global];
      const int nN1_subdomain = nodeNumbering_global2subdomain[nN1_global];
      mesh.subdomainp->edgeNodesArray[2*i+0]=nN0_subdomain;
      mesh.subdomainp->edgeNodesArray[2*i+1]=nN1_subdomain;
      edgeNumbering_subdomain2global[i] = ig;
	
    }

  //now build rest of subdomain mesh connectivity information etc
  mesh.subdomainp->px = mesh.px;
  mesh.subdomainp->py = mesh.py;
  mesh.subdomainp->pz = mesh.pz;

  if (mesh.subdomainp->px != 0)
    {
      //constructElementBoundaryElementsArray_tetrahedron(*mesh.subdomainp);
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_tetrahedron(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_NURBS(*mesh.subdomainp);
      allocateGeometricInfo_NURBS(*mesh.subdomainp);
      computeGeometricInfo_NURBS(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 2)
    {
      //constructElementBoundaryElementsArray_edge(*mesh.subdomainp);
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_edge(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_edge(*mesh.subdomainp);
      allocateGeometricInfo_edge(*mesh.subdomainp);
      computeGeometricInfo_edge(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 3)
    {
      //constructElementBoundaryElementsArray_triangle(*mesh.subdomainp);
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_triangle(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_triangle(*mesh.subdomainp);
      allocateGeometricInfo_triangle(*mesh.subdomainp);
      computeGeometricInfo_triangle(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 4)
    {
      //constructElementBoundaryElementsArray_tetrahedron(*mesh.subdomainp);
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_tetrahedron(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_tetrahedron(*mesh.subdomainp);
      allocateGeometricInfo_tetrahedron(*mesh.subdomainp);
      computeGeometricInfo_tetrahedron(*mesh.subdomainp);
    }

  else if (mesh.subdomainp->nNodes_element == 8)
    {
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_hexahedron(*mesh.subdomainp);
      allocateGeometricInfo_hexahedron(*mesh.subdomainp);
      computeGeometricInfo_hexahedron(*mesh.subdomainp);
    }


  if (mesh.elementBoundaryMaterialTypes != NULL)
    {
      assert(mesh.elementBoundariesArray != NULL);
      assert(mesh.subdomainp->elementBoundariesArray != NULL);
      for (int eN=0;eN<mesh.subdomainp->nElements_global;eN++)
        {
          int eN_global_new = elementNumbering_subdomain2global[eN];
          int eN_global_old = elementNumbering_global_new2old[eN_global_new];
          for (int ebN_element = 0; ebN_element < mesh.nElementBoundaries_element; ebN_element++)
            {
              int ebN_global_old = mesh.elementBoundariesArray[eN_global_old*mesh.nElementBoundaries_element+ebN_element];
              int ebN_subdomain = mesh.subdomainp->elementBoundariesArray[eN*mesh.nElementBoundaries_element+ebN_element];
              mesh.subdomainp->elementBoundaryMaterialTypes[ebN_subdomain] = mesh.elementBoundaryMaterialTypes[ebN_global_old];
            }
        }
    }
  //transfer information about owned nodes and elements to mesh
  if (mesh.nodeOffsets_subdomain_owned) 
    delete [] mesh.nodeOffsets_subdomain_owned;
  if (mesh.elementOffsets_subdomain_owned)
    delete [] mesh.elementOffsets_subdomain_owned;
  if (mesh.elementBoundaryOffsets_subdomain_owned)
    delete [] mesh.elementBoundaryOffsets_subdomain_owned;
  if (mesh.edgeOffsets_subdomain_owned)
    delete [] mesh.edgeOffsets_subdomain_owned;
  mesh.nodeOffsets_subdomain_owned    = new int[size+1];
  mesh.elementOffsets_subdomain_owned = new int[size+1];
  mesh.elementBoundaryOffsets_subdomain_owned = new int[size+1];
  mesh.edgeOffsets_subdomain_owned = new int[size+1];
  for (int sdN = 0; sdN < size+1; sdN++)
    {
      mesh.nodeOffsets_subdomain_owned[sdN]    = nodeOffsets_new[sdN];
      mesh.elementOffsets_subdomain_owned[sdN] = elementOffsets_new[sdN];
      mesh.elementBoundaryOffsets_subdomain_owned[sdN] = elementBoundaryOffsets_new[sdN];
      mesh.edgeOffsets_subdomain_owned[sdN] = edgeOffsets_new[sdN];
    }
  if (mesh.nodeNumbering_subdomain2global)
    delete [] mesh.nodeNumbering_subdomain2global;
  mesh.nodeNumbering_subdomain2global = new int[mesh.subdomainp->nNodes_global];
  for (int nN = 0; nN < mesh.subdomainp->nNodes_global; nN++)
    mesh.nodeNumbering_subdomain2global[nN] = nodeNumbering_subdomain2global[nN];
  if (mesh.elementNumbering_subdomain2global)
    delete [] mesh.elementNumbering_subdomain2global;
  mesh.elementNumbering_subdomain2global = new int[mesh.subdomainp->nElements_global];
  for (int eN = 0; eN < mesh.subdomainp->nElements_global; eN++)
    mesh.elementNumbering_subdomain2global[eN] = elementNumbering_subdomain2global[eN];
  //
  if (mesh.elementBoundaryNumbering_subdomain2global)
    delete [] mesh.elementBoundaryNumbering_subdomain2global;
  mesh.elementBoundaryNumbering_subdomain2global = new int[mesh.subdomainp->nElementBoundaries_global];
  for (int ebN = 0; ebN < mesh.subdomainp->nElementBoundaries_global; ebN++)
    mesh.elementBoundaryNumbering_subdomain2global[ebN] = elementBoundaryNumbering_subdomain2global[ebN];
  //
  if (mesh.edgeNumbering_subdomain2global)
    delete [] mesh.edgeNumbering_subdomain2global;
  mesh.edgeNumbering_subdomain2global = new int[mesh.subdomainp->nEdges_global];
  for (int i=0; i< mesh.subdomainp->nEdges_global; i++)
    mesh.edgeNumbering_subdomain2global[i] = edgeNumbering_subdomain2global[i];

  //cleanup
  ISRestoreIndices(nodeNumberingIS_global_old2new,&nodeNumbering_global_old2new);

  ISDestroy(&nodePartitioningIS_new);
  ISDestroy(&nodeNumberingIS_subdomain_old2new);
  ISDestroy(&nodeNumberingIS_global_old2new);
  
  ISRestoreIndices(elementNumberingIS_global_new2old,&elementNumbering_global_new2old);
  
  ISDestroy(&elementNumberingIS_subdomain_new2old);
  ISDestroy(&elementNumberingIS_global_new2old);

  ISRestoreIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);
  
  ISDestroy(&elementBoundaryNumberingIS_subdomain_new2old);
  ISDestroy(&elementBoundaryNumberingIS_global_new2old);

  ISRestoreIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);
  
  ISDestroy(&edgeNumberingIS_subdomain_new2old);
  ISDestroy(&edgeNumberingIS_global_new2old);

  return 0;
}

int enforceMemoryLimit(int rank, double max_rss_gb,const char* msg)
{
  double current, current_global,gb(1.0e-9);
  PetscBarrier(NULL);
  current = double(getCurrentRSS())*gb;
  PetscBarrier(NULL);
  current_global=0.0;
  MPI_Allreduce(&current,&current_global,1,MPI_DOUBLE,MPI_MAX,PROTEUS_COMM_WORLD);
  if (current > max_rss_gb)
    {
      std::cout<<"Raising PETSC_ERR_MEM, Memory usage  on rank "<<rank<<'\t'<<current<<"GB"<<'\t'<<"limit "<<max_rss_gb<<std::endl;
      SETERRABORT(PROTEUS_COMM_WORLD,PETSC_ERR_MEM,"Exceeded Proteus memory limit");
    }
  if (rank ==  0)
    std::cout<<msg<<std::endl
	     <<"Max memory usage per core "<<current_global<<"GB"<<std::endl;
  return 0;
}

int partitionNodesFromTetgenFiles(const char* filebase, int indexBase, Mesh& newMesh, int nNodes_overlap)
{
  using namespace std;
  PetscErrorCode ierr;
  PetscMPIInt size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  PetscLogStage partitioning_stage;
  PetscLogStageRegister("Mesh Partition",&partitioning_stage);
  PetscLogStagePush(partitioning_stage);

  /***********************************************************************
    partition domain based on the nodes without reading in the global mesh.

    1. Partition nodes in the default partition of contiguous chunks with
       input ordering

    2. Determine nodal connectivity on local processor

    3. Generate new nodal partition using PETSc interface

    4. Collect subdomain elements: any that contain owned nodes; tag ownership

    4a,b. Collect element boundaries and edges on the subdomain elements and tag ownership

    5. Generate global element, element boundary, and edge  numbering for new subdomain ownership

    6. Create overlap (ghost) information for nodes and elements 

    7. March through additional layers of overlap if requested,

    8. Build subdomain meshes in new numbering
  ***********************************************************************/
  //
  //0. Set up tetgen files. Note, tetgen should have been run with -feen to get alll the faces and elements
  //
  bool failed = false;
  const int simplexDim = 4;
  const int vertexDim  = 3;
  using namespace IOutils;
  std::string vertexFileName  = std::string(filebase) + ".node" ;
  std::string elementFileName = std::string(filebase) + ".ele" ;
  std::string elementBoundaryFileName  = std::string(filebase) + ".face" ;
  std::string edgeFileName  = std::string(filebase) + ".edge" ;
  
  //
  //1. Build default nodal partition
  //
  //compute offsets to build subdomain2global ordering for nodes
  //in default partitioning of the domain into subdomains. All
  //we need is the number of nodes in the global mesh
  //
  //read nodes for tetgen format
  //first just read the number of nodes and whether or not there are node tags
  int read_elements_event;
  PetscLogEventRegister("Read eles",0,&read_elements_event);
  PetscLogEventBegin(read_elements_event,0,0,0,0);
  std::ifstream vertexFile(vertexFileName.c_str());  
  if (!vertexFile.good())
    {
      std::cerr<<"cannot open Tetgen node file "
               <<vertexFileName<<std::endl;
      failed = true;
      return failed;
    }
  int hasVertexMarkers(0),hasVertexAttributes(0),nSpace(3),nNodes_global;
  //read number of vertices and whether node flags are provided
  vertexFile >> eatcomments >> nNodes_global >> nSpace >> hasVertexAttributes >> hasVertexMarkers >> eatline ;
  assert(nNodes_global > 0);
  assert(nSpace == 3);
  newMesh.nNodes_global = nNodes_global;
  newMesh.nNodes_element = simplexDim;
  newMesh.nNodes_elementBoundary = simplexDim-1;
  newMesh.nElementBoundaries_element = simplexDim;
  if (hasVertexAttributes > 0)
    {
      std::cerr<<"WARNING Tetgen nodes hasAttributes= "<<hasVertexAttributes
               <<" > 0 will treat first value as integer id for boundary!!"<<std::endl;
      hasVertexMarkers = 1;
    }
  //don't need to read anymore from nodes for now

  //offsets provide the lower and upper bounds for the global numbering
  //first we just partition the nodes approximately equally ignoring connectivity
  valarray<int> nodeOffsets_old(size+1);
  nodeOffsets_old[0] = 0;
  for (int sdN=0; sdN < size; sdN++)
    {
      nodeOffsets_old[sdN+1] = nodeOffsets_old[sdN] + 
        int(nNodes_global)/size + (int(nNodes_global)%size > sdN);
    }
  int nNodes_subdomain_old = nodeOffsets_old[rank+1] - nodeOffsets_old[rank];
      
  //
  //2. Determine nodal connectivity (nodeStarArray) for nodes on subdomain
  //
  //connectivty commes from the topology (elements) file. We just grab the elements
  //that contain currently owned nodes
  std::ifstream elementFile(elementFileName.c_str());
  if (!elementFile.good())
    {
      std::cerr<<"cannot open Tetgen element file "
               <<elementFileName<<std::endl;
      failed = true;
      return failed;
    }
  //read elements
  int nNodesPerSimplex(simplexDim),hasElementMarkers = 0,nElements_global;
  elementFile >> eatcomments >> nElements_global >> nNodesPerSimplex >> hasElementMarkers >> eatline;
  assert(nElements_global > 0);
  assert(nNodesPerSimplex == simplexDim);
  newMesh.nElements_global = nElements_global;
  vector<int> element_nodes_old(4);
  vector<set<int> > nodeStar(nNodes_subdomain_old);
  map<int,vector<int> > elements_old;//elementNodesMap_old
  for (int ie = 0; ie < nElements_global; ie++)
    {
      int ne, nv;
      elementFile >> eatcomments >> ne;
      ne -= indexBase;
      assert(0 <= ne && ne < nElements_global && elementFile.good());
      for (int iv = 0; iv < simplexDim; iv++)
        {
          elementFile >> nv; 
          nv -= indexBase;
          assert(0 <= nv && nv < nNodes_global);
          element_nodes_old[iv] = nv;
        }
      //for each node on the element
      for (int iv = 0; iv < simplexDim; iv++)
        {
          //check if the node is owned by the subdomain
          int nN_star = element_nodes_old[iv];
          bool inSubdomain=false;
          if (nN_star >= nodeOffsets_old[rank] && nN_star < nodeOffsets_old[rank+1])
            {
              //this node is owned by the subdomain so 
              inSubdomain = true;
              for (int jv = 0; jv < simplexDim; jv++)
                {
                  if (iv != jv)
                    {
                      int nN_star_subdomain = nN_star-nodeOffsets_old[rank];
                      nodeStar[nN_star_subdomain].insert(element_nodes_old[jv]);
                    }
                }
            }
          if (inSubdomain)
            elements_old[ie] = element_nodes_old;
        }
      elementFile >> eatline;
    }//end ie
  elementFile.close();
  PetscLogEventEnd(read_elements_event,0,0,0,0);
  int repartition_nodes_event;
  PetscLogEventRegister("Repart nodes",0,&repartition_nodes_event);
  PetscLogEventBegin(repartition_nodes_event,0,0,0,0);
  //done reading element file for first time; will need to read again after node partitioning
  //build compact data structure for nodeStar
  valarray<int> nodeStarOffsets(nNodes_subdomain_old+1);
  nodeStarOffsets[0] = 0;
  for (int nN=1;nN<nNodes_subdomain_old+1;nN++)
    nodeStarOffsets[nN] = nodeStarOffsets[nN-1] + nodeStar[nN-1].size();
  valarray<int> nodeStarArray(nodeStarOffsets[nNodes_subdomain_old]);
  for (int nN=0,offset=0;nN<nNodes_subdomain_old;nN++)
    for (set<int>::iterator nN_star=nodeStar[nN].begin();nN_star!=nodeStar[nN].end();nN_star++,offset++)
      nodeStarArray[offset] = *nN_star;
  //find maximum number of nodes in any star
  int max_nNodeNeighbors_node=0;
  for (int nN=0;nN<nNodes_subdomain_old;nN++)
    max_nNodeNeighbors_node=max(max_nNodeNeighbors_node,nodeStarOffsets[nN+1]-nodeStarOffsets[nN]);
  //build connectivity data structures for PETSc
  PetscBool isInitialized;
  PetscInitialized(&isInitialized);
  PetscInt *nodeNeighborsOffsets_subdomain,*nodeNeighbors_subdomain,*weights_subdomain,*vertex_weights_subdomain;
  PetscReal *partition_weights;
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain_old+1),&nodeNeighborsOffsets_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain_old*max_nNodeNeighbors_node),&nodeNeighbors_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain_old*max_nNodeNeighbors_node),&weights_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nNodes_subdomain_old),&vertex_weights_subdomain);
  PetscMalloc(sizeof(PetscReal)*(size),&partition_weights);
  for (int sd=0;sd<size;sd++)
    partition_weights[sd] = 1.0/double(size);
  nodeNeighborsOffsets_subdomain[0] = 0;
  //I think we can simplify this now that nodeStarArray is local to the subdomain, could just use nodeStar instead of nodeStarArray
  for (int nN = 0,offset=0; nN < nNodes_subdomain_old; nN++)
    {
      for (int offset_subdomain = nodeStarOffsets[nN]; 
           offset_subdomain < nodeStarOffsets[nN+1]; 
           offset_subdomain++)
        {
          nodeNeighbors_subdomain[offset++] = nodeStarArray[offset_subdomain];
        }
      nodeNeighborsOffsets_subdomain[nN+1]=offset;
      sort(&nodeNeighbors_subdomain[nodeNeighborsOffsets_subdomain[nN]],&nodeNeighbors_subdomain[nodeNeighborsOffsets_subdomain[nN+1]]);
      //weight nodes by size of star
      int weight= (nodeNeighborsOffsets_subdomain[nN+1] - nodeNeighborsOffsets_subdomain[nN]);
      vertex_weights_subdomain[nN] = weight;
      for (int k=nodeNeighborsOffsets_subdomain[nN];k<nodeNeighborsOffsets_subdomain[nN+1];k++)
        weights_subdomain[k] = weight;
    }
  //
  //3. Generate new nodal partition using PETSc interface
  //
  Mat petscAdjacency;
  int nNodes_subdomain_max=0;
  MPI_Allreduce(&nNodes_subdomain_old,
		&nNodes_subdomain_max,
		1,
		MPI_INT,
		MPI_MAX,
		PROTEUS_COMM_WORLD);
  if (rank ==  0)
    std::cout<<"Max nNodes_subdomain "<<nNodes_subdomain_max<<" nNodes_global "<<nNodes_global<<std::endl;
  ierr = MatCreateMPIAdj(PROTEUS_COMM_WORLD,
			 nNodes_subdomain_old, 
			 nNodes_global,
			 nodeNeighborsOffsets_subdomain, 
			 nodeNeighbors_subdomain,
			 weights_subdomain,
			 &petscAdjacency);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  //const double max_rss_gb(0.75*3.25);//half max mem per  core  on topaz
  const double max_rss_gb(0.75*990.0);//half max mem per  core  on topaz
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done allocating MPIAdj");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  MatPartitioning petscPartition;
  ierr = MatPartitioningCreate(PROTEUS_COMM_WORLD,&petscPartition);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MatPartitioningSetAdjacency(petscPartition,petscAdjacency);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MatPartitioningSetFromOptions(petscPartition);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MatPartitioningSetVertexWeights(petscPartition,vertex_weights_subdomain);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MatPartitioningSetPartitionWeights(petscPartition,partition_weights);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  //get petsc index set that has the new subdomain number for each node
  IS nodePartitioningIS_new;
  ierr = MatPartitioningApply(petscPartition,&nodePartitioningIS_new);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  ierr = MatPartitioningDestroy(&petscPartition);CHKERRABORT(PROTEUS_COMM_WORLD, ierr); //gets petscAdjacency too I believe
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done applying partition");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);

  //determine the number of nodes per subdomain in new partitioning
  valarray<int> nNodes_subdomain_new(size);
  ISPartitioningCount(nodePartitioningIS_new,size,&nNodes_subdomain_new[0]);

  //need new offsets for subdomain to global numbering
  valarray<int> nodeOffsets_new(size+1);
  nodeOffsets_new[0] = 0;
  for (int sdN = 0; sdN < size; sdN++)
    {
      nodeOffsets_new[sdN+1] = nodeOffsets_new[sdN] + nNodes_subdomain_new[sdN];
    }

  //get the new node numbers for nodes on this subdomain
  IS nodeNumberingIS_subdomain_old2new;
  ISPartitioningToNumbering(nodePartitioningIS_new,&nodeNumberingIS_subdomain_old2new);
  //
  //try out of core
  //
  /* 
   * Set up file access property list with parallel I/O access
   */
  MPI_Info info  = MPI_INFO_NULL;
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, PROTEUS_COMM_WORLD, info);
  
  /*
   * Create a new file collectively and release property list identifier.
   */
  const char* H5FILE_NAME("mappings.h5");
  hid_t file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
   
  
  /*
   * Create the dataspace for the dataset.
   */
  hsize_t     dimsf[1];
  dimsf[0] = nNodes_global;
#define RANK   1
  hid_t filespace = H5Screate_simple(RANK, dimsf, NULL); 
  
  /*
   * Create the dataset with default properties and close filespace.
   */
  hid_t dset_id = H5Dcreate(file_id, "nodeNumbering_old2new", H5T_NATIVE_INT, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);
  
  /* 
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */
  hsize_t	count[1];	          /* hyperslab selection parameters */
  hsize_t	offset[1];
  count[0] = nNodes_subdomain_old;
  offset[0] = nodeOffsets_old[rank];
  hid_t memspace = H5Screate_simple(RANK, count, NULL);

  /*
   * Select hyperslab in the file.
   */
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

  /*
   * Initialize data buffer 
   */
  // data = (int *) malloc(sizeof(int)*count[0]*count[1]);
  // for (i=0; i < count[0]*count[1]; i++) {
  //   data[i] = mpi_rank + 10;
  // }
  const PetscInt* data;
  ISGetIndices(nodeNumberingIS_subdomain_old2new, &data);
  
  /*
   * Create property list for collective dataset write.
   */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  
  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                    plist_id, data);
  //free(data);
  ISRestoreIndices(nodeNumberingIS_subdomain_old2new, &data);
  /*
   * Close/release resources.
   */
  H5Dclose(dset_id);
  //
  //end try out of core
  //
  //collect new node numbers for whole mesh so that subdomain reordering and renumbering
  //can be done easily
  
  IS nodeNumberingIS_global_old2new;
  ISAllGather(nodeNumberingIS_subdomain_old2new,&nodeNumberingIS_global_old2new);
  const PetscInt * nodeNumbering_global_old2new;//needs restore call
  ISGetIndices(nodeNumberingIS_global_old2new,&nodeNumbering_global_old2new);
  //
  //test out of core
  //
  if (rank == 0)
    {
      hid_t       dataset_id;  /* identifiers */
      herr_t      status;
      int         dset_data[nNodes_global];
      
      /* Open an existing file. */
      //file_id = H5Fopen("mappings.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
      
      /* Open an existing dataset. */
      dataset_id = H5Dopen2(file_id, "/nodeNumbering_old2new", H5P_DEFAULT);
      
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                       dset_data);

      /* Close the dataset. */
      status = H5Dclose(dataset_id);
      
      for (int i=0;i<nNodes_global;i++)
        assert(nodeNumbering_global_old2new[i] == dset_data[i]);
      std::cout<<"==================out of core old2new is correct!===================="<<std::endl;
    }
  //
  //end test out of core
  //
  //reverse mapping for node numbers too
  //cek hack, not needed
  /*
  valarray<int> nodeNumbering_global_new2old(nNodes_global);
  for (int nN = 0; nN < nNodes_global; nN++)
    nodeNumbering_global_new2old[nodeNumbering_global_old2new[nN]] = nN;
  */
  PetscLogEventEnd(repartition_nodes_event,0,0,0,0);
  int receive_element_mask_event;
  PetscLogEventRegister("Recv. ele mask",0,&receive_element_mask_event);
  PetscLogEventBegin(receive_element_mask_event,0,0,0,0);
  //
  //4. To build subdomain meshes, go through and collect elements containing
  //   the locally owned nodes. Assign processor ownership of elements 
  //  
  PetscLogEventEnd(receive_element_mask_event,0,0,0,0);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done with masks");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  int build_subdomains_reread_elements_event;
  PetscLogEventRegister("Reread eles",0,&build_subdomains_reread_elements_event);
  PetscLogEventBegin(build_subdomains_reread_elements_event,0,0,0,0);
  //
  //mark the unmarked elements on this subdomain and store element numbers (in old numbering)
  //
  //We need to find the elements in the support of nodes in the new
  //partitioning, so we'll need to re-read the elements file to get
  //the the elements for the nodes in the new partitioning. We will be collecting OLD element numbers.
  std::ifstream elementFile2(elementFileName.c_str());

  if (!elementFile2.good())
    {
      std::cerr<<"cannot open Tetgen elements file"
               <<elementFileName<<std::endl;
      failed = true;
      return failed;
    }
  elementFile2 >> eatcomments >> nElements_global >> nNodesPerSimplex >> hasElementMarkers >> eatline;
  assert(nElements_global > 0);
  assert(nNodesPerSimplex == simplexDim); 
  set<int> elements_subdomain_owned;
  vector<int> element_nodes_new(4);
  int element_nodes_new_array[4];
  vector<set<int> > nodeElementsStar(nNodes_subdomain_new[rank]);
  vector<set<int> > nodeStarNew(nNodes_subdomain_new[rank]);
  map<int,vector<int> > elementNodesArrayMap;
  map<int,long int> elementMaterialTypesMap;
  map<NodeTuple<3>,ElementNeighbors> elementBoundaryElementsMap;
  map<NodeTuple<2>,set<pair<int,int> > > edgeElementsMap;
  //note any element index containers are in the old element numbering
  for (int ie = 0; ie < nElements_global; ie++)
    {
      int ne, nv, elementId(0);
      long double elementId_double;
      elementFile2 >> eatcomments >> ne;
      ne -= indexBase;
      assert(0 <= ne && ne < nElements_global && elementFile.good());
      for (int iv = 0; iv < simplexDim; iv++)
        {
          elementFile2 >> nv ; 
          nv -= indexBase;
          assert(0 <= nv && nv < nNodes_global);
          element_nodes_old[iv] = nv;
          element_nodes_new[iv] = nodeNumbering_global_old2new[nv];
          element_nodes_new_array[iv] = element_nodes_new[iv];
        }
      NodeTuple<4> nodeTuple(element_nodes_new_array);      
      for (int iv = 0; iv < simplexDim; iv++)
        {
          int nN_star_new = element_nodes_new[iv];
          bool inSubdomain=false;
          if (nN_star_new >= nodeOffsets_new[rank] && nN_star_new < nodeOffsets_new[rank+1])
            {
              inSubdomain = true;
              //add all the element boundaries of this element
              for (int ebN=0;ebN < 4 ; ebN++)
                {
                  int nodes[3] = { element_nodes_new[(ebN+1) % 4], 
                                   element_nodes_new[(ebN+2) % 4], 
                                   element_nodes_new[(ebN+3) % 4]};
                  NodeTuple<3> nodeTuple(nodes);
                  if(elementBoundaryElementsMap.find(nodeTuple) != elementBoundaryElementsMap.end())
                    {
                      if (elementBoundaryElementsMap[nodeTuple].right == -1 && ne != elementBoundaryElementsMap[nodeTuple].left)
                        {
                          elementBoundaryElementsMap[nodeTuple].right=ne;
                          elementBoundaryElementsMap[nodeTuple].right_ebN_element=ebN;
                        }
                    }
                  else
                    {
                      elementBoundaryElementsMap[nodeTuple] = ElementNeighbors(ne,ebN);
                    }
                }
              //add all the edges of this element
              for (int nNL=0,edN=0;nNL < 4 ; nNL++)
                for(int nNR=nNL+1;nNR < 4;nNR++,edN++)
                  {
                    int nodes[2] = { element_nodes_new[nNL], 
                                     element_nodes_new[nNR]}; 
                    NodeTuple<2> nodeTuple(nodes);
                    edgeElementsMap[nodeTuple].insert(pair<int,int>(ne,edN));
                  }
              //add all the nodes to the node star
              int nN_star_new_subdomain = nN_star_new - nodeOffsets_new[rank];
              nodeElementsStar[nN_star_new_subdomain].insert(ne);
              for (int jv = 0; jv < simplexDim; jv++)
                {
                  if (iv != jv)
                    {
                      int nN_point_new = element_nodes_new[jv];
                      nodeStarNew[nN_star_new_subdomain].insert(nN_point_new);
                    }
                }
            }
          if (inSubdomain)
            {
              elementNodesArrayMap[ne] = element_nodes_new;
            }
        }
      if (elementNodesArrayMap.find(ne) != elementNodesArrayMap.end())//this element contains a node owned by this subdomain
        {
          if (nodeTuple.nodes[1] >= nodeOffsets_new[rank] && nodeTuple.nodes[1] < nodeOffsets_new[rank+1])
            elements_subdomain_owned.insert(ne);
          if (hasElementMarkers > 0)
            {
              elementFile2 >> elementId_double;
              elementId = static_cast<long int>(elementId_double);
              elementMaterialTypesMap[ne] = elementId;
            }
        }
      elementFile2 >> eatline;
    }
  elementFile2.close();
  int nElements_owned_subdomain(elements_subdomain_owned.size()),
    nElements_owned_new=0;
  MPI_Allreduce(&nElements_owned_subdomain,&nElements_owned_new,1,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  assert(nElements_owned_new == nElements_global);
  PetscLogEventEnd(build_subdomains_reread_elements_event,0,0,0,0);
  int build_subdomains_send_marked_elements_event;
  PetscLogEventRegister("Mark/send eles",0,&build_subdomains_send_marked_elements_event);
  PetscLogEventBegin(build_subdomains_send_marked_elements_event,0,0,0,0);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done marking elements");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  //
  //done with the element file
  //
  //construct compact nodeElementsArray
  valarray<int> nodeElementOffsets(nNodes_subdomain_new[rank]+1);
  nodeElementOffsets[0] = 0;
  for (int nN = 0; nN < nNodes_subdomain_new[rank]; nN++)
    nodeElementOffsets[nN+1] = nodeElementOffsets[nN]+nodeElementsStar[nN].size();
  valarray<int> nodeElementsArray(nodeElementOffsets[nNodes_subdomain_new[rank]]);
  for (int nN=0,offset=0; nN < nNodes_subdomain_new[rank]; nN++)
    {
      for (set<int>::iterator eN_star = nodeElementsStar[nN].begin(); eN_star != nodeElementsStar[nN].end();
           eN_star++,offset++)
        {
          nodeElementsArray[offset] = *eN_star;
        }
    }
  //construct compact nodeStarArray
  valarray<int> nodeStarOffsetsNew(nNodes_subdomain_new[rank]+1);
  nodeStarOffsetsNew[0] = 0;
  for (int nN=1;nN<nNodes_subdomain_new[rank]+1;nN++)
    nodeStarOffsetsNew[nN] = nodeStarOffsetsNew[nN-1] + nodeStarNew[nN-1].size();
  valarray<int> nodeStarArrayNew(nodeStarOffsetsNew[nNodes_subdomain_new[rank]]);
  for (int nN=0,offset=0;nN<nNodes_subdomain_new[rank];nN++)
    {
      for (set<int>::iterator nN_star=nodeStarNew[nN].begin();nN_star!=nodeStarNew[nN].end();nN_star++,offset++)
        {
          nodeStarArrayNew[offset] = *nN_star;
        }
    }
  PetscLogEventEnd(build_subdomains_send_marked_elements_event,0,0,0,0);
  int build_subdomains_global_numbering_elements_event;
  PetscLogEventRegister("Global ele nmbr",0,&build_subdomains_global_numbering_elements_event);
  PetscLogEventBegin(build_subdomains_global_numbering_elements_event,0,0,0,0);
  //
  //5. Generate global element numbering corresponding to new subdomain ownership
  //
  valarray<int> nElements_subdomain_new(size),
    elementOffsets_new(size+1);
  for (int sdN = 0; sdN < size; sdN++)
    {
      if (sdN == rank)
        nElements_subdomain_new[sdN] = int(elements_subdomain_owned.size());
      else
        nElements_subdomain_new[sdN] = 0;
    }
  valarray<int> nElements_subdomain_new_send = nElements_subdomain_new;
  MPI_Allreduce(&nElements_subdomain_new_send[0],&nElements_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  //construct new offsets for elements
  elementOffsets_new[0] = 0;
  for (int sdN = 0; sdN < size; sdN++)
    elementOffsets_new[sdN+1] = elementOffsets_new[sdN] + nElements_subdomain_new[sdN];  
  //map to old element numbering
  valarray<int> elementNumbering_subdomain_new2old(elements_subdomain_owned.size());
  set<int>::iterator eN_ownedp = elements_subdomain_owned.begin();
  for (int eN = 0; eN < int(elements_subdomain_owned.size()); eN++,eN_ownedp++)
    {
      elementNumbering_subdomain_new2old[eN] = *eN_ownedp;
    }
  //use Petsc IS to get global new2old numbering
  IS elementNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,elements_subdomain_owned.size(),&elementNumbering_subdomain_new2old[0],PETSC_COPY_VALUES,
                  &elementNumberingIS_subdomain_new2old);
  IS elementNumberingIS_global_new2old;
  ISAllGather(elementNumberingIS_subdomain_new2old,&elementNumberingIS_global_new2old);
  
  const PetscInt *elementNumbering_global_new2old;//needs to be restored
  ISGetIndices(elementNumberingIS_global_new2old,&elementNumbering_global_new2old);
  //construct reverse mapping
  valarray<int> elementNumbering_global_old2new(nElements_global);
  for (int eN = 0; eN < nElements_global; eN++)
    {
      elementNumbering_global_old2new[elementNumbering_global_new2old[eN]] = eN;
    }
  PetscLogEventEnd(build_subdomains_global_numbering_elements_event,0,0,0,0);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done allocating element numbering new2old/old2new");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  int build_subdomains_faces_event;
  PetscLogEventRegister("Subd faces",0,&build_subdomains_faces_event);
  PetscLogEventBegin(build_subdomains_faces_event,0,0,0,0);
  //
  //4b,5b. repeat process to build global face (elementBoundary) numbering
  //
  //first read element boundaries to create nodeElementBoundariesArray
  //for all element boundaries on this subdomain, which we'll use to
  //grab element boundaries from the bit array

  std::ifstream elementBoundaryFile(elementBoundaryFileName.c_str());
  
  if (!elementBoundaryFile.good())
    {
      std::cerr<<"cannot open Tetgen face file "
               <<elementBoundaryFileName<<std::endl;
      failed = true;
      return failed;
    }

  bool hasElementBoundaryMarkers = false;
  int nElementBoundaries_global;
  int ihasElementBoundaryMarkers(0);
  elementBoundaryFile >> eatcomments >> nElementBoundaries_global >> ihasElementBoundaryMarkers >> eatline ;
  assert(nElementBoundaries_global > 0);
  if (ihasElementBoundaryMarkers > 0)
    {
      hasElementBoundaryMarkers = true;
    }
  newMesh.nElementBoundaries_global = nElementBoundaries_global;
  //note, these will be in the new element numbering
  set<int> elementBoundaries_subdomain_owned;
  vector<set<int> > nodeElementBoundariesStar(nNodes_subdomain_new[rank]);
  map<int,int> elementBoundaryMaterialTypesMap;
  map<int,vector<int> > elementBoundariesMap;
  set<int> supportedElementBoundaries;
  for (int ieb = 0; ieb < nElementBoundaries_global; ieb++)
    {
      int neb,nn0,nn1,nn2; int ebId(0);
      elementBoundaryFile >> eatcomments >> neb >> nn0 >> nn1 >> nn2;
      if (ihasElementBoundaryMarkers > 0)
        elementBoundaryFile >> ebId;
      neb -= indexBase;
      nn0 -= indexBase;
      nn1 -= indexBase;
      nn2 -= indexBase;
      assert(0 <= neb && neb < nElementBoundaries_global && elementBoundaryFile.good());
      //grab the element boundaries for the node if the node is owned by the subdomain
      //this will miss the element boundaries on the "outside boundary" of the star, which will grab later
      int nn0_new = nodeNumbering_global_old2new[nn0];
      if (nn0_new >= nodeOffsets_new[rank] && nn0_new < nodeOffsets_new[rank+1])
        {
          nodeElementBoundariesStar[nn0_new-nodeOffsets_new[rank]].insert(neb);
          supportedElementBoundaries.insert(neb);
        }
      int nn1_new = nodeNumbering_global_old2new[nn1];
      if (nn1_new >= nodeOffsets_new[rank] && nn1_new < nodeOffsets_new[rank+1])
        {
          nodeElementBoundariesStar[nn1_new-nodeOffsets_new[rank]].insert(neb);
          supportedElementBoundaries.insert(neb);
        }
      int nn2_new = nodeNumbering_global_old2new[nn2];
      if (nn2_new >= nodeOffsets_new[rank] && nn2_new < nodeOffsets_new[rank+1])
        {
          nodeElementBoundariesStar[nn2_new-nodeOffsets_new[rank]].insert(neb);
          supportedElementBoundaries.insert(neb);
        }
      int nodes[3] = {nn0_new,nn1_new,nn2_new};
      NodeTuple<3> nodeTuple(nodes);      
      elementBoundaryFile >> eatline;
      if (elementBoundaryElementsMap.find(nodeTuple) != elementBoundaryElementsMap.end())//this element boundary is on an element in the subdomain
        {
          if (nodeTuple.nodes[1] >= nodeOffsets_new[rank] && nodeTuple.nodes[1] < nodeOffsets_new[rank+1])
            elementBoundaries_subdomain_owned.insert(neb);
          if (ihasElementBoundaryMarkers > 0)
            elementBoundaryMaterialTypesMap[neb]=ebId;
          int eN_left = elementNumbering_global_old2new[elementBoundaryElementsMap[nodeTuple].left];
          if (elementBoundariesMap.find(eN_left) != elementBoundariesMap.end())
            {
              elementBoundariesMap[eN_left][elementBoundaryElementsMap[nodeTuple].left_ebN_element] = neb;
            }
          else
            {
              //initialize
              vector<int> elementBoundaries_element(4,-1);
              elementBoundariesMap[eN_left] = elementBoundaries_element;
              //assign
              elementBoundariesMap[eN_left][elementBoundaryElementsMap[nodeTuple].left_ebN_element] = neb;
            }
          if (elementBoundaryElementsMap[nodeTuple].right >= 0)
            {
              int eN_right = elementNumbering_global_old2new[elementBoundaryElementsMap[nodeTuple].right];
              if (elementBoundariesMap.find(eN_right) != elementBoundariesMap.end())
                {
                  elementBoundariesMap[eN_right][elementBoundaryElementsMap[nodeTuple].right_ebN_element] = neb;
                }
              else
                {
                  //initialize
                  vector<int> elementBoundaries_element(4,-1);
                  elementBoundariesMap[eN_right] = elementBoundaries_element;
                  //assign
                  elementBoundariesMap[eN_right][elementBoundaryElementsMap[nodeTuple].right_ebN_element] = neb;
                }
            }
        }
    }
  //done reading element boundaries
  elementBoundaryFile.close();
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done reading element boundaries");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  int nElementBoundaries_owned_subdomain=elementBoundaries_subdomain_owned.size(),
    nElementBoundaries_owned_new=0;
  MPI_Allreduce(&nElementBoundaries_owned_subdomain,&nElementBoundaries_owned_new,1,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  assert(nElementBoundaries_owned_new == nElementBoundaries_global);

  //now get the element boundaries on the outside of the star
  for (map<int,vector<int> >::iterator elementBoundariesp=elementBoundariesMap.begin();
       elementBoundariesp!=elementBoundariesMap.end();
       elementBoundariesp++)
    {
      //loop over the nodes of this element for the owned nodes
      for (int iv=0;iv<4;iv++)
        {
          //the elementNodesArrayMap is in the old element numbering while the elementBoundariesMap is in the new element numbering
          int nN_global = elementNodesArrayMap[elementNumbering_global_new2old[elementBoundariesp->first]][iv];
          if (nN_global >= nodeOffsets_new[rank] && nN_global < nodeOffsets_new[rank+1])
            {
              //add all the faces to this node star
              for(int eb=0;eb<4;eb++)
                {
                  nodeElementBoundariesStar[nN_global-nodeOffsets_new[rank]].insert(elementBoundariesp->second[eb]);
                }
            }
        }
    }
  //build compact structures for nodeElementBoundariesArray
  valarray<int> nodeElementBoundaryOffsets(nNodes_subdomain_new[rank]+1);
  nodeElementBoundaryOffsets[0] = 0;
  for (int nN = 0; nN < nNodes_subdomain_new[rank]; nN++)
    nodeElementBoundaryOffsets[nN+1] = nodeElementBoundaryOffsets[nN]+nodeElementBoundariesStar[nN].size();
  valarray<int> nodeElementBoundariesArray(nodeElementBoundaryOffsets[nNodes_subdomain_new[rank]]);
  for (int nN=0,offset=0; nN < nNodes_subdomain_new[rank]; nN++)
    {
      for (set<int>::iterator ebN_star = nodeElementBoundariesStar[nN].begin(); ebN_star != nodeElementBoundariesStar[nN].end();
           ebN_star++,offset++)
        {
          nodeElementBoundariesArray[offset] = *ebN_star;
        }
    }

  //get the number of elementBoundaries owned on each processor
  valarray<int> nElementBoundaries_subdomain_new(size),
    elementBoundaryOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nElementBoundaries_subdomain_new[sdN] = elementBoundaries_subdomain_owned.size();
    else
      nElementBoundaries_subdomain_new[sdN] = 0;
  valarray<int> nElementBoundaries_subdomain_new_send=nElementBoundaries_subdomain_new;
  MPI_Allreduce(&nElementBoundaries_subdomain_new_send[0],&nElementBoundaries_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  elementBoundaryOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    elementBoundaryOffsets_new[sdN+1] = elementBoundaryOffsets_new[sdN]+nElementBoundaries_subdomain_new[sdN];
  //
  //Now as with elements and nodes build a global face numbering
  //resetting the face-based information is a little different since much of this is currently built below based
  //on the element and node information
  //
  valarray<int> elementBoundaryNumbering_new2old(elementBoundaries_subdomain_owned.size());
  set<int>::iterator ebN_ownedp=elementBoundaries_subdomain_owned.begin();
  for (int ebN=0;ebN<int(elementBoundaries_subdomain_owned.size());ebN++)
    {
      elementBoundaryNumbering_new2old[ebN] = *ebN_ownedp++;
    }
  IS elementBoundaryNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,elementBoundaries_subdomain_owned.size(),&elementBoundaryNumbering_new2old[0],PETSC_COPY_VALUES,&elementBoundaryNumberingIS_subdomain_new2old);
  IS elementBoundaryNumberingIS_global_new2old;
  ISAllGather(elementBoundaryNumberingIS_subdomain_new2old,&elementBoundaryNumberingIS_global_new2old);
  const PetscInt *elementBoundaryNumbering_global_new2old;
  valarray<int> elementBoundaryNumbering_global_old2new(newMesh.nElementBoundaries_global);
  ISGetIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Allocating elementBoudnary old2new/new2old");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  for (int ebN=0;ebN<newMesh.nElementBoundaries_global;ebN++)
    {
      elementBoundaryNumbering_global_old2new[elementBoundaryNumbering_global_new2old[ebN]] = ebN;
    }
  ISRestoreIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);
  ISDestroy(&elementBoundaryNumberingIS_subdomain_new2old);
  ISDestroy(&elementBoundaryNumberingIS_global_new2old);
  PetscLogEventEnd(build_subdomains_faces_event,0,0,0,0);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done allocating elementBoudnary old2new/new2old");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  int build_subdomains_edges_event;
  PetscLogEventRegister("Subd edges",0,&build_subdomains_edges_event);
  PetscLogEventBegin(build_subdomains_edges_event,0,0,0,0);
  
  //
  //4c,5c. Repeate the process for edges
  //
  std::ifstream edgeFile(edgeFileName.c_str());
  
  if (!edgeFile.good())
    {
      std::cerr<<"cannot open Tetgen edge file"
               <<edgeFileName<<std::endl;
      failed = true;
      return failed;
    }

  bool hasEdgeMarkers = false;
  int nEdges_global;
  int ihasEdgeMarkers(0);
  edgeFile >> eatcomments >> nEdges_global >> eatline;// edge file doesn currently contain markers >> ihasEdgeMarkers >> eatline ;
  assert(nEdges_global > 0);
  if (ihasEdgeMarkers > 0)
    {
      hasEdgeMarkers = true;
    }
  newMesh.nEdges_global = nEdges_global;
  set<int> edges_subdomain_owned;
  vector<set<int> > nodeEdgesStar(nNodes_subdomain_new[rank]);
  map<int,int> edgeMaterialTypesMap;
  map<int,vector<int> > elementEdgesMap;
  map<int,pair<int,int> > edgeNodesMap;
  set<int> supportedEdges;
  for (int ied = 0; ied < nEdges_global; ied++)
    {
      int ned,nn0,nn1; int edId(0);
      edgeFile >> eatcomments >> ned >> nn0 >> nn1;
      if (ihasEdgeMarkers > 0)
        edgeFile >> edId;
      ned -= indexBase;
      nn0 -= indexBase;
      nn1 -= indexBase;
      assert(0 <= ned && ned < nEdges_global && edgeFile.good());
      int nn0_new = nodeNumbering_global_old2new[nn0];
      if (nn0_new >= nodeOffsets_new[rank] && nn0_new < nodeOffsets_new[rank+1])
        {
          nodeEdgesStar.at(nn0_new-nodeOffsets_new[rank]).insert(ned);
          supportedEdges.insert(ned);
        }
      int nn1_new = nodeNumbering_global_old2new[nn1];
      if (nn1_new >= nodeOffsets_new[rank] && nn1_new < nodeOffsets_new[rank+1])
        {
          nodeEdgesStar.at(nn1_new-nodeOffsets_new[rank]).insert(ned);
          supportedEdges.insert(ned);
        }
      int nodes[2] = {nn0_new,nn1_new};
      NodeTuple<2> nodeTuple(nodes);      
      edgeFile >> eatline;
      if (edgeElementsMap.find(nodeTuple) != edgeElementsMap.end())//this edge is on an element in the subdomain
        {
          if (nodeTuple.nodes[0] >= nodeOffsets_new[rank] && nodeTuple.nodes[0] < nodeOffsets_new[rank+1])
            edges_subdomain_owned.insert(ned);
          //pick up all the edges on the subdomain and store their nodes
          edgeNodesMap[ned].first = nodeTuple.nodes[0];
          edgeNodesMap[ned].second = nodeTuple.nodes[1];
          if (ihasEdgeMarkers > 0)
            edgeMaterialTypesMap[ned]=edId;
          for (set<pair<int,int> >::iterator elementp=edgeElementsMap[nodeTuple].begin();
               elementp != edgeElementsMap[nodeTuple].end();
               elementp++)
            {
              int eN = elementNumbering_global_old2new[elementp->first];
              if (elementEdgesMap.find(eN) != elementEdgesMap.end())
                {
                  elementEdgesMap[eN][elementp->second] = ned;
                }
              else
                {
                  std::vector<int> init(6,-1);
                  elementEdgesMap[eN] = init;
                  elementEdgesMap[eN][elementp->second] = ned;
                }
            }
        }
    }//end iv
  edgeFile.close();
  int nEdges_owned_subdomain=edges_subdomain_owned.size(),
    nEdges_owned_new=0;
  MPI_Allreduce(&nEdges_owned_subdomain,&nEdges_owned_new,1,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  assert(nEdges_owned_new == nEdges_global);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done reading edges");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  //done with edge file
  //
  //just as with faces, we need to add edges along outer boundaries of star
  //not sure if we need to collect nodeEdges star above anymore, since we're doing this
  //
  for (map<int,vector<int> >::iterator edgesp=elementEdgesMap.begin();
       edgesp!=elementEdgesMap.end();
       edgesp++)
    {
      //loop over the nodes of this element for the owned nodes
      for (int iv=0;iv<4;iv++)
        {
          //the elementNodesArrayMap is in the old elemetn numbering while the elementEdgesMap is in the new element numbering
          int nN_global = elementNodesArrayMap[elementNumbering_global_new2old[edgesp->first]][iv];
          if (nN_global >= nodeOffsets_new[rank] && nN_global < nodeOffsets_new[rank+1])
            {
              //add all the edges to this node star
              for(int ed=0;ed<6;ed++)
                {
                  nodeEdgesStar.at(nN_global-nodeOffsets_new[rank]).insert(edgesp->second[ed]);
                }
            }
        }
    }
  //build compact data structures for nodeEdgesArray
  valarray<int> nodeEdgeOffsets(nNodes_subdomain_new[rank]+1);
  nodeEdgeOffsets[0] = 0;
  for (int nN = 0; nN < nNodes_subdomain_new[rank]; nN++)
    nodeEdgeOffsets[nN+1] = nodeEdgeOffsets[nN]+nodeEdgesStar.at(nN).size();
  valarray<int> nodeEdgesArray(nodeEdgeOffsets[nNodes_subdomain_new[rank]]);
  for (int nN=0,offset=0; nN < nNodes_subdomain_new[rank]; nN++)
    {
      for (set<int>::iterator edN_star = nodeEdgesStar.at(nN).begin(); 
           edN_star != nodeEdgesStar.at(nN).end();
           edN_star++,offset++)
        {
          nodeEdgesArray[offset] = *edN_star;
        }
    }
  newMesh.subdomainp->nEdges_global = edgeNodesMap.size();
  //get the number of edges on each processor
  valarray<int> nEdges_subdomain_new(size),
    edgeOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nEdges_subdomain_new[sdN] = edges_subdomain_owned.size();
    else
      nEdges_subdomain_new[sdN] = 0;
  valarray<int> nEdges_subdomain_new_send=nEdges_subdomain_new;
  MPI_Allreduce(&nEdges_subdomain_new_send[0],&nEdges_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  //
  //construct new offsets for owned edges
  //
  edgeOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    edgeOffsets_new[sdN+1] = edgeOffsets_new[sdN]+nEdges_subdomain_new[sdN];
  //
  //Now as with elementBoundaries, build a global face numbering
  //resetting the edge based information is a little different since much of this is currently built below based
  //on the element and node information
  //
  valarray<int> edgeNumbering_new2old(edges_subdomain_owned.size());
  set<int>::iterator edN_ownedp=edges_subdomain_owned.begin();
  for (int edN=0;edN<int(edges_subdomain_owned.size());edN++,edN_ownedp++)
    {
      edgeNumbering_new2old[edN] = *edN_ownedp;
    }
  IS edgeNumberingIS_subdomain_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,edges_subdomain_owned.size(),&edgeNumbering_new2old[0],PETSC_COPY_VALUES,&edgeNumberingIS_subdomain_new2old);
  IS edgeNumberingIS_global_new2old;
  ISAllGather(edgeNumberingIS_subdomain_new2old,&edgeNumberingIS_global_new2old);
  const PetscInt *edgeNumbering_global_new2old;
  valarray<int> edgeNumbering_global_old2new(newMesh.nEdges_global);
  ISGetIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Setting edgeNumering old2new/new2old");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  for (int edN=0;edN<newMesh.nEdges_global;edN++)
    {
      edgeNumbering_global_old2new[edgeNumbering_global_new2old[edN]] = edN;
    }
  ISRestoreIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);
  ISDestroy(&edgeNumberingIS_subdomain_new2old);
  ISDestroy(&edgeNumberingIS_global_new2old);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done allocating edgeNumering old2new/new2old");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  //
  //6. Figure out what is in the node stars but not locally owned, create ghost information
  //
  
  set<int> elements_overlap,nodes_overlap,elementBoundaries_overlap,edges_overlap;
  for (int nN = 0; nN < nNodes_subdomain_new[rank]; nN++)
    {
      //nodes
      for (int offset = nodeStarOffsetsNew[nN];offset<nodeStarOffsetsNew[nN+1];offset++)
        {
          int nN_point_global = nodeStarArrayNew[offset];
          bool offproc = nN_point_global < nodeOffsets_new[rank] || nN_point_global >= nodeOffsets_new[rank+1];
          if (offproc)
            nodes_overlap.insert(nN_point_global);
        }
      //elements
      for (int eN_star_offset = nodeElementOffsets[nN]; 
           eN_star_offset < nodeElementOffsets[nN+1]; eN_star_offset++)
        {
          int eN_star_old = nodeElementsArray[eN_star_offset];
          int eN_star_new = elementNumbering_global_old2new[eN_star_old];
          bool offproc = eN_star_new >= elementOffsets_new[rank+1] || eN_star_new < elementOffsets_new[rank];
          if (offproc)
            elements_overlap.insert(eN_star_new);
        }
      //element boundaries
      for (int ebN_star_offset = nodeElementBoundaryOffsets[nN]; 
           ebN_star_offset < nodeElementBoundaryOffsets[nN+1]; ebN_star_offset++)
        {
          int ebN_star_old = nodeElementBoundariesArray[ebN_star_offset];
          int ebN_star_new = elementBoundaryNumbering_global_old2new[ebN_star_old];
          bool offproc = ebN_star_new >= elementBoundaryOffsets_new[rank+1] || ebN_star_new < elementBoundaryOffsets_new[rank];
          if (offproc)
            elementBoundaries_overlap.insert(ebN_star_new);
        }
      //edges in overlap
      for (int edN_star_offset = nodeEdgeOffsets[nN]; 
           edN_star_offset < nodeEdgeOffsets[nN+1]; edN_star_offset++)
        {
          int edN_star_old = nodeEdgesArray[edN_star_offset];
          int edN_star_new = edgeNumbering_global_old2new[edN_star_old];
          bool offproc = edN_star_new >= edgeOffsets_new[rank+1] || edN_star_new < edgeOffsets_new[rank];
          if (offproc)
            edges_overlap.insert(edN_star_new);
        }
    }//nodes on this processor
  elementNumbering_global_old2new.resize(0);
  //cek debugging, edge overlap seems to be messed up. Check global node tuples of edges vs global edge numbers
  assert(edges_overlap.size() + nEdges_subdomain_new[rank] == edgeNodesMap.size());
  //
  //enumerate the overlap
  //
  int nN_subdomain = nNodes_subdomain_new[rank];
  map<int,int> nodes_overlap_global2subdomainMap;
  for (set<int>::iterator nN_globalp=nodes_overlap.begin();nN_globalp != nodes_overlap.end(); nN_globalp++,nN_subdomain++)
    nodes_overlap_global2subdomainMap[*nN_globalp] = nN_subdomain;

  // map<int,int> elements_overlap_global2subdomainMap;
  // for (set<int>::iterator eN_globalp=elements_overlap.begin();eN_globalp != elements_overlap.end(); eN_globalp++,eN_subdomain++)
  //   elements_overlap_global2subdomainMap[*eN_globalp] = eN_subdomain;

  // int ebN_subdomain = nElementBoundaries_subdomain_new[rank];
  // map<int,int> elementBoundaries_overlap_global2subdomainMap;
  // for (set<int>::iterator ebN_globalp=elementBoundaries_overlap.begin();ebN_globalp != elementBoundaries_overlap.end(); ebN_globalp++,ebN_subdomain++)
  //   elementBoundaries_overlap_global2subdomainMap[*ebN_globalp] = ebN_subdomain;

  // int edN_subdomain = nEdges_subdomain_new[rank];
  // map<int,int> edges_overlap_global2subdomainMap;
  // for (set<int>::iterator edN_globalp=edges_overlap.begin();edN_globalp != edges_overlap.end(); edN_globalp++,edN_subdomain++)
  //   edges_overlap_global2subdomainMap[*edN_globalp] = edN_subdomain;
  //
  //7. add any addtional overlap, skip for now
  //

  PetscLogEventEnd(build_subdomains_edges_event,0,0,0,0);
  int build_subdomains_renumber_event;
  PetscLogEventRegister("Subd's renumber",0,&build_subdomains_renumber_event);
  PetscLogEventBegin(build_subdomains_renumber_event,0,0,0,0);
  //
  //8. Build subdomain meshes in new numbering, assumes memory not allocated in subdomain mesh
  //   
  if(rank==0){
    std::cerr<<"USER WARNING: In order to avoid a segmentation fault, you need to have supplied the 'f' flag to the triangleOptions input."<<std::endl;
    std::cerr<<"USER WARNING: In order to avoid an edge assertion error, you need to have supplied the 'ee' flag to the triangleOptions input."<<std::endl; 
  }

  if (newMesh.subdomainp == NULL)
    newMesh.subdomainp = new Mesh();
  newMesh.subdomainp->nElements_global = nElements_subdomain_new[rank] + elements_overlap.size();
  newMesh.subdomainp->nNodes_global    = nNodes_subdomain_new[rank] + nodes_overlap.size();
  newMesh.subdomainp->nElementBoundaries_global = nElementBoundaries_subdomain_new[rank]+elementBoundaries_overlap.size();
  assert(int(edges_subdomain_owned.size()+edges_overlap.size()) == newMesh.subdomainp->nEdges_global);
  //newMesh.subdomainp->nEdges_global   = edges_subdomain_owned.size()+edges_overlap.size();
  newMesh.subdomainp->nNodes_element   = newMesh.nNodes_element;
  newMesh.subdomainp->nNodes_elementBoundary = newMesh.nNodes_elementBoundary;
  newMesh.subdomainp->nElementBoundaries_element = newMesh.nElementBoundaries_element;
  //subdomain 2 global mappings (including ghost info)
  valarray<int> nodeNumbering_subdomain2global(newMesh.subdomainp->nNodes_global);
  valarray<int> elementNumbering_subdomain2global(newMesh.subdomainp->nElements_global);
  valarray<int> elementBoundaryNumbering_subdomain2global(newMesh.subdomainp->nElementBoundaries_global);
  valarray<int> edgeNumbering_subdomain2global(newMesh.subdomainp->nEdges_global);
  map<int,int> nodeNumbering_global2subdomainMap;
  map<int,int> elementBoundaryNumbering_global2subdomainMap;
  map<int,int> edgeNumbering_global2subdomainMap;
  newMesh.subdomainp->nodeArray = new double[newMesh.subdomainp->nNodes_global*3];
  newMesh.subdomainp->nodeMaterialTypes = new int[newMesh.subdomainp->nNodes_global];
  //
  //now finally finish reading node coordinates and node flags
  //
  for (int iv = 0; iv < nNodes_global; iv++)
    {
      int nv; double x,y,z; int nodeId(0);
      vertexFile >> eatcomments >> nv >> x >> y >> z;
      if (hasVertexMarkers > 0)
        vertexFile >> nodeId;
      nv -= indexBase;
      assert(0 <= nv && nv < nNodes_global && vertexFile.good());
      int nN_global_new = nodeNumbering_global_old2new[nv];
      //local
      if (nN_global_new >= nodeOffsets_new[rank] && nN_global_new < nodeOffsets_new[rank+1])
        {
          int nv_subdomain_new = nN_global_new - nodeOffsets_new[rank];
          nodeNumbering_subdomain2global[nv_subdomain_new] = nN_global_new;
          nodeNumbering_global2subdomainMap[nN_global_new] = nv_subdomain_new;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 0] = x;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 1] = y;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 2] = z;
          if (hasVertexMarkers > 0)
            newMesh.subdomainp->nodeMaterialTypes[nv_subdomain_new] = nodeId;
        }
      //overlap
      if (nodes_overlap.count(nN_global_new) == 1)
        {
          int nv_subdomain_new = nodes_overlap_global2subdomainMap[nN_global_new];
          nodeNumbering_subdomain2global[nv_subdomain_new] = nN_global_new;
          nodeNumbering_global2subdomainMap[nN_global_new] = nv_subdomain_new;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 0] = x;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 1] = y;
          newMesh.subdomainp->nodeArray[vertexDim*nv_subdomain_new + 2] = z;
          if (hasVertexMarkers > 0)
            newMesh.subdomainp->nodeMaterialTypes[nv_subdomain_new] = nodeId;
        }
      vertexFile >> eatline;
    }//end iv
  vertexFile.close();
  ISRestoreIndices(nodeNumberingIS_global_old2new,&nodeNumbering_global_old2new);
  ISDestroy(&nodePartitioningIS_new);
  ISDestroy(&nodeNumberingIS_subdomain_old2new);
  ISDestroy(&nodeNumberingIS_global_old2new);
  //done with vertex file (and all file reads at this point)
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done reading vertices");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);

  newMesh.subdomainp->elementNodesArray = new int[newMesh.subdomainp->nElements_global*newMesh.subdomainp->nNodes_element];
  newMesh.subdomainp->elementMaterialTypes = new int[newMesh.subdomainp->nElements_global];
  //
  //elements
  //
  //locally owned
  //
  for (int eN = 0; eN < nElements_subdomain_new[rank]; eN++)
    {
      int eN_global_new = elementOffsets_new[rank] + eN;
      int eN_global_old = elementNumbering_global_new2old[eN_global_new];
      elementNumbering_subdomain2global[eN] = eN_global_new;
      newMesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypesMap[eN_global_old];
      for (int nN =  0; nN < newMesh.subdomainp->nNodes_element; nN++)
        {
          int nN_global_new = elementNodesArrayMap[eN_global_old][nN];
          int nN_subdomain  = nodeNumbering_global2subdomainMap[nN_global_new];
          newMesh.subdomainp->elementNodesArray[eN*newMesh.subdomainp->nNodes_element + nN]= nN_subdomain;
        }
    }
  //
  //ghost
  //
  set<int>::iterator eN_p = elements_overlap.begin();
  for (int eN = nElements_subdomain_new[rank]; eN < nElements_subdomain_new[rank] + int(elements_overlap.size()); eN++,eN_p++)
    {
      int eN_global_new = *eN_p;
      int eN_global_old = elementNumbering_global_new2old[eN_global_new];
      elementNumbering_subdomain2global[eN] = eN_global_new;
      newMesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypesMap[eN_global_old];
      for (int nN =  0; nN < newMesh.subdomainp->nNodes_element; nN++)
        {
          int nN_global_new = elementNodesArrayMap[eN_global_old][nN];
          int nN_subdomain  = nodeNumbering_global2subdomainMap[nN_global_new];
          newMesh.subdomainp->elementNodesArray[eN*newMesh.subdomainp->nNodes_element + nN]= nN_subdomain;
        }
    }
  ISRestoreIndices(elementNumberingIS_global_new2old,&elementNumbering_global_new2old);
  ISDestroy(&elementNumberingIS_subdomain_new2old);
  ISDestroy(&elementNumberingIS_global_new2old);
  //
  //element boundaries
  //
  //owned
  //
  for (int ebN=0; ebN < nElementBoundaries_subdomain_new[rank]; ebN++)
    {
      int ebN_global = ebN + elementBoundaryOffsets_new[rank];
      elementBoundaryNumbering_subdomain2global[ebN]=ebN_global;
      elementBoundaryNumbering_global2subdomainMap[ebN_global] = ebN;
    }
  //
  //ghost
  //
  set<int>::iterator ebN_p = elementBoundaries_overlap.begin();
  for(int ebN=nElementBoundaries_subdomain_new[rank];ebN < nElementBoundaries_subdomain_new[rank] + int(elementBoundaries_overlap.size()); ebN++,ebN_p++)
    {
      int ebN_global = *ebN_p;
      elementBoundaryNumbering_subdomain2global[ebN] = ebN_global;
      elementBoundaryNumbering_global2subdomainMap[ebN_global] = ebN;
    }
  //
  //need elementBoundariesArray to assign consistent numbering on subdomain
  //
  //local
  //
  newMesh.subdomainp->elementBoundariesArray = 
    new int[newMesh.subdomainp->nElements_global*newMesh.subdomainp->nElementBoundaries_element];
  for (int eN=0;eN<nElements_subdomain_new[rank];eN++)
    {
      int eN_global = eN+elementOffsets_new[rank];
      for (int ebN=0;ebN<newMesh.subdomainp->nElementBoundaries_element;ebN++)
        {
          newMesh.subdomainp->elementBoundariesArray[eN*newMesh.subdomainp->nElementBoundaries_element+ebN] =
            elementBoundaryNumbering_global2subdomainMap[elementBoundaryNumbering_global_old2new[elementBoundariesMap[eN_global][ebN]]];
        }
    }
  //
  //ghost elements
  //
  eN_p = elements_overlap.begin();
  for (int eN = nElements_subdomain_new[rank]; eN < nElements_subdomain_new[rank] + int(elements_overlap.size()); eN++,eN_p++)
    {
      int eN_global_new = *eN_p;
      for (int ebN=0;ebN<newMesh.subdomainp->nElementBoundaries_element;ebN++)
        {
          newMesh.subdomainp->elementBoundariesArray[eN*newMesh.subdomainp->nElementBoundaries_element+ebN] =
            elementBoundaryNumbering_global2subdomainMap[elementBoundaryNumbering_global_old2new[elementBoundariesMap[eN_global_new][ebN]]];
        }
    }      
  //
  //edges
  //
  //local
  //
  for (int edN=0; edN < nEdges_subdomain_new[rank]; edN++)
    {
      int edN_global = edN + edgeOffsets_new[rank];
      edgeNumbering_subdomain2global[edN]=edN_global;
      edgeNumbering_global2subdomainMap[edN_global] = edN;
    }
  //
  //ghost
  //
  set<int>::iterator edN_p = edges_overlap.begin();
  for(int edN=nEdges_subdomain_new[rank];edN < nEdges_subdomain_new[rank] + int(edges_overlap.size()); edN++,edN_p++)
    {
      int edN_global = *edN_p;
      edgeNumbering_subdomain2global[edN] = edN_global;
      edgeNumbering_global2subdomainMap[edN_global] = edN;
    }
  //
  //now build edgeNodes array in new numberings
  //
  newMesh.subdomainp->edgeNodesArray = new int[newMesh.subdomainp->nEdges_global*2];
  for (int i=0;i<newMesh.subdomainp->nEdges_global*2;i++)
    newMesh.subdomainp->edgeNodesArray[i] = -1;
  for (map<int,pair<int,int> >::iterator edgep=edgeNodesMap.begin();
       edgep!=edgeNodesMap.end();
       edgep++)
    {
      int edN_global_old = edgep->first;
      int edN_global_new = edgeNumbering_global_old2new[edN_global_old];
      assert(edgeNumbering_global2subdomainMap.find(edN_global_new) != edgeNumbering_global2subdomainMap.end());
      int edN_subdomain  = edgeNumbering_global2subdomainMap[edN_global_new];
      newMesh.subdomainp->edgeNodesArray[edN_subdomain*2+0] = nodeNumbering_global2subdomainMap[edgep->second.first];
      newMesh.subdomainp->edgeNodesArray[edN_subdomain*2+1] = nodeNumbering_global2subdomainMap[edgep->second.second];
    }
  edgeNumbering_global_old2new.resize(0);
  //
  //end edges
  //

  //now build rest of subdomain mesh connectivity information etc
  bool callOld=true;
  if(callOld)
    constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_tetrahedron(*newMesh.subdomainp);
  else
    {
      //const int DEFAULT_ELEMENT_MATERIAL=0;
      const int DEFAULT_NODE_MATERIAL=-1;
      const int INTERIOR_NODE_MATERIAL=0;
      const int EXTERIOR_NODE_MATERIAL=1;
      const int INTERIOR_ELEMENT_BOUNDARY_MATERIAL=0;
      const int EXTERIOR_ELEMENT_BOUNDARY_MATERIAL=1;

      newMesh.subdomainp->nNodes_elementBoundary = 3;
      newMesh.subdomainp->nElementBoundaries_element = 4;
      assert(newMesh.subdomainp->elementBoundariesArray);
      using namespace std;
      //double start,stop;
      map<NodeTuple<3>,
          ElementNeighbors> elementBoundaryElements;
      map<NodeTuple<3>,
          int> elementBoundaryIds;
      //start=CurrentTime();
      //cout<<"Extracting boundary elements"<<endl;
      for(int eN=0;eN<newMesh.subdomainp->nElements_global;eN++)
        for(int ebN=0;ebN<newMesh.subdomainp->nElementBoundaries_element;ebN++)
          {
            register int ebN_global = newMesh.subdomainp->elementBoundariesArray[eN*newMesh.subdomainp->nElementBoundaries_element+ebN];
            register int nodes[3];
            nodes[0] = newMesh.subdomainp->elementNodesArray[eN*4+((ebN+1)%4)];
            nodes[1] = newMesh.subdomainp->elementNodesArray[eN*4+((ebN+2)%4)];
            nodes[2] = newMesh.subdomainp->elementNodesArray[eN*4+((ebN+3)%4)];
            NodeTuple<3> ebt(nodes);
            if(elementBoundaryElements.find(ebt) != elementBoundaryElements.end())
              {
                elementBoundaryElements[ebt].right=eN;
                elementBoundaryElements[ebt].right_ebN_element=ebN;
                assert(elementBoundaryIds[ebt] == ebN_global);
              }
            else
              {
                elementBoundaryElements.insert(elementBoundaryElements.end(),make_pair(ebt,ElementNeighbors(eN,ebN)));
                elementBoundaryIds.insert(elementBoundaryIds.end(),make_pair(ebt,ebN_global));
              }
          }
      //stop = CurrentTime();
      //cout<<"Elapsed time for building element boundary elements map= "<<(stop-start)<<"s"<<endl;
      newMesh.subdomainp->nElementBoundaries_global = elementBoundaryElements.size();
      //cout<<"nElementBoundaries_global = "<<newMesh.subdomainp->nElementBoundaries_global<<endl;

      //cout<<"Allocating Arrays"<<endl;
      //start = CurrentTime();
      set<int> interiorElementBoundaries,exteriorElementBoundaries;
      newMesh.subdomainp->elementBoundaryNodesArray =  new int[newMesh.subdomainp->nElementBoundaries_global*newMesh.subdomainp->nNodes_elementBoundary];
      newMesh.subdomainp->elementBoundaryElementsArray = new int[newMesh.subdomainp->nElementBoundaries_global*2];
      newMesh.subdomainp->elementBoundaryLocalElementBoundariesArray = new int[newMesh.subdomainp->nElementBoundaries_global*2];
      newMesh.subdomainp->elementNeighborsArray = new int[newMesh.subdomainp->nElements_global*newMesh.subdomainp->nElementBoundaries_element];
      //stop = CurrentTime();
      //cout<<"Elapsed time for allocating arrays = "<<(stop-start)<<"s"<<endl;

      //cout<<"Generating elementBoundaryElementsArray and elementBoundaryNodesArray"<<endl;
      //start = CurrentTime();
      for(map<NodeTuple<3>,ElementNeighbors>::iterator eb=elementBoundaryElements.begin();
          eb != elementBoundaryElements.end();
          eb++)
        {
          int ebN = elementBoundaryIds[eb->first];
          newMesh.subdomainp->elementBoundaryNodesArray[ebN*3 + 0] = eb->first.nodes[0];
          newMesh.subdomainp->elementBoundaryNodesArray[ebN*3 + 1] = eb->first.nodes[1];
          newMesh.subdomainp->elementBoundaryNodesArray[ebN*3 + 2] = eb->first.nodes[2];

          newMesh.subdomainp->elementBoundaryElementsArray[ebN*2 + 0] = eb->second.left;
          newMesh.subdomainp->elementBoundaryLocalElementBoundariesArray[ebN*2 + 0] = eb->second.left_ebN_element;
          newMesh.subdomainp->elementBoundaryElementsArray[ebN*2 + 1] = eb->second.right;
          newMesh.subdomainp->elementBoundaryLocalElementBoundariesArray[ebN*2 + 1] = eb->second.right_ebN_element;
          newMesh.subdomainp->elementNeighborsArray[eb->second.left*newMesh.subdomainp->nElementBoundaries_element + eb->second.left_ebN_element] = eb->second.right; 
          if(eb->second.right != -1)
            {
              interiorElementBoundaries.insert(ebN);
              newMesh.subdomainp->elementNeighborsArray[eb->second.right*newMesh.subdomainp->nElementBoundaries_element + eb->second.right_ebN_element] = eb->second.left; 
            }
          else
            exteriorElementBoundaries.insert(ebN);          
          assert(newMesh.subdomainp->elementBoundariesArray[eb->second.left*newMesh.subdomainp->nElementBoundaries_element + eb->second.left_ebN_element] == ebN);
          if (eb->second.right != -1)
            {
              assert(newMesh.subdomainp->elementBoundariesArray[eb->second.right*newMesh.subdomainp->nElementBoundaries_element + eb->second.right_ebN_element] == ebN);
            }
        }
      newMesh.subdomainp->nInteriorElementBoundaries_global = interiorElementBoundaries.size();
      newMesh.subdomainp->interiorElementBoundariesArray = new int[newMesh.subdomainp->nInteriorElementBoundaries_global];
      newMesh.subdomainp->nExteriorElementBoundaries_global = exteriorElementBoundaries.size();
      newMesh.subdomainp->exteriorElementBoundariesArray = new int[newMesh.subdomainp->nExteriorElementBoundaries_global];
      int ebNI=0,ebNE=0;
      for (set<int>::iterator ebN=interiorElementBoundaries.begin();ebN != interiorElementBoundaries.end(); ebN++,ebNI++)
        newMesh.subdomainp->interiorElementBoundariesArray[ebNI] = *ebN;
      for (set<int>::iterator ebN=exteriorElementBoundaries.begin();ebN != exteriorElementBoundaries.end(); ebN++,ebNE++)
        newMesh.subdomainp->exteriorElementBoundariesArray[ebNE] = *ebN;
      set<NodeTuple<2> > edges;
      for (int eN=0;eN<newMesh.subdomainp->nElements_global;eN++)
        {
          int nodes[2];
          for (int nN_L=0;nN_L<newMesh.subdomainp->nNodes_element;nN_L++)
            for (int nN_R=nN_L+1;nN_R<newMesh.subdomainp->nNodes_element;nN_R++)
              {
                nodes[0] = newMesh.subdomainp->elementNodesArray[eN*4+nN_L];
                nodes[1] = newMesh.subdomainp->elementNodesArray[eN*4+nN_R];
                edges.insert(NodeTuple<2>(nodes));
              }
        }
      assert(newMesh.subdomainp->nEdges_global == int(edges.size()));
      vector<set<int> > nodeStar(newMesh.subdomainp->nNodes_global);
      for (int edgeN=0;edgeN<newMesh.subdomainp->nEdges_global;edgeN++)
        {
          nodeStar[newMesh.subdomainp->edgeNodesArray[edgeN*2+0]].insert(newMesh.subdomainp->edgeNodesArray[edgeN*2+1]);
          nodeStar[newMesh.subdomainp->edgeNodesArray[edgeN*2+1]].insert(newMesh.subdomainp->edgeNodesArray[edgeN*2+0]);
        }
      newMesh.subdomainp->nodeStarOffsets = new int[newMesh.subdomainp->nNodes_global+1];
      newMesh.subdomainp->nodeStarOffsets[0] = 0;
      for (int nN=1;nN<newMesh.subdomainp->nNodes_global+1;nN++)
        newMesh.subdomainp->nodeStarOffsets[nN] = newMesh.subdomainp->nodeStarOffsets[nN-1] + nodeStar[nN-1].size();
      newMesh.subdomainp->nodeStarArray = new int[newMesh.subdomainp->nodeStarOffsets[newMesh.subdomainp->nNodes_global]];
      for (int nN=0,offset=0;nN<newMesh.subdomainp->nNodes_global;nN++)
        for (set<int>::iterator nN_star=nodeStar[nN].begin();nN_star!=nodeStar[nN].end();nN_star++,offset++)
          newMesh.subdomainp->nodeStarArray[offset] = *nN_star;
      //stop = CurrentTime();
      newMesh.subdomainp->max_nNodeNeighbors_node=0;
      for (int nN=0;nN<newMesh.subdomainp->nNodes_global;nN++)
        newMesh.subdomainp->max_nNodeNeighbors_node=max(newMesh.subdomainp->max_nNodeNeighbors_node,newMesh.subdomainp->nodeStarOffsets[nN+1]-newMesh.subdomainp->nodeStarOffsets[nN]);
      //mwf repeat for node-->elements arrays
      vector<set<int> > nodeElementsStar(newMesh.subdomainp->nNodes_global);
      for (int eN = 0; eN < newMesh.subdomainp->nElements_global; eN++)
        {
          for (int nN = 0; nN < newMesh.subdomainp->nNodes_element; nN++)
            nodeElementsStar[newMesh.subdomainp->elementNodesArray[eN*newMesh.subdomainp->nNodes_element+nN]].insert(eN);
        }
      newMesh.subdomainp->nodeElementOffsets = new int[newMesh.subdomainp->nNodes_global+1];
      newMesh.subdomainp->nodeElementOffsets[0] = 0;
      for (int nN = 0; nN < newMesh.subdomainp->nNodes_global; nN++)
        newMesh.subdomainp->nodeElementOffsets[nN+1] = newMesh.subdomainp->nodeElementOffsets[nN]+nodeElementsStar[nN].size();
      newMesh.subdomainp->nodeElementsArray  = new int[newMesh.subdomainp->nodeElementOffsets[newMesh.subdomainp->nNodes_global]];
      for (int nN=0,offset=0; nN < newMesh.subdomainp->nNodes_global; nN++)
        {
          for (set<int>::iterator eN_star = nodeElementsStar[nN].begin(); eN_star != nodeElementsStar[nN].end();
               eN_star++,offset++)
            {
              newMesh.subdomainp->nodeElementsArray[offset] = *eN_star;
            }
        }
      //mwf end node-->elements construction
      newMesh.subdomainp->elementBoundaryMaterialTypes = new int[newMesh.subdomainp->nElementBoundaries_global];
      //if nodeMaterial is DEFAULT, go ahead and set to interior or exterior
      //depending on which boundary node belongs to. 
      //If node on at least one exterior boundary then it's exterior
      for (int ebNE = 0; ebNE < newMesh.subdomainp->nExteriorElementBoundaries_global; ebNE++)
        {
          int ebN = newMesh.subdomainp->exteriorElementBoundariesArray[ebNE];
          newMesh.subdomainp->elementBoundaryMaterialTypes[ebN] = EXTERIOR_ELEMENT_BOUNDARY_MATERIAL;
          for (int nN_local = 0; nN_local < newMesh.subdomainp->nNodes_elementBoundary; nN_local++)
            {
              int nN = newMesh.subdomainp->elementBoundaryNodesArray[ebN*newMesh.subdomainp->nNodes_elementBoundary+nN_local];
              if (newMesh.subdomainp->nodeMaterialTypes[nN] == DEFAULT_NODE_MATERIAL)
                newMesh.subdomainp->nodeMaterialTypes[nN] = EXTERIOR_NODE_MATERIAL;
            }
        }
      for (int ebNI = 0; ebNI < newMesh.subdomainp->nInteriorElementBoundaries_global; ebNI++)
        {
          int ebN = newMesh.subdomainp->interiorElementBoundariesArray[ebNI];
          newMesh.subdomainp->elementBoundaryMaterialTypes[ebN] = INTERIOR_ELEMENT_BOUNDARY_MATERIAL;
          for (int nN_local = 0; nN_local < newMesh.subdomainp->nNodes_elementBoundary; nN_local++)
            {
              int nN = newMesh.subdomainp->elementBoundaryNodesArray[ebN*newMesh.subdomainp->nNodes_elementBoundary+nN_local];
              if (newMesh.subdomainp->nodeMaterialTypes[nN] == DEFAULT_NODE_MATERIAL)
                newMesh.subdomainp->nodeMaterialTypes[nN] = INTERIOR_NODE_MATERIAL;
            }
        }
      //cout<<"Elapsed time for populating arrays = "<<(stop-start)<<"s"<<endl;
    }
  //build local geometric info
  allocateGeometricInfo_tetrahedron(*newMesh.subdomainp);
  computeGeometricInfo_tetrahedron(*newMesh.subdomainp);

  if (hasElementBoundaryMarkers)
    {
      assert(newMesh.subdomainp->elementBoundariesArray != NULL);
      for (map<int,int>::iterator ebmp = elementBoundaryMaterialTypesMap.begin(); ebmp != elementBoundaryMaterialTypesMap.end();ebmp++)
        {
          int ebN_global_new = elementBoundaryNumbering_global_old2new[ebmp->first];
          assert(elementBoundaryNumbering_global2subdomainMap.find(ebN_global_new) != elementBoundaryNumbering_global2subdomainMap.end());
          int ebN_subdomain  = elementBoundaryNumbering_global2subdomainMap[ebN_global_new];
          newMesh.subdomainp->elementBoundaryMaterialTypes[ebN_subdomain] = ebmp->second;
        }
      if (!hasVertexMarkers)
        for (int ebNE = 0; ebNE < newMesh.subdomainp->nExteriorElementBoundaries_global; ebNE++)
          {
            int ebN = newMesh.subdomainp->exteriorElementBoundariesArray[ebNE];
            for (int nN_local = 0; nN_local < newMesh.subdomainp->nNodes_elementBoundary; nN_local++)
              {
                int nN = newMesh.subdomainp->elementBoundaryNodesArray[ebN*newMesh.subdomainp->nNodes_elementBoundary+nN_local];
                newMesh.subdomainp->nodeMaterialTypes[nN] = newMesh.subdomainp->elementBoundaryMaterialTypes[ebN];
              }
          }
    }
  elementBoundaryNumbering_global_old2new.resize(0);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done with material types");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  PetscLogEventEnd(build_subdomains_renumber_event,0,0,0,0);
  int build_subdomains_cleanup_event;
  PetscLogEventRegister("Cleanup",0,&build_subdomains_cleanup_event);
  PetscLogEventBegin(build_subdomains_cleanup_event,0,0,0,0);
  //transfer information about owned nodes and elements to mesh
  if (newMesh.nodeOffsets_subdomain_owned) 
    delete [] newMesh.nodeOffsets_subdomain_owned;
  if (newMesh.elementOffsets_subdomain_owned)
    delete [] newMesh.elementOffsets_subdomain_owned;
  if (newMesh.elementBoundaryOffsets_subdomain_owned)
    delete [] newMesh.elementBoundaryOffsets_subdomain_owned;
  if (newMesh.edgeOffsets_subdomain_owned)
    delete [] newMesh.edgeOffsets_subdomain_owned;
  newMesh.nodeOffsets_subdomain_owned    = new int[size+1];
  newMesh.elementOffsets_subdomain_owned = new int[size+1];
  newMesh.elementBoundaryOffsets_subdomain_owned = new int[size+1];
  newMesh.edgeOffsets_subdomain_owned = new int[size+1];
  for (int sdN = 0; sdN < size+1; sdN++)
    {
      newMesh.nodeOffsets_subdomain_owned[sdN]    = nodeOffsets_new[sdN];
      newMesh.elementOffsets_subdomain_owned[sdN] = elementOffsets_new[sdN];
      newMesh.elementBoundaryOffsets_subdomain_owned[sdN] = elementBoundaryOffsets_new[sdN];
      newMesh.edgeOffsets_subdomain_owned[sdN] = edgeOffsets_new[sdN];
    }
  if (newMesh.nodeNumbering_subdomain2global)
    delete [] newMesh.nodeNumbering_subdomain2global;
  newMesh.nodeNumbering_subdomain2global = new int[newMesh.subdomainp->nNodes_global];
  for (int nN = 0; nN < newMesh.subdomainp->nNodes_global; nN++)
    newMesh.nodeNumbering_subdomain2global[nN] = nodeNumbering_subdomain2global[nN];
  if (newMesh.elementNumbering_subdomain2global)
    delete [] newMesh.elementNumbering_subdomain2global;
  newMesh.elementNumbering_subdomain2global = new int[newMesh.subdomainp->nElements_global];
  for (int eN = 0; eN < newMesh.subdomainp->nElements_global; eN++)
    newMesh.elementNumbering_subdomain2global[eN] = elementNumbering_subdomain2global[eN];
  //
  if (newMesh.elementBoundaryNumbering_subdomain2global)
    delete [] newMesh.elementBoundaryNumbering_subdomain2global;
  newMesh.elementBoundaryNumbering_subdomain2global = new int[newMesh.subdomainp->nElementBoundaries_global];
  for (int ebN = 0; ebN < newMesh.subdomainp->nElementBoundaries_global; ebN++)
    newMesh.elementBoundaryNumbering_subdomain2global[ebN] = elementBoundaryNumbering_subdomain2global[ebN];
  //
  if (newMesh.edgeNumbering_subdomain2global)
    delete [] newMesh.edgeNumbering_subdomain2global;
  newMesh.edgeNumbering_subdomain2global = new int[newMesh.subdomainp->nEdges_global];
  for (int i=0; i< newMesh.subdomainp->nEdges_global; i++)
    newMesh.edgeNumbering_subdomain2global[i] = edgeNumbering_subdomain2global[i];
  //cleanup
  /* out of core*/
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);
  /* out of core */
  PetscLogEventEnd(build_subdomains_cleanup_event,0,0,0,0);
  PetscLogStagePop();
  PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
  ierr = enforceMemoryLimit(rank, max_rss_gb,"Done with partitioning!");CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  return 0;
}

//todo add overlap for element based partitions
int partitionElements(Mesh& mesh, int nElements_overlap)
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  //Contents 
  //
  //1. Partition the elements in the "default" partition (contiguous
  //chunks in given ordering) 
  //
  //2. Partition the elementNeighbors based on this partition
  //
  //3. Pass to Parmetis to build a better partition of the elements
  //
  //4. Tag a subset of the nodes and faces on the subdomain elements as owned
  //using a mark and pass approach.** 
  //
  //5. Extract the nodes and faces in the
  //overlapping elements.** 
  //
  //6. Build the subdomain mesh from the
  //subdomain elements 
  //
  //**To be more general we could get all the support (i.e. faces
  //and edges) and partitiong them, but the main reason for
  //partitioning is to keep track of a global numbering for degrees
  //of freedom that live on each type of geometric entity. We only
  //have node and element based DOFs so I just rebuild the other
  //information once we have elements and nodes partitioned.
  //
  // \todo check that I restore all data that PETSc expects to have
  // back, add PETSc error checking macros
  //
  //1. Build default partitioning
  //
  //get offsets so we can calculate the processor to global mapping
  //for elements in the old (default) partitioning
  valarray<int> elementOffsets_old(size+1);
  elementOffsets_old[0] = 0;
  for(int sdN=0;sdN<size;sdN++)
    elementOffsets_old[sdN+1] = elementOffsets_old[sdN] + 
      int(mesh.nElements_global)/size + 
      (int(mesh.nElements_global)%size > sdN);
    
  //2. Extract subdomain element adjacency information could read
  //only the required portion from a file
  int nElements_subdomain = (elementOffsets_old[rank+1] - elementOffsets_old[rank]);
  PetscInt *elementNeighborsOffsets_subdomain,*elementNeighbors_subdomain,*weights_subdomain;
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain+1),&elementNeighborsOffsets_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain*mesh.nElementBoundaries_element),&elementNeighbors_subdomain);
  PetscMalloc(sizeof(PetscInt)*(nElements_subdomain*mesh.nElementBoundaries_element),&weights_subdomain);
  //this wastes a little space
  elementNeighborsOffsets_subdomain[0] = 0;
  for (int eN=0,offset=0; eN < nElements_subdomain; eN++)
    {
      int eN_global = elementOffsets_old[rank] + eN;
      int offsetStart=offset;
      for (int ebN=0; ebN< mesh.nElementBoundaries_element;ebN++)
        {
          int eN_neighbor_global = mesh.elementNeighborsArray[eN_global*mesh.nElementBoundaries_element + ebN];
          if (eN_neighbor_global >= 0 )
            elementNeighbors_subdomain[offset++] = eN_neighbor_global;
        }
      elementNeighborsOffsets_subdomain[eN+1]=offset;
      sort(&elementNeighbors_subdomain[offsetStart],&elementNeighbors_subdomain[offset]);
      int weight = (elementNeighborsOffsets_subdomain[eN+1] - elementNeighborsOffsets_subdomain[eN]);
      for (int k=elementNeighborsOffsets_subdomain[eN];k<elementNeighborsOffsets_subdomain[eN+1];k++)
        weights_subdomain[k] = weight;
    }
  //3. Generate the  new partitiong using PETSc, this is done in parallel using parmetis
  Mat petscAdjacency;
  //     MatCreateMPIAdj(PROTEUS_COMM_WORLD,
  //                     nElements_subdomain, mesh.nElements_global, 
  //                     &elementNeighborsOffsets_subdomain[0], &elementNeighbors_subdomain[0], 
  //                     &weights_subdomain[0],//PETSC_NULL, 
  //                     &petscAdjacency);
  ierr = MatCreateMPIAdj(PROTEUS_COMM_WORLD,
			 nElements_subdomain, 
			 mesh.nElements_global,
			 elementNeighborsOffsets_subdomain, 
			 elementNeighbors_subdomain,
			 PETSC_NULL,//weights_subdomain,
			 &petscAdjacency);CHKERRABORT(PROTEUS_COMM_WORLD, ierr);
  MatPartitioning petscPartition;
  MatPartitioningCreate(PROTEUS_COMM_WORLD,&petscPartition);
  MatPartitioningSetAdjacency(petscPartition,petscAdjacency);
  MatPartitioningSetFromOptions(petscPartition);
  
  //get a petsc index set that has the new submdomain number for each element
  IS elementPartitioningIS_new;
  MatPartitioningApply(petscPartition,&elementPartitioningIS_new); 
  MatPartitioningDestroy(&petscPartition);
  //MatDestroy(&petscAdjacency);
  //ISView(elementPartitioningIS_new,PETSC_VIEWER_STDOUT_WORLD);
    
  //experiment with metis
  //mwf set some defaults and not call if size == 1 since metis crashes
  //cek commenting out for now
  //int etype=1,edgecut=0,base=0;
  //epart assign everything to processor zero by default
  //valarray<int> epart(0,mesh.nElements_global),npart(mesh.nNodes_global);
  //if (size > 1)
  //  METIS_PartMeshNodal(&mesh.nElements_global,&mesh.nNodes_global,mesh.elementNodesArray,&etype,&base,&size,&edgecut,&epart[0],&npart[0]);
  //ISCreateGeneralWithArray(PETSC_COMM_SELF,mesh.nElements_global,&epart[0],&elementPartitioningIS_new);
  //write mesh to view with showme
  //     std::ofstream nodeout("mesh.node"),eleout("mesh.ele"),partout("mesh.part");
  //     eleout<<mesh.nElements_global<<" 3 0"<<std::endl;
  //     partout<<mesh.nElements_global<<"\t"<<size<<std::endl;
  //     for (int eN=0;eN<mesh.nElements_global;eN++)
  //       {
  // 	partout<<(eN+1)<<"\t"<<(1+epart[eN])<<std::endl;
  // 	eleout<<(eN+1)<<"\t"<<(1+mesh.elementNodesArray[eN*3+0])
  // 	      <<"\t"<<(1+mesh.elementNodesArray[eN*3+1])
  // 	      <<"\t"<<(1+mesh.elementNodesArray[eN*3+2])
  // 	      <<std::endl;
  //       }
  //     nodeout<<mesh.nNodes_global<<" 2 0 0"<<std::endl;
  //     for (int nN=0;nN<mesh.nNodes_global;nN++)
  //       {
  // 	nodeout<<(nN+1)<<"\t"<<mesh.nodeArray[nN*3+0]<<"\t"<<mesh.nodeArray[nN*3+1]<<std::endl;
  //       }
  //     eleout.close();
  //     partout.close();
  //count the new number of elements on each subdomain
  valarray<int> nElements_subdomain_new(size);
  ISPartitioningCount(elementPartitioningIS_new,size,&nElements_subdomain_new[0]);

  //get the new offsets for the subdomain to global numbering
  valarray<int> elementOffsets_new(size+1);
  elementOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    elementOffsets_new[sdN+1] = elementOffsets_new[sdN] + nElements_subdomain_new[sdN];

  //get the new element numbers for the elements on this  subdomain
  IS elementNumberingIS_subdomain_old2new;
  ISPartitioningToNumbering(elementPartitioningIS_new,&elementNumberingIS_subdomain_old2new);

  //now get the new element numbers for the whole mesh so that we
  //can just read this processors elements, reorder, and renumber**
  //
  //**We could do this in parallel by scattering all the element
  //information
  IS elementNumberingIS_global_old2new;
  ISAllGather(elementNumberingIS_subdomain_old2new,&elementNumberingIS_global_old2new);
  const PetscInt *elementNumbering_global_old2new;
  ISGetIndices(elementNumberingIS_global_old2new,&elementNumbering_global_old2new);
  valarray<int> elementNumbering_global_new2old(mesh.nElements_global);
  for(int eN=0;eN<mesh.nElements_global;eN++)
    elementNumbering_global_new2old[elementNumbering_global_old2new[eN]] = eN;

  //Sort element based arrays, maybe I don't need to do this, maybe
  //I just need to start writing into the subdomain mesh here and
  //preserve subdomain2old and subdomain2global mappings
  valarray<int> elementNodesArray_new(mesh.nElements_global*mesh.nNodes_element),
    elementNeighborsArray_new(mesh.nElements_global*mesh.nElementBoundaries_element),
    elementMaterialTypes_new(mesh.nElements_global),
    elementBoundariesArray_new(mesh.nElements_global*mesh.nElementBoundaries_element),
    elementBoundaryElementsArray_new(mesh.nElementBoundaries_global*2),
    elementBoundaryLocalElementBoundariesArray_new(mesh.nElementBoundaries_global*2),
    elementBoundaryNodesArray_new(mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary),
    edgeNodesArray_new(mesh.nEdges_global*2);
  valarray<int> elementBoundaryMaterialTypes_new(mesh.nElementBoundaries_global);
  for (int eN=0;eN<mesh.nElements_global;eN++)
    {
      for (int nN=0;nN<mesh.nNodes_element;nN++)
        elementNodesArray_new[eN*mesh.nNodes_element + nN] = 
          mesh.elementNodesArray[elementNumbering_global_new2old[eN]*mesh.nNodes_element+nN];
      for (int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
        {
          elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN] = 
            mesh.elementNeighborsArray[elementNumbering_global_new2old[eN]*mesh.nElementBoundaries_element+ebN];
          //need new elements --> old element boundary numbers for now
          elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN] =
            mesh.elementBoundariesArray[elementNumbering_global_new2old[eN]*mesh.nElementBoundaries_element+ebN];
        }
      elementMaterialTypes_new[eN] = mesh.elementMaterialTypes[elementNumbering_global_new2old[eN]];
	
    }
  //renumber references to element numbers
  for (int eN=0;eN<mesh.nElements_global;eN++)
    {
      for (int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
        {
          int eN_ebN = elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN];
          if (eN_ebN >= 0)
            elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN] = 
              elementNumbering_global_old2new[eN_ebN];
        }
    }
  //4. now we need to build a new node ordering with better data locality for C0 finite elements
  //otherwise we could just grab the nodes on the subdomain and not worry about ownership
  //in the long run it wouldn't be bad to do a global repartition of faces and edges for mixed hybrid
  //and non-conforming finite elements
  MPI_Status status;
  PetscBT nodeMask; 
  PetscBTCreate(mesh.nNodes_global,&nodeMask);
  if (rank > 0) 
    {
      MPI_Recv(nodeMask,PetscBTLength(mesh.nNodes_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status);
    }
  //mark the unmarked nodes on this subdomain and store the node numbers
  set<int> nodes_subdomain_owned;
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    for(int nN=0;nN<mesh.nNodes_element;nN++)
      {
        int nN_global = elementNodesArray_new[eN*mesh.nNodes_element+nN];
        if (!PetscBTLookupSet(nodeMask,nN_global))
          nodes_subdomain_owned.insert(nN_global);
      }
  //ship off the mask
  if (rank < size-1)
    MPI_Send(nodeMask,PetscBTLength(mesh.nNodes_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&nodeMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy"<<endl;
  //get the number of nodes on each processor
  valarray<int> nNodes_subdomain_new(size),
    nodeOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nNodes_subdomain_new[sdN] = nodes_subdomain_owned.size();
    else
      nNodes_subdomain_new[sdN] = 0;
  valarray<int> nNodes_subdomain_new_send=nNodes_subdomain_new;
  MPI_Allreduce(&nNodes_subdomain_new_send[0],&nNodes_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  nodeOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    nodeOffsets_new[sdN+1] = nodeOffsets_new[sdN]+nNodes_subdomain_new[sdN];
      
  assert(nodeOffsets_new[size]==mesh.nNodes_global);
      
  //Now as with elements build a global node numbering, sort node
  //based information, and renumber references to node numbers
  valarray<int> nodeNumbering_new2old(nodes_subdomain_owned.size());
  set<int>::iterator nN_ownedp=nodes_subdomain_owned.begin();
  for (int nN=0;nN<int(nodes_subdomain_owned.size());nN++)
    {
      nodeNumbering_new2old[nN] = *nN_ownedp++;
    }
  IS nodeNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nodes_subdomain_owned.size(),&nodeNumbering_new2old[0],PETSC_COPY_VALUES,&nodeNumberingIS_new2old);
  IS nodeNumberingIS_global_new2old;
  ISAllGather(nodeNumberingIS_new2old,&nodeNumberingIS_global_new2old);
  const PetscInt *nodeNumbering_global_new2old;
  valarray<int> nodeNumbering_old2new_global(mesh.nNodes_global);
  ISGetIndices(nodeNumberingIS_global_new2old,&nodeNumbering_global_new2old);
  for (int nN=0;nN<mesh.nNodes_global;nN++)
    {
      nodeNumbering_old2new_global[nodeNumbering_global_new2old[nN]] = nN;
    }
  for (int eN=0;eN < mesh.nElements_global; eN++)
    {
      int nN_old;
      for (int nN=0;nN < mesh.nNodes_element; nN++)
        {
          nN_old = elementNodesArray_new[eN*mesh.nNodes_element+nN];
          elementNodesArray_new[eN*mesh.nNodes_element+nN] = nodeNumbering_old2new_global[nN_old];
        }
    }
  //mwf postpone until after element boundary renumbering
  //     for (int i=0;i<mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary;i++)
  //       {
  //         int nN_old = mesh.elementBoundaryNodesArray[i];
  //         elementBoundaryNodesArray_new[i] = nodeNumbering_old2new_global[nN_old];
  //       }
  for (int i=0;i<mesh.nEdges_global*2;i++)
    {
      int nN_old = mesh.edgeNodesArray[i];
      edgeNodesArray_new[i] = nodeNumbering_old2new_global[nN_old];
    }
  valarray<int> nodeStarArray_new(mesh.nodeStarOffsets[mesh.nNodes_global]);
  for (int i=0;i<mesh.nodeStarOffsets[mesh.nNodes_global];i++)
    {
      int nN_old = mesh.nodeStarArray[i];
      nodeStarArray_new[i] = nodeNumbering_old2new_global[nN_old];
    }
  valarray<double> nodeArray_new(mesh.nNodes_global*3);
  valarray<int>  nodeMaterialTypes_new(mesh.nNodes_global);
  for (int nN=0;nN<mesh.nNodes_global;nN++)
    {
      int nN_new = nodeNumbering_old2new_global[nN];
      nodeArray_new[nN_new*3+0] = mesh.nodeArray[nN*3+0];
      nodeArray_new[nN_new*3+1] = mesh.nodeArray[nN*3+1];
      nodeArray_new[nN_new*3+2] = mesh.nodeArray[nN*3+2];
      nodeMaterialTypes_new[nN_new] = mesh.nodeMaterialTypes[nN];
    }

  //4b. repeat process to build global face numbering
  MPI_Status status_elementBoundaries;
  PetscBT elementBoundaryMask; 
  PetscBTCreate(mesh.nElementBoundaries_global,&elementBoundaryMask);
  if (rank > 0) 
    {
      MPI_Recv(elementBoundaryMask,PetscBTLength(mesh.nElementBoundaries_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status_elementBoundaries);
    }
  //mark the unmarked faces on this subdomain and store the global face numbers
  set<int> elementBoundaries_subdomain_owned;
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    for(int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
      {
        int ebN_global = elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN];
        if(!PetscBTLookup(elementBoundaryMask,ebN_global))
          {
            bool notFound=true;
            for(int nN=0;nN<mesh.nNodes_elementBoundary;nN++)
              {
                int nN_global_old = mesh.elementBoundaryNodesArray[ebN_global*mesh.nNodes_elementBoundary+nN];
                if (nodes_subdomain_owned.count(nN_global_old) > 0)
                  {
                    notFound=false;
                  }
              }
            if (notFound)
              {
                //std::cout<<"=========================Face has no owned nodes"<<std::endl;
                PetscBTSet(elementBoundaryMask,ebN_global);
                elementBoundaries_subdomain_owned.insert(ebN_global);
              }
            else
              {
                PetscBTSet(elementBoundaryMask,ebN_global);
                elementBoundaries_subdomain_owned.insert(ebN_global);
              }
          }
      }
  //std::cout<<"Done marking element boundares "<<elementBoundaries_subdomain_owned.size()<<std::endl;
  //ship off the mask
  if (rank < size-1)
    MPI_Send(elementBoundaryMask,PetscBTLength(mesh.nElementBoundaries_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&elementBoundaryMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy for elementBoundaries"<<endl;
  //get the number of elementBoundaries on each processor
  valarray<int> nElementBoundaries_subdomain_new(size),
    elementBoundaryOffsets_new(size+1);
  for (int sdN=0;sdN<size;sdN++)
    if (sdN == rank)
      nElementBoundaries_subdomain_new[sdN] = elementBoundaries_subdomain_owned.size();
    else
      nElementBoundaries_subdomain_new[sdN] = 0;
  valarray<int> nElementBoundaries_subdomain_new_send=nElementBoundaries_subdomain_new;
  MPI_Allreduce(&nElementBoundaries_subdomain_new_send[0],&nElementBoundaries_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  elementBoundaryOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    elementBoundaryOffsets_new[sdN+1] = elementBoundaryOffsets_new[sdN]+nElementBoundaries_subdomain_new[sdN];
  //Now as with elements and nodes build a global face numbering
  //resetting the face based information is a little different since much of this is currently built below based
  //on the element and node information
  //
  valarray<int> elementBoundaryNumbering_new2old(elementBoundaries_subdomain_owned.size());
  set<int>::iterator ebN_ownedp=elementBoundaries_subdomain_owned.begin();
  for (int ebN=0;ebN<int(elementBoundaries_subdomain_owned.size());ebN++)
    {
      elementBoundaryNumbering_new2old[ebN] = *ebN_ownedp++;
    }
  IS elementBoundaryNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,elementBoundaries_subdomain_owned.size(),&elementBoundaryNumbering_new2old[0],PETSC_COPY_VALUES,&elementBoundaryNumberingIS_new2old);
  IS elementBoundaryNumberingIS_global_new2old;
  ISAllGather(elementBoundaryNumberingIS_new2old,&elementBoundaryNumberingIS_global_new2old);
  const PetscInt *elementBoundaryNumbering_global_new2old;
  valarray<int> elementBoundaryNumbering_old2new_global(mesh.nElementBoundaries_global);
  ISGetIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);
  for (int ebN=0;ebN<mesh.nElementBoundaries_global;ebN++)
    {
      elementBoundaryNumbering_old2new_global[elementBoundaryNumbering_global_new2old[ebN]] = ebN;
    }
  for (int eN=0;eN < mesh.nElements_global; eN++)
    {
      int ebN_old;
      for (int ebN=0;ebN < mesh.nElementBoundaries_element; ebN++)
        {
          ebN_old = elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN];
          elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN] = elementBoundaryNumbering_old2new_global[ebN_old];
        }
    }
  //redo
  for (int ebN=0;ebN<mesh.nElementBoundaries_global;ebN++)
    {
      int ebN_old=elementBoundaryNumbering_global_new2old[ebN];
      for (int nN=0; nN<mesh.nNodes_elementBoundary; nN++)
        {
          int nN_old = mesh.elementBoundaryNodesArray[ebN_old*mesh.nNodes_elementBoundary+nN];
          elementBoundaryNodesArray_new[ebN*mesh.nNodes_elementBoundary+nN]=nodeNumbering_old2new_global[nN_old];
        }
      int eN_L_old = mesh.elementBoundaryElementsArray[ebN_old*2+0],
        eN_R_old = mesh.elementBoundaryElementsArray[ebN_old*2+1];
      elementBoundaryElementsArray_new[ebN*2+0] = elementNumbering_global_old2new[eN_L_old];
      if(eN_R_old >= 0)
        elementBoundaryElementsArray_new[ebN*2+1] = elementNumbering_global_old2new[eN_R_old];

      elementBoundaryMaterialTypes_new[ebN] = mesh.elementBoundaryMaterialTypes[ebN_old];
    }
  //     //mwf debug check constistency
  //     for (int ebN=0; ebN<mesh.nElementBoundaries_global;ebN++)
  //       {
  // 	int eN_left=elementBoundaryElementsArray_new[ebN*2+0];
  // 	int eN_right=elementBoundaryElementsArray_new[ebN*2+1];
  // 	assert(eN_left>=0);
  // 	bool found_ebN_left=false;
  // 	for (int ebN_element=0; ebN_element<mesh.nElementBoundaries_element; ebN_element++)
  // 	  {
  // 	    if (ebN == elementBoundariesArray_new[eN_left*mesh.nElementBoundaries_element+ebN_element])
  // 	      {
  // 		assert(ebN_element==elementBoundaryLocalElementBoundariesArray_new[ebN*2+0]);
  // 		found_ebN_left=true;
  // 	      }
  // 	  }
  // 	assert(found_ebN_left);
  // 	if (eN_right>=0)
  // 	  {
  // 	    bool found_ebN_right=false;
  // 	    for (int ebN_element=0; ebN_element<mesh.nElementBoundaries_element; ebN_element++)
  // 	      {
  // 		if (ebN == elementBoundariesArray_new[eN_right*mesh.nElementBoundaries_element+ebN_element])
  // 		  {
  // 		    assert(ebN_element==elementBoundaryLocalElementBoundariesArray_new[ebN*2+1]);
  // 		    found_ebN_right=true;
  // 		  }
  // 	      }
  // 	    assert(found_ebN_right);
  // 	  }
  //       }

  //do not renumber interior and exterior element boundary arrays yet

  //write partitioned mesh to view with "showme"
  //     std::ofstream nodeout("mesh.node"),eleout("mesh.ele"),partout("mesh.part");
  //     eleout<<mesh.nElements_global<<" 3 0"<<std::endl;
  //     partout<<mesh.nElements_global<<"\t"<<size<<std::endl;
  //     for (int eN=0;eN<mesh.nElements_global;eN++)
  //       {
  // 	partout<<(eN+1)<<"\t"<<(1+epart[elementNumbering_global_new2old[eN]])<<std::endl;
  // 	//partout<<(eN+1)<<"\t"<<(1+epart[eN])<<std::endl;
  // 	eleout<<(eN+1)<<"\t"<<(1+elementNodesArray_new[eN*3+0])
  // 	      <<"\t"<<(1+elementNodesArray_new[eN*3+1])
  // 	      <<"\t"<<(1+elementNodesArray_new[eN*3+2])
  // 	      <<std::endl;
  //       }
  //     nodeout<<mesh.nNodes_global<<" 2 0 0"<<std::endl;
  //     for (int nN=0;nN<mesh.nNodes_global;nN++)
  //       {
  // 	nodeout<<(nN+1)<<"\t"<<nodeArray_new[nN*3+0]<<"\t"<<nodeArray_new[nN*3+1]<<std::endl;
  //       }
  //     eleout.close();
  //     partout.close();

  //4c. Build global edge numbering as well
  // ownership is determined by the edges on owned elements, then
  // who owns the left (0) node  of the edge
  MPI_Status status_edges;
  PetscBT edgesMask; 
  PetscBTCreate(mesh.nEdges_global,&edgesMask);
  if (rank > 0) 
    {
      MPI_Recv(edgesMask,PetscBTLength(mesh.nEdges_global),MPI_CHAR,rank-1,0,PROTEUS_COMM_WORLD,&status_edges);
    }
  //mark the unmarked faces on this subdomain and store the global face numbers
  map<NodeTuple<2>, int> nodesEdgeMap_global; //new global node numbers --> original edge numbering
  set<int> edges_subdomain_owned;
  for (int ig = 0; ig < mesh.nEdges_global; ig++)
    {
      int nodes[2];
      nodes[0] = edgeNodesArray_new[2*ig];
      nodes[1] = edgeNodesArray_new[2*ig+1];
      NodeTuple<2> et(nodes);
      nodesEdgeMap_global[et] = ig;
    }
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    {
      for (int nN0=0;nN0<mesh.nNodes_element;nN0++)//assume all nodes in element may be edges
        for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
          {
            int nodes[2];
            nodes[0] = elementNodesArray_new[eN*mesh.nNodes_element+nN0];
            nodes[1] = elementNodesArray_new[eN*mesh.nNodes_element+nN1];
            NodeTuple<2> et(nodes);
            if (nodesEdgeMap_global.find(et) != nodesEdgeMap_global.end())
              {
                int edge_global = nodesEdgeMap_global[et];
                if(!PetscBTLookup(edgesMask,edge_global))
                  {
                    PetscBTSet(edgesMask,edge_global);
                    edges_subdomain_owned.insert(nodesEdgeMap_global[et]);
                  }
              }
          }
    }
  if (rank < size-1)
    MPI_Send(edgesMask,PetscBTLength(mesh.nEdges_global),MPI_CHAR,rank+1,0,PROTEUS_COMM_WORLD);
  ierr = PetscBTDestroy(&edgesMask);
  if (ierr)
    cerr<<"Error in PetscBTDestroy for edges"<<endl;

  valarray<int> nEdges_subdomain_new(size),
    edgeOffsets_new(size+1);
  
  for (int sdN=0; sdN < size; sdN++)
    {
      if (sdN == rank)
        nEdges_subdomain_new[sdN] = edges_subdomain_owned.size();
      else
        nEdges_subdomain_new[sdN] = 0;
    }
  //collect ownership info
  valarray<int> nEdges_subdomain_new_send=nEdges_subdomain_new;
  MPI_Allreduce(&nEdges_subdomain_new_send[0],&nEdges_subdomain_new[0],size,MPI_INT,MPI_SUM,PROTEUS_COMM_WORLD);
  edgeOffsets_new[0] = 0;
  for (int sdN=0;sdN<size;sdN++)
    edgeOffsets_new[sdN+1] = edgeOffsets_new[sdN]+nEdges_subdomain_new[sdN];
    
  //build new petsc numbering and global maps from old2new and new2old 
  valarray<int> edgeNumbering_new2old(edges_subdomain_owned.size());
  set<int>::iterator edges_ownedp = edges_subdomain_owned.begin();
  for (int i=0; i < int(edges_subdomain_owned.size());i++)
    edgeNumbering_new2old[i] = *edges_ownedp++;
    
  IS edgeNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,edges_subdomain_owned.size(),&edgeNumbering_new2old[0],PETSC_COPY_VALUES,&edgeNumberingIS_new2old);
  IS edgeNumberingIS_global_new2old;
  ISAllGather(edgeNumberingIS_new2old,&edgeNumberingIS_global_new2old);
  const PetscInt *edgeNumbering_global_new2old;

  
  valarray<int> edgeNumbering_old2new_global(mesh.nEdges_global);
  ISGetIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);
  for (int ig=0;ig<mesh.nEdges_global;ig++)
    {
      edgeNumbering_old2new_global[edgeNumbering_global_new2old[ig]] = ig;
    }
  
  //create  array with (new edge) --> (new node 0, new node 1)
  //and map from (new node 0, new node 1) --> (new global edge)
  valarray<int> edgeNodesArray_newNodesAndEdges(2*mesh.nEdges_global);
  map<NodeTuple<2>,int > nodesEdgeMap_global_new;
  for (int ig = 0; ig < mesh.nEdges_global; ig++)
    {
      int nodes[2];
      const int edge_old = edgeNumbering_global_new2old[ig];
      nodes[0] = edgeNodesArray_new[edge_old*2+0];
      nodes[1] = edgeNodesArray_new[edge_old*2+1];
      NodeTuple<2> et(nodes);
      edgeNodesArray_newNodesAndEdges[ig*2+0] = edgeNodesArray_new[edge_old*2+0];
      edgeNodesArray_newNodesAndEdges[ig*2+1] = edgeNodesArray_new[edge_old*2+1];
      assert(nodesEdgeMap_global_new.find(et) == nodesEdgeMap_global_new.end());
      nodesEdgeMap_global_new[et] = ig;
    }

  //5. At this point we have new, renumbered and sorted the global element, and node based information, and
  //we have it all on each processor so we can add overlap
  set<int> elements_overlap,nodes_overlap,elementBoundaries_overlap,edges_overlap;
  for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
    {
      for (int nN=0;nN<mesh.nNodes_element;nN++)
        {
          int nN_global = elementNodesArray_new[eN*mesh.nNodes_element+nN];
          if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
            nodes_overlap.insert(nN_global);
        }
      for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
        {
          int ebN_global=elementBoundariesArray_new[eN*mesh.nElementBoundaries_element+ebN];
          if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
            elementBoundaries_overlap.insert(ebN_global);
        }
      for (int nN0=0;nN0<mesh.nNodes_element;nN0++)//assume all nodes in element may be edges
        for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
          {
            int nodes[2];
            nodes[0] = elementNodesArray_new[eN*mesh.nNodes_element+nN0];
            nodes[1] = elementNodesArray_new[eN*mesh.nNodes_element+nN1];
            NodeTuple<2> et(nodes);
            const int nN0_global = elementNodesArray_new[eN*mesh.nNodes_element+nN0];
            const int nN1_global = elementNodesArray_new[eN*mesh.nNodes_element+nN1];
            if (nN0_global < nodeOffsets_new[rank] || nN0_global >= nodeOffsets_new[rank+1])
              assert(nodes_overlap.find(nN0_global) != nodes_overlap.end());
            if (nN1_global < nodeOffsets_new[rank] || nN1_global >= nodeOffsets_new[rank+1])
              assert(nodes_overlap.find(nN1_global) != nodes_overlap.end());
            if (nodesEdgeMap_global_new.find(et) != nodesEdgeMap_global_new.end())
              {
                const int edge_global = nodesEdgeMap_global_new[et];
                if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                  edges_overlap.insert(edge_global);
              }
          }//edges
    }

  if (nElements_overlap > 0)
    {
      //get all elements in the node stars and their nodes
      // for (int nN_new=nodeOffsets_new[rank]; nN_new < nodeOffsets_new[rank+1]; nN_new++)
      //   {
      //     int nN = nodeNumbering_global_new2old[nN_new];
      //     for (int offset =mesh.nodeElementOffsets[nN];offset<mesh.nodeElementOffsets[nN+1];offset++)
      //       {
      //         int eN = mesh.nodeElementsArray[offset];
      //         int eN_new = elementNumbering_global_old2new[eN];
      //         if (eN_new < elementOffsets_new[rank] or eN_new >= elementOffsets_new[rank+1])
      //           {
      //             elements_overlap.insert(eN_new);
      //             for (int nN_element=0;nN_element<mesh.nNodes_element;nN_element++)
      //               {
      //                 int nN_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN_element];
      //                 if (nN_global < nodeOffsets_new[rank] or nN_global >= nodeOffsets_new[rank+1])
      //                   nodes_overlap.insert(nN_global);
      //               }
      //           }
      //       }
      //   }
      //get all elements in the node stars of owned nodes and those elements' nodes and faces and edges
      for (set<int>::iterator nN=nodes_subdomain_owned.begin();nN != nodes_subdomain_owned.end();nN++)
        {
          for (int offset =mesh.nodeElementOffsets[*nN];offset<mesh.nodeElementOffsets[(*nN)+1];offset++)
            {
              int eN = mesh.nodeElementsArray[offset];
              int eN_new = elementNumbering_global_old2new[eN];
              if (eN_new < elementOffsets_new[rank] or eN_new >= elementOffsets_new[rank+1])
                {
                  elements_overlap.insert(eN_new);
                  for (int nN_element=0;nN_element<mesh.nNodes_element;nN_element++)
                    {
                      int nN_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN_element];
                      if (nN_global < nodeOffsets_new[rank] or nN_global >= nodeOffsets_new[rank+1])
                        nodes_overlap.insert(nN_global);
                    }
                  for (int ebN=0; ebN<mesh.nElementBoundaries_element;ebN++)
                    {
                      int ebN_global=elementBoundariesArray_new[eN_new*mesh.nElementBoundaries_element+ebN];
                      if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                        elementBoundaries_overlap.insert(ebN_global);
                    }
                  //edges too
                  for (int nN0=0;nN0<mesh.nNodes_element;nN0++)
                    for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
                      {
                        int nodes[2];
                        nodes[0] = elementNodesArray_new[eN_new*mesh.nNodes_element+nN0];
                        nodes[1] = elementNodesArray_new[eN_new*mesh.nNodes_element+nN1];
                        NodeTuple<2> et(nodes);
                        const int nN0_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN0];
                        const int nN1_global = elementNodesArray_new[eN_new*mesh.nNodes_element+nN1];
                        if (nodesEdgeMap_global_new.find(et) != nodesEdgeMap_global_new.end())
                          {
                            const int edge_global = nodesEdgeMap_global_new[et];
                            if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                              edges_overlap.insert(edge_global);
                          }
                      }//edges
                }
            }
        }
      //get all the element neighbors
      for(int eN=elementOffsets_new[rank];eN<elementOffsets_new[rank+1];eN++)
        for(int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
          {
            int eN_ebN = elementNeighborsArray_new[eN*mesh.nElementBoundaries_element+ebN];
            if (eN_ebN >= 0 && 
                (eN_ebN < elementOffsets_new[rank] || eN_ebN >= elementOffsets_new[rank+1]))
              {
                elements_overlap.insert(eN_ebN);
                for (int nN=0;nN<mesh.nNodes_element;nN++)
                  {
                    int nN_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN];
                    if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
                      nodes_overlap.insert(nN_global);
                  }
                for (int ebN_element=0; ebN_element<mesh.nElementBoundaries_element;ebN_element++)
                  {
                    int ebN_global=elementBoundariesArray_new[eN_ebN*mesh.nElementBoundaries_element+ebN_element];
                    if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                      elementBoundaries_overlap.insert(ebN_global);
                  }
                //edges too
                for (int nN0=0;nN0<mesh.nNodes_element;nN0++)
                  for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
                    {
                      int nodes[2];
                      nodes[0] = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN0];
                      nodes[1] = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN1];
                      NodeTuple<2> et(nodes);
                      const int nN0_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN0];
                      const int nN1_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN1];
                      bool foundEdge = false;
                      const int edge_global = nodesEdgeMap_global_new[et];
                      if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                        edges_overlap.insert(edge_global);
                    }//edges
              }
          }
    }
  for (int layer=1;layer<nElements_overlap;layer++)
    {
      for (set<int>::iterator eN_p=elements_overlap.begin();eN_p != elements_overlap.end();eN_p++)
        {
          int eN_global = *eN_p;
          for(int ebN=0;ebN<mesh.nElementBoundaries_element;ebN++)
            {
              int eN_ebN = elementNeighborsArray_new[eN_global*mesh.nElementBoundaries_element+ebN];
              if (eN_ebN >= 0 &&
                  (eN_ebN < elementOffsets_new[rank] || eN_ebN >= elementOffsets_new[rank+1]))
                {
                  elements_overlap.insert(eN_ebN);
                  for (int nN=0;nN<mesh.nNodes_element;nN++)
                    {
                      int nN_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN];
                      if (nN_global < nodeOffsets_new[rank] || nN_global >= nodeOffsets_new[rank+1])
                        nodes_overlap.insert(nN_global);
                    }
                  for (int ebN_element=0; ebN_element<mesh.nElementBoundaries_element;ebN_element++)
                    {
                      int ebN_global=elementBoundariesArray_new[eN_ebN*mesh.nElementBoundaries_element+ebN_element];
                      if (ebN_global < elementBoundaryOffsets_new[rank] || ebN_global >= elementBoundaryOffsets_new[rank+1])
                        elementBoundaries_overlap.insert(ebN_global);
                    }
                  //edges too
                  for (int nN0=0;nN0<mesh.nNodes_element;nN0++)
                    for (int nN1=nN0+1; nN1<mesh.nNodes_element;nN1++)
                      {
                        int nodes[2];
                        nodes[0] = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN0];
                        nodes[1] = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN1];
                        NodeTuple<2> et(nodes);
                        const int nN0_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN0];
                        const int nN1_global = elementNodesArray_new[eN_ebN*mesh.nNodes_element+nN1];
                        bool foundEdge = false;
                        if (nodesEdgeMap_global_new.find(et) != nodesEdgeMap_global_new.end())
                          {
                            const int edge_global = nodesEdgeMap_global_new[et];
                            if (edge_global < edgeOffsets_new[rank] || edge_global >= edgeOffsets_new[rank+1])
                              edges_overlap.insert(edge_global);
                          }
                      }//edges
                }
            }
        }            
    }
  //
  //6. Now build subdomain mesh
  //
  //set what we know
  if (mesh.subdomainp == NULL)
    mesh.subdomainp = new Mesh();
  mesh.subdomainp->nElements_global = nElements_subdomain_new[rank] + elements_overlap.size();
  mesh.subdomainp->nNodes_global = nNodes_subdomain_new[rank] + nodes_overlap.size();
  //better be true ... check below
  mesh.subdomainp->nElementBoundaries_global = nElementBoundaries_subdomain_new[rank]+elementBoundaries_overlap.size();
  mesh.subdomainp->nEdges_global = nEdges_subdomain_new[rank] + edges_overlap.size();

  mesh.subdomainp->nNodes_element = mesh.nNodes_element;
  mesh.subdomainp->nNodes_elementBoundary = mesh.nNodes_elementBoundary;
  mesh.subdomainp->nElementBoundaries_element = mesh.nElementBoundaries_element;
  //load the elements, nodes. and elementBoundaries
  valarray<int> nodeNumbering_subdomain2global(mesh.subdomainp->nNodes_global);
  map<int,int> nodeNumbering_global2subdomain;
  mesh.subdomainp->nodeArray = new double[mesh.subdomainp->nNodes_global*3];
  mesh.subdomainp->nodeMaterialTypes = new int[mesh.subdomainp->nNodes_global];
  for(int nN=0;nN<nNodes_subdomain_new[rank];nN++)
    {
      int nN_global = nN + nodeOffsets_new[rank];
      nodeNumbering_subdomain2global[nN] = nN_global;
      nodeNumbering_global2subdomain[nN_global] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = nodeArray_new[nN_global*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = nodeArray_new[nN_global*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = nodeArray_new[nN_global*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN]= nodeMaterialTypes_new[nN_global];
    }
  //note: sets in C++ are sorted so the overlap is laid out in
  //contiguous chunks corresponding to the partitions
  set<int>::iterator nN_p=nodes_overlap.begin();
  for(int nN=nNodes_subdomain_new[rank];nN < nNodes_subdomain_new[rank] + int(nodes_overlap.size()); nN++)
    {
      int nN_global = *nN_p++;
      nodeNumbering_subdomain2global[nN] = nN_global;
      nodeNumbering_global2subdomain[nN_global] = nN;
      mesh.subdomainp->nodeArray[nN*3+0] = nodeArray_new[nN_global*3+0];
      mesh.subdomainp->nodeArray[nN*3+1] = nodeArray_new[nN_global*3+1];
      mesh.subdomainp->nodeArray[nN*3+2] = nodeArray_new[nN_global*3+2];
      mesh.subdomainp->nodeMaterialTypes[nN]= nodeMaterialTypes_new[nN_global];
    }

  //see what we must/can set for element boundaries here
  valarray<int> elementBoundaryNumbering_subdomain2global(mesh.subdomainp->nElementBoundaries_global);
  map<int,int> elementBoundaryNumbering_global2subdomain;
  for (int ebN=0; ebN < nElementBoundaries_subdomain_new[rank]; ebN++)
    {
      int ebN_global = ebN + elementBoundaryOffsets_new[rank];
      elementBoundaryNumbering_subdomain2global[ebN]=ebN_global;
      elementBoundaryNumbering_global2subdomain[ebN_global] = ebN;
    }
  set<int>::iterator ebN_p = elementBoundaries_overlap.begin();
  for(int ebN=nElementBoundaries_subdomain_new[rank];ebN < nElementBoundaries_subdomain_new[rank] + int(elementBoundaries_overlap.size()); ebN++)
    {
      int ebN_global = *ebN_p++;
      elementBoundaryNumbering_subdomain2global[ebN] = ebN_global;
      elementBoundaryNumbering_global2subdomain[ebN_global] = ebN;
    }


  mesh.subdomainp->elementNodesArray = new int[mesh.subdomainp->nElements_global*mesh.subdomainp->nNodes_element];
  mesh.subdomainp->elementMaterialTypes = new int[mesh.subdomainp->nElements_global];
  //try to use elementBoundariesArray to set unique element boundary id below
  mesh.subdomainp->elementBoundariesArray = 
    new int[mesh.subdomainp->nElements_global*mesh.subdomainp->nElementBoundaries_element];
  valarray<int> elementNumbering_subdomain2global(mesh.subdomainp->nElements_global);
  for (int eN=0;eN<nElements_subdomain_new[rank];eN++)
    {
      int eN_global = eN+elementOffsets_new[rank];
      elementNumbering_subdomain2global[eN] = eN_global;
      mesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypes_new[eN_global];
      for (int nN=0;nN<mesh.subdomainp->nNodes_element;nN++)
        {
          mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+nN] = 
            nodeNumbering_global2subdomain[elementNodesArray_new[eN_global*mesh.nNodes_element + nN]];
        }
      for (int ebN=0;ebN<mesh.subdomainp->nElementBoundaries_element;ebN++)
        mesh.subdomainp->elementBoundariesArray[eN*mesh.subdomainp->nElementBoundaries_element+ebN] =
          elementBoundaryNumbering_global2subdomain[elementBoundariesArray_new[eN_global*mesh.nElementBoundaries_element + ebN]];
    }
  set<int>::iterator eN_p=elements_overlap.begin();
  for(int eN=nElements_subdomain_new[rank];eN < nElements_subdomain_new[rank]+int(elements_overlap.size());eN++)
    {
      int eN_global = *eN_p++;
      mesh.subdomainp->elementMaterialTypes[eN] = elementMaterialTypes_new[eN_global];
      elementNumbering_subdomain2global[eN] = eN_global;
      for (int nN=0;nN<mesh.subdomainp->nNodes_element;nN++)
        {
          mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+nN] = 
            nodeNumbering_global2subdomain[elementNodesArray_new[eN_global*mesh.nNodes_element + nN]];
        }
      for (int ebN=0;ebN<mesh.subdomainp->nElementBoundaries_element;ebN++)
        mesh.subdomainp->elementBoundariesArray[eN*mesh.subdomainp->nElementBoundaries_element+ebN] =
          elementBoundaryNumbering_global2subdomain[elementBoundariesArray_new[eN_global*mesh.nElementBoundaries_element + ebN]];
    }

  //set edge information, main piece necessary for setting subdomain data structures consistently is a map
  //from (nN0,nN1) --> edge_subdomain
  valarray<int> edgeNumbering_subdomain2global(mesh.subdomainp->nEdges_global);
  mesh.subdomainp->edgeNodesArray = new int[mesh.subdomainp->nEdges_global*2];
  for (int i=0; i < nEdges_subdomain_new[rank]; i++)
    {
      const int ig = i+edgeOffsets_new[rank];
      const int nN0_global = edgeNodesArray_newNodesAndEdges[ig*2+0];
      const int nN1_global = edgeNodesArray_newNodesAndEdges[ig*2+1];
      //mwf todo double check can always count on having nodes on this processor
      assert(nodeNumbering_global2subdomain.find(nN0_global) != nodeNumbering_global2subdomain.end());
      assert(nodeNumbering_global2subdomain.find(nN1_global) != nodeNumbering_global2subdomain.end());
      const int nN0_subdomain = nodeNumbering_global2subdomain[nN0_global];
      const int nN1_subdomain = nodeNumbering_global2subdomain[nN1_global];
      mesh.subdomainp->edgeNodesArray[2*i+0]=nN0_subdomain;
      mesh.subdomainp->edgeNodesArray[2*i+1]=nN1_subdomain;
      edgeNumbering_subdomain2global[i] = ig;
    }
  set<int>::iterator edge_p = edges_overlap.begin();
  for (int i=nEdges_subdomain_new[rank]; i < nEdges_subdomain_new[rank] + int(edges_overlap.size()); i++)
    {
      const int ig =*edge_p++;
      const int nN0_global = edgeNodesArray_newNodesAndEdges[ig*2+0];
      const int nN1_global = edgeNodesArray_newNodesAndEdges[ig*2+1];
      //mwf todo make sure always have nodes for the edge on this processor 
      const int nN0_subdomain = nodeNumbering_global2subdomain[nN0_global];
      const int nN1_subdomain = nodeNumbering_global2subdomain[nN1_global];
      mesh.subdomainp->edgeNodesArray[2*i+0]=nN0_subdomain;
      mesh.subdomainp->edgeNodesArray[2*i+1]=nN1_subdomain;
      edgeNumbering_subdomain2global[i] = ig;
    }
  //now fill in rest of boundary information, etc
  mesh.subdomainp->px = mesh.px;
  mesh.subdomainp->py = mesh.py;
  mesh.subdomainp->pz = mesh.pz; 

  if (mesh.px != 0)    
    {  
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_tetrahedron(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_NURBS(*mesh.subdomainp);
      allocateGeometricInfo_NURBS(*mesh.subdomainp);
      computeGeometricInfo_NURBS(*mesh.subdomainp);

    }   
  else if (mesh.subdomainp->nNodes_element == 2)
    {
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_edge(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_edge(*mesh.subdomainp);
      allocateGeometricInfo_edge(*mesh.subdomainp);
      computeGeometricInfo_edge(*mesh.subdomainp);
    }

  else if (mesh.subdomainp->nNodes_element == 3)
    {
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_triangle(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_triangle(*mesh.subdomainp);
      allocateGeometricInfo_triangle(*mesh.subdomainp);
      computeGeometricInfo_triangle(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 4)
    {
      //constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_tetrahedron(*mesh.subdomainp);
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_tetrahedron(*mesh.subdomainp);
      allocateGeometricInfo_tetrahedron(*mesh.subdomainp);
      computeGeometricInfo_tetrahedron(*mesh.subdomainp);
    }
  else if (mesh.subdomainp->nNodes_element == 8)
    {
      constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_hexahedron(*mesh.subdomainp);
      allocateGeometricInfo_hexahedron(*mesh.subdomainp);
      computeGeometricInfo_hexahedron(*mesh.subdomainp);
    }
 
  if (mesh.elementBoundaryMaterialTypes != NULL)
    {
      assert(mesh.elementBoundariesArray != NULL);
      //todo need to check that local element boundary numbering for 
      //element boundaries stays the same 
      for (int eN=0;eN<mesh.subdomainp->nElements_global;eN++)
        {
          int eN_global_new = elementNumbering_subdomain2global[eN];
          int eN_global_old = elementNumbering_global_new2old[eN_global_new];
          for (int ebN_element = 0; ebN_element < mesh.nElementBoundaries_element; ebN_element++)
            {
              int ebN_global_old = mesh.elementBoundariesArray[eN_global_old*mesh.nElementBoundaries_element+ebN_element];
              int ebN_subdomain = mesh.subdomainp->elementBoundariesArray[eN*mesh.nElementBoundaries_element+ebN_element];
              mesh.subdomainp->elementBoundaryMaterialTypes[ebN_subdomain] = mesh.elementBoundaryMaterialTypes[ebN_global_old];
            }
        }
    }
  //now we've got the old mesh in the old ordering and the subdomain mesh in the new ordering
  //the first chunk of nodes and elements are the owned elements so we need to know how many of those there are
  //and the offset of the first one so we can compute subdomain2global
  mesh.nodeOffsets_subdomain_owned = new int[size+1];
  mesh.elementOffsets_subdomain_owned = new int[size+1];
  mesh.elementBoundaryOffsets_subdomain_owned = new int[size+1];
  mesh.edgeOffsets_subdomain_owned = new int[size+1];
  for (int sdN=0;sdN<size+1;sdN++)
    {
      mesh.nodeOffsets_subdomain_owned[sdN] = nodeOffsets_new[sdN];
      mesh.elementOffsets_subdomain_owned[sdN] = elementOffsets_new[sdN];
      mesh.elementBoundaryOffsets_subdomain_owned[sdN] = elementBoundaryOffsets_new[sdN];
      mesh.edgeOffsets_subdomain_owned[sdN] = edgeOffsets_new[sdN];
    }
  //we also need the subdomain 2 new global mappings
  mesh.nodeNumbering_subdomain2global = new int[mesh.subdomainp->nNodes_global];
  for (int nN=0;nN<mesh.subdomainp->nNodes_global;nN++)
    mesh.nodeNumbering_subdomain2global[nN] = nodeNumbering_subdomain2global[nN];
  mesh.elementNumbering_subdomain2global = new int[mesh.subdomainp->nElements_global];
  for (int eN=0;eN<mesh.subdomainp->nElements_global;eN++)
    mesh.elementNumbering_subdomain2global[eN] = elementNumbering_subdomain2global[eN];
  mesh.elementBoundaryNumbering_subdomain2global = new int[mesh.subdomainp->nElementBoundaries_global];
  for (int ebN=0;ebN<mesh.subdomainp->nElementBoundaries_global;ebN++)
    mesh.elementBoundaryNumbering_subdomain2global[ebN] = elementBoundaryNumbering_subdomain2global[ebN];
  mesh.edgeNumbering_subdomain2global = new int[mesh.subdomainp->nEdges_global];
  for (int i=0; i< mesh.subdomainp->nEdges_global; i++)
    mesh.edgeNumbering_subdomain2global[i] = edgeNumbering_subdomain2global[i];

  //
  //go ahead and renumber global mesh
  //
  //std::cout<<"nodeArray"<<std::endl;
  for (int i=0;i<mesh.nNodes_global*3;i++)
    mesh.nodeArray[i] = nodeArray_new[i];
  //std::cout<<"elementNodesArray"<<std::endl;
  for (int i=0;i<mesh.nElements_global*mesh.nNodes_element;i++)
    mesh.elementNodesArray[i] = elementNodesArray_new[i];
  //std::cout<<"elementBoundaryNodesArray"<<std::endl;
  for (int i=0;i<mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary;i++)
    mesh.elementBoundaryNodesArray[i] = elementBoundaryNodesArray_new[i];
  //std::cout<<"edgeNodesArray"<<std::endl;
  for (int i=0;i<mesh.nEdges_global*2;i++)
    mesh.edgeNodesArray[i] = edgeNodesArray_newNodesAndEdges[i];
  //std::cout<<"nodeStarArray"<<std::endl;
  for (int i=0;i<mesh.nodeStarOffsets[mesh.nNodes_global];i++)
    mesh.nodeStarArray[i] = nodeStarArray_new[i];
  //std::cout<<"elementNeighborsArray"<<std::endl;
  for (int i=0; i<mesh.nElements_global*mesh.nElementBoundaries_element; i++)
    mesh.elementNeighborsArray[i] = elementNeighborsArray_new[i];
  //std::cout<<"elementBoundariesArray"<<std::endl;
  for (int i=0; i<mesh.nElements_global*mesh.nElementBoundaries_element; i++)
    mesh.elementBoundariesArray[i] = elementBoundariesArray_new[i];
  //std::cout<<"elementBoundaryElementsArray"<<std::endl;
  for (int i=0; i<mesh.nElementBoundaries_global*2; i++)
    mesh.elementBoundaryElementsArray[i] = elementBoundaryElementsArray_new[i];
  //std::cout<<"elementBoundaryLocalElementBoundariesArray"<<std::endl;
  for (int i=0; i<mesh.nElementBoundaries_global*2; i++)
    mesh.elementBoundaryLocalElementBoundariesArray[i] = elementBoundaryLocalElementBoundariesArray_new[i];
  //
  //need material properties too
  //
  if (mesh.elementMaterialTypes != NULL)
    {
      //std::cout<<"elementMaterialTypes"<<std::endl;
      for (int i=0; i < mesh.nElements_global; i++)
        mesh.elementMaterialTypes[i] = elementMaterialTypes_new[i];
    }
  if (mesh.elementBoundaryMaterialTypes != NULL)
    {
      //std::cout<<"elementBoundaryMaterialTypes"<<std::endl;
      for (int i=0; i < mesh.nElementBoundaries_global; i++)
        mesh.elementBoundaryMaterialTypes[i] = elementBoundaryMaterialTypes_new[i];
    }
   
  ISRestoreIndices(elementNumberingIS_global_old2new,&elementNumbering_global_old2new);

  ISDestroy(&elementPartitioningIS_new);
  ISDestroy(&elementNumberingIS_subdomain_old2new);
  ISDestroy(&elementNumberingIS_global_old2new);

  ISRestoreIndices(nodeNumberingIS_global_new2old,&nodeNumbering_global_new2old);

  ISDestroy(&nodeNumberingIS_new2old);
  ISDestroy(&nodeNumberingIS_global_new2old);

  ISRestoreIndices(elementBoundaryNumberingIS_global_new2old,&elementBoundaryNumbering_global_new2old);

  ISDestroy(&elementBoundaryNumberingIS_new2old);
  ISDestroy(&elementBoundaryNumberingIS_global_new2old);

  ISRestoreIndices(edgeNumberingIS_global_new2old,&edgeNumbering_global_new2old);

  ISDestroy(&edgeNumberingIS_new2old);
  ISDestroy(&edgeNumberingIS_global_new2old);

  return 0;
}

int buildQuadraticSubdomain2GlobalMappings_1d(Mesh& mesh, 
                                              const int *elementOffsets_subdomain_owned,
                                              const int *nodeOffsets_subdomain_owned,
                                              const int *elementNumbering_subdomain2global,
                                              const int *nodeNumbering_subdomain2global,
                                              int& nDOF_all_processes,//total number of dofs in whole domain
                                              int& nDOF_subdomain,//total number of dofs in sub-domain
                                              int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                              int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                              int * subdomain_l2g, //local to global dof mapping on subdomain
                                              int* subdomain2global,//subdomain dof to global (parallel) numbering
                                              double * lagrangeNodesArray)//location of nodes corresponding to dofs
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);


  //In 1d the quadratic dofs can be associated with nodes and elements
  //assuming have ownership info and consistent local/global mappings for nodes, elements
  //build a mapping from local quadratic dofs to global quadratic dofs for petsc
  //assume a processor owns a dof if it owns that node or element
  assert(elementOffsets_subdomain_owned);
  assert(nodeOffsets_subdomain_owned);
  assert(elementNumbering_subdomain2global);
  assert(nodeNumbering_subdomain2global);
  assert(mesh.subdomainp);

  const int nNodes_owned    = nodeOffsets_subdomain_owned[rank+1]-nodeOffsets_subdomain_owned[rank];
  const int nElements_owned = elementOffsets_subdomain_owned[rank+1]-elementOffsets_subdomain_owned[rank];
  
  //start with a logical global ordering of dofs as
  //[global nodes, global elements]
  //want to create global numbering
  //[nodes proc 0,elements proc 0,nodes proc 1, elements proc 1,...]
  valarray<int> quadraticNumbering_new2old(nNodes_owned + nElements_owned);
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global = nodeNumbering_subdomain2global[nN];
      quadraticNumbering_new2old[nN] = dof_global;
    }
  for (int eN = 0; eN < nElements_owned; eN++)
    {
      int dof_global = mesh.nNodes_global + elementNumbering_subdomain2global[eN];
      quadraticNumbering_new2old[nNodes_owned + eN] = dof_global;
    }

  //build an index set for new numbering
  IS quadraticNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nNodes_owned + nElements_owned,&quadraticNumbering_new2old[0],PETSC_COPY_VALUES,&quadraticNumberingIS_new2old);
  IS quadraticNumberingIS_global_new2old;
  ISAllGather(quadraticNumberingIS_new2old,&quadraticNumberingIS_global_new2old);
  //get old 2 new mapping for dofs
  const PetscInt *quadraticNumbering_global_new2old;
  valarray<int> quadraticNumbering_old2new_global(mesh.nNodes_global + mesh.nElements_global);
  ISGetIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  for (int id = 0; id < mesh.nNodes_global + mesh.nElements_global; id++)
    quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]] = id;
  //   //mwf debug
  //   for (int id = 0; id < mesh.nNodes_global + mesh.nElements_global; id++)
  //     {
  //       std::cout<<" rank= "<<rank<<" build c0p2 mappings new2old["<<id<<"]= "<<quadraticNumbering_global_new2old[id]
  // 	       <<" old2new["<<quadraticNumbering_global_new2old[id]<<"]= "<<quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]]
  // 	       <<std::endl;
  //     }
  assert(offsets_subdomain_owned);
  assert(subdomain2global);
  assert(subdomain_l2g);
  assert(lagrangeNodesArray);
  nDOF_all_processes = mesh.nNodes_global+mesh.nElements_global;
  nDOF_subdomain = mesh.subdomainp->nNodes_global+mesh.subdomainp->nElements_global;
  max_dof_neighbors = 2*mesh.max_nNodeNeighbors_node;

  for (int sdN=0; sdN < size+1; sdN++)
    offsets_subdomain_owned[sdN] = nodeOffsets_subdomain_owned[sdN]+elementOffsets_subdomain_owned[sdN];

  //loop through owned and ghost dofs build subdomain mapping by 
  //going from old --> new
  int localOffset = 0;
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nNodes_owned;
  for (int eN = 0; eN < nElements_owned; eN++)
    {
      int dof_global_old = mesh.nNodes_global + elementNumbering_subdomain2global[eN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + eN] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" eN= "<<eN<<" local dof= "<<localOffset+eN<<" nElements_owned= "<<nElements_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nElements_owned;
  for (int nN = nNodes_owned; nN < mesh.subdomainp->nNodes_global; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN -nNodes_owned] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN-nNodes_owned<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += mesh.subdomainp->nNodes_global - nNodes_owned;
  for (int eN = nElements_owned; eN < mesh.subdomainp->nElements_global; eN++)
    {
      int dof_global_old = mesh.nNodes_global + elementNumbering_subdomain2global[eN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + eN-nElements_owned] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" eN= "<<eN<<" local dof= "<<localOffset+eN-nElements_owned<<" nElements_owned= "<<nElements_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  //setup local to global mapping on the subdomain for finite elements
  const int nDOF_element = 3;
  const int ghostNodeOffset = nNodes_owned + nElements_owned;
  const int ghostElementOffset = ghostNodeOffset + mesh.subdomainp->nNodes_global-nNodes_owned;
 
  //for lagrange nodes
  int nN_global_subdomain[2];
  for (int eN=0; eN < mesh.subdomainp->nElements_global; eN++)
    {
      for (int nN=0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          nN_global_subdomain[nN] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          if (nN_global_subdomain[nN] < nNodes_owned)
            subdomain_l2g[eN*nDOF_element+nN] = nN_global_subdomain[nN];
          else
            subdomain_l2g[eN*nDOF_element+nN] = nN_global_subdomain[nN]-nNodes_owned + ghostNodeOffset;
          //vertex dof
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+nN]*3+eI] = mesh.subdomainp->nodeArray[3*nN_global_subdomain[nN]+eI];
        }
      if (eN < nElements_owned)
        subdomain_l2g[eN*nDOF_element+2] = nNodes_owned + eN;
      else
        subdomain_l2g[eN*nDOF_element+2] = eN - nElements_owned + ghostElementOffset; 
      //vertex dof
      for (int eI=0; eI < 3; eI++)
        lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+2]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[3*nN_global_subdomain[0]+eI]+mesh.subdomainp->nodeArray[3*nN_global_subdomain[1]+eI]);
    }

  ISRestoreIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  ISDestroy(&quadraticNumberingIS_new2old);
  ISDestroy(&quadraticNumberingIS_global_new2old);

  return 0;
}

int buildQuadraticSubdomain2GlobalMappings_2d(Mesh& mesh, 
                                              const int *elementBoundaryOffsets_subdomain_owned,
                                              const int *nodeOffsets_subdomain_owned,
                                              const int *elementBoundaryNumbering_subdomain2global,
                                              const int *nodeNumbering_subdomain2global,
                                              int& nDOF_all_processes,//total number of dofs in whole domain
                                              int& nDOF_subdomain,//total number of dofs in sub-domain
                                              int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                              int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                              int *subdomain_l2g, //local to global dof mapping on subdomain
                                              int *subdomain2global,//subdomain dof to global (parallel) numbering
                                              double * lagrangeNodesArray)//location of nodes corresponding to dofs
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  //In 2d the quadratic dofs can be associated with nodes and element Boundaries
  //assuming have ownership info and consistent local/global mappings for nodes, elementBoundaries
  //build a mapping from local quadratic dofs to global quadratic dofs for petsc
  //assume a processor owns a dof if it owns that node or element
  assert(elementBoundaryOffsets_subdomain_owned);
  assert(nodeOffsets_subdomain_owned);
  assert(elementBoundaryNumbering_subdomain2global);
  assert(nodeNumbering_subdomain2global);
  assert(mesh.subdomainp);

  const int nNodes_owned    = nodeOffsets_subdomain_owned[rank+1]-nodeOffsets_subdomain_owned[rank];
  const int nElementBoundaries_owned = elementBoundaryOffsets_subdomain_owned[rank+1]-elementBoundaryOffsets_subdomain_owned[rank];
  const int nDOFs_owned = nNodes_owned + nElementBoundaries_owned;
  const int nDOFs_global= mesh.nNodes_global + mesh.nElementBoundaries_global;
  //start with a logical global ordering of dofs as
  //[global nodes, global elementBoundaries]
  //want to create global numbering
  //[nodes proc 0,elementBoundaries proc 0,nodes proc 1, elementBoundaries proc 1,...]
  valarray<int> quadraticNumbering_new2old(nDOFs_owned);
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global = nodeNumbering_subdomain2global[nN];
      quadraticNumbering_new2old[nN] = dof_global;
    }
  for (int ebN = 0; ebN < nElementBoundaries_owned; ebN++)
    {
      int dof_global = mesh.nNodes_global + elementBoundaryNumbering_subdomain2global[ebN];
      quadraticNumbering_new2old[nNodes_owned + ebN] = dof_global;
    }

  //build an index set for new numbering
  IS quadraticNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nDOFs_owned,&quadraticNumbering_new2old[0],PETSC_COPY_VALUES,&quadraticNumberingIS_new2old);
  IS quadraticNumberingIS_global_new2old;
  ISAllGather(quadraticNumberingIS_new2old,&quadraticNumberingIS_global_new2old);
  //get old 2 new mapping for dofs
  const PetscInt *quadraticNumbering_global_new2old;
  valarray<int> quadraticNumbering_old2new_global(nDOFs_global);
  ISGetIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  for (int id = 0; id < nDOFs_global; id++)
    quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]] = id;
  //mwf debug
  //for (int id = 0; id < nDOFs_global; id++)
  //  {
  //    std::cout<<" rank= "<<rank<<" build 2d c0p2 mappings new2old["<<id<<"]= "<<quadraticNumbering_global_new2old[id]
  //	       <<" old2new["<<quadraticNumbering_global_new2old[id]<<"]= "<<quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]]
  //	       <<std::endl;
  //  }
  assert(offsets_subdomain_owned);
  assert(subdomain2global);
  assert(subdomain_l2g);
  assert(lagrangeNodesArray);
  for (int sdN=0; sdN < size+1; sdN++)
    offsets_subdomain_owned[sdN] = nodeOffsets_subdomain_owned[sdN]+elementBoundaryOffsets_subdomain_owned[sdN];
  nDOF_all_processes = mesh.nNodes_global+mesh.nElementBoundaries_global;
  nDOF_subdomain = mesh.subdomainp->nNodes_global+mesh.subdomainp->nElementBoundaries_global;
  max_dof_neighbors = 2*mesh.max_nNodeNeighbors_node;

  //loop through owned and ghost dofs build subdomain mapping by 
  //going from old --> new
  int localOffset = 0;
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nNodes_owned;
  for (int ebN = 0; ebN < nElementBoundaries_owned; ebN++)
    {
      int dof_global_old = mesh.nNodes_global + elementBoundaryNumbering_subdomain2global[ebN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + ebN] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" ebN= "<<ebN<<" local dof= "<<localOffset+ebN<<" nElementBoundaries_owned= "<<nElementBoundaries_owned
      //	       <<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nElementBoundaries_owned;
  for (int nN = nNodes_owned; nN < mesh.subdomainp->nNodes_global; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN -nNodes_owned] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN-nNodes_owned<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += mesh.subdomainp->nNodes_global - nNodes_owned;
  for (int ebN = nElementBoundaries_owned; ebN < mesh.subdomainp->nElementBoundaries_global; ebN++)
    {
      int dof_global_old = mesh.nNodes_global + elementBoundaryNumbering_subdomain2global[ebN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + ebN-nElementBoundaries_owned] = dof_global_new;
      //mwf debug
      //std::cout<<" rank= "<<rank<<" ebN= "<<ebN<<" local dof= "<<localOffset+ebN-nElementBoundaries_owned<<" nElementBoundaries_owned= "<<nElementBoundaries_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  //setup local to global mapping on the subdomain for finite elements
  const int nDOF_element = 6;
  const int ghostNodeOffset = nNodes_owned + nElementBoundaries_owned;
  const int ghostElementBoundaryOffset = ghostNodeOffset + mesh.subdomainp->nNodes_global-nNodes_owned;
 
  //for lagrange nodes
  int nN_global_subdomain[3];
  for (int eN=0; eN < mesh.subdomainp->nElements_global; eN++)
    {
      for (int nN=0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          nN_global_subdomain[nN] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          if (nN_global_subdomain[nN] < nNodes_owned)
            subdomain_l2g[eN*nDOF_element+nN] = nN_global_subdomain[nN];
          else
            subdomain_l2g[eN*nDOF_element+nN] = nN_global_subdomain[nN]-nNodes_owned + ghostNodeOffset;
          //vertex dof
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+nN]*3+eI] = mesh.subdomainp->nodeArray[3*nN_global_subdomain[nN]+eI];
        }
      for (int ebN=0; ebN < mesh.subdomainp->nElementBoundaries_element; ebN++)
        {
          //take into account numbering of edges according to
          //vertex they are across from 
          int ebN_global = mesh.subdomainp->elementBoundariesArray[eN*mesh.subdomainp->nElementBoundaries_element+((ebN+2)%3)];
          if (ebN_global < nElementBoundaries_owned)
            subdomain_l2g[eN*nDOF_element+3+ebN] = nNodes_owned + ebN_global;
          else
            subdomain_l2g[eN*nDOF_element+3+ebN] = ebN_global - nElementBoundaries_owned + ghostElementBoundaryOffset; 
          //center of edge dof
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+3+ebN]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[3*nN_global_subdomain[(ebN+0)%3]+eI]+mesh.subdomainp->nodeArray[3*nN_global_subdomain[(ebN+1)%3]+eI]);
        }
    }//eN

  ISRestoreIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  ISDestroy(&quadraticNumberingIS_new2old);
  ISDestroy(&quadraticNumberingIS_global_new2old);

  return 0;
}


int buildQuadraticSubdomain2GlobalMappings_3d(Mesh& mesh, 
                                              const int *edgeOffsets_subdomain_owned,
                                              const int *nodeOffsets_subdomain_owned,
                                              const int *edgeNumbering_subdomain2global,
                                              const int *nodeNumbering_subdomain2global,
                                              int& nDOF_all_processes,//total number of dofs in whole domain
                                              int& nDOF_subdomain,//total number of dofs in sub-domain
                                              int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                              int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                              int *subdomain_l2g, //local to global dof mapping on subdomain
                                              int *subdomain2global,//subdomain dof to global (parallel) numbering
                                              double * lagrangeNodesArray)//location of nodes corresponding to dofs
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  //In 2d the quadratic dofs can be associated with nodes and element Boundaries
  //assuming have ownership info and consistent local/global mappings for nodes, elementBoundaries
  //build a mapping from local quadratic dofs to global quadratic dofs for petsc
  //assume a processor owns a dof if it owns that node or element
  assert(edgeOffsets_subdomain_owned);
  assert(nodeOffsets_subdomain_owned);
  assert(edgeNumbering_subdomain2global);
  assert(nodeNumbering_subdomain2global);
  assert(mesh.subdomainp);

  const int nNodes_owned = nodeOffsets_subdomain_owned[rank+1]-nodeOffsets_subdomain_owned[rank];
  const int nEdges_owned = edgeOffsets_subdomain_owned[rank+1]-edgeOffsets_subdomain_owned[rank];
  const int nDOFs_owned = nNodes_owned + nEdges_owned;
  const int nDOFs_global= mesh.nNodes_global + mesh.nEdges_global;
  //start with a logical global ordering of dofs as
  //[global nodes, global edges]
  //want to create global numbering
  //[nodes proc 0,edges proc 0,nodes proc 1, edges proc 1,...]
  valarray<int> quadraticNumbering_new2old(nDOFs_owned);
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global = nodeNumbering_subdomain2global[nN];
      quadraticNumbering_new2old[nN] = dof_global;
    }
  for (int i = 0; i < nEdges_owned; i++)
    {
      int dof_global = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      quadraticNumbering_new2old[nNodes_owned + i] = dof_global;
    }

  //build an index set for new numbering
  IS quadraticNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nDOFs_owned,&quadraticNumbering_new2old[0],PETSC_COPY_VALUES,&quadraticNumberingIS_new2old);
  IS quadraticNumberingIS_global_new2old;
  ISAllGather(quadraticNumberingIS_new2old,&quadraticNumberingIS_global_new2old);
  //get old 2 new mapping for dofs
  const PetscInt *quadraticNumbering_global_new2old;
  valarray<int> quadraticNumbering_old2new_global(nDOFs_global);
  ISGetIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  for (int id = 0; id < nDOFs_global; id++)
    quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]] = id;
  //mwf debug
  //   if (rank == 0)
  //     {
  //       for (int id = 0; id < nDOFs_global; id++)
  // 	{
  // 	  std::cout<<" rank= "<<rank<<" build 2d c0p2 mappings new2old["<<id<<"]= "<<quadraticNumbering_global_new2old[id]
  // 		   <<" old2new["<<quadraticNumbering_global_new2old[id]<<"]= "<<quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]]
  // 		   <<std::endl;
  // 	}
  //     }
  assert(offsets_subdomain_owned);
  assert(subdomain2global);
  assert(subdomain_l2g);
  assert(lagrangeNodesArray);
  for (int sdN=0; sdN < size+1; sdN++)
    offsets_subdomain_owned[sdN] = nodeOffsets_subdomain_owned[sdN]+edgeOffsets_subdomain_owned[sdN];
  nDOF_all_processes = mesh.nNodes_global+mesh.nEdges_global;
  nDOF_subdomain = mesh.subdomainp->nNodes_global+mesh.subdomainp->nEdges_global;
  max_dof_neighbors = 2*mesh.max_nNodeNeighbors_node;

  //loop through owned and ghost dofs build subdomain mapping by 
  //going from old --> new
  int localOffset = 0;
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN] = dof_global_new;
      //mwf debug
      //      if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nNodes_owned;
  for (int i = 0; i < nEdges_owned; i++)
    {
      int dof_global_old = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	  std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i<<" nEdges_owned= "<<nEdges_owned
      // 		   <<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += nEdges_owned;
  for (int nN = nNodes_owned; nN < mesh.subdomainp->nNodes_global; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN -nNodes_owned] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN-nNodes_owned<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  localOffset += mesh.subdomainp->nNodes_global - nNodes_owned;
  for (int i = nEdges_owned; i < mesh.subdomainp->nEdges_global; i++)
    {
      int dof_global_old = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i-nEdges_owned] = dof_global_new;
      //       //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i-nEdges_owned<<" nEdges_owned= "<<nEdges_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
  //setup local to global mapping on the subdomain for finite elements
  const int nDOF_element = 10;
  const int ghostNodeOffset = nNodes_owned + nEdges_owned;
  const int ghostEdgeOffset = ghostNodeOffset + mesh.subdomainp->nNodes_global-nNodes_owned;
  //need mapping from nodes to edge to setup element based relationship
  map<NodeTuple<2>, int> nodesEdgeMap_subdomain;
  for (int i=0; i < mesh.subdomainp->nEdges_global; i++)
    {
      int nodes[2];
      nodes[0] = mesh.subdomainp->edgeNodesArray[2*i];
      nodes[1] = mesh.subdomainp->edgeNodesArray[2*i+1];
      NodeTuple<2> et(nodes);
      const int nN0 = mesh.subdomainp->edgeNodesArray[2*i];
      const int nN1 = mesh.subdomainp->edgeNodesArray[2*i+1];
      nodesEdgeMap_subdomain[et] = i;
    }
  for (int eN=0; eN < mesh.subdomainp->nElements_global; eN++)
    {
      int local_offset = 0;
      for (int nN=0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          //node_i --> unique subdomain node  0 <= i <= 3
          const int nN_global_subdomain = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (nN_global_subdomain < nNodes_owned)
            {
              subdomain_l2g[eN*nDOF_element + local_offset + nN] = nN_global_subdomain;
            }
          else
            {
              subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostNodeOffset + nN_global_subdomain - nNodes_owned;
            }
          //unique subdomain id --> unique cross processor id
          //subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = nodeNumbering_subdomain2global[nN_global_subdomain];
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" node= "<<nN<<" sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element + local_offset +nN]*3+eI] = mesh.subdomainp->nodeArray[nN_global_subdomain*3+eI];
        }
      local_offset += mesh.subdomainp->nNodes_element;
      //(node_i,node_{i+1}) --> unique subdomain edge, 0 <= i < 3
      for (int nN = 0; nN < mesh.subdomainp->nNodes_element-1; nN++)
        {
          const int nN_neig = (nN+1)%mesh.subdomainp->nNodes_element;
          const int nN_global_subdomain     = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          const int nN_neig_global_subdomain= mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN_neig];
          //unique edge subdomain id and global id
          int edge_subdomain = -1,edge_global=-1; 
          //see if edge is (nN,nN_neig) or vice versa
          int nodes[2];
          nodes[0] = nN_global_subdomain;
          nodes[1] = nN_neig_global_subdomain;
          NodeTuple<2> et(nodes);
          edge_subdomain = nodesEdgeMap_subdomain[et];
          edge_global    = edgeNumbering_subdomain2global[edge_subdomain];
          assert(edge_subdomain >= 0 && edge_global >= 0);
	  
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (edge_subdomain < nEdges_owned)
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = nNodes_owned + edge_subdomain;
          else
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostEdgeOffset + edge_subdomain - nEdges_owned;
          //unique subdomain id --> unique cross processor id
          //set above could do here subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global + edge_global;
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }  
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset+nN]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[nN_global_subdomain*3+eI]+
                                                                                           mesh.subdomainp->nodeArray[nN_neig_global_subdomain*3+eI]);
        }
      local_offset +=  mesh.subdomainp->nNodes_element-1;
      //(node_i,node_{i+2}) --> unique subdomain edge, 0 <= i < 2
      for (int nN = 0; nN < mesh.subdomainp->nNodes_element-2; nN++)
        {
          const int nN_neig = (nN+2)%mesh.subdomainp->nNodes_element;
          const int nN_global_subdomain     = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          const int nN_neig_global_subdomain= mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN_neig];
          //unique edge subdomain id and global id
          int edge_subdomain = -1, edge_global = -1;
          //see if edge is (nN,nN_neig) or vice versa
          int nodes[2];
          nodes[0] = nN_global_subdomain;
          nodes[1] = nN_neig_global_subdomain;
          NodeTuple<2> et(nodes);
          edge_subdomain = nodesEdgeMap_subdomain[et];
          edge_global =    edgeNumbering_subdomain2global[edge_subdomain];
          assert(edge_subdomain >= 0 && edge_global >= 0);
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (edge_subdomain < nEdges_owned)
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = nNodes_owned +  edge_subdomain;
          else
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostEdgeOffset + edge_subdomain - nEdges_owned;
          //unique subdomain id --> unique cross processor id
          //subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global + edge_global;
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset+nN]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[nN_global_subdomain*3+eI]+
                                                                                           mesh.subdomainp->nodeArray[nN_neig_global_subdomain*3+eI]);
        }
      local_offset += mesh.subdomainp->nNodes_element-2;
      //(node_i,node_{i+3}) --> unique subdomain edge, 0 = i
      for (int nN = 0; nN < mesh.subdomainp->nNodes_element-3; nN++)
        {
          const int nN_neig = (nN+3)%mesh.subdomainp->nNodes_element;
          const int nN_global_subdomain     = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          const int nN_neig_global_subdomain= mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN_neig];
          //unique edge subdomain id and global id
          int edge_subdomain = -1, edge_global = -1;
          //see if edge is (nN,nN_neig) or vice versa
          int nodes[2];
          nodes[0] = nN_global_subdomain;
          nodes[1] = nN_neig_global_subdomain;
          NodeTuple<2> et(nodes);
          edge_subdomain = nodesEdgeMap_subdomain[et];
          edge_global    = edgeNumbering_subdomain2global[edge_subdomain];
          assert(edge_subdomain >= 0 && edge_global >= 0);
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (edge_subdomain < nEdges_owned)
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = nNodes_owned +  edge_subdomain;
          else
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostEdgeOffset + edge_subdomain - nEdges_owned;
          //unique subdomain id --> unique cross processor id
          //subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global +  edge_global;
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset+nN]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[nN_global_subdomain*3+eI]+
                                                                                           mesh.subdomainp->nodeArray[nN_neig_global_subdomain*3+eI]);
        }
      
    }//eN

  ISRestoreIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  ISDestroy(&quadraticNumberingIS_new2old);
  ISDestroy(&quadraticNumberingIS_global_new2old);

  return 0;
}

int buildQuadraticCubeSubdomain2GlobalMappings_3d(Mesh& mesh, 
                                                  const int *edgeOffsets_subdomain_owned,
                                                  const int *nodeOffsets_subdomain_owned,
                                                  const int *edgeNumbering_subdomain2global,
                                                  const int *nodeNumbering_subdomain2global,
                                                  int& nDOF_all_processes,//total number of dofs in whole domain
                                                  int& nDOF_subdomain,//total number of dofs in sub-domain
                                                  int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                  int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                  int *subdomain_l2g, //local to global dof mapping on subdomain
                                                  int *subdomain2global,//subdomain dof to global (parallel) numbering
                                                  double * lagrangeNodesArray)//location of nodes corresponding to dofs
{
  using namespace std;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return -1;
  }

  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  //In 2d the quadratic dofs can be associated with nodes and element Boundaries
  //assuming have ownership info and consistent local/global mappings for nodes, elementBoundaries
  //build a mapping from local quadratic dofs to global quadratic dofs for petsc
  //assume a processor owns a dof if it owns that node or element
  assert(edgeOffsets_subdomain_owned);
  assert(nodeOffsets_subdomain_owned);
  assert(edgeNumbering_subdomain2global);
  assert(nodeNumbering_subdomain2global);
  assert(mesh.subdomainp);

  const int *elementBoundaryOffsets_subdomain_owned=mesh.elementBoundaryOffsets_subdomain_owned; 
  const int *elementOffsets_subdomain_owned=mesh.elementOffsets_subdomain_owned; 
  const int *elementBoundaryNumbering_subdomain2global=mesh.elementBoundaryNumbering_subdomain2global;
  const int *elementNumbering_subdomain2global=mesh.elementNumbering_subdomain2global;

  const int nNodes_owned      = nodeOffsets_subdomain_owned[rank+1]-nodeOffsets_subdomain_owned[rank];
  const int nEdges_owned      = edgeOffsets_subdomain_owned[rank+1]-edgeOffsets_subdomain_owned[rank];
  const int nBoundaries_owned = elementBoundaryOffsets_subdomain_owned[rank+1]-elementBoundaryOffsets_subdomain_owned[rank];
  const int nElements_owned   = elementOffsets_subdomain_owned[rank+1]-elementOffsets_subdomain_owned[rank];

  const int nDOFs_owned = nNodes_owned + nEdges_owned + nBoundaries_owned + nElements_owned;
  const int nDOFs_global= mesh.nNodes_global + mesh.nEdges_global + mesh.nElementBoundaries_global + mesh.nElements_global;
  //start with a logical global ordering of dofs as
  //[global nodes, global edges]
  //want to create global numbering
  //[nodes proc 0,edges proc 0,nodes proc 1, edges proc 1,...]
  valarray<int> quadraticNumbering_new2old(nDOFs_owned);
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global = nodeNumbering_subdomain2global[nN];
      quadraticNumbering_new2old[nN] = dof_global;
    }
  for (int i = 0; i < nEdges_owned; i++)
    {
      int dof_global = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      quadraticNumbering_new2old[nNodes_owned + i] = dof_global;
    }

  //-------------------
  for (int i = 0; i < nBoundaries_owned; i++)
    {
      int dof_global =  mesh.nNodes_global + mesh.nEdges_global + elementBoundaryNumbering_subdomain2global[i];
      quadraticNumbering_new2old[nNodes_owned + nEdges_owned + i] = dof_global;
    }
  for (int i = 0; i < nElements_owned; i++)
    {
      int dof_global = mesh.nNodes_global + mesh.nEdges_global + mesh.nElementBoundaries_global  + elementNumbering_subdomain2global[i];
      quadraticNumbering_new2old[nNodes_owned + nEdges_owned + nBoundaries_owned + i] = dof_global;
    }  
  //-------------------


  //build an index set for new numbering
  IS quadraticNumberingIS_new2old;
  ISCreateGeneral(PROTEUS_COMM_WORLD,nDOFs_owned,&quadraticNumbering_new2old[0],PETSC_COPY_VALUES,&quadraticNumberingIS_new2old);
  IS quadraticNumberingIS_global_new2old;
  ISAllGather(quadraticNumberingIS_new2old,&quadraticNumberingIS_global_new2old);
  //get old 2 new mapping for dofs
  const PetscInt *quadraticNumbering_global_new2old;
  valarray<int> quadraticNumbering_old2new_global(nDOFs_global);
  ISGetIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  for (int id = 0; id < nDOFs_global; id++)
    quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]] = id;
  //mwf debug
  //   if (rank == 0)
  //     {
  //       for (int id = 0; id < nDOFs_global; id++)
  // 	{
  // 	  std::cout<<" rank= "<<rank<<" build 2d c0p2 mappings new2old["<<id<<"]= "<<quadraticNumbering_global_new2old[id]
  // 		   <<" old2new["<<quadraticNumbering_global_new2old[id]<<"]= "<<quadraticNumbering_old2new_global[quadraticNumbering_global_new2old[id]]
  // 		   <<std::endl;
  // 	}
  //     }
  assert(offsets_subdomain_owned);
  assert(subdomain2global);
  assert(subdomain_l2g);
  assert(lagrangeNodesArray);
  for (int sdN=0; sdN < size+1; sdN++)
    offsets_subdomain_owned[sdN] = nodeOffsets_subdomain_owned[sdN]+edgeOffsets_subdomain_owned[sdN]+elementBoundaryOffsets_subdomain_owned[sdN]+elementOffsets_subdomain_owned[sdN];
  nDOF_all_processes = mesh.nNodes_global+mesh.nEdges_global+ mesh.nElementBoundaries_global+mesh.nElements_global;
  nDOF_subdomain = mesh.subdomainp->nNodes_global+mesh.subdomainp->nEdges_global+mesh.subdomainp->nElementBoundaries_global+mesh.subdomainp->nElements_global;
  max_dof_neighbors = 2*mesh.max_nNodeNeighbors_node;//8*27?

  //loop through owned and ghost dofs build subdomain mapping by 
  //going from old --> new
  int localOffset = 0;
  for (int nN = 0; nN < nNodes_owned; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN] = dof_global_new;
      //mwf debug
      //      if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
    
  localOffset += nNodes_owned;
  for (int i = 0; i < nEdges_owned; i++)
    {
      int dof_global_old = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	  std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i<<" nEdges_owned= "<<nEdges_owned
      // 		   <<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }


  localOffset += nEdges_owned;
  for (int i = 0; i < nBoundaries_owned; i++)
    {
      int dof_global_old = mesh.nNodes_global + mesh.nEdges_global +  elementBoundaryNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	  std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i<<" nEdges_owned= "<<nEdges_owned
      // 		   <<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }

  localOffset += nBoundaries_owned;
  for (int i = 0; i < nElements_owned; i++)
    {
      int dof_global_old = mesh.nNodes_global + mesh.nEdges_global +  mesh.nElementBoundaries_global + elementNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	  std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i<<" nEdges_owned= "<<nEdges_owned
      // 		   <<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }

  localOffset += nElements_owned;
  for (int nN = nNodes_owned; nN < mesh.subdomainp->nNodes_global; nN++)
    {
      int dof_global_old = nodeNumbering_subdomain2global[nN];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + nN - nNodes_owned] = dof_global_new;
      //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" nN= "<<nN<<" local dof= "<<localOffset+nN-nNodes_owned<<" nNodes_owned= "<<nNodes_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
    
  localOffset += mesh.subdomainp->nNodes_global - nNodes_owned;
  for (int i = nEdges_owned; i < mesh.subdomainp->nEdges_global; i++)
    {
      int dof_global_old = mesh.nNodes_global + edgeNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i - nEdges_owned] = dof_global_new;
      //       //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i-nEdges_owned<<" nEdges_owned= "<<nEdges_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
    
  localOffset += mesh.subdomainp->nEdges_global - nEdges_owned;
  for (int i = nBoundaries_owned; i < mesh.subdomainp->nElementBoundaries_global; i++)
    {
      int dof_global_old = mesh.nNodes_global + mesh.nEdges_global + elementBoundaryNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i - nBoundaries_owned] = dof_global_new;
      //       //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i-nEdges_owned<<" nEdges_owned= "<<nEdges_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }
    
  localOffset += mesh.subdomainp->nElementBoundaries_global - nBoundaries_owned;
  for (int i = nElements_owned; i < mesh.subdomainp->nElements_global; i++)
    {
      int dof_global_old = mesh.nNodes_global + mesh.nEdges_global + mesh.nElementBoundaries_global + elementNumbering_subdomain2global[i];
      int dof_global_new = quadraticNumbering_old2new_global[dof_global_old];
      subdomain2global[localOffset + i - nElements_owned] = dof_global_new;
      //       //mwf debug
      //       if (rank == 0)
      // 	std::cout<<" rank= "<<rank<<" i= "<<i<<" local dof= "<<localOffset+i-nEdges_owned<<" nEdges_owned= "<<nEdges_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
    }

  //setup local to global mapping on the subdomain for finite elements
  const int nDOF_element = 27;
  const int ghostNodeOffset = nNodes_owned + nEdges_owned + nBoundaries_owned + nElements_owned;
  const int ghostEdgeOffset = ghostNodeOffset + mesh.subdomainp->nNodes_global-nNodes_owned;
  const int ghostBoundaryOffset = ghostEdgeOffset + mesh.subdomainp->nEdges_global-nEdges_owned;
  const int ghostElementOffset = ghostBoundaryOffset + mesh.subdomainp->nElementBoundaries_global-nBoundaries_owned;

  int ledge[12][2] = {{0,1},{1,2},{2,3},{3,0},
                      {0,4},{1,5},{2,6},{3,7},
                      {4,5},{5,6},{6,7},{7,4}};

  int nEdges_element = 12;

  int lface[6][4] = {{0,1,2,3},
                     {0,1,5,4},
                     {1,2,6,5},
                     {2,3,7,6},
                     {3,0,4,7},
                     {4,5,6,7}};
   
  assert(mesh.subdomainp->nElementBoundaries_element == 6);
  //need mapping from nodes to edge to setup element based relationship
  map<NodeTuple<2>, int> nodesEdgeMap_subdomain;
  for (int i=0; i < mesh.subdomainp->nEdges_global; i++)
    {
      register int nodes[2];
      nodes[0] = mesh.subdomainp->edgeNodesArray[2*i+0];
      nodes[1] = mesh.subdomainp->edgeNodesArray[2*i+1];
      NodeTuple<2> nt(nodes);
      nodesEdgeMap_subdomain[nt] = i;
    }
   
  map<NodeTuple<4>, int> nodesBoundaryMap_subdomain;
  for (int i=0; i < mesh.subdomainp->nElementBoundaries_global; i++)
    {
       
      register int nodes[4];
      nodes[0] = mesh.subdomainp->elementBoundaryNodesArray[i*4+0];
      nodes[1] = mesh.subdomainp->elementBoundaryNodesArray[i*4+1];
      nodes[2] = mesh.subdomainp->elementBoundaryNodesArray[i*4+2];
      nodes[3] = mesh.subdomainp->elementBoundaryNodesArray[i*4+3];
      NodeTuple<4> ebt(nodes);
      assert(nodesBoundaryMap_subdomain.find(ebt) == nodesBoundaryMap_subdomain.end());
      nodesBoundaryMap_subdomain[ebt] = i;
    }
   
  for (int eN=0; eN < mesh.subdomainp->nElements_global; eN++)
    {
      int local_offset = 0;
      for (int nN=0; nN < mesh.subdomainp->nNodes_element; nN++)
        {
          //node_i --> unique subdomain node  0 <= i <= 3
          const int nN_global_subdomain = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + nN];
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (nN_global_subdomain < nNodes_owned)
            {
              subdomain_l2g[eN*nDOF_element + local_offset + nN] = nN_global_subdomain;
            }
          else
            {
              subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostNodeOffset + nN_global_subdomain - nNodes_owned;
            }
          //unique subdomain id --> unique cross processor id
          //subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = nodeNumbering_subdomain2global[nN_global_subdomain];
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" node= "<<nN<<" sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element + local_offset +nN]*3+eI] = mesh.subdomainp->nodeArray[nN_global_subdomain*3+eI];
        }
      local_offset += mesh.subdomainp->nNodes_element;
      for (int nN = 0; nN < nEdges_element; nN++)
        {
          register int nodes[2];
          nodes[0] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + ledge[nN][0]];
          nodes[1] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + ledge[nN][1]];
          NodeTuple<2> nt(nodes);
          int edge_subdomain = nodesEdgeMap_subdomain[nt];
          int edge_global    = edgeNumbering_subdomain2global[edge_subdomain];
          assert(edge_subdomain >= 0 && edge_global >= 0);
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (edge_subdomain < nEdges_owned)
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = nNodes_owned + edge_subdomain;
          else
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostEdgeOffset + edge_subdomain - nEdges_owned;
          //unique subdomain id --> unique cross processor id
          //set above could do here subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global + edge_global;
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }  
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset+nN]*3+eI] = 0.5*(mesh.subdomainp->nodeArray[nodes[0]*3+eI]+
                                                                                           mesh.subdomainp->nodeArray[nodes[1]*3+eI]);
        }
       
      local_offset += nEdges_element;
      for (int nN = 0; nN <  mesh.subdomainp->nElementBoundaries_element; nN++)
        {
           
          register int nodes[4];
          nodes[0] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+lface[nN][0]];
          nodes[1] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+lface[nN][1]];
          nodes[2] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+lface[nN][2]];      
          nodes[3] = mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element+lface[nN][3]];
          NodeTuple<4> ebt(nodes);
           
           
          int boundary_subdomain = nodesBoundaryMap_subdomain[ebt];
          int boundary_global    = elementBoundaryNumbering_subdomain2global[boundary_subdomain];
           
          assert(boundary_subdomain >= 0 && boundary_global >= 0);
           
          //assign unique subdomain id based on |owned nodes|owned edges|ghost nodes|ghost edges| 
          if (boundary_subdomain < nBoundaries_owned)
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = nNodes_owned + nEdges_owned + boundary_subdomain;
          else
            subdomain_l2g[eN*nDOF_element + local_offset + nN] = ghostBoundaryOffset + boundary_subdomain - nBoundaries_owned;
          //unique subdomain id --> unique cross processor id
          //set above could do here subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global + edge_global;
          //mwf debug
          // 	  if (rank == 0)
          // 	    {
          // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
          // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
          // 	    }  
          for (int eI=0; eI < 3; eI++)
            lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset+nN]*3+eI] = 0.25*(mesh.subdomainp->nodeArray[nodes[0]*3+eI]+
                                                                                            mesh.subdomainp->nodeArray[nodes[1]*3+eI]+
                                                                                            mesh.subdomainp->nodeArray[nodes[2]*3+eI]+
                                                                                            mesh.subdomainp->nodeArray[nodes[3]*3+eI]);
        }
       
      local_offset += mesh.subdomainp->nElementBoundaries_element;
       
      // Interior node!!!
       
      if (eN < nElements_owned)
        subdomain_l2g[eN*nDOF_element + local_offset] = nNodes_owned + nEdges_owned + nBoundaries_owned + eN;
      else
        subdomain_l2g[eN*nDOF_element + local_offset] = ghostElementOffset + eN - nElements_owned;
      //unique subdomain id --> unique cross processor id
      //set above could do here subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]] = mesh.nNodes_global + edge_global;
      //mwf debug
      // 	  if (rank == 0)
      // 	    {
      // 	      std::cout<<"build loc2glob c0p2 3d eN= "<<eN<<" edge("<<nN<<","<<nN_neig<<") sub_dof= "<< subdomain_l2g[eN*nDOF_element + local_offset + nN]
      // 		       <<" glob_dof= "<<subdomain2global[subdomain_l2g[eN*nDOF_element + local_offset + nN]]<<std::endl;
      // 	    }  
      for (int eI=0; eI < 3; eI++)
        lagrangeNodesArray[subdomain_l2g[eN*nDOF_element+local_offset]*3+eI] = 0.125*(mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 0]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 1]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 2]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 3]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 4]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 5]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 6]*3+eI]+
                                                                                      mesh.subdomainp->nodeArray[mesh.subdomainp->elementNodesArray[eN*mesh.subdomainp->nNodes_element + 7]*3+eI]);
       
       
    }//eN
   
  ISRestoreIndices(quadraticNumberingIS_global_new2old,&quadraticNumbering_global_new2old);
  ISDestroy(&quadraticNumberingIS_new2old);
  ISDestroy(&quadraticNumberingIS_global_new2old);
   
  return 0;
}





int buildDiscontinuousGalerkinSubdomain2GlobalMappings(Mesh& mesh, 
                                                       const int *elementOffsets_subdomain_owned,
                                                       const int *elementNumbering_subdomain2global,
                                                       int nDOF_element,
                                                       int& nDOF_all_processes,//total number of dofs in whole domain
                                                       int& nDOF_subdomain,//total number of dofs in sub-domain
                                                       int& max_dof_neighbors,//maximum number of neighbors for connectivity of dofs
                                                       int *offsets_subdomain_owned, //starting point of local dofs on each processor (nProcs+1)
                                                       int * subdomain_l2g, //local to global dof mapping on subdomain
                                                       int* subdomain2global)//subdomain dof to global (parallel) numbering

{
  using namespace std;
  int ierr,size,rank;
  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  if (!ensure_comm()) {
    return -1;
  }

  //DG dofs stored element wise
  //[...,eN_0,eN_1,..,eN_ndof_local,...]
  //assume a processor owns a dof if it owns that element
  assert(elementOffsets_subdomain_owned);
  assert(elementNumbering_subdomain2global);
  assert(mesh.subdomainp);

  const int nElements_owned = elementOffsets_subdomain_owned[rank+1]-elementOffsets_subdomain_owned[rank];
  
  assert(offsets_subdomain_owned);
  assert(subdomain2global);
  assert(subdomain_l2g);
  nDOF_all_processes = mesh.nElements_global*nDOF_element;
  nDOF_subdomain = mesh.subdomainp->nElements_global*nDOF_element;
  max_dof_neighbors = mesh.nElementBoundaries_element*nDOF_element;

  for (int sdN=0; sdN < size+1; sdN++)
    offsets_subdomain_owned[sdN] = elementOffsets_subdomain_owned[sdN]*nDOF_element;

  //loop through owned and ghost dofs build subdomain mapping 
  for (int eN = 0; eN < mesh.subdomainp->nElements_global; eN++)
    {
      for (int i = 0; i < nDOF_element; i++)
        {
          subdomain2global[eN*nDOF_element + i] = elementNumbering_subdomain2global[eN]*nDOF_element + i;
          //mwf debug
          //std::cout<<" rank= "<<rank<<" eN= "<<eN<<" local dof= "<<localOffset+eN<<" nElements_owned= "<<nElements_owned<<" dof_old= "<<dof_global_old<<" new= "<<dof_global_new<<std::endl;
          subdomain_l2g[eN*nDOF_element+i] = eN*nDOF_element + i;
        }
    }

  return 0;
}

static PyObject* flcbdfWrappersGlobalSum(PyObject* self, PyObject* args)
{
  using namespace std;
  double value,value_new;

  if (!ensure_comm()) {
    return NULL;
  }

  if (!PyArg_ParseTuple(args,
                        "d",
                        &value))
    return NULL;
  MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_SUM,PROTEUS_COMM_WORLD);
  return Py_BuildValue("d",value_new);
}

static PyObject* flcbdfWrappersGlobalMax(PyObject* self, PyObject* args)
{
  using namespace std;
  double value,value_new;

  if (!ensure_comm()) {
    return NULL;
  }

  if (!PyArg_ParseTuple(args,
                        "d",
                        &value))
    return NULL;

  MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_MAX,PROTEUS_COMM_WORLD);
  return Py_BuildValue("d",value_new);
}
static PyObject* flcbdfWrappersGlobalMin(PyObject* self, PyObject* args)
{
  using namespace std;
  double value,value_new;

  if (!ensure_comm()) {
    return NULL;
  }

  if (!PyArg_ParseTuple(args,
                        "d",
                        &value))
    return NULL;
  MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_MIN,PROTEUS_COMM_WORLD);
  return Py_BuildValue("d",value_new);
}


static PyObject* flcbdfWrappersPartitionElements(PyObject* self,
                                                 PyObject* args)
{
  using namespace std;
  int nLayersOfOverlap;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *nodeOffsets_subdomain_owned,
    *nodeNumbering_subdomain2global,
    *elementBoundaryOffsets_subdomain_owned,
    *elementBoundaryNumbering_subdomain2global,
    *edgeOffsets_subdomain_owned,
    *edgeNumbering_subdomain2global;
  if (!PyArg_ParseTuple(args,
                        "iOO",
                        &nLayersOfOverlap,
                        &cmesh,
                        &subdomain_cmesh))
    return NULL;
  MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
  int ierr,size,rank;

  if (!ensure_comm()) {
    return NULL;
  }


  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
  partitionElements(MESH(cmesh),nLayersOfOverlap);
 

  //   Vec u2;
  //   int n = MESH(cmesh).nodeOffsets_subdomain_owned[rank+1] - MESH(cmesh).nodeOffsets_subdomain_owned[rank],
  //     N = MESH(cmesh).nNodes_global,
  //     nghost = MESH(cmesh).subdomainp->nNodes_global - n;
  //   valarray<int> ghost(nghost);
  //   for (int ii=0;ii<nghost;ii++)
  //     ghost[ii] = MESH(cmesh).nodeNumbering_subdomain2global[n+ii];
  //   VecCreateGhost(PROTEUS_COMM_WORLD,
  //                  n,
  //                  N,
  //                  nghost,
  //                  &ghost[0],
  //                  &u2);
  //   VecDestroy(&u2);

  int dims[1];
  //build handles to python arrays
  dims[0] = size+1;
  elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).elementOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElements_global;
  elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                              dims,
                                                              PyArray_INT,
                                                              (char*)MESH(cmesh).elementNumbering_subdomain2global);
  dims[0] = size+1;
  nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).nodeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nNodes_global;
  nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).nodeNumbering_subdomain2global);
  dims[0] = size+1;
  elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                   dims,
                                                                   PyArray_INT,
                                                                   (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
  elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                      dims,
                                                                      PyArray_INT,
                                                                      (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
  dims[0] = size+1;
  edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).edgeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nEdges_global;
  edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).edgeNumbering_subdomain2global);
  return Py_BuildValue("OOOOOOOO",
                       elementOffsets_subdomain_owned,
                       elementNumbering_subdomain2global,
                       nodeOffsets_subdomain_owned,
                       nodeNumbering_subdomain2global,
                       elementBoundaryOffsets_subdomain_owned,
                       elementBoundaryNumbering_subdomain2global,
                       edgeOffsets_subdomain_owned,
                       edgeNumbering_subdomain2global);
}
static PyObject* flcbdfWrappersPartitionNodes(PyObject* self,
                                              PyObject* args)
{
  using namespace std;
  int nLayersOfOverlap;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *nodeOffsets_subdomain_owned,
    *nodeNumbering_subdomain2global,
    *elementBoundaryOffsets_subdomain_owned,
    *elementBoundaryNumbering_subdomain2global,
    *edgeOffsets_subdomain_owned,
    *edgeNumbering_subdomain2global;
  if (!PyArg_ParseTuple(args,
                        "iOO",
                        &nLayersOfOverlap,
                        &cmesh,
                        &subdomain_cmesh))
    return NULL;

  if (!ensure_comm()) {
    return NULL;
  }

  MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
  int ierr,size,rank;
  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
  partitionNodes(MESH(cmesh),nLayersOfOverlap);
 

  //   Vec u2;
  //   int n = MESH(cmesh).nodeOffsets_subdomain_owned[rank+1] - MESH(cmesh).nodeOffsets_subdomain_owned[rank],
  //     N = MESH(cmesh).nNodes_global,
  //     nghost = MESH(cmesh).subdomainp->nNodes_global - n;
  //   valarray<int> ghost(nghost);
  //   for (int ii=0;ii<nghost;ii++)
  //     ghost[ii] = MESH(cmesh).nodeNumbering_subdomain2global[n+ii];
  //   VecCreateGhost(PROTEUS_COMM_WORLD,
  //                  n,
  //                  N,
  //                  nghost,
  //                  &ghost[0],
  //                  &u2);
  //   VecDestroy(&u2);

  int dims[1];
  //build handles to python arrays
  dims[0] = size+1;
  elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).elementOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElements_global;
  elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                              dims,
                                                              PyArray_INT,
                                                              (char*)MESH(cmesh).elementNumbering_subdomain2global);
  dims[0] = size+1;
  nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).nodeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nNodes_global;
  nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).nodeNumbering_subdomain2global);
  dims[0] = size+1;
  elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                   dims,
                                                                   PyArray_INT,
                                                                   (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
  elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                      dims,
                                                                      PyArray_INT,
                                                                      (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
  dims[0] = size+1;
  edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).edgeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nEdges_global;
  edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).edgeNumbering_subdomain2global);

  return Py_BuildValue("OOOOOOOO",
                       elementOffsets_subdomain_owned,
                       elementNumbering_subdomain2global,
                       nodeOffsets_subdomain_owned,
                       nodeNumbering_subdomain2global,
                       elementBoundaryOffsets_subdomain_owned,
                       elementBoundaryNumbering_subdomain2global,
                       edgeOffsets_subdomain_owned,
                       edgeNumbering_subdomain2global);
}

static PyObject* flcbdfWrappersConvertPUMIPartitionToPython(PyObject* self,
							     PyObject* args)
{
  using namespace std;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *elementNumbering_global2original,
    *nodeOffsets_subdomain_owned,
    *nodeNumbering_subdomain2global,
    *nodeNumbering_global2original,
    *elementBoundaryOffsets_subdomain_owned,
    *elementBoundaryNumbering_subdomain2global,
    *elementBoundaryNumbering_global2original,
    *edgeOffsets_subdomain_owned,
    *edgeNumbering_subdomain2global,
    *edgeNumbering_global2original;
  if (!PyArg_ParseTuple(args,
                        "OO",
                        &cmesh,
                        &subdomain_cmesh))
    return NULL;
  MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
  int ierr,size,rank;
  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

  int dims[1];
  //build handles to python arrays
  dims[0] = size+1;
  elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).elementOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElements_global;
  elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                              dims,
                                                              PyArray_INT,
                                                              (char*)MESH(cmesh).elementNumbering_subdomain2global);
  dims[0] = size+1;
  nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).nodeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nNodes_global;
  nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).nodeNumbering_subdomain2global);
  dims[0] = size+1;
  elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
								   dims,
								   PyArray_INT,
								   (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
  elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
								      dims,
								      PyArray_INT,
								      (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
  dims[0] = size+1;
  edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
							dims,
							PyArray_INT,
							(char*)MESH(cmesh).edgeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nEdges_global;
  edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
							   dims,
							   PyArray_INT,
							   (char*)MESH(cmesh).edgeNumbering_subdomain2global);
  return Py_BuildValue("OOOOOOOO",
                       elementOffsets_subdomain_owned,
                       elementNumbering_subdomain2global,
                       nodeOffsets_subdomain_owned,
                       nodeNumbering_subdomain2global,
                       elementBoundaryOffsets_subdomain_owned,
                       elementBoundaryNumbering_subdomain2global,
                       edgeOffsets_subdomain_owned,
                       edgeNumbering_subdomain2global);

}                      
static PyObject* flcbdfWrappersPartitionNodesFromTetgenFiles(PyObject* self,
                                                             PyObject* args)
{
  using namespace std;
  int nLayersOfOverlap, indexBase;
  char* filebase;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *nodeOffsets_subdomain_owned,
    *nodeNumbering_subdomain2global,
    *elementBoundaryOffsets_subdomain_owned,
    *elementBoundaryNumbering_subdomain2global,
    *edgeOffsets_subdomain_owned,
    *edgeNumbering_subdomain2global;
  if (!PyArg_ParseTuple(args,
                        "siiOO",
                        &filebase,
                        &indexBase,
                        &nLayersOfOverlap,
                        &cmesh,
                        &subdomain_cmesh))
    return NULL;
  MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;

  if (!ensure_comm()) {
    return NULL;
  }

  int ierr,size,rank;
  ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
  partitionNodesFromTetgenFiles(filebase,indexBase,MESH(cmesh),nLayersOfOverlap);

  //   Vec u2;
  //   int n = MESH(cmesh).nodeOffsets_subdomain_owned[rank+1] - MESH(cmesh).nodeOffsets_subdomain_owned[rank],
  //     N = MESH(cmesh).nNodes_global,
  //     nghost = MESH(cmesh).subdomainp->nNodes_global - n;
  //   valarray<int> ghost(nghost);
  //   for (int ii=0;ii<nghost;ii++)
  //     ghost[ii] = MESH(cmesh).nodeNumbering_subdomain2global[n+ii];
  //   VecCreateGhost(PROTEUS_COMM_WORLD,
  //                  n,
  //                  N,
  //                  nghost,
  //                  &ghost[0],
  //                  &u2);
  //   VecDestroy(&u2);

  int dims[1];
  //build handles to python arrays
  dims[0] = size+1;
  elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).elementOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElements_global;
  elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                              dims,
                                                              PyArray_INT,
                                                              (char*)MESH(cmesh).elementNumbering_subdomain2global);
  dims[0] = size+1;
  nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).nodeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nNodes_global;
  nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).nodeNumbering_subdomain2global);
  dims[0] = size+1;
  elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                   dims,
                                                                   PyArray_INT,
                                                                   (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
  elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
								      dims,
								      PyArray_INT,
								      (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
  dims[0] = size+1;
  edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                        dims,
                                                        PyArray_INT,
                                                        (char*)MESH(cmesh).edgeOffsets_subdomain_owned);
  
  dims[0] = MESH(cmesh).subdomainp->nEdges_global;
  edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
							   dims,
							   PyArray_INT,
							   (char*)MESH(cmesh).edgeNumbering_subdomain2global);
  return Py_BuildValue("OOOOOOOO",
                       elementOffsets_subdomain_owned,
                       elementNumbering_subdomain2global,
                       nodeOffsets_subdomain_owned,
                       nodeNumbering_subdomain2global,
		       elementBoundaryOffsets_subdomain_owned,
                       elementBoundaryNumbering_subdomain2global,
		       edgeOffsets_subdomain_owned,
                       edgeNumbering_subdomain2global);
}


static PyObject* flcbdfWrappersBuildQuadraticLocal2GlobalMappings(PyObject* self,
                                                                  PyObject* args)
{
  using namespace std;
  int nSpace;
  int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *nodeOffsets_subdomain_owned,
    *elementBoundaryOffsets_subdomain_owned,
    *edgeOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *nodeNumbering_subdomain2global,
    *elementBoundaryNumbering_subdomain2global,
    *edgeNumbering_subdomain2global,
    *quadratic_dof_offsets_subdomain_owned,
    *quadraticNumbering_subdomain2global,
    *quadratic_subdomain_l2g,
    *quadratic_lagrangeNodes;
  if (!PyArg_ParseTuple(args,
                        "iOOOOOOOOOOOOOO",
                        &nSpace,
                        &cmesh,
                        &subdomain_cmesh,
                        &elementOffsets_subdomain_owned,
                        &nodeOffsets_subdomain_owned,
                        &elementBoundaryOffsets_subdomain_owned,
                        &edgeOffsets_subdomain_owned,
                        &elementNumbering_subdomain2global,
                        &nodeNumbering_subdomain2global,
                        &elementBoundaryNumbering_subdomain2global,
                        &edgeNumbering_subdomain2global,
                        &quadratic_dof_offsets_subdomain_owned,
                        &quadratic_subdomain_l2g,
                        &quadraticNumbering_subdomain2global,
                        &quadratic_lagrangeNodes))
    return NULL;
  //MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  if (nSpace == 1)
    {
      buildQuadraticSubdomain2GlobalMappings_1d(MESH(cmesh),
                                                IDATA(elementOffsets_subdomain_owned),
                                                IDATA(nodeOffsets_subdomain_owned),
                                                IDATA(elementNumbering_subdomain2global),
                                                IDATA(nodeNumbering_subdomain2global),
                                                nDOF_all_processes,
                                                nDOF_subdomain,
                                                max_dof_neighbors,
                                                IDATA(quadratic_dof_offsets_subdomain_owned),
                                                IDATA(quadratic_subdomain_l2g),
                                                IDATA(quadraticNumbering_subdomain2global),
                                                DDATA(quadratic_lagrangeNodes));
    }
  else if (nSpace == 2)
    {
      buildQuadraticSubdomain2GlobalMappings_2d(MESH(cmesh),
                                                IDATA(elementBoundaryOffsets_subdomain_owned),
                                                IDATA(nodeOffsets_subdomain_owned),
                                                IDATA(elementBoundaryNumbering_subdomain2global),
                                                IDATA(nodeNumbering_subdomain2global),
                                                nDOF_all_processes,
                                                nDOF_subdomain,
                                                max_dof_neighbors,
                                                IDATA(quadratic_dof_offsets_subdomain_owned),
                                                IDATA(quadratic_subdomain_l2g),
                                                IDATA(quadraticNumbering_subdomain2global),
                                                DDATA(quadratic_lagrangeNodes));

    }
  else
    {
      buildQuadraticSubdomain2GlobalMappings_3d(MESH(cmesh),
                                                IDATA(edgeOffsets_subdomain_owned),
                                                IDATA(nodeOffsets_subdomain_owned),
                                                IDATA(edgeNumbering_subdomain2global),
                                                IDATA(nodeNumbering_subdomain2global),
                                                nDOF_all_processes,
                                                nDOF_subdomain,
                                                max_dof_neighbors,
                                                IDATA(quadratic_dof_offsets_subdomain_owned),
                                                IDATA(quadratic_subdomain_l2g),
                                                IDATA(quadraticNumbering_subdomain2global),
                                                DDATA(quadratic_lagrangeNodes));
    }
  

  return Py_BuildValue("iii",
                       nDOF_all_processes,
                       nDOF_subdomain,
                       max_dof_neighbors);

}


static PyObject* flcbdfWrappersBuildQuadraticCubeLocal2GlobalMappings(PyObject* self,
                                                                      PyObject* args)
{
  using namespace std;
  int nSpace;
  int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *nodeOffsets_subdomain_owned,
    *elementBoundaryOffsets_subdomain_owned,
    *edgeOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *nodeNumbering_subdomain2global,
    *elementBoundaryNumbering_subdomain2global,
    *edgeNumbering_subdomain2global,
    *quadratic_dof_offsets_subdomain_owned,
    *quadraticNumbering_subdomain2global,
    *quadratic_subdomain_l2g,
    *quadratic_lagrangeNodes;
  if (!PyArg_ParseTuple(args,
                        "iOOOOOOOOOOOOOO",
                        &nSpace,
                        &cmesh,
                        &subdomain_cmesh,
                        &elementOffsets_subdomain_owned,
                        &nodeOffsets_subdomain_owned,
                        &elementBoundaryOffsets_subdomain_owned,
                        &edgeOffsets_subdomain_owned,
                        &elementNumbering_subdomain2global,
                        &nodeNumbering_subdomain2global,
                        &elementBoundaryNumbering_subdomain2global,
                        &edgeNumbering_subdomain2global,
                        &quadratic_dof_offsets_subdomain_owned,
                        &quadratic_subdomain_l2g,
                        &quadraticNumbering_subdomain2global,
                        &quadratic_lagrangeNodes))
    return NULL;
  //MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
  if (nSpace == 1)
    {
      /*buildQuadraticCubeSubdomain2GlobalMappings_1d(MESH(cmesh),
        IDATA(elementOffsets_subdomain_owned),
        IDATA(nodeOffsets_subdomain_owned),
        IDATA(elementNumbering_subdomain2global),
        IDATA(nodeNumbering_subdomain2global),
        nDOF_all_processes,
        nDOF_subdomain,
        max_dof_neighbors,
        IDATA(quadratic_dof_offsets_subdomain_owned),
        IDATA(quadratic_subdomain_l2g),
        IDATA(quadraticNumbering_subdomain2global),
        DDATA(quadratic_lagrangeNodes));*/
      std::cout<<"buildQuadraticCubeSubdomain2GlobalMappings_1d not implemented!!"<<std::endl;
    }
  else if (nSpace == 2)
    {
      /* buildQuadraticCubeSubdomain2GlobalMappings_2d(MESH(cmesh),
         IDATA(elementBoundaryOffsets_subdomain_owned),
         IDATA(nodeOffsets_subdomain_owned),
         IDATA(elementBoundaryNumbering_subdomain2global),
         IDATA(nodeNumbering_subdomain2global),
         nDOF_all_processes,
         nDOF_subdomain,
         max_dof_neighbors,
         IDATA(quadratic_dof_offsets_subdomain_owned),
         IDATA(quadratic_subdomain_l2g),
         IDATA(quadraticNumbering_subdomain2global),
         DDATA(quadratic_lagrangeNodes));*/
      std::cout<<"buildQuadraticCubeSubdomain2GlobalMappings_2d not implemented!!"<<std::endl;
    }
  else
    {
      buildQuadraticCubeSubdomain2GlobalMappings_3d(MESH(cmesh),
                                                    IDATA(edgeOffsets_subdomain_owned),
                                                    IDATA(nodeOffsets_subdomain_owned),
                                                    IDATA(edgeNumbering_subdomain2global),
                                                    IDATA(nodeNumbering_subdomain2global),
                                                    nDOF_all_processes,
                                                    nDOF_subdomain,
                                                    max_dof_neighbors,
                                                    IDATA(quadratic_dof_offsets_subdomain_owned),
                                                    IDATA(quadratic_subdomain_l2g),
                                                    IDATA(quadraticNumbering_subdomain2global),
                                                    DDATA(quadratic_lagrangeNodes));
    }
  

  return Py_BuildValue("iii",
                       nDOF_all_processes,
                       nDOF_subdomain,
                       max_dof_neighbors);

}




static PyObject* flcbdfWrappersBuildDiscontinuousGalerkinLocal2GlobalMappings(PyObject* self,
                                                                              PyObject* args)
{
  using namespace std;
  int nDOF_element;
  int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
  PyObject *cmesh,*subdomain_cmesh,
    *elementOffsets_subdomain_owned,
    *elementNumbering_subdomain2global,
    *dg_dof_offsets_subdomain_owned,
    *dgNumbering_subdomain2global,
    *dg_subdomain_l2g;
  if (!PyArg_ParseTuple(args,
                        "iOOOOOOO",
                        &nDOF_element,
                        &cmesh,
                        &subdomain_cmesh,
                        &elementOffsets_subdomain_owned,
                        &elementNumbering_subdomain2global,
                        &dg_dof_offsets_subdomain_owned,
                        &dg_subdomain_l2g,
                        &dgNumbering_subdomain2global))

    return NULL;
  buildDiscontinuousGalerkinSubdomain2GlobalMappings(MESH(cmesh),
                                                     IDATA(elementOffsets_subdomain_owned),
                                                     IDATA(elementNumbering_subdomain2global),
                                                     nDOF_element,
                                                     nDOF_all_processes,
                                                     nDOF_subdomain,
                                                     max_dof_neighbors,
                                                     IDATA(dg_dof_offsets_subdomain_owned),
                                                     IDATA(dg_subdomain_l2g),
                                                     IDATA(dgNumbering_subdomain2global));
 

  return Py_BuildValue("iii",
                       nDOF_all_processes,
                       nDOF_subdomain,
                       max_dof_neighbors);

}

static PyMethodDef DaetkPetscSys_methods[] = {
  {"barrier",
   (PyCFunction)DaetkPetscSys_barrier,
   METH_VARARGS,
   "set an MPI barrier"},
  {"isMaster",
   (PyCFunction)DaetkPetscSys_isMaster,
   METH_VARARGS,
   "return whether this is the master process"},
  {"isInitialized",
   (PyCFunction)DaetkPetscSys_isInitialized,
   METH_VARARGS,
   "return whether the MPI communicator is initialized"},
  {"beginSequential",
   (PyCFunction)DaetkPetscSys_beginSequential,
   METH_VARARGS,
   "begin a sequential section of code"},
  {"endSequential",
   (PyCFunction)DaetkPetscSys_endSequential,
   METH_VARARGS,
   "end a sequential section of code"},
  {"catchError",
   (PyCFunction)DaetkPetscSys_catchError,
   METH_VARARGS,
   "catch an error on any processes"},
  {"rank",
   (PyCFunction)DaetkPetscSys_rank,
   METH_VARARGS,
   "get the process rank"},
  {"size",
   (PyCFunction)DaetkPetscSys_size,
   METH_VARARGS,
   "get the MPI communicator size"},
  {NULL,NULL}
};

static PyObject*
DaetkPetscSys_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  DaetkPetscSys *self;
  self = (DaetkPetscSys *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
DaetkPetscSys_init(DaetkPetscSys *self, PyObject *args, PyObject *kwds)
{
  int argc;
  char **argv;
  PyObject *sys_argv;
  /*mwf add for petsc database?*/
  int isInitialized;
  char *petscDatabaseFilename(0);
  if(!PyArg_ParseTuple(args,
                       "iO|s",
                       &isInitialized,&sys_argv,&petscDatabaseFilename))
    return -1;
  argc = PyList_Size(sys_argv);
  argv = new char* [argc+1];//don't  know why needs one past end cek
  for (int i=0;i<argc;i++)
    {
      argv[i] = PyString_AsString(PyList_GetItem(sys_argv,i));
    }
  argv[argc] = new char[1];
  argv[argc][0] = '\0';
  //cek need to think more about how to handle petsc4py and daetk. 
  //if (isInitialized)
  //Daetk::Petsc::Sys::initialized=true;
  if (petscDatabaseFilename)
    self->petscSys = new Daetk::Petsc::Sys(argc,argv,(char*)("Initializing petsc for Proteus, with options database\n"),
                                           petscDatabaseFilename);
  else
    self->petscSys = new Daetk::Petsc::Sys(argc,argv,(char*)("Initializing petsc for Proteus\n"));
  PROTEUS_COMM_WORLD = Daetk::Petsc::cc::PETSC_COMM_WORLD;
  delete [] argv;
  return 0;
}

static  void
DaetkPetscSys_dealloc(DaetkPetscSys* self)
{
  delete self->petscSys;
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject DaetkPetscSysType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "flcbdfWrappers.DaetkPetscSys",             /*tp_name*/
  sizeof(DaetkPetscSys), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)DaetkPetscSys_dealloc,                         /*tp_dealloc*/
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
  "DaetkPetscSys objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  DaetkPetscSys_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)DaetkPetscSys_init,      /* tp_init */
  0,                         /* tp_alloc */
  DaetkPetscSys_new,                 /* tp_new */
};

static PyObject* 
DaetkPetscSys_get(PyObject* self, PyObject* args)
{
  DaetkPetscSys* daetkPetscSys;
  daetkPetscSys = PyObject_NEW(DaetkPetscSys, &DaetkPetscSysType);
  daetkPetscSys->petscSys = new Daetk::Petsc::Sys();
  return (PyObject*)daetkPetscSys;
}


static PyMethodDef flcbdfWrappersMethods[] = {
  { "globalSum",
    flcbdfWrappersGlobalSum,
    METH_VARARGS, 
    "sum the value over all subdomains(processes)"},
  { "globalMax",
    flcbdfWrappersGlobalMax,
    METH_VARARGS, 
    "take the max of the value over all subdomains(processes)"},
  { "globalMin",
    flcbdfWrappersGlobalMin,
    METH_VARARGS, 
    "take the max of the value over all subdomains(processes)"},
  { "partitionElements",
    flcbdfWrappersPartitionElements,
    METH_VARARGS, 
    "partition the mesh using an element-based partitioning"},
  { "partitionNodes",
    flcbdfWrappersPartitionNodes,
    METH_VARARGS, 
    "partition the mesh using a node-based partitioning"},
  {"convertPUMIPartitionToPython", 
    flcbdfWrappersConvertPUMIPartitionToPython,
    METH_VARARGS,
    "Convert C structures to python for PUMI partitioned mesh"},
  { "partitionNodesFromTetgenFiles",
    flcbdfWrappersPartitionNodesFromTetgenFiles,
    METH_VARARGS, 
    "partition the mesh using a node-based partitioning"},
  { "buildQuadraticLocal2GlobalMappings",
    flcbdfWrappersBuildQuadraticLocal2GlobalMappings,
    METH_VARARGS, 
    "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping"},
  { "buildQuadraticCubeLocal2GlobalMappings",
    flcbdfWrappersBuildQuadraticCubeLocal2GlobalMappings,
    METH_VARARGS, 
    "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping for cubes"},
  { "buildDiscontinuousGalerkinLocal2GlobalMappings",
    flcbdfWrappersBuildDiscontinuousGalerkinLocal2GlobalMappings,
    METH_VARARGS, 
    "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping"},
  { "getDaetkPetscSys",
    DaetkPetscSys_get,
    METH_VARARGS, 
    "get the communicator  object"},
  { NULL,NULL,0,NULL}
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initflcbdfWrappers(void) 
{
  PyObject *m,*d,*c_api_object;
  static void* PyFLCBDFWrappers_API[1];
  if (PyType_Ready(&FLCBDF_integratorType) < 0)
    return;
  if (PyType_Ready(&DaetkPetscSysType) < 0)
    return;
  if (PyType_Ready(&ParVecType) < 0)
    return;
  if (PyType_Ready(&ParMatType) < 0)
    return;
  if (PyType_Ready(&CKSPType) < 0)
    return;
  m = Py_InitModule3("flcbdfWrappers", 
                     flcbdfWrappersMethods,
                     "flcbdf wrappers module");
  d = PyModule_GetDict(m);
  import_array();
  Py_INCREF(&FLCBDF_integratorType);
  PyModule_AddObject(m, "FLCBDF_integrator", (PyObject *)&FLCBDF_integratorType);
  Py_INCREF(&DaetkPetscSysType);
  PyModule_AddObject(m, "DaetkPetscSys", (PyObject *)&DaetkPetscSysType);
  Py_INCREF(&CKSPType);
  PyModule_AddObject(m, "KSP", (PyObject *)&CKSPType);
  Py_INCREF(&ParVecType);
  PyModule_AddObject(m, "ParVec", (PyObject *)&ParVecType);
  Py_INCREF(&ParMatType);
  PyModule_AddObject(m, "ParMat", (PyObject *)&ParMatType);

  // ensure PETSc, then DAETK, are initialized
  // PETSc first, via the proteus.Comm module
  //this wasn't working anyway... 
  //PyRun_SimpleString("proteus.Comm.init()");

  // DAETK, PETSc is initialized so only initialize subsystem
  Daetk::Petsc::Sys::initialized=true;
  int ignore1;
  char** ignore2;
  DAETK_SYS = new Daetk::Petsc::Sys(ignore1, ignore2, PETSC_NULL, PETSC_NULL);

  // Set up default Proteus communicator
  PROTEUS_COMM_WORLD = PETSC_COMM_WORLD;

  PyFLCBDFWrappers_API[0] = (void*)(&PROTEUS_COMM_WORLD);
  c_api_object = PyCObject_FromVoidPtr((void*)PyFLCBDFWrappers_API,NULL);
  PyModule_AddObject(m,"_C_API",c_api_object);
}
}
/** @} */
