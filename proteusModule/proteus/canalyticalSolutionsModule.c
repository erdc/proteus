#include "Python.h"
#include "numpy/arrayobject.h"
#include "analyticalSolutions.h"
/** 
    \file analyticalSolutions.h
    \defgroup canalyticalSolutions canalyticalSolutions
    \brief Python interface to analyticalSolutions library 
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

/* wrapper for C library function, only needs name changed*/
static PyObject* canalyticalSolutions_PoiseuillePipeFlow_P(PyObject* self,
                                                        PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  PoiseuillePipeFlow_P(IDATA(iwork),
                     DDATA(rwork),
                     nPoints,
                     t,
                     DDATA(x),
                     DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_PoiseuillePipeFlow(PyObject* self,
                                                        PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  PoiseuillePipeFlow(IDATA(iwork),
                     DDATA(rwork),
                     nPoints,
                     t,
                     DDATA(x),
                     DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_PlanePoiseuilleFlow_u(PyObject* self,
                                                            PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  PlanePoiseuilleFlow_u(IDATA(iwork),
                        DDATA(rwork),
                        nPoints,
                        t,
                        DDATA(x),
                        DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_PlaneCouetteFlow_u(PyObject* self,
                                                         PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  PlaneCouetteFlow_u(IDATA(iwork),
                     DDATA(rwork),
                     nPoints,
                     t,
                     DDATA(x),
                     DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_STflowSphere_P(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  STflowSphere_P(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_STflowSphere_Vz(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  STflowSphere_Vz(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* canalyticalSolutions_STflowSphere_Vy(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  STflowSphere_Vy(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_STflowSphere_Vx(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  STflowSphere_Vx(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_diffusionSin1D(PyObject* self, 
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  diffusionSin1D(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* canalyticalSolutions_diffusionSin2D(PyObject* self, 
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
    nPoints *= SHAPE(x)[i];
  diffusionSin2D(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* canalyticalSolutions_diffusionSin3D(PyObject* self, 
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
    nPoints *= SHAPE(x)[i];
  diffusionSin3D(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* canalyticalSolutions_diffusionSin1D_r(PyObject* self,
                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  diffusionSin1D_r(IDATA(iwork),
                   DDATA(rwork),
                   nPoints,
                   t,
                   DDATA(x),
                   DDATA(u),
                   DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_diffusionSin2D_r(PyObject* self,
                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
    nPoints *= SHAPE(x)[i];
  diffusionSin2D_r(IDATA(iwork),
                   DDATA(rwork),
                   nPoints,
                   t,
                   DDATA(x),
                   DDATA(u),
                   DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* canalyticalSolutions_diffusionSin3D_r(PyObject* self,
                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
    nPoints *= SHAPE(x)[i];
  diffusionSin3D_r(IDATA(iwork),
                   DDATA(rwork),
                   nPoints,
                   t,
                   DDATA(x),
                   DDATA(u),
                   DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* canalyticalSolutions_NonlinearDAE(PyObject* self,
                                                  PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearDAE(IDATA(iwork),
               DDATA(rwork),
               nPoints,
               t,
               DDATA(x),
               DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* canalyticalSolutions_NonlinearDAE_f(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearDAE_f(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_SteadyState(PyObject* self,
                                                            PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_SteadyState(IDATA(iwork),
                       DDATA(rwork),
                       nPoints,
                       t,
                       DDATA(x),
                       DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_NonlinearAD_SteadyState(PyObject* self,
                                                               PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearAD_SteadyState(IDATA(iwork),
                            DDATA(rwork),
                            nPoints,
                            t,
                            DDATA(x),
                            DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine(PyObject* self,
                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_du(PyObject* self,
                                                       PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_du(IDATA(iwork),
                    DDATA(rwork),
                    nPoints,
                    t,
                    DDATA(x),
                    DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_r(PyObject* self,
                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_r(IDATA(iwork),
                   DDATA(rwork),
                   nPoints,
                   t,
                   DDATA(x),
                   DDATA(u),
                   DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_dr(PyObject* self,
                                                       PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_dr(IDATA(iwork),
                    DDATA(rwork),
                    nPoints,
                    t,
                    DDATA(x),
                    DDATA(u),
                    DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_advectiveVelocity(PyObject* self,
                                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_advectiveVelocity(IDATA(iwork),
                                   DDATA(rwork),
                                   nPoints,
                                   t,
                                   DDATA(x),
                                   DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_diffusiveVelocity(PyObject* self,
                                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_diffusiveVelocity(IDATA(iwork),
                                   DDATA(rwork),
                                   nPoints,
                                   t,
                                   DDATA(x),
                                   DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Sine_totalVelocity(PyObject* self,
                                                                  PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Sine_totalVelocity(IDATA(iwork),
                               DDATA(rwork),
                               nPoints,
                               t,
                               DDATA(x),
                               DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp1D(PyObject* self,
                                                           PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp1D(IDATA(iwork),
                        DDATA(rwork),
                        nPoints,
                        t,
                        DDATA(x),
                        DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp2D(PyObject* self,
                                                           PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp2D(IDATA(iwork),
                        DDATA(rwork),
                        nPoints,
                        t,
                        DDATA(x),
                        DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp3D(PyObject* self,
                                                           PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp3D(IDATA(iwork),
                 DDATA(rwork),
                 nPoints,
                 t,
                 DDATA(x),
                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp1D_r(PyObject* self,
                                                             PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp1D_r(IDATA(iwork),
                          DDATA(rwork),
                          nPoints,
                          t,
                          DDATA(x),
                          DDATA(u),
                          DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp2D_r(PyObject* self,
                                                             PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp2D_r(IDATA(iwork),
                          DDATA(rwork),
                          nPoints,
                          t,
                          DDATA(x),
                          DDATA(u),
                          DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp3D_r(PyObject* self,
                                                             PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp3D_r(IDATA(iwork),
                          DDATA(rwork),
                          nPoints,
                          t,
                          DDATA(x),
                          DDATA(u),
                          DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_poissonsEquationExp3D_dr(PyObject* self,
                                                              PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  poissonsEquationExp3D_dr(IDATA(iwork),
                           DDATA(rwork),
                           nPoints,
                           t,
                           DDATA(x),
                           DDATA(u),
                           DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_DiracIC(PyObject* self,
                                                      PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_DiracIC(IDATA(iwork),
                   DDATA(rwork),
                   nPoints,
                   t,
                   DDATA(x),
                   DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_DiracIC_du(PyObject* self,
                                                         PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_DiracIC_du(IDATA(iwork),
                      DDATA(rwork),
                      nPoints,
                      t,
                      DDATA(x),
                      DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_DiracIC_advectiveVelocity(PyObject* self,
                                                                        PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_DiracIC_advectiveVelocity(IDATA(iwork),
                                     DDATA(rwork),
                                     nPoints,
                                     t,
                                     DDATA(x),
                                     DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_DiracIC_diffusiveVelocity(PyObject* self,
                                                                        PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_DiracIC_diffusiveVelocity(IDATA(iwork),
                                     DDATA(rwork),
                                     nPoints,
                                     t,
                                     DDATA(x),
                                     DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearAD_DiracIC_totalVelocity(PyObject* self,
                                                                    PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearAD_DiracIC_totalVelocity(IDATA(iwork),
                                 DDATA(rwork),
                                 nPoints,
                                 t,
                                 DDATA(x),
                                 DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Decay_DiracIC(PyObject* self,
                                                             PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Decay_DiracIC(IDATA(iwork),
                          DDATA(rwork),
                          nPoints,
                          t,
                          DDATA(x),
                          DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Decay_DiracIC_r(PyObject* self,
                                                               PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Decay_DiracIC_r(IDATA(iwork),
                            DDATA(rwork),
                            nPoints,
                            t,
                            DDATA(x),
                            DDATA(u),
                            DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_LinearADR_Decay_DiracIC_dr(PyObject* self,
                                                                PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  LinearADR_Decay_DiracIC_dr(IDATA(iwork),
                             DDATA(rwork),
                             nPoints,
                             t,
                             DDATA(x),
                             DDATA(u),
                             DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* canalyticalSolutions_NonlinearADR_Decay_DiracIC(PyObject* self,
                                                                PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearADR_Decay_DiracIC(IDATA(iwork),
                             DDATA(rwork),
                             nPoints,
                             t,
                             DDATA(x),
                             DDATA(u));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* canalyticalSolutions_NonlinearADR_Decay_DiracIC_dr(PyObject* self,
                                                                   PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearADR_Decay_DiracIC_dr(IDATA(iwork),
                                DDATA(rwork),
                                nPoints,
                                t,
                                DDATA(x),
                                DDATA(u),
                                DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* canalyticalSolutions_NonlinearADR_Decay_DiracIC_r(PyObject* self,
                                                                  PyObject* args)
{
  /* declare arguments */
  int i,nPoints=1;
  double t;
  PyObject *iwork,*rwork,*x,*u, *r;
  /* extract arguments using PyArg_ParseTuple from python.h */
  /* O = python object, d = double, i = int, c = char */
  if(!PyArg_ParseTuple(args,"OOdOO",
                       &iwork,
                       &rwork,
                       &t,
                       &x,
                       &u,
                       &r))
    return NULL;
  /* now call actual C function that does the work*/
  /* use DDATA macro to convert python array to double* */
  /* use IDATA macro to convert python int array to int* */
  /* use SHAPE to get shape array from python array */
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  NonlinearADR_Decay_DiracIC_r(IDATA(iwork),
                               DDATA(rwork),
                               nPoints,
                               t,
                               DDATA(x),
                               DDATA(u),
                               DDATA(r));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef canalyticalSolutionsMethods[] = {
  { "PoiseuillePipeFlow_P",
    canalyticalSolutions_PoiseuillePipeFlow_P,
    METH_VARARGS,
    "need to add doc"},
  { "PoiseuillePipeFlow",
    canalyticalSolutions_PoiseuillePipeFlow,
    METH_VARARGS,
    "need to add doc"},
  { "PlanePoiseuilleFlow_u",
    canalyticalSolutions_PlanePoiseuilleFlow_u,
    METH_VARARGS,
    "need to add doc"},
  { "PlaneCouetteFlow_u",
    canalyticalSolutions_PlaneCouetteFlow_u,
    METH_VARARGS,
    "need to add doc"},
  { "STflowSphere_P",
    canalyticalSolutions_STflowSphere_P,
    METH_VARARGS,
    "need to add doc"},
  { "STflowSphere_Vz",
    canalyticalSolutions_STflowSphere_Vz,
    METH_VARARGS,
    "need to add doc"},
  { "STflowSphere_Vy",
    canalyticalSolutions_STflowSphere_Vy,
    METH_VARARGS,
    "need to add doc"},
  { "STflowSphere_Vx",
    canalyticalSolutions_STflowSphere_Vx,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin1D",
    canalyticalSolutions_diffusionSin1D,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin2D",
    canalyticalSolutions_diffusionSin2D,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin3D",
    canalyticalSolutions_diffusionSin3D,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin1D_r",
    canalyticalSolutions_diffusionSin1D_r,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin2D_r",
    canalyticalSolutions_diffusionSin2D_r,
    METH_VARARGS,
    "need to add doc"},
  { "diffusionSin3D_r",
    canalyticalSolutions_diffusionSin3D_r,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearDAE",
    canalyticalSolutions_NonlinearDAE,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearDAE_f",
    canalyticalSolutions_NonlinearDAE_f,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_SteadyState",
    canalyticalSolutions_LinearAD_SteadyState,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearAD_SteadyState",
    canalyticalSolutions_NonlinearAD_SteadyState,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine",
    canalyticalSolutions_LinearADR_Sine,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_du",
    canalyticalSolutions_LinearADR_Sine_du,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_r",
    canalyticalSolutions_LinearADR_Sine_r,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_dr",
    canalyticalSolutions_LinearADR_Sine_dr,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_advectiveVelocity",
    canalyticalSolutions_LinearADR_Sine_advectiveVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_diffusiveVelocity",
    canalyticalSolutions_LinearADR_Sine_diffusiveVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Sine_totalVelocity",
    canalyticalSolutions_LinearADR_Sine_totalVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp1D",
    canalyticalSolutions_poissonsEquationExp1D,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp2D",
    canalyticalSolutions_poissonsEquationExp2D,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp3D",
    canalyticalSolutions_poissonsEquationExp3D,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp1D_r",
    canalyticalSolutions_poissonsEquationExp1D_r,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp2D_r",
    canalyticalSolutions_poissonsEquationExp2D_r,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp3D_r",
    canalyticalSolutions_poissonsEquationExp3D_r,
    METH_VARARGS,
    "need to add doc"},
  { "poissonsEquationExp3D_dr",
    canalyticalSolutions_poissonsEquationExp3D_dr,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_DiracIC",
    canalyticalSolutions_LinearAD_DiracIC,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_DiracIC_du",
    canalyticalSolutions_LinearAD_DiracIC_du,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_DiracIC_advectiveVelocity",
    canalyticalSolutions_LinearAD_DiracIC_advectiveVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_DiracIC_diffusiveVelocity",
    canalyticalSolutions_LinearAD_DiracIC_diffusiveVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "LinearAD_DiracIC_totalVelocity",
    canalyticalSolutions_LinearAD_DiracIC_totalVelocity,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Decay_DiracIC",
    canalyticalSolutions_LinearADR_Decay_DiracIC,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Decay_DiracIC_r",
    canalyticalSolutions_LinearADR_Decay_DiracIC_r,
    METH_VARARGS,
    "need to add doc"},
  { "LinearADR_Decay_DiracIC_dr",
    canalyticalSolutions_LinearADR_Decay_DiracIC_dr,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearADR_Decay_DiracIC",
    canalyticalSolutions_NonlinearADR_Decay_DiracIC,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearADR_Decay_DiracIC_r",
    canalyticalSolutions_NonlinearADR_Decay_DiracIC_r,
    METH_VARARGS,
    "need to add doc"},
  { "NonlinearADR_Decay_DiracIC_dr",
    canalyticalSolutions_NonlinearADR_Decay_DiracIC_dr,
    METH_VARARGS,
    "need to add doc"},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcanalyticalSolutions(void)
{
  PyObject *m,*d;
  m = Py_InitModule("canalyticalSolutions", canalyticalSolutionsMethods);
  d = PyModule_GetDict(m);
  import_array();
}

/** @} */
