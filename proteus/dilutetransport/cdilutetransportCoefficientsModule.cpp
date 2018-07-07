#include "Python.h"
#include "numpy/arrayobject.h"
#include "dilutetransportCoefficients.h"
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject* cdilutetransportCoefficientsConMassFluidEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double K,grav;
  PyObject *u,*m,*dm,*f,*df,*phi,*dphi,*a,*r,*x,*mass_frac;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOO",
                       &K,
                       &grav,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &phi,
                       &dphi,
                       &a,
                       &r,
                       &x,
                       &mass_frac))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  ConMassFluidEvaluate(nPoints,
                       SHAPE(f)[ND(f)-1],
                       K,
                       grav,
                       DDATA(u),
                       DDATA(m),
                       DDATA(dm),
                       DDATA(f),
                       DDATA(df),
                       DDATA(phi),
                       DDATA(dphi),
                       DDATA(a),
                       DDATA(r),
                       DDATA(x),
                       DDATA(mass_frac));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cdilutetransportCoefficientsDiluteDispersionEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double poro,diff,alpha_L;
  PyObject *u,*m,*dm,*f,*df,*a,*da,*velocity;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOO",
                       &poro,
                       &diff,
                       &alpha_L,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a,
                       &da,
                       &velocity))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  DiluteDispersionEvaluate(nPoints,
                       SHAPE(f)[ND(f)-1],
                       poro,
                       diff,
                       alpha_L,
                       DDATA(u),
                       DDATA(m),
                       DDATA(dm),
                       DDATA(f),
                       DDATA(df),
                       DDATA(a),
                       DDATA(da),
                       DDATA(velocity));
  Py_INCREF(Py_None); 
  return Py_None;
}





static PyMethodDef cdilutetransportCoefficientsMethods[] = {
  { "ConMassFluidEvaluate", 
      cdilutetransportCoefficientsConMassFluidEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients of Darcy's Law"}, 
  { "DiluteDispersionEvaluate", 
      cdilutetransportCoefficientsDiluteDispersionEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients for Dilute Diffusion"}, 
};


PyMODINIT_FUNC initcdilutetransportCoefficients(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cdilutetransportCoefficients", cdilutetransportCoefficientsMethods);
  d = PyModule_GetDict(m);
  import_array();
}

