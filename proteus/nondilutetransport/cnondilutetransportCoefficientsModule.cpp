#include "Python.h"
#include "numpy/arrayobject.h"
#include "nondilutetransportCoefficients.h"
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject* cnondilutetransportCoefficientsConMassFluidEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double K,grav,dt,poro,L;
  PyObject *u,*m,*dm,*f,*df,*phi,*dphi,*a,*da,*r,*x,*mass_frac,*mass_frac_old;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOO",
                       &K,
                       &grav,
                       &dt,
                       &poro,
                       &L,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &phi,
                       &dphi,
                       &a,
                       &da,
                       &r,
                       &x,
                       &mass_frac,
                       &mass_frac_old))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  ConMassFluidEvaluate(nPoints,
                       SHAPE(f)[ND(f)-1],
                       K,
                       grav,
                       dt,
                       poro,
                       L,
                       DDATA(u),
                       DDATA(m),
                       DDATA(dm),
                       DDATA(f),
                       DDATA(df),
                       DDATA(phi),
                       DDATA(dphi),
                       DDATA(a),
                       DDATA(da),
                       DDATA(r),
                       DDATA(x),
                       DDATA(mass_frac),
                       DDATA(mass_frac_old));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cnondilutetransportCoefficientsNonDiluteDispersionEvaluate(PyObject* self, 
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
  NonDiluteDispersionEvaluate(nPoints,
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


static PyObject* cnondilutetransportCoefficientsNonDiluteEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double poro,diff,alpha_L,R,theta,MW_a,MW_b,beta1,beta2;
  PyObject *w,*act,*m,*dm,*f,*df,*r,*dr0,*dr1,*a00,*da000,*a01,*da010,*da011,*velocity,*pressure,*grad_w,*grad_act,*phi0,*dphi00,*dphi01,*phi1,*dphi10,*dphi11,*x;
  if(!PyArg_ParseTuple(args,"dddddddddOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &poro,
                       &diff,
                       &alpha_L,
                       &R,
                       &theta,
                       &MW_a,
                       &MW_b,
                       &beta1,
                       &beta2,
                       &w,
                       &act,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &r,
                       &dr0,
                       &dr1,
                       &a00,
                       &da000,
                       &a01,
                       &da010,
                       &da011,
                       &velocity,
                       &pressure,
                       &grad_w,
                       &grad_act,
                       &phi0,
                       &dphi00,
                       &dphi01,
                       &phi1,
                       &dphi10,
                       &dphi11,
                       &x))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  NonDiluteEvaluate(nPoints,
                     SHAPE(f)[ND(f)-1],
                     poro,
                     diff,
                     alpha_L,
                     R,
                     theta,
                     MW_a,
                     MW_b,
                     beta1,
                     beta2,
                     DDATA(w),
                     DDATA(act),
                     DDATA(m),
                     DDATA(dm),
                     DDATA(f),
                     DDATA(df),
                     DDATA(r),
                     DDATA(dr0),
                     DDATA(dr1),
                     DDATA(a00),
                     DDATA(da000),
                     DDATA(a01),
                     DDATA(da010),
                     DDATA(da011),
                     DDATA(velocity),
                     DDATA(pressure),
                     DDATA(grad_w),
                     DDATA(grad_act),
                     DDATA(phi0),
                     DDATA(dphi00),
                     DDATA(dphi01),
                     DDATA(phi1),
                     DDATA(dphi10),
                     DDATA(dphi11),
                     DDATA(x));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cnondilutetransportCoefficientsNonDilutePhiEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double poro,diff,alpha_L,R,theta,MW_a,MW_b,beta1,beta2;
  PyObject *w,*act,*phi0,*dphi00,*dphi01,*phi1,*dphi10,*dphi11,*f;
  if(!PyArg_ParseTuple(args,"dddddddddOOOOOOOOO",
                       &poro,
                       &diff,
                       &alpha_L,
                       &R,
                       &theta,
                       &MW_a,
                       &MW_b,
                       &beta1,
                       &beta2,
                       &w,
                       &act,
                       &phi0,
                       &dphi00,
                       &dphi01,
                       &phi1,
                       &dphi10,
                       &dphi11,
                       &f))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  NonDilutePhiEvaluate(nPoints,
                     SHAPE(f)[ND(f)-1],
                     poro,
                     diff,
                     alpha_L,
                     R,
                     theta,
                     MW_a,
                     MW_b,
                     beta1,
                     beta2,
                     DDATA(w),
                     DDATA(act),
                     DDATA(phi0),
                     DDATA(dphi00),
                     DDATA(dphi01),
                     DDATA(phi1),
                     DDATA(dphi10),
                     DDATA(dphi11),
                     DDATA(f));
  Py_INCREF(Py_None); 
  return Py_None;
}




static PyObject* cnondilutetransportCoefficientsAdvectionEvaluate(PyObject* self, 
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
  AdvectionEvaluate(nPoints,
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







static PyObject* cnondilutetransportCoefficientsNonDiluteTESTEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double poro,diff,alpha_L,R,theta,MW_a,MW_b,beta1,beta2;
  PyObject *w,*m,*dm,*f,*df,*a00,*da000,*velocity,*pressure,*grad_w,*x;
  if(!PyArg_ParseTuple(args,"dddddddddOOOOOOOOOOO",
                       &poro,
                       &diff,
                       &alpha_L,
                       &R,
                       &theta,
                       &MW_a,
                       &MW_b,
                       &beta1,
                       &beta2,
                       &w,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a00,
                       &da000,
                       &velocity,
                       &pressure,
                       &grad_w,
                       &x))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  NonDiluteTESTEvaluate(nPoints,
                     SHAPE(f)[ND(f)-1],
                     poro,
                     diff,
                     alpha_L,
                     R,
                     theta,
                     MW_a,
                     MW_b,
                     beta1,
                     beta2,
                     DDATA(w),
                     DDATA(m),
                     DDATA(dm),
                     DDATA(f),
                     DDATA(df),
                     DDATA(a00),
                     DDATA(da000),
                     DDATA(velocity),
                     DDATA(pressure),
                     DDATA(grad_w),
                     DDATA(x));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cnondilutetransportCoefficientsNonDilutePhiTESTEvaluate(PyObject* self, 
                                                               PyObject* args)
{
  int i,nPoints=1;
  double poro,diff,alpha_L,R,theta,MW_a,MW_b,beta1,beta2;
  PyObject *w,*act,*phi0,*dphi00,*dphi01,*phi1,*dphi10,*dphi11,*f;
  if(!PyArg_ParseTuple(args,"dddddddddOOOOOOOOO",
                       &poro,
                       &diff,
                       &alpha_L,
                       &R,
                       &theta,
                       &MW_a,
                       &MW_b,
                       &beta1,
                       &beta2,
                       &w,
                       &act,
                       &phi0,
                       &dphi00,
                       &dphi01,
                       &phi1,
                       &dphi10,
                       &dphi11,
                       &f))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  NonDilutePhiTESTEvaluate(nPoints,
                     SHAPE(f)[ND(f)-1],
                     poro,
                     diff,
                     alpha_L,
                     R,
                     theta,
                     MW_a,
                     MW_b,
                     beta1,
                     beta2,
                     DDATA(w),
                     DDATA(act),
                     DDATA(phi0),
                     DDATA(dphi00),
                     DDATA(dphi01),
                     DDATA(phi1),
                     DDATA(dphi10),
                     DDATA(dphi11),
                     DDATA(f));
  Py_INCREF(Py_None); 
  return Py_None;
}

















static PyMethodDef cnondilutetransportCoefficientsMethods[] = {
  { "ConMassFluidEvaluate", 
      cnondilutetransportCoefficientsConMassFluidEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients of Darcy's Law NonDilute"}, 
  { "NonDiluteDispersionEvaluate", 
      cnondilutetransportCoefficientsNonDiluteDispersionEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients for NonDilute Diffusion"}, 
  { "NonDiluteEvaluate", 
      cnondilutetransportCoefficientsNonDiluteEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients for NonDilute Transport"}, 
  { "NonDilutePhiEvaluate", 
      cnondilutetransportCoefficientsNonDilutePhiEvaluate,
      METH_VARARGS, 
      "Evaluate the phi coefficients for NonDilute Transport"}, 
  { "AdvectionEvaluate", 
      cnondilutetransportCoefficientsAdvectionEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients for Dilute Advective Transport"}, 
  { "NonDiluteTESTEvaluate", 
      cnondilutetransportCoefficientsNonDiluteTESTEvaluate,
      METH_VARARGS, 
      "Evaluate the coefficients for NonDilute Transport"}, 
  { "NonDilutePhiTESTEvaluate", 
      cnondilutetransportCoefficientsNonDilutePhiTESTEvaluate,
      METH_VARARGS, 
      "Evaluate the phi coefficients for NonDilute Transport"}, 
};


PyMODINIT_FUNC initcnondilutetransportCoefficients(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cnondilutetransportCoefficients", cnondilutetransportCoefficientsMethods);
  d = PyModule_GetDict(m);
  import_array();
}

