#include "Python.h"
#include "numpy/oldnumeric.h"
#include "transportCoefficients.h"

#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define LIDATA(p) ((int *) (((PyArrayObject *)p)->data))


static PyObject* transportCoefficientsLinearADR_ConstatCoefficientsEvaluate(PyObject* self, 
                                                                            PyObject* args)
{
  int nPoints,nSpace;
  double M,C,t;
  PyObject *B,*A,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"iidOOddOOOOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &M,
                       &A,
                       &B,
                       &C,
                       &t,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a,
                       &da,
                       &phi,
                       &dphi,
                       &r,
                       &dr))
    return NULL;
  linearADR_ConstantCoefficientsEvaluate(nPoints,
                                         nSpace,
                                         M,
                                         DDATA(A),
                                         DDATA(B),
                                         C,
                                         t,
                                         DDATA(x),
                                         DDATA(u),
                                         DDATA(m),
                                         DDATA(dm),
                                         DDATA(f),
                                         DDATA(df),
                                         DDATA(a),
                                         DDATA(da),
                                         DDATA(phi),
                                         DDATA(dphi),
                                         DDATA(r),
                                         DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsNonlinearADR_pqrstEvaluate(PyObject* self, 
                                                                 PyObject* args)
{
  int nPoints,nSpace;
  double M,C,t,p_pow,q_pow,r_pow,s_pow,t_pow;
  PyObject *B,*A,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"iidOOdddddddOOOOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &M,
                       &A,
                       &B,
                       &C,
                       &p_pow,
                       &q_pow,
                       &r_pow,
                       &s_pow,
                       &t_pow,
                       &t,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a,
                       &da,
                       &phi,
                       &dphi,
                       &r,
                       &dr))
    return NULL;
  nonlinearADR_pqrstEvaluate(nPoints,
                             nSpace,
                             M,
                             DDATA(A),
                             DDATA(B),
                             C,
                             p_pow,
                             q_pow,
                             r_pow,
                             s_pow,
                             t_pow,
                             t,
                             DDATA(x),
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df),
                             DDATA(a),
                             DDATA(da),
                             DDATA(phi),
                             DDATA(dphi),
                             DDATA(r),
                             DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsNonlinearADR_pqrstDualEvaluate(PyObject* self, 
                                                                     PyObject* args)
{
  int nPoints,nSpace;
  double M,C,t,p1,q1,r1,s1,t1,p2,q2,r2,s2,t2;
  PyObject *B,*A,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"iidOOddddddddddddOOOOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &M,
                       &A,
                       &B,
                       &C,
                       &p1,
                       &q1,
                       &r1,
                       &s1,
                       &t1,
                       &p2,
                       &q2,
                       &r2,
                       &s2,
                       &t2,
                       &t,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a,
                       &da,
                       &phi,
                       &dphi,
                       &r,
                       &dr))
    return NULL;
  nonlinearADR_pqrstDualEvaluate(nPoints,
                                 nSpace,
                                 M,
                                 DDATA(A),
                                 DDATA(B),
                                 C,
                                 p1,
                                 q1,
                                 r1,
                                 s1,
                                 t1,
                                 p2,
                                 q2,
                                 r2,
                                 s2,
                                 t2,
                                 t,
                                 DDATA(x),
                                 DDATA(u),
                                 DDATA(m),
                                 DDATA(dm),
                                 DDATA(f),
                                 DDATA(df),
                                 DDATA(a),
                                 DDATA(da),
                                 DDATA(phi),
                                 DDATA(dphi),
                                 DDATA(r),
                                 DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsUnitSquareRotationEvaluate(PyObject* self, 
                                                                 PyObject* args)
{
  int nPoints;
  int nSpace;
  PyObject *x,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"iiOOOOOO",
                       &nPoints,
                       &nSpace,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  unitSquareRotationEvaluate(nPoints,
                             nSpace,
                             DDATA(x),
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsRotatingPulseVelEvaluate(PyObject* self,
                                                               PyObject* args)
{
  int nPoints;
  int nSpace;
  double self_a;
  PyObject *x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"iidOOOOOOOOOO",
                        &nPoints,
                        &nSpace,
                        &self_a,
                        &x,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a,
                        &da,
                        &phi,
                        &dphi))
    return NULL;
  rotatingPulseVelEvaluate(nPoints,
                           nSpace,
                           self_a,
                           DDATA(x),
                           DDATA(u),
                           DDATA(m),
                           DDATA(dm),
                           DDATA(f),
                           DDATA(df),
                           DDATA(a),
                           DDATA(da),
                           DDATA(phi),
                           DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsDisRotatingPulseVelEvaluate(PyObject* self,
                                                                  PyObject* args)
{
  int nPoints;
  int nSpace;
  double self_a;
  PyObject *x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"iidOOOOOOOOOO",
                        &nPoints,
                        &nSpace,
                        &self_a,
                        &x,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a,
                        &da,
                        &phi,
                        &dphi))
    return NULL;
  disRotatingPulseVelEvaluate(nPoints,
                              nSpace,
                              self_a,
                              DDATA(x),
                              DDATA(u),
                              DDATA(m),
                              DDATA(dm),
                              DDATA(f),
                              DDATA(df),
                              DDATA(a),
                              DDATA(da),
                              DDATA(phi),
                              DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsDisVelEvaluate(PyObject* self,
                                                     PyObject* args)
{
  int nPoints;
  int nSpace;
  double self_a;
  PyObject *x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"iidOOOOOOOOOO",
                        &nPoints,
                        &nSpace,
                        &self_a,
                        &x,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a,
                        &da,
                        &phi,
                        &dphi))
    return NULL;
  disVelEvaluate(nPoints,
                 nSpace,
                 self_a,
                 DDATA(x),
                 DDATA(u),
                 DDATA(m),
                 DDATA(dm),
                 DDATA(f),
                 DDATA(df),
                 DDATA(a),
                 DDATA(da),
                 DDATA(phi),
                 DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsBurgersDiagonalVelEvaluate(PyObject* self,
                                                                 PyObject* args)
{
  int nPoints;
  int nSpace;
  double self_a;
  PyObject *u,*m,*dm,*f,*df,*a,*da,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"iidOOOOOOOOO",
                        &nPoints,
                        &nSpace,
                        &self_a,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a,
                        &da,
                        &phi,
                        &dphi))
    return NULL;
  burgersDiagonalVelEvaluate(nPoints,
                             nSpace,
                             self_a,
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df),
                             DDATA(a),
                             DDATA(da),
                             DDATA(phi),
                             DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophasePotentialFlowEvaluate(PyObject* self, 
                                                                    PyObject* args)
{
  int nPoints,nSpace;
  double t;
  PyObject *M,*B,*Bcon,*A,*C,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"iiOOOOOOOOOOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &M,
                       &A,
                       &B,
                       &Bcon,
                       &C,
                       &t,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a,
                       &da,
                       &phi,
                       &dphi,
                       &r,
                       &dr))
    return NULL;
  twophasePotentialFlowEvaluate(nPoints,
                                nSpace,
                                DDATA(M),
                                DDATA(A),
                                DDATA(B),
                                DDATA(Bcon),
                                DDATA(C),
                                t,
                                DDATA(x),
                                DDATA(u),
                                DDATA(m),
                                DDATA(dm),
                                DDATA(f),
                                DDATA(df),
                                DDATA(a),
                                DDATA(da),
                                DDATA(phi),
                                DDATA(dphi),
                                DDATA(r),
                                DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophasePotentialFlowUpdateFreeSurface(PyObject* self, 
                                                                             PyObject* args)
{
  int nPoints,nSpace;
  double eps,M1,M2,C1,C2;
  PyObject *u_levelSet,*M,*B,*Bcon,*A1,*A2,*A,*C,*B1,*B2,*Bcon1,*Bcon2;
  if(!PyArg_ParseTuple(args,"iidOddOOOOOOOOOOddO",
                       &nPoints,
                       &nSpace,
                       &eps,
                       &u_levelSet,
                       &M1,&M2,&M,
                       &A1,&A2,&A,
                       &B1,&B2,&B,
                       &Bcon1,&Bcon2,&Bcon,
                       &C1,&C2,&C))
    return NULL;
  twophasePotentialFlowUpdateFreeSurface(nPoints,
                                         nSpace,
                                         eps,
                                         DDATA(u_levelSet),
                                         M1,M2,DDATA(M),
                                         DDATA(A1),DDATA(A2),DDATA(A),
                                         DDATA(B1),DDATA(B2),DDATA(B),
                                         DDATA(Bcon1),DDATA(Bcon2),DDATA(Bcon),
                                         C1,C2,DDATA(C));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophaseLevelSetCoefficientsUpdateVelocity(PyObject* self, 
                                                                                 PyObject* args)
{
  int nPoints,nSpace;
  double  v_scale;
  PyObject *vIn,*vOut;
  if(!PyArg_ParseTuple(args,"iidOO",
                       &nPoints,
                       &nSpace,
                       &v_scale,
                       &vIn,
                       &vOut))
    return NULL;
  twophaseLevelSetCoefficientsUpdateVelocity(nPoints,
                                             nSpace,
                                             v_scale,
                                             DDATA(vIn),
                                             DDATA(vOut));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophaseLevelSetCoefficientsEvaluate(PyObject* self, 
                                                                           PyObject* args)
{
  int nPoints,nSpace;
  double  t;
  PyObject *B,*x,*u,*grad_u,*m,*dm,*h,*dh,*rh;
  if(!PyArg_ParseTuple(args,"iiOdOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &B,
                       &t,
                       &x,
                       &u,
                       &grad_u,
                       &m,&dm,
                       &h,&dh,
                       &rh))
    return NULL;
  twophaseLevelSetCoefficientsEvaluate(nPoints,
                                       nSpace,
                                       DDATA(B),
                                       t,
                                       DDATA(x),
                                       DDATA(u),
                                       DDATA(grad_u),
                                       DDATA(m),
                                       DDATA(dm),
                                       DDATA(h),
                                       DDATA(dh),
                                       DDATA(rh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophaseLevelSetCoefficientsEvaluateCI(PyObject* self, 
                                                                             PyObject* args)
{
  int nPoints,nSpace;
  double  t;
  PyObject *B,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"iiOdOOOOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &B,
                       &t,
                       &x,
                       &u,
                       &m,&dm,
                       &f,&df,
                       &a,&da,
                       &phi,&dphi,
                       &r,&dr))
    return NULL;
  twophaseLevelSetCoefficientsEvaluateCI(nPoints,
                                         nSpace,
                                         DDATA(B),
                                         t,
                                         DDATA(x),
                                         DDATA(u),
                                         DDATA(m),
                                         DDATA(dm),
                                         DDATA(f),
                                         DDATA(df),
                                         DDATA(a),
                                         DDATA(da),
                                         DDATA(phi),
                                         DDATA(dphi),
                                         DDATA(r),
                                         DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophaseSignedDistanceCoefficientsUpdateSignFunction(PyObject* self, 
                                                                                           PyObject* args)
{
  int nPoints;
  double  smearingFactor;
  PyObject *u,*S;
  if(!PyArg_ParseTuple(args,"idOO",
                       &nPoints,
                       &smearingFactor,
                       &u,
                       &S))
    return NULL;
  twophaseSignedDistanceCoefficientsUpdateSignFunction(nPoints,
                                                       smearingFactor,
                                                       DDATA(u),
                                                       DDATA(S));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsTwophaseSignedDistanceCoefficientsEvaluate(PyObject* self, 
                                                                                 PyObject* args)
{
  int nPoints,nSpace;
  PyObject *S,*u,*grad_u,*m,*dm,*h,*dh,*rh;
  if(!PyArg_ParseTuple(args,"iiOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &S,
                       &u,
                       &grad_u,
                       &m,&dm,
                       &h,&dh,
                       &rh))
    return NULL;
  twophaseSignedDistanceCoefficientsEvaluate(nPoints,
                                             nSpace,
                                             DDATA(S),
                                             DDATA(u),
                                             DDATA(grad_u),
                                             DDATA(m),
                                             DDATA(dm),
                                             DDATA(h),
                                             DDATA(dh),
                                             DDATA(rh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHomEvaluate(PyObject* self, 
                                                                                         PyObject* args)
{
  int nPoints,nSpace;
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df;
  if(!PyArg_ParseTuple(args,"iidOddddddOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
                       &m,
                       &thetaR,
                       &thetaSR,
                       &KWs,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  conservativeHeadRichardsMualemVanGenuchtenHomEvaluate(nPoints,
                                                        nSpace,
                                                        rho,
                                                        DDATA(gravity),
                                                        alpha,
                                                        n,
                                                        m,
                                                        thetaR,
                                                        thetaSR,
                                                        KWs,
                                                        DDATA(u),
                                                        DDATA(mass),
                                                        DDATA(dmass),
                                                        DDATA(f),
                                                        DDATA(df),
                                                        DDATA(a),
                                                        DDATA(da));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* transportCoefficientsConservativeSatRichardsMualemVanGenuchtenHomEvaluate(PyObject* self, 
                                                                                           PyObject* args)
{
  int nPoints,nSpace;
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*phi,*dphi;
  if(!PyArg_ParseTuple(args,"iidOddddddOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
                       &m,
                       &thetaR,
                       &thetaSR,
                       &KWs,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da,
                       &phi,
                       &dphi))
    return NULL;
  conservativeSatRichardsMualemVanGenuchtenHomEvaluate(nPoints,
                                                       nSpace,
                                                       rho,
                                                       DDATA(gravity),
                                                       alpha,
                                                       n,
                                                       m,
                                                       thetaR,
                                                       thetaSR,
                                                       KWs,
                                                       DDATA(u),
                                                       DDATA(mass),
                                                       DDATA(dmass),
                                                       DDATA(f),
                                                       DDATA(df),
                                                       DDATA(a),
                                                       DDATA(da),
                                                       DDATA(phi),
                                                       DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef transportCoefficientsMethods[] = {
  { "linearADR_ConstantCoefficientsEvaluate", 
    transportCoefficientsLinearADR_ConstatCoefficientsEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar transport equation for the linear, constant coefficients case"}, 
  { "nonlinearADR_pqrstEvaluate", 
    transportCoefficientsNonlinearADR_pqrstEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar transport equation for simple monomial coefficients"}, 
  { "nonlinearADR_pqrstDualEvaluate", 
    transportCoefficientsNonlinearADR_pqrstDualEvaluate, 
    METH_VARARGS, 
    "Coefficients for doubly degenerate pqrst equation"},
  { "unitSquareRotationEvaluate", 
    transportCoefficientsUnitSquareRotationEvaluate, 
    METH_VARARGS, 
    "Coefficients for rotating velocity field on the unit square"}, 
  { "rotatingPulseVelEvaluate", 
    transportCoefficientsRotatingPulseVelEvaluate, 
    METH_VARARGS, 
    "coefficients advection-diffusion with a  rotating velocity field"}, 
  { "disRotatingPulseVelEvaluate", 
    transportCoefficientsDisRotatingPulseVelEvaluate, 
    METH_VARARGS, 
    "transport coefficients for a discontinuous rotating pulse"}, 
  { "disVelEvaluate", 
    transportCoefficientsDisVelEvaluate, 
    METH_VARARGS, 
    "transport coefficients for a discontinuous velocity field"}, 
  { "burgersDiagonalVelEvaluate", 
    transportCoefficientsBurgersDiagonalVelEvaluate, 
    METH_VARARGS, 
    "Burgers equations coefficients for a diagonal velocity field"}, 
  { "twophasePotentialFlowEvaluate", 
    transportCoefficientsTwophasePotentialFlowEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar transport equation for spatially variable coefficients"}, 
  { "twophasePotentialFlowUpdateFreeSurface", 
    transportCoefficientsTwophasePotentialFlowUpdateFreeSurface, 
    METH_VARARGS, 
    "Update  the coefficients of a scalar transport equation based on the location of a free surface"}, 
  { "twophaseLevelSetCoefficientsUpdateVelocity", 
    transportCoefficientsTwophaseLevelSetCoefficientsUpdateVelocity, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseLevelSetCoefficientsEvaluate", 
    transportCoefficientsTwophaseLevelSetCoefficientsEvaluate, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseLevelSetCoefficientsEvaluateCI", 
    transportCoefficientsTwophaseLevelSetCoefficientsEvaluateCI, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseSignedDistanceUpdateSignFunction", 
    transportCoefficientsTwophaseSignedDistanceCoefficientsUpdateSignFunction, 
    METH_VARARGS, 
    "Update  the level set sign function"}, 
  { "twophaseSignedDistanceCoefficientsEvaluate", 
    transportCoefficientsTwophaseSignedDistanceCoefficientsEvaluate,
    METH_VARARGS, 
    "Update  the eikonal equation coefficients"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHomEvaluate", 
    transportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeSatRichardsMualemVanGenuchtenHomEvaluate", 
    transportCoefficientsConservativeSatRichardsMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { NULL,NULL,0,NULL}
};


PyMODINIT_FUNC inittransportCoefficients(void)
{
  PyObject *m,*d;
  m = Py_InitModule("transportCoefficients", transportCoefficientsMethods);
  d = PyModule_GetDict(m);
  import_array();
}
