#include "Python.h"
#include "numpy/oldnumeric.h"
#include "femIntegrals.h"
#include "superluWrappersModule.h"

#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define LIDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

static PyObject* femIntegralsUpdateMass(PyObject* self, 
                                        PyObject* args)
{
  PyObject* m;
  PyObject* w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &m,
                       &w,
                       &residual)) 
    return NULL;
  updateMass(SHAPE(w)[0],
             SHAPE(w)[1],
             SHAPE(w)[2],
             DDATA(m),
             DDATA(w),
             DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateMassJacobian(PyObject* self, 
                                                PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* dm;
  PyObject* v_x_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &dm,
                       &v_x_w,
                       &jacobian)) 
    return NULL;
  updateMassJacobian(nElements_global,
                     nQuadraturePoints_element,
                     nDOF_element,
                     DDATA(dm),
                     DDATA(v_x_w),
                     DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateAdvection(PyObject* self, 
                                             PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* f;
  PyObject* w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &f,
                       &w,
                       &residual))
    return NULL;
  updateAdvection(nElements_global,
                  nQuadraturePoints_element,
                  nDOF_element,
                  nSpace,
                  DDATA(f),
                  DDATA(w),
                  DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateAdvectionJacobian(PyObject* self, 
                                                     PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* df;
  PyObject* v_x_grad_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &df,
                       &v_x_grad_w,
                       &jacobian))
    return NULL;
  updateAdvectionJacobian(nElements_global,
                          nQuadraturePoints_element,
                          nDOF_element,
                          nSpace,
                          DDATA(df),
                          DDATA(v_x_grad_w),
                          DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateHamiltonian(PyObject* self, 
                                               PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* h;
  PyObject* w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &h,
                       &w,
                       &residual))
    return NULL;
  updateHamiltonian(nElements_global,
                    nQuadraturePoints_element,
                    nDOF_element,
                    DDATA(h),
                    DDATA(w),
                    DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateHamiltonianJacobian(PyObject* self, 
                                                       PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace_global;
  PyObject* dh;
  PyObject* grad_v_x_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace_global,
                       &dh,
                       &grad_v_x_w,
                       &jacobian))
    return NULL;
  updateHamiltonianJacobian(nElements_global,
                            nQuadraturePoints_element,
                            nDOF_element,
                            nSpace_global,
                            DDATA(dh),
                            DDATA(grad_v_x_w),
                            DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* femIntegralsUpdateDiffusion(PyObject* self, 
                                             PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* a;
  PyObject* grad_phi_x_grad_w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &a,
                       &grad_phi_x_grad_w,
                       &residual))
    return NULL;
  updateDiffusion(nElements_global,
                  nQuadraturePoints_element,
                  nDOF_element,
                  nSpace,
                  DDATA(a),
                  DDATA(grad_phi_x_grad_w),
                  DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* femIntegralsUpdateDiffusionJacobian(PyObject* self, 
                                                     PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* l2g;
  PyObject* a;
  PyObject* da;
  PyObject* grad_phi_x_grad_w;
  PyObject* dphi;
  PyObject* v;
  PyObject* grad_v_x_grad_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &l2g,
                       &a,
                       &da,
                       &grad_phi_x_grad_w,
                       &dphi,
                       &v,
                       &grad_v_x_grad_w,
                       &jacobian))
    return NULL;
  updateDiffusionJacobian(nElements_global,
                          nQuadraturePoints_element,
                          nDOF_element,
                          nSpace,
                          LIDATA(l2g),
                          DDATA(a),
                          DDATA(da),
                          DDATA(grad_phi_x_grad_w),
                          DDATA(dphi),
                          DDATA(v),
                          DDATA(grad_v_x_grad_w),
                          DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateReaction(PyObject* self, 
                           PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* r;
  PyObject* w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &r,
                       &w,
                       &residual)) 
    return NULL;
  updateReaction(nElements_global,
                 nQuadraturePoints_element,
                 nDOF_element,
                 DDATA(r),
                 DDATA(w),
                 DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateReactionJacobian(PyObject* self, 
                                   PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* dr;
  PyObject* v_x_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &dr,
                       &v_x_w,
                       &jacobian)) 
    return NULL;
  updateReactionJacobian(nElements_global,
                         nQuadraturePoints_element,
                         nDOF_element,
                         DDATA(dr),
                         DDATA(v_x_w),
                         DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateStabilization(PyObject* self, 
                                PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* tau;
  PyObject* pdeResidual;
  PyObject* LstarW;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &tau,
                       &pdeResidual,
                       &LstarW,
                       &residual)) 
    return NULL;
  updateStabilization(nElements_global,
                      nQuadraturePoints_element,
                      nDOF_element,
                      DDATA(tau),
                      DDATA(pdeResidual),
                      DDATA(LstarW),
                      DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateStabilizationJacobian(PyObject* self, 
                                        PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject* tau;
  PyObject* dpdeResidual;
  PyObject* LstarW;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &tau,
                       &dpdeResidual,
                       &LstarW,
                       &jacobian)) 
    return NULL;
  updateStabilizationJacobian(nElements_global,
                              nQuadraturePoints_element,
                              nDOF_element,
                              DDATA(tau),
                              DDATA(dpdeResidual),
                              DDATA(LstarW),
                              DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateShockCapturing(PyObject* self, 
                                 PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* numDiff;
  PyObject* grad_u_x_grad_w;
  PyObject* residual;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &numDiff,
                       &grad_u_x_grad_w,
                       &residual)) 
    return NULL;
  updateShockCapturing(nElements_global,
                       nQuadraturePoints_element,
                       nDOF_element,
                       nSpace,
                       DDATA(numDiff),
                       DDATA(grad_u_x_grad_w),
                       DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsUpdateShockCapturingJacobian(PyObject* self, 
                                         PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject* numDiff;
  PyObject* grad_v_x_grad_w;
  PyObject* jacobian;
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &numDiff,
                       &grad_v_x_grad_w,
                       &jacobian)) 
    return NULL;
  updateShockCapturingJacobian(nElements_global,
                               nQuadraturePoints_element,
                               nDOF_element,
                               nSpace,
                               DDATA(numDiff),
                               DDATA(grad_v_x_grad_w),
                               DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsCalculateDiv_f(PyObject* self, 
                           PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *df,*grad_u,*grad_v,*div_f,*ddiv_f;
  
  if(!PyArg_ParseTuple(args,"iiiiOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &df,
                       &grad_u,
                       &grad_v,
                       &div_f,
                       &ddiv_f))
    return NULL;
  calculateDiv_f(nElements_global,
                 nQuadraturePoints_element,
                 nDOF_element,
                 nSpace,
                 DDATA(df),
                 DDATA(grad_u),
                 DDATA(grad_v),
                 DDATA(div_f),
                 DDATA(ddiv_f));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateDiv_a(PyObject* self, 
                           PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *l2g,*da,*dphi,*grad_phi,*grad_u,*grad_v,*div_a,*ddiv_a;
  
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &l2g,
                       &da,
                       &dphi,
                       &grad_phi,
                       &grad_u,
                       &grad_v,
                       &div_a,
                       &ddiv_a))
    return NULL;
  calculateDiv_a(nElements_global,
                 nQuadraturePoints_element,
                 nDOF_element,
                 nSpace,
                 LIDATA(l2g),
                 DDATA(da),
                 DDATA(dphi),
                 DDATA(grad_phi),
                 DDATA(grad_u),
                 DDATA(grad_v),
                 DDATA(div_a),
                 DDATA(ddiv_a));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject*
femIntegralsCalculateScalarScalarProduct(PyObject* self, 
                                         PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  PyObject *sL,*sR,*sResult;

  if(!PyArg_ParseTuple(args,"iiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &sL,
                       &sR,
                       &sResult))
    return NULL;
  calculateScalarScalarProduct(nElements_global,
                               nQuadraturePoints_element,
                               DDATA(sL),
                               DDATA(sR),
                               DDATA(sResult));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateVectorScalarProduct(PyObject* self, 
                                         PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace;
  PyObject *v,*s,*vResult;

  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &v,
                       &s,
                       &vResult)) return NULL;
  calculateVectorScalarProduct(nElements_global,
                               nQuadraturePoints_element,
                               nSpace,
                               DDATA(v),
                               DDATA(s),
                               DDATA(vResult));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateTensorScalarProduct(PyObject* self, 
                                         PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace;
  PyObject *t,*s,*tResult;

  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &t,
                       &s,
                       &tResult)) return NULL;
  calculateTensorScalarProduct(nElements_global,
                               nQuadraturePoints_element,
                               nSpace,
                               DDATA(t),
                               DDATA(s),
                               DDATA(tResult));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsCalculateAdjointADR(PyObject* self, 
                                PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *w,*grad_w,*df,*da,*grad_phi,*dr,*LstarW;

  if(!PyArg_ParseTuple(args,"iiiiOOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &w,
                       &grad_w,
                       &df,
                       &da,
                       &grad_phi,
                       &dr,
                       &LstarW))
    return NULL;
  calculateAdjointADR(nElements_global,
                      nQuadraturePoints_element,
                      nDOF_element,
                      nSpace,
                      DDATA(w),
                      DDATA(grad_w),
                      DDATA(df),
                      DDATA(da),
                      DDATA(grad_phi),
                      DDATA(dr),
                      DDATA(LstarW));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsCalculatePDEResidualADR(PyObject* self, 
                                    PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  PyObject *div_f,*div_a,*r,*mt,*pdeResidual;
  if(!PyArg_ParseTuple(args,"iiOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &div_f,
                       &div_a,
                       &r,
                       &mt,
                       &pdeResidual)) return NULL;
  calculatePDEResidualADR(nElements_global,
                          nQuadraturePoints_element,
                          DDATA(div_f),
                          DDATA(div_a),
                          DDATA(r),
                          DDATA(mt),
                          DDATA(pdeResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculatePDEResidualJacobianADR(PyObject* self, 
                                            PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject *ddiv_f,*ddiv_a,*dr,*dmt,*v,*dpdeResidual;

  if(!PyArg_ParseTuple(args,"iiiOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &ddiv_f,
                       &ddiv_a,
                       &dr,
                       &dmt,
                       &v,
                       &dpdeResidual))
    return NULL;
  calculatePDEResidualJacobianADR(nElements_global,
                                  nQuadraturePoints_element,
                                  nDOF_element,
                                  DDATA(ddiv_f),
                                  DDATA(ddiv_a),
                                  DDATA(dr),
                                  DDATA(dmt),
                                  DDATA(v),
                                  DDATA(dpdeResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculatePDEResidualHJ(PyObject* self, 
                                   PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace_global;
  PyObject *dh,*grad_u,*rh,*mt,*pdeResidual;
  if(!PyArg_ParseTuple(args,"iiiOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace_global,
                       &dh,
                       &grad_u,
                       &rh,
                       &mt,
                       &pdeResidual)) return NULL;
  calculatePDEResidualHJ(nElements_global,
                         nQuadraturePoints_element,
                         nSpace_global,
                         DDATA(dh),
                         DDATA(grad_u),
                         DDATA(rh),
                         DDATA(mt),
                         DDATA(pdeResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculatePDEResidualJacobianHJ(PyObject* self, 
                                           PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace_global;
  PyObject *dh,*grad_v,*dmt,*v,*dpdeResidual;

  if(!PyArg_ParseTuple(args,"iiiiOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace_global,
                       &dh,
                       &grad_v,
                       &dmt,
                       &v,
                       &dpdeResidual))
    return NULL;
  calculatePDEResidualJacobianHJ(nElements_global,
                                 nQuadraturePoints_element,
                                 nDOF_element,
                                 nSpace_global,
                                 DDATA(dh),
                                 DDATA(grad_v),
                                 DDATA(dmt),
                                 DDATA(v),
                                 DDATA(dpdeResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateStabilizationADR(PyObject* self, 
                                      PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace;
  char stabilizationType;
  PyObject *elementDiameter,*df,*a,*da,*grad_phi,*dphi,*dr,*dmt,*pe,*cfl,*tau;

  if(!PyArg_ParseTuple(args,"iiicOOOOOOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &stabilizationType,
                       &elementDiameter,
                       &df,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &dr,
                       &dmt,
                       &pe,
                       &cfl,
                       &tau)) 
    return NULL;
  calculateStabilizationADR(nElements_global,
                            nQuadraturePoints_element,
                            nSpace,
                            stabilizationType,
                            DDATA(elementDiameter),
                            DDATA(df),
                            DDATA(a),
                            DDATA(da),
                            DDATA(grad_phi),
                            DDATA(dphi),
                            DDATA(dr),
                            DDATA(dmt),
                            DDATA(pe),
                            DDATA(cfl),
                            DDATA(tau));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateDimensionlessNumbersADR(PyObject* self, 
                                             PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace;
  PyObject *elementDiameter,*df,*a,*dphi,*dr,*dmt,*pe,*cfl;

  if(!PyArg_ParseTuple(args,"iiiOOOOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &elementDiameter,
                       &df,
                       &a,
                       &dphi,
                       &dr,
                       &dmt,
                       &pe,
                       &cfl)) 
    return NULL;
  calculateDimensionlessNumbersADR(nElements_global,
                                   nQuadraturePoints_element,
                                   nSpace,
                                   DDATA(elementDiameter),
                                   DDATA(df),
                                   DDATA(a),
                                   DDATA(dphi),
                                   DDATA(dr),
                                   DDATA(dmt),
                                   DDATA(pe),
                                   DDATA(cfl));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateShockCapturingADR(PyObject* self, 
                                       PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace;
  char shockCapturingType;
  double shockCapturingDiffusion;

  PyObject *elementDiameter,*pdeResidual,*mt,*grad_u,*numDiff;
  
  if(!PyArg_ParseTuple(args,"iiicdOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &shockCapturingType,
                       &shockCapturingDiffusion,
                       &elementDiameter,
                       &pdeResidual,
                       &mt,
                       &grad_u,
                       &numDiff))
    return NULL;
  calculateShockCapturingADR(nElements_global,
                             nQuadraturePoints_element,
                             nSpace,
                             shockCapturingType,
                             shockCapturingDiffusion,
                             DDATA(elementDiameter),
                             DDATA(pdeResidual),
                             DDATA(mt),
                             DDATA(grad_u),
                             DDATA(numDiff));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateAdjointHJ(PyObject* self, 
                               PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace_global;
  PyObject *w,*grad_w,*dh,*rh,*LstarW;

  if(!PyArg_ParseTuple(args,"iiiiOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace_global,
                       &w,
                       &grad_w,
                       &dh,
                       &rh,
                       &LstarW))
    return NULL;
  calculateAdjointHJ(nElements_global,
                     nQuadraturePoints_element,
                     nDOF_element,
                     nSpace_global,
                     DDATA(w),
                     DDATA(grad_w),
                     DDATA(dh),
                     DDATA(rh),
                     DDATA(LstarW));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
femIntegralsCalculateStabilizationHJ(PyObject* self, 
                                     PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace_global;
  char stabilizationType;
  PyObject *elementDiameter,*dh,*dmt,*cfl,*tau;

  if(!PyArg_ParseTuple(args,"iiicOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace_global,
                       &stabilizationType,
                       &elementDiameter,
                       &dh,
                       &dmt,
                       &cfl,
                       &tau)) 
    return NULL;
  calculateStabilizationHJ(nElements_global,
                           nQuadraturePoints_element,
                           nSpace_global,
                           stabilizationType,
                           DDATA(elementDiameter),
                           DDATA(dh),
                           DDATA(dmt),
                           DDATA(cfl),
                           DDATA(tau));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateShockCapturingHJ(PyObject* self, 
                                      PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nSpace_global;
  char shockCapturingType;
  double shockCapturingDiffusion;
  PyObject *elementDiameter,*pdeResidual,*mt,*dh,*numDiff;
  
  if(!PyArg_ParseTuple(args,"iiicdOOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace_global,
                       &shockCapturingType,
                       &shockCapturingDiffusion,
                       &elementDiameter,
                       &pdeResidual,
                       &mt,
                       &dh,
                       &numDiff))
    return NULL;
  calculateShockCapturingHJ(nElements_global,
                            nQuadraturePoints_element,
                            nSpace_global,
                            shockCapturingType,
                            shockCapturingDiffusion,
                            DDATA(elementDiameter),
                            DDATA(pdeResidual),
                            DDATA(mt),
                            DDATA(dh),
                            DDATA(numDiff));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateShape_x_Shape(PyObject* self, 
                                   PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  PyObject *v,*w,*v_x_w;
  
  if(!PyArg_ParseTuple(args,"iiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &v,
                       &w,
                       &v_x_w))
    return NULL;
  calculateShape_x_Shape(nElements_global,
                         nQuadraturePoints_element,
                         nDOF_element,
                         DDATA(v),
                         DDATA(w),
                         DDATA(v_x_w));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateShape_x_GradShape(PyObject* self, 
                                       PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *v,*grad_w,*v_x_grad_w;
  
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &v,
                       &grad_w,
                       &v_x_grad_w))
    return NULL;
  calculateShape_x_GradShape(nElements_global,
                             nQuadraturePoints_element,
                             nDOF_element,
                             nSpace,
                             DDATA(v),
                             DDATA(grad_w),
                             DDATA(v_x_grad_w));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateGradShape_x_Shape(PyObject* self, 
                                       PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *grad_v,*w,*grad_v_x_w;
  
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &grad_v,
                       &w,
                       &grad_v_x_w))
    return NULL;
  calculateGradShape_x_Shape(nElements_global,
                             nQuadraturePoints_element,
                             nDOF_element,
                             nSpace,
                             DDATA(grad_v),
                             DDATA(w),
                             DDATA(grad_v_x_w));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateGradShape_x_GradShape(PyObject* self, 
                                           PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nSpace;
  PyObject *grad_v,*grad_w,*grad_v_x_grad_w;
  
  if(!PyArg_ParseTuple(args,"iiiiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nSpace,
                       &grad_v,
                       &grad_w,
                       &grad_v_x_grad_w))
    return NULL;
  calculateGradShape_x_GradShape(nElements_global,
                                 nQuadraturePoints_element,
                                 nDOF_element,
                                 nSpace,
                                 DDATA(grad_v),
                                 DDATA(grad_w),
                                 DDATA(grad_v_x_grad_w));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateIntegrationWeights(PyObject* self, 
                                        PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  PyObject *abs_det_J,*referenceWeights,*weights;
  
  if(!PyArg_ParseTuple(args,"iiOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &abs_det_J,
                       &referenceWeights,
                       &weights))
    return NULL;
  calculateIntegrationWeights(nElements_global,
                              nQuadraturePoints_element,
                              DDATA(abs_det_J),
                              DDATA(referenceWeights),
                              DDATA(weights));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFiniteElementFunctionValues(PyObject* self, 
                                                 PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nComponents;
  PyObject *l2g,*dof,*v,*values;
  
  if(!PyArg_ParseTuple(args,"iiiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nComponents,
                       &l2g,
                       &dof,
                       &v,
                       &values))
    return NULL;
  calculateFiniteElementFunctionValues(nElements_global,
                                       nQuadraturePoints_element,
                                       nDOF_element,
                                       nComponents,
                                       LIDATA(l2g),
                                       DDATA(dof),
                                       DDATA(v),
                                       DDATA(values));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFiniteElementFunctionGradientValues(PyObject* self, 
                                                         PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nComponents;
  int nSpace;
  PyObject *l2g,*dof,*grad_v,*grad_u;
  
  if(!PyArg_ParseTuple(args,"iiiiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nComponents,
                       &nSpace,
                       &l2g,
                       &dof,
                       &grad_v,
                       &grad_u))
    return NULL;
  calculateFiniteElementFunctionGradientValues(nElements_global,
                                               nQuadraturePoints_element,
                                               nDOF_element,
                                               nComponents,
                                               nSpace,
                                               LIDATA(l2g),
                                               DDATA(dof),
                                               DDATA(grad_v),
                                               DDATA(grad_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFiniteElementFunctionGradientTensorValues(PyObject* self, 
                                                               PyObject* args)
{
  int nElements_global;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nComponents;
  int nSpace;
  PyObject *l2g,
    *dof,
    *grad_v_grad_w_q,
    *grad_u_grad_w_q;

  if(!PyArg_ParseTuple(args,"iiiiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nComponents,
                       &nSpace,
                       &l2g,
                       &dof,
                       &grad_v_grad_w_q,
                       &grad_u_grad_w_q)) 
    return NULL;
  calculateFiniteElementFunctionGradientTensorValues(nElements_global,
                                                     nQuadraturePoints_element,
                                                     nDOF_element,
                                                     nComponents,
                                                     nSpace,
                                                     LIDATA(l2g),
                                                     DDATA(dof),
                                                     DDATA(grad_v_grad_w_q),
                                                     DDATA(grad_u_grad_w_q));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFiniteElementFunctionValuesTrace(PyObject* self, 
                                                      PyObject* args)
{
  int nElements_global;
  int nElementBoundaries_element;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nComponents;
  PyObject *l2g,*dof,*v,*u;
  
  if(!PyArg_ParseTuple(args,"iiiiiOOOO",
                       &nElements_global,
                       &nElementBoundaries_element,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nComponents,
                       &l2g,
                       &dof,
                       &v,
                       &u))
    return NULL;
  calculateFiniteElementFunctionValuesTrace(nElements_global,
                                            nQuadraturePoints_element,
                                            nElementBoundaries_element,
                                            nDOF_element,
                                            nComponents,
                                            LIDATA(l2g),
                                            DDATA(dof),
                                            DDATA(v),
                                            DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFiniteElementFunctionGradientValuesTrace(PyObject* self, 
                                                              PyObject* args)
{
  int nElements_global;
  int nElementBoundaries_element;
  int nQuadraturePoints_element;
  int nDOF_element;
  int nComponents;
  int nSpace;
  PyObject *l2g,*dof,*grad_v,*grad_u;
  
  if(!PyArg_ParseTuple(args,"iiiiiiOOOO",
                       &nElements_global,
                       &nElementBoundaries_element,
                       &nQuadraturePoints_element,
                       &nDOF_element,
                       &nComponents,
                       &nSpace,
                       &l2g,
                       &dof,
                       &grad_v,
                       &grad_u))
    return NULL;
  calculateFiniteElementFunctionGradientValuesTrace(nElements_global,
                                                    nElementBoundaries_element,
                                                    nQuadraturePoints_element,
                                                    nDOF_element,
                                                    nComponents,
                                                    nSpace,
                                                    LIDATA(l2g),
                                                    DDATA(dof),
                                                    DDATA(grad_v),
                                                    DDATA(grad_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateGlobalResidualFromElementResidual(PyObject* self, 
                                                    PyObject* args)
{
  int nElements_global;
  int nDOF_element;
  PyObject *nFreeDOF_element,*freeLocal,*freeGlobal,*elementResidual,*globalResidual;
  
  if(!PyArg_ParseTuple(args,"iiOOOOO",
                       &nElements_global,
                       &nDOF_element,
                       &nFreeDOF_element,
                       &freeLocal,
                       &freeGlobal,
                       &elementResidual,
                       &globalResidual))
    return NULL;
  updateGlobalResidualFromElementResidual(nElements_global,
                                          nDOF_element,
                                          LIDATA(nFreeDOF_element),
                                          LIDATA(freeLocal),
                                          LIDATA(freeGlobal),
                                          DDATA(elementResidual),
                                          DDATA(globalResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateGlobalJacobianFromElementJacobian_CSR(PyObject* self, 
                                                        PyObject* args)
{
  int nElements_global;
  int nDOF_element;
  PyObject *nFreeDOF_element,*csrRowIndeces,*csrColumnOffsets,*freeLocal,*elementJacobian,*globalJacobianCSR;
  
  if(!PyArg_ParseTuple(args,"iiOOOOOO",
                       &nElements_global,
                       &nDOF_element,
                       &nFreeDOF_element,
                       &freeLocal,
                       &csrRowIndeces,
                       &csrColumnOffsets,
                       &elementJacobian,
                       &globalJacobianCSR))
    return NULL;
  updateGlobalJacobianFromElementJacobian_CSR(nElements_global,
                                              nDOF_element,
                                              LIDATA(nFreeDOF_element),
                                              LIDATA(freeLocal),
                                              LIDATA(csrRowIndeces),
                                              LIDATA(csrColumnOffsets),
                                              DDATA(elementJacobian),
                                              CSRVAL(globalJacobianCSR));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateGlobalJacobianFromElementJacobian_dense(PyObject* self, 
                                                        PyObject* args)
{
  int nElements_global;
  int nDOF_element;
  int nFreeDOF_global;
  PyObject *nFreeDOF_element,*freeLocal,*freeGlobal,*elementJacobian,*jacobian;
  
  if(!PyArg_ParseTuple(args,"iiiOOOOO",
                       &nElements_global,
                       &nDOF_element,
		       &nFreeDOF_global,
                       &nFreeDOF_element,
                       &freeLocal,
		       &freeGlobal,
                       &elementJacobian,
                       &jacobian))
    return NULL;
  updateGlobalJacobianFromElementJacobian_dense(nElements_global,
						nDOF_element,
						nFreeDOF_global,
						LIDATA(nFreeDOF_element),
						LIDATA(freeLocal),
						LIDATA(freeGlobal),
						DDATA(elementJacobian),
						DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateFlowVelocity(PyObject* self, 
                                  PyObject* args)
{
  int nElements_global,nQuadraturePoints_element,nSpace;
  PyObject *f,*a,*grad_phi,*v;
  
  if(!PyArg_ParseTuple(args,"iiiOOOO",
                       &nElements_global,
                       &nQuadraturePoints_element,
                       &nSpace,
                       &f,
                       &a,
                       &grad_phi,
                       &v))
    return NULL;
  calculateFlowVelocity(nElements_global,
                        nQuadraturePoints_element,
                        nSpace,
                        DDATA(f),
                        DDATA(a),
                        DDATA(grad_phi),
                        DDATA(v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateAddJacobian_CSR(PyObject* self, 
                                  PyObject* args)
{
  int jacIndex;
  double  val;
  PyObject *jac;
  
  if(!PyArg_ParseTuple(args,"idO",
                       &jacIndex,
                       &val,
                       &jac))
    return NULL;
  updateAddJacobian_CSR(jacIndex,val,CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsZeroJacobian_CSR(PyObject* self, 
                             PyObject* args)
{
  int nNonzeros;
  PyObject *jac;
  if(!PyArg_ParseTuple(args,"iO",
                       &nNonzeros,
                       &jac))
    return NULL;
  zeroJacobian_CSR(nNonzeros,CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(PyObject* self, 
                                                                            PyObject* args)
{
  int nInteriorElementBoundaries_global,nDOF_element,nQuadraturePoints_elementBoundary,nElementBoundaries_element;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*nFreeDOF_element,*freeLocal,*csrRowIndeces,*csrColumnOffsets_eb,*elementBoundaryFluxJacobian,*w,*jac;  
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nDOF_element,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element,
                       &freeLocal,
                       &csrRowIndeces,
                       &csrColumnOffsets_eb,
                       &elementBoundaryFluxJacobian,
                       &w,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(nInteriorElementBoundaries_global,
                                                                  nDOF_element,
                                                                  nQuadraturePoints_elementBoundary,
                                                                  nElementBoundaries_element,
                                                                  LIDATA(interiorElementBoundaries),
                                                                  LIDATA(elementBoundaryElements),
                                                                  LIDATA(elementBoundaryLocalElementBoundaries),
                                                                  LIDATA(nFreeDOF_element),
                                                                  LIDATA(freeLocal),
                                                                  LIDATA(csrRowIndeces),
                                                                  LIDATA(csrColumnOffsets_eb),
                                                                  DDATA(elementBoundaryFluxJacobian),
                                                                  DDATA(w),
                                                                  CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject*
femIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(PyObject* self, 
                                                                            PyObject* args)
{
  int nExteriorElementBoundaries_global,nDOF_element,nQuadraturePoints_elementBoundary,nElementBoundaries_element;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*nFreeDOF_element,*freeLocal,*csrRowIndeces,*csrColumnOffsets_eb,*elementBoundaryFluxJacobian,*w,*jac;  
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nDOF_element,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element,
                       &freeLocal,
                       &csrRowIndeces,
                       &csrColumnOffsets_eb,
                       &elementBoundaryFluxJacobian,
                       &w,
                       &jac))
    return NULL;
  updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(nExteriorElementBoundaries_global,
                                                                  nDOF_element,
                                                                  nQuadraturePoints_elementBoundary,
                                                                  nElementBoundaries_element,
                                                                  LIDATA(exteriorElementBoundaries),
                                                                  LIDATA(elementBoundaryElements),
                                                                  LIDATA(elementBoundaryLocalElementBoundaries),
                                                                  LIDATA(nFreeDOF_element),
                                                                  LIDATA(freeLocal),
                                                                  LIDATA(csrRowIndeces),
                                                                  LIDATA(csrColumnOffsets_eb),
                                                                  DDATA(elementBoundaryFluxJacobian),
                                                                  DDATA(w),
                                                                  CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject*
femIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(PyObject* self, 
									      PyObject* args)
{
  int nInteriorElementBoundaries_global,nDOF_element,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nFreeDOF_global;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*nFreeDOF_element,*freeLocal,*freeGlobal,*elementBoundaryFluxJacobian,*w,*jac;  
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nDOF_element,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
		       &nFreeDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element,
                       &freeLocal,
                       &freeGlobal,
		       &elementBoundaryFluxJacobian,
                       &w,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(nInteriorElementBoundaries_global,
								    nDOF_element,
								    nQuadraturePoints_elementBoundary,
								    nElementBoundaries_element,
								    nFreeDOF_global,
								    LIDATA(interiorElementBoundaries),
								    LIDATA(elementBoundaryElements),
								    LIDATA(elementBoundaryLocalElementBoundaries),
								    LIDATA(nFreeDOF_element),
								    LIDATA(freeLocal),
								    LIDATA(freeGlobal),
								    DDATA(elementBoundaryFluxJacobian),
								    DDATA(w),
								    DDATA(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject*
femIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(PyObject* self, 
                                                                            PyObject* args)
{
  int nExteriorElementBoundaries_global,nDOF_element,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nFreeDOF_global;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*nFreeDOF_element,*freeLocal,*freeGlobal,*elementBoundaryFluxJacobian,*w,*jac;  
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nDOF_element,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
		       &nFreeDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element,
                       &freeLocal,
                       &freeGlobal,
		       &elementBoundaryFluxJacobian,
                       &w,
                       &jac))
    return NULL;
  updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(nExteriorElementBoundaries_global,
								    nDOF_element,
								    nQuadraturePoints_elementBoundary,
								    nElementBoundaries_element,
								    nFreeDOF_global,
								    LIDATA(exteriorElementBoundaries),
								    LIDATA(elementBoundaryElements),
								    LIDATA(elementBoundaryLocalElementBoundaries),
								    LIDATA(nFreeDOF_element),
								    LIDATA(freeLocal),
								    LIDATA(freeGlobal),
								    DDATA(elementBoundaryFluxJacobian),
								    DDATA(w),
								    DDATA(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateInteriorElementBoundaryFlux(PyObject* self, 
                                              PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*flux,*w,*residual;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  updateInteriorElementBoundaryFlux(nInteriorElementBoundaries_global,
                                    nQuadraturePoints_elementBoundary,
                                    nElementBoundaries_element,
                                    nDOF_element,
                                    LIDATA(interiorElementBoundaries),
                                    LIDATA(elementBoundaryElements),
                                    LIDATA(elementBoundaryLocalElementBoundaries),
                                    DDATA(flux),
                                    DDATA(w),
                                    DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateExteriorElementBoundaryFlux(PyObject* self, 
                                              PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element;
  PyObject *exterioElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*flux,*w,*residual;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &exterioElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  updateExteriorElementBoundaryFlux(nExteriorElementBoundaries_global,
                                    nQuadraturePoints_elementBoundary,
                                    nElementBoundaries_element,
                                    nDOF_element,
                                    LIDATA(exterioElementBoundaries),
                                    LIDATA(elementBoundaryElements),
                                    LIDATA(elementBoundaryLocalElementBoundaries),
                                    DDATA(flux),
                                    DDATA(w),
                                    DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateInteriorNumericalAdvectiveFlux(PyObject* self, 
                                                    PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*f,*df,*dx_f,*flux;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &f,
                       &df,
                       &dx_f,
                       &flux))
    return NULL;
  calculateInteriorNumericalAdvectiveFlux(nInteriorElementBoundaries_global,
                                          nQuadraturePoints_elementBoundary,
                                          nElementBoundaries_element,
                                          nSpace,
                                          LIDATA(interiorElementBoundaries),
                                          LIDATA(elementBoundaryElements),
                                          LIDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(f),
                                          DDATA(df),
                                          DDATA(dx_f),
                                          DDATA(flux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateInteriorNumericalAdvectiveFluxJacobian(PyObject* self, 
                                                         PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element,nSpace;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*df,*v,*dx_f,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &nSpace,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &df,
		       &v,
                       &dx_f,
                       &fluxJacobian))
    return NULL;
  updateInteriorNumericalAdvectiveFluxJacobian(nInteriorElementBoundaries_global,
                                               nQuadraturePoints_elementBoundary,
                                               nElementBoundaries_element,
                                               nDOF_element,
                                               nSpace,
                                               LIDATA(interiorElementBoundaries),
                                               LIDATA(elementBoundaryElements),
                                               LIDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(df),
					       DDATA(v),
                                               DDATA(dx_f),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateExteriorNumericalAdvectiveFlux(PyObject* self, 
                                                    PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*f,*df,*dx_f,*flux,*inflowBoundary,*inflowFlux;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowBoundary,
                       &inflowFlux,
                       &n,
                       &f,
                       &df,
                       &dx_f,
                       &flux))
    return NULL;
  calculateExteriorNumericalAdvectiveFlux(nExteriorElementBoundaries_global,
                                          nQuadraturePoints_elementBoundary,
                                          nElementBoundaries_element,
                                          nSpace,
                                          LIDATA(exteriorElementBoundaries),
                                          LIDATA(elementBoundaryElements),
                                          LIDATA(elementBoundaryLocalElementBoundaries),
                                          LIDATA(inflowBoundary),
                                          DDATA(inflowFlux),
                                          DDATA(n),
                                          DDATA(f),
                                          DDATA(df),
                                          DDATA(dx_f),
                                          DDATA(flux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsSetInflowFlux(PyObject* self, 
                                                    PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary;
  PyObject *exteriorElementBoundaries,*flux,*inflowFlux;
  if(!PyArg_ParseTuple(args,"iiOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &exteriorElementBoundaries,
                       &flux,
                       &inflowFlux))
    return NULL;
  setInflowFlux(nExteriorElementBoundaries_global,
                nQuadraturePoints_elementBoundary,
                LIDATA(exteriorElementBoundaries),
                DDATA(flux),
                DDATA(inflowFlux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateExteriorNumericalAdvectiveFluxJacobian(PyObject* self, 
                                                         PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element,nSpace;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*df,*v,*dx_f,*fluxJacobian,*inflowBoundary;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &nSpace,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowBoundary,
                       &n,
                       &df,
		       &v,
                       &dx_f,
                       &fluxJacobian))
    return NULL;
  updateExteriorNumericalAdvectiveFluxJacobian(nExteriorElementBoundaries_global,
                                               nQuadraturePoints_elementBoundary,
                                               nElementBoundaries_element,
                                               nDOF_element,
                                               nSpace,
                                               LIDATA(exteriorElementBoundaries),
                                               LIDATA(elementBoundaryElements),
                                               LIDATA(elementBoundaryLocalElementBoundaries),
                                               LIDATA(inflowBoundary),
                                               DDATA(n),
                                               DDATA(df),
					       DDATA(v),
                                               DDATA(dx_f),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateInteriorNumericalDiffusiveFlux(PyObject* self, 
                                                    PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*dx_a,*flux;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &dx_a,
                       &flux))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux(nInteriorElementBoundaries_global,
                                          nQuadraturePoints_elementBoundary,
                                          nElementBoundaries_element,
                                          nSpace,
                                          LIDATA(interiorElementBoundaries),
                                          LIDATA(elementBoundaryElements),
                                          LIDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(a),
                                          DDATA(grad_phi),
                                          DDATA(u),
                                          DDATA(penalty),
                                          DDATA(dx_a),
                                          DDATA(flux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateInteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                         PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element,nSpace;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*dx_a,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &nSpace,
                       &l2g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &dx_a,
                       &fluxJacobian))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian(nInteriorElementBoundaries_global,
                                               nQuadraturePoints_elementBoundary,
                                               nElementBoundaries_element,
                                               nDOF_element,
                                               nSpace,
                                               LIDATA(l2g),
                                               LIDATA(interiorElementBoundaries),
                                               LIDATA(elementBoundaryElements),
                                               LIDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(a),
                                               DDATA(da),
                                               DDATA(grad_phi),
                                               DDATA(dphi),
                                               DDATA(v),
                                               DDATA(grad_v),
                                               DDATA(penalty),
                                               DDATA(dx_a),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateExteriorNumericalDiffusiveFlux(PyObject* self, 
                                                    PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*dx_a,*flux;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &dx_a,
                       &flux))
    return NULL;
  calculateExteriorNumericalDiffusiveFlux(nExteriorElementBoundaries_global,
                                          nQuadraturePoints_elementBoundary,
                                          nElementBoundaries_element,
                                          nSpace,
                                          LIDATA(exteriorElementBoundaries),
                                          LIDATA(elementBoundaryElements),
                                          LIDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(a),
                                          DDATA(grad_phi),
                                          DDATA(u),
                                          DDATA(penalty),
                                          DDATA(dx_a),
                                          DDATA(flux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsUpdateExteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                         PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nDOF_element,nSpace;
  PyObject *exterioElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*dx_a,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nDOF_element,
                       &nSpace,
                       &l2g,
                       &exterioElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &dx_a,
                       &fluxJacobian))
    return NULL;
  updateExteriorNumericalDiffusiveFluxJacobian(nExteriorElementBoundaries_global,
                                               nQuadraturePoints_elementBoundary,
                                               nElementBoundaries_element,
                                               nDOF_element,
                                               nSpace,
                                               LIDATA(l2g),
                                               LIDATA(exterioElementBoundaries),
                                               LIDATA(elementBoundaryElements),
                                               LIDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(a),
                                               DDATA(da),
                                               DDATA(grad_phi),
                                               DDATA(dphi),
                                               DDATA(v),
                                               DDATA(grad_v),
                                               DDATA(penalty),
                                               DDATA(dx_a),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateInteriorElementBoundaryVelocities(PyObject* self, 
						       PyObject* args)
{
  int nInteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*m,*a,*grad_phi,*f,*vAverage,*vJump,*mAverage,*mJump;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOOO",
                       &nInteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &m,
                       &a,
		       &grad_phi,
		       &f,
		       &vAverage,
		       &vJump,
		       &mAverage,
		       &mJump))
    return NULL;
  calculateInteriorElementBoundaryVelocities(nInteriorElementBoundaries_global,
					     nQuadraturePoints_elementBoundary,
					     nElementBoundaries_element,
					     nSpace,
					     LIDATA(interiorElementBoundaries),
					     LIDATA(elementBoundaryElements),
					     LIDATA(elementBoundaryLocalElementBoundaries),
					     DDATA(m),
					     DDATA(a),
					     DDATA(grad_phi),
					     DDATA(f),
					     DDATA(vAverage),
					     DDATA(vJump),
					     DDATA(mAverage),
					     DDATA(mJump));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateExteriorElementBoundaryVelocities(PyObject* self, 
						       PyObject* args)
{
  int nExteriorElementBoundaries_global,nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*m,*a,*grad_phi,*f,*vAverage,*vJump,*mAverage,*mJump;
  if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOOO",
                       &nExteriorElementBoundaries_global,
                       &nQuadraturePoints_elementBoundary,
                       &nElementBoundaries_element,
                       &nSpace,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &m,
                       &a,
		       &grad_phi,
		       &f,
		       &vAverage,
		       &vJump,
		       &mAverage,
		       &mJump))
    return NULL;
  calculateExteriorElementBoundaryVelocities(nExteriorElementBoundaries_global,
					     nQuadraturePoints_elementBoundary,
					     nElementBoundaries_element,
					     nSpace,
					     LIDATA(exteriorElementBoundaries),
					     LIDATA(elementBoundaryElements),
					     LIDATA(elementBoundaryLocalElementBoundaries),
					     DDATA(m),
					     DDATA(a),
					     DDATA(grad_phi),
					     DDATA(f),
					     DDATA(vAverage),
					     DDATA(vJump),
					     DDATA(mAverage),
					     DDATA(mJump));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateConservationResidualPWL(PyObject* self, 
					     PyObject* args)
{
  int nElements_global;
  int nInteriorElementBoundaries_global;
  int nExteriorElementBoundaries_global;
  int nQuadraturePoints_elementBoundary;
  int nElementBoundaries_element;
  int nNodes_element;
  int nSpace;
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nodeStarOffset;
  PyObject* nElements_node;
  PyObject* elementResidual;
  PyObject* vAverage;
  PyObject* starU;
  PyObject* w;
  PyObject* n;
  PyObject* dx;
  PyObject* conservationResidual;
  PyObject* starConservationResidual;
  PyObject* vConservative;
  PyObject* vConservative_element;
  if(!PyArg_ParseTuple(args,"iiiiiiiOOOOOOOOOOOOOOOOOOO",
		       &nElements_global,
		       &nInteriorElementBoundaries_global,
		       &nExteriorElementBoundaries_global,
		       &nQuadraturePoints_elementBoundary,
		       &nElementBoundaries_element,
		       &nNodes_element,
		       &nSpace,
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nodeStarOffset,
		       &nElements_node,
		       &elementResidual,
		       &vAverage,
		       &starU,
		       &w,
		       &n,
		       &dx,
		       &conservationResidual,
		       &starConservationResidual,
		       &vConservative,
                       &vConservative_element))
    
    return NULL;
  calculateConservationResidualPWL(nElements_global,
				   nInteriorElementBoundaries_global,
				   nExteriorElementBoundaries_global,
				   nQuadraturePoints_elementBoundary,
				   nElementBoundaries_element,
				   nNodes_element,
				   nSpace,
				   LIDATA(interiorElementBoundaries),
				   LIDATA(exteriorElementBoundaries),
				   LIDATA(elementBoundaryElements),
				   LIDATA(elementBoundaryLocalElementBoundaries),
				   LIDATA(elementNodes),
				   LIDATA(nodeStarElements),
				   LIDATA(nodeStarElementNeighbors),
				   LIDATA(nodeStarOffset),
				   LIDATA(nElements_node),
				   DDATA(elementResidual),
				   DDATA(vAverage),
				   DDATA(starU),
				   DDATA(w),
				   DDATA(n),
				   DDATA(dx),
				   DDATA(conservationResidual),
				   DDATA(starConservationResidual),
				   DDATA(vConservative),
                                   DDATA(vConservative_element));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateConservationJacobianPWL(PyObject* self, 
					     PyObject* args)
{
  int nNodes_global;
  int nNodes_internal;
  int nElements_global;
  int nInteriorElementBoundaries_global;
  int nExteriorElementBoundaries_global;
  int nQuadraturePoints_elementBoundary;
  int nElementBoundaries_element;
  int nNodes_element;
  int nSpace;
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nodeStarOffset;
  PyObject* nodeStarJacobianOffset;
  PyObject* nElements_node;
  PyObject* internalNodes;
  PyObject* starJacobian;
  PyObject* w;
  PyObject* n;
  PyObject* dx;
  if(!PyArg_ParseTuple(args,"iiiiiiiiiOOOOOOOOOOOOOOO",
                       &nNodes_global,
                       &nNodes_internal,
		       &nElements_global,
		       &nInteriorElementBoundaries_global,
		       &nExteriorElementBoundaries_global,
		       &nQuadraturePoints_elementBoundary,
		       &nElementBoundaries_element,
		       &nNodes_element,
		       &nSpace,
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nodeStarOffset,
		       &nodeStarJacobianOffset,
		       &nElements_node,
                       &internalNodes,
		       &w,
		       &n,
		       &dx,
		       &starJacobian))
    
    return NULL;
  calculateConservationJacobianPWL(nNodes_global,
                                   nNodes_internal,
                                   nElements_global,
				   nInteriorElementBoundaries_global,
				   nExteriorElementBoundaries_global,
				   nQuadraturePoints_elementBoundary,
				   nElementBoundaries_element,
				   nNodes_element,
				   nSpace,
				   LIDATA(interiorElementBoundaries),
				   LIDATA(exteriorElementBoundaries),
				   LIDATA(elementBoundaryElements),
				   LIDATA(elementBoundaryLocalElementBoundaries),
				   LIDATA(elementNodes),
				   LIDATA(nodeStarElements),
				   LIDATA(nodeStarElementNeighbors),
				   LIDATA(nodeStarOffset),
				   LIDATA(nodeStarJacobianOffset),
				   LIDATA(nElements_node),
                                   LIDATA(internalNodes),
				   DDATA(w),
				   DDATA(n),
				   DDATA(dx),
				   DDATA(starJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsCalculateConservationFluxPWL(PyObject* self, 
                                         PyObject* args)
{
  int nNodes_global,nNodes_internal;
  PyObject *nElements_node,
    *nodeStarOffsets,
    *nodeStarJacobianOffsets,
    *internalNodes,
    *starR,
    *starJ,
    *starU;
  if(!PyArg_ParseTuple(args,"iiOOOOOOO",
                       &nNodes_global,
                       &nNodes_internal,
                       &nElements_node,
                       &nodeStarOffsets,
                       &nodeStarJacobianOffsets,
                       &internalNodes,
                       &starR,
                       &starJ,
                       &starU))    
    return NULL;
  calculateConservationFluxPWL(nNodes_global,
                               nNodes_internal,
                               LIDATA(nElements_node),
                               LIDATA(nodeStarOffsets),
                               LIDATA(nodeStarJacobianOffsets),
                               LIDATA(internalNodes),
                               DDATA(starR),
                               DDATA(starJ),
                               DDATA(starU));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
femIntegralsWriteDOF(PyObject* self, 
                     PyObject* args)
{
  int j,nDOF,nComponents_DOF,component;
  const  char* format;
  double* dofptr;
  FILE* fileptr;
  PyObject *dof,*file;
  if(!PyArg_ParseTuple(args,"iiisOO",
                       &nDOF,
                       &nComponents_DOF,
                       &component,
                       &format,
                       &dof,
                       &file))
    return NULL;
  dofptr = DDATA(dof);
  fileptr = PyFile_AsFile(file);
  for(j=0;j<nDOF;j++)
    fprintf(fileptr,format,dofptr[j*nComponents_DOF+component]);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef femIntegralsMethods[] = {
  { "updateMass", 
    femIntegralsUpdateMass, 
    METH_VARARGS, 
    "Update residual with the mass integral"}, 
  { "updateMassJacobian", 
    femIntegralsUpdateMassJacobian, 
    METH_VARARGS, 
    "Update jacobian with the jacobian of the mass integral"}, 
  { "updateAdvection", 
    femIntegralsUpdateAdvection, 
    METH_VARARGS, 
    "Update residual with the advection integral value"}, 
  { "updateAdvectionJacobian", 
    femIntegralsUpdateAdvectionJacobian, 
    METH_VARARGS, 
    "Update jacobian with the jacobian of the advection integral"}, 
  { "updateHamiltonian", 
    femIntegralsUpdateHamiltonian, 
    METH_VARARGS, 
    "Update residual with the Hamiltonian integral value"}, 
  { "updateHamiltonianJacobian", 
    femIntegralsUpdateHamiltonianJacobian, 
    METH_VARARGS, 
    "Update jacobian with the jacobian of the Hamiltonian integral"}, 
  { "updateDiffusion", 
    femIntegralsUpdateDiffusion, 
    METH_VARARGS, 
    "Update residual with the diffusion integral"}, 
  { "updateDiffusionJacobian", 
    femIntegralsUpdateDiffusionJacobian, 
    METH_VARARGS, 
    "Update Jacobian with the Jacobian of the diffusion integral"}, 
  { "updateReaction", 
    femIntegralsUpdateReaction, 
    METH_VARARGS, 
    "Update residual with the reaction integral"}, 
  { "updateReactionJacobian", 
    femIntegralsUpdateReactionJacobian, 
    METH_VARARGS, 
    "Update Jacobian with the Jacobian of the reaction integral"}, 
  { "updateStabilization", 
    femIntegralsUpdateStabilization, 
    METH_VARARGS, 
    "Update residual with the stabilization integral"}, 
  { "updateStabilizationJacobian", 
    femIntegralsUpdateStabilizationJacobian, 
    METH_VARARGS, 
    "Update Jacobian with the Jacobian of the stabilization integral"}, 
  { "updateShockCapturing", 
    femIntegralsUpdateShockCapturing, 
    METH_VARARGS, 
    "Update residual with the shock capturing integral"}, 
  { "updateShockCapturingJacobian", 
    femIntegralsUpdateShockCapturingJacobian, 
    METH_VARARGS, 
    "Update Jacobian with the Jacobian of the shock capturing integral"}, 
  { "calculateDiv_f", 
    femIntegralsCalculateDiv_f, 
    METH_VARARGS, 
    "The divergence of f and its derivative w.r.t u"},
  { "calculateDiv_a", 
    femIntegralsCalculateDiv_a, 
    METH_VARARGS, 
    "The divergence of f and its derivative w.r.t u"},
  { "calculateScalarScalarProduct", 
    femIntegralsCalculateScalarScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a scalar"},
  { "calculateVectorScalarProduct", 
    femIntegralsCalculateVectorScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a vector"},
  { "calculateTensorScalarProduct", 
    femIntegralsCalculateTensorScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a tensor"},
  { "calculateAdjointADR", 
    femIntegralsCalculateAdjointADR, 
    METH_VARARGS, 
    "Calculate the action of the linearized adjoint of an Advection-Diffusion-Reaction operator on a test functions"},
  { "calculatePDEResidualADR", 
    femIntegralsCalculatePDEResidualADR, 
    METH_VARARGS, 
    "The PDE residual of an Advection-Diffusion-Reaction equation over the element"},
  { "calculatePDEResidualJacobianADR", 
    femIntegralsCalculatePDEResidualJacobianADR, 
    METH_VARARGS, 
    "The Jacobian of the PDE residual of an Advection-Diffusion-Reaction equation over the element"},
  { "calculateStabilizationADR", 
    femIntegralsCalculateStabilizationADR, 
    METH_VARARGS, 
    "The stabilization parameter for an Advection-Diffusion-Reaction equation"},
  { "calculateDimensionlessNumbersADR", 
    femIntegralsCalculateDimensionlessNumbersADR, 
    METH_VARARGS, 
    "The Courant and Peclet numbers"},
  { "calculateShockCapturingADR", 
    femIntegralsCalculateShockCapturingADR, 
    METH_VARARGS, 
    "The shock capturing diffusion"},
  { "calculateAdjointHJ", 
    femIntegralsCalculateAdjointHJ, 
    METH_VARARGS, 
    "Calculate the action of the linearized adjoint of an Hamilton-Jacobi operator on a test functions"},
  { "calculatePDEResidualHJ", 
    femIntegralsCalculatePDEResidualHJ, 
    METH_VARARGS, 
    "The PDE residual of a Hamilton-Jacobi equation over the element"},
  { "calculatePDEResidualJacobianHJ", 
    femIntegralsCalculatePDEResidualJacobianHJ, 
    METH_VARARGS, 
    "The Jacobian of the PDE residual of a Hamilton-Jacobi equation over the element"},
  { "calculateStabilizationHJ", 
    femIntegralsCalculateStabilizationHJ, 
    METH_VARARGS, 
    "The stabilization parameter for a Hamilton-Jacobi equation"},
  { "calculateShockCapturingHJ", 
    femIntegralsCalculateShockCapturingHJ, 
    METH_VARARGS, 
    "The shock capturing diffusion for a Hamilton-Jacobi equation"},
  { "calculateShape_x_Shape", 
    femIntegralsCalculateShape_x_Shape,
    METH_VARARGS, 
    "The tensor product of  shape  functions"},
  { "calculateShape_x_GradShape", 
    femIntegralsCalculateShape_x_GradShape,
    METH_VARARGS, 
    "The tensor product of  shape and shape  gradient functions"},
  { "calculateGradShape_x_Shape", 
    femIntegralsCalculateGradShape_x_Shape,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape functions"},
  { "calculateGradShape_x_GradShape", 
    femIntegralsCalculateGradShape_x_GradShape,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape gradient functions"},
  { "calculateIntegrationWeights", 
    femIntegralsCalculateIntegrationWeights,
    METH_VARARGS, 
    "The physical space integration weights"},
  { "calculateFiniteElementFunctionValues", 
    femIntegralsCalculateFiniteElementFunctionValues,
    METH_VARARGS, 
    "Calculate the values of the finite element function at quadrature points"},
  { "calculateFiniteElementFunctionGradientValues", 
    femIntegralsCalculateFiniteElementFunctionGradientValues,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points"},
  { "calculateFiniteElementFunctionValuesTrace", 
    femIntegralsCalculateFiniteElementFunctionValuesTrace,
    METH_VARARGS, 
    "Calculate the values of the finite element function at quadrature points on the element boundary"},
  { "calculateFiniteElementFunctionGradientValuesTrace", 
    femIntegralsCalculateFiniteElementFunctionGradientValuesTrace,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points on the element boundary"},
  { "calculateFiniteElementFunctionGradientTensorValues", 
    femIntegralsCalculateFiniteElementFunctionGradientTensorValues, 
    METH_VARARGS, 
    "Calculate tensor product  of finite element function gradients with test function gradients"},
  { "updateGlobalResidualFromElementResidual", 
    femIntegralsUpdateGlobalResidualFromElementResidual,
    METH_VARARGS, 
    "load  the global residual"},
  { "updateGlobalJacobianFromElementJacobian_CSR", 
    femIntegralsUpdateGlobalJacobianFromElementJacobian_CSR,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "updateGlobalJacobianFromElementJacobian_dense", 
    femIntegralsUpdateGlobalJacobianFromElementJacobian_dense,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "calculateFlowVelocity", 
    femIntegralsCalculateFlowVelocity,
    METH_VARARGS, 
    "calculate the total flow velocity"},
  { "updateAddJacobian_CSR", 
    femIntegralsUpdateAddJacobian_CSR,
    METH_VARARGS, 
    "add a jacobian entry"},
  { "zeroJacobian_CSR", 
    femIntegralsZeroJacobian_CSR,
    METH_VARARGS, 
    "zero a jacobian entry"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR", 
    femIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR", 
    femIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense", 
    femIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense", 
    femIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  {"updateInteriorElementBoundaryFlux",
   femIntegralsUpdateInteriorElementBoundaryFlux,
   METH_VARARGS,
   "update the residual with the element boundary fluxes on the interior element boundaries"},
  {"updateExteriorElementBoundaryFlux",
   femIntegralsUpdateExteriorElementBoundaryFlux,
   METH_VARARGS,
   "update the residual with the element boundary fluxes on the exterior element boundaries"},
  {"calculateInteriorNumericalAdvectiveFlux",
   femIntegralsCalculateInteriorNumericalAdvectiveFlux,
   METH_VARARGS,
   "calculate the numerical advective flux on interior element boundaries"},
  {"updateInteriorNumericalAdvectiveFluxJacobian",
   femIntegralsUpdateInteriorNumericalAdvectiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the advective flux on interior element boundaries"},
  {"calculateExteriorNumericalAdvectiveFlux",
   femIntegralsCalculateExteriorNumericalAdvectiveFlux,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"setInflowFlux",
   femIntegralsSetInflowFlux,
   METH_VARARGS,
   "set the inflow flux to the current flux"},
  {"updateExteriorNumericalAdvectiveFluxJacobian",
   femIntegralsUpdateExteriorNumericalAdvectiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the advective flux on exterior element boundaries"},
  {"calculateInteriorNumericalDiffusiveFlux",
   femIntegralsCalculateInteriorNumericalDiffusiveFlux,
   METH_VARARGS,
   "calculate the numerical diffusive flux on interior element boundaries"},
  {"updateInteriorNumericalDiffusiveFluxJacobian",
   femIntegralsUpdateInteriorNumericalDiffusiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the diffusive flux on interior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFlux",
   femIntegralsCalculateExteriorNumericalDiffusiveFlux,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries"},
  {"updateExteriorNumericalDiffusiveFluxJacobian",
   femIntegralsUpdateExteriorNumericalDiffusiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries"},
  {"calculateInteriorElementBoundaryVelocities",
   femIntegralsCalculateInteriorElementBoundaryVelocities,
   METH_VARARGS,
   "calcualte the averages and jumps of mass and velocity on the interior element boundaries"},
  {"calculateExteriorElementBoundaryVelocities",
   femIntegralsCalculateExteriorElementBoundaryVelocities,
   METH_VARARGS,
   "calcualte the averages and jumps of mass and velocity on the interior element boundaries"},
  {"calculateConservationResidualPWL",
   femIntegralsCalculateConservationResidualPWL,
   METH_VARARGS,
   "calculate the conservation residuals for star shaped subdomains"},
  {"calculateConservationJacobianPWL",
   femIntegralsCalculateConservationJacobianPWL,
   METH_VARARGS,
   "calculate the (LU-factorized) conservation jacobians for  star shaped subdomains"},
  {"calculateConservationFluxPWL",
   femIntegralsCalculateConservationFluxPWL,
   METH_VARARGS,
   "solve the mass conservation system and calculate the conservative flux and residual"},
  {"writeDOF",
   femIntegralsWriteDOF,
   METH_VARARGS,
   "write a component of the DOF to a file using the given format"},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initfemIntegrals(void)
{
  PyObject *m,*d;
  m = Py_InitModule("femIntegrals", femIntegralsMethods);
  d = PyModule_GetDict(m);
  import_array();
}
