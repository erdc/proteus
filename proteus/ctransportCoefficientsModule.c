#include "Python.h"
#include "numpy/arrayobject.h"
#include "transportCoefficients.h"
/** \file ctransportCoefficientsModule.c
    \defgroup ctransportCoefficients ctransportCoefficients
    \brief Python interface to transportCoefficients library
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject*
ctransportCoefficientsApplyContactLineSlip(PyObject* self, 
                                           PyObject* args)
{
  double eps;
  PyObject *isDOFBoundary,*phi,*advectiveFlux,*diffusiveFlux;
  if(!PyArg_ParseTuple(args,"dOOOO",
                       &eps,
                       &isDOFBoundary,
                       &phi,
                       &advectiveFlux,
                       &diffusiveFlux))
    return NULL;
  applyContactLineSlip(SHAPE(phi)[0],
                       SHAPE(phi)[1],
                       eps,
                       IDATA(isDOFBoundary),
                       DDATA(phi),
                       DDATA(advectiveFlux),
                       DDATA(diffusiveFlux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
ctransportCoefficientsApplyContactLineSlipJacobian(PyObject* self, 
                                           PyObject* args)
{
  double eps;
  PyObject *isDOFBoundary,*phi,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"dOOO",
                       &eps,
                       &isDOFBoundary,
                       &phi,
                       &fluxJacobian))
    return NULL;
  applyContactLineSlipJacobian(SHAPE(fluxJacobian)[0],
                               SHAPE(fluxJacobian)[1],
                               SHAPE(fluxJacobian)[2],
                               eps,
                               IDATA(isDOFBoundary),
                               DDATA(phi),
                               DDATA(fluxJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsLinearADR_ConstatCoefficientsEvaluate(PyObject* self, 
                                                                            PyObject* args)
{
  int i,nPoints=1;
  double M,C,t;
  PyObject *B,*A,*x,*u,*m,*dm,*f,*df,*a,*r,*dr;
  if(!PyArg_ParseTuple(args,"dOOddOOOOOOOOO",
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
                       &r,
                       &dr))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  linearADR_ConstantCoefficientsEvaluate(nPoints,
                                         SHAPE(f)[ND(f)-1],
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
                                         DDATA(r),
                                         DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsGroundwaterTransportCoefficientsEvaluate(PyObject* self, 
                                                                                PyObject* args)
{
  int i,nPoints=1;
  double omega,d,alpha_L,alpha_T;
  PyObject *v,*u,*m,*dm,*f,*df,*a;
  if(!PyArg_ParseTuple(args,"ddddOOOOOOO",
                       &omega,
                       &d,
                       &alpha_L,
                       &alpha_T,
                       &v,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  groundwaterTransportCoefficientsEvaluate(nPoints,
                                           SHAPE(f)[ND(f)-1],
                                           omega,
                                           d,
                                           alpha_L,
                                           alpha_T,
                                           DDATA(v),
                                           DDATA(u),
                                           DDATA(m),
                                           DDATA(dm),
                                           DDATA(f),
                                           DDATA(df),
                                           DDATA(a));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsGroundwaterBiodegradation01EvaluateFC(PyObject* self, 
									     PyObject* args)
{
  int i,nPoints=1;
  double omega,d_c,d_e,alpha_L,alpha_T,Kox_max,Kox_C,Kox_E,Kox_X,Yield,k_d;
  PyObject *v, 
    *c_c,
    *c_e,
    *c_x,
    *m_c,
    *dm_c,
    *m_e,
    *dm_e,
    *m_x,
    *dm_x,
    *f_c,
    *df_c,
    *f_e,
    *df_e,
    *a_c,
    *a_e,
    *r_c,
    *dr_c_dc,
    *dr_c_de,
    *dr_c_dx,
    *r_e,
    *dr_e_dc,
    *dr_e_de,
    *dr_e_dx,
    *r_x,
    *dr_x_dc,
    *dr_x_de,
    *dr_x_dx;
  if(!PyArg_ParseTuple(args,"dddddddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &omega,
                       &d_c,
                       &d_e,
                       &alpha_L,
                       &alpha_T,
		       &Kox_max,
		       &Kox_C,
		       &Kox_E,
		       &Kox_X,
		       &Yield,
		       &k_d,
                       &v,
                       &c_c,
                       &c_e,
                       &c_x,
		       &m_c,
		       &dm_c,
		       &m_e,
		       &dm_e,
		       &m_x,
		       &dm_x,
		       &f_c,
		       &df_c,
		       &f_e,
		       &df_e,
		       &a_c,
		       &a_e,
		       &r_c,
		       &dr_c_dc,
		       &dr_c_de,
		       &dr_c_dx,
		       &r_e,
		       &dr_e_dc,
		       &dr_e_de,
		       &dr_e_dx,
		       &r_x,
		       &dr_x_dc,
		       &dr_x_de,
		       &dr_x_dx))
    return NULL;
  for(i=0;i<ND(f_c)-1;i++)
      nPoints *= SHAPE(f_c)[i];
  groundwaterBiodegradation01EvaluateFC(nPoints,
					SHAPE(f_c)[ND(f_c)-1],
					omega,
					d_c,
					d_e,
					alpha_L,
					alpha_T,
					Kox_max,
					Kox_C,
					Kox_E,
					Kox_X,
					Yield,
					k_d,
					DDATA(v),
					DDATA(c_c),
					DDATA(c_e),
					DDATA(c_x),
					DDATA(m_c),
					DDATA(dm_c),
					DDATA(m_e),
					DDATA(dm_e),
					DDATA(m_x),
					DDATA(dm_x),
					DDATA(f_c),
					DDATA(df_c),
					DDATA(f_e),
					DDATA(df_e),
					DDATA(a_c),
					DDATA(a_e),
					DDATA(r_c),
					DDATA(dr_c_dc),
					DDATA(dr_c_de),
					DDATA(dr_c_dx),
					DDATA(r_e),
					DDATA(dr_e_dc),
					DDATA(dr_e_de),
					DDATA(dr_e_dx),
					DDATA(r_x),
					DDATA(dr_x_dc),
					DDATA(dr_x_de),
					DDATA(dr_x_dx));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsGroundwaterBryantDawsonIonExEvaluateFC(PyObject* self, 
									      PyObject* args)
{
  int i,nPoints=1;
  double omega,d_m,d_h,alpha_L,alpha_T,K_m,K_h,K_w,Z_tot;
  PyObject *v, 
    *c_m,
    *c_h,
    *m_m,
    *dm_m_m,
    *dm_m_h,
    *m_h,
    *dm_h_m,
    *dm_h_h,
    *f_m,
    *df_m,
    *f_h,
    *df_h,
    *a_m,
    *a_h,
    *phi_h,
    *dphi_h,
    *r_m,
    *dr_m_dm,
    *dr_m_dh,
    *r_h,
    *dr_h_dm,
    *dr_h_dh;

  if(!PyArg_ParseTuple(args,"dddddddddOOOOOOOOOOOOOOOOOOOOOOO",
		       &v, 
		       &c_m,
		       &c_h,
		       &m_m,
		       &dm_m_m,
		       &dm_m_h,
		       &m_h,
		       &dm_h_m,
		       &dm_h_h,
		       &f_m,
		       &df_m,
		       &f_h,
		       &df_h,
		       &a_m,
		       &a_h,
		       &phi_h,
		       &dphi_h,
		       &r_m,
		       &dr_m_dm,
		       &dr_m_dh,
		       &r_h,
		       &dr_h_dm,
		       &dr_h_dh))
    return NULL;
  for(i=0;i<ND(f_m)-1;i++)
      nPoints *= SHAPE(f_m)[i];
  groundwaterBryantDawsonIonExEvaluateFC(nPoints,
					 SHAPE(f_m)[ND(f_m)-1],
					 omega,
					 d_m,
					 d_h,
					 alpha_L,
					 alpha_T,
					 K_m,
					 K_h,
					 K_w,
					 Z_tot,
					 DDATA(v),
					 DDATA(c_m),
					 DDATA(c_h),
					 DDATA(m_m),
					 DDATA(dm_m_m),
					 DDATA(dm_m_h),
					 DDATA(m_h),
					 DDATA(dm_h_m),
					 DDATA(dm_h_h),
					 DDATA(f_m),
					 DDATA(df_m),
					 DDATA(f_h),
					 DDATA(df_h),
					 DDATA(a_m),
					 DDATA(a_h),
					 DDATA(phi_h),
					 DDATA(dphi_h),
					 DDATA(r_m),
					 DDATA(dr_m_dm),
					 DDATA(dr_m_dh),
					 DDATA(r_h),
					 DDATA(dr_h_dm),
					 DDATA(dr_h_dh));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsGroundwaterTransportCoefficientsEvaluate_hetMat(PyObject* self, 
										       PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double d;
  PyObject *materialTypes,*omega,*alpha_L,*alpha_T,*v,*u,*m,*dm,*f,*df,*a;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOO",
		       &d,
		       &materialTypes,
                       &omega,
                       &alpha_L,
                       &alpha_T,
                       &v,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a))
    return NULL;
  assert(ND(u) == 2 || ND(u) == 3);
  if (ND(u) == 2)
    {
      nSimplex=SHAPE(u)[0]; nPointsPerSimplex=SHAPE(u)[1];
    }
  else
    {
      nSimplex=SHAPE(u)[0]*SHAPE(u)[1]; nPointsPerSimplex=SHAPE(u)[2];
    }
  groundwaterTransportCoefficientsEvaluate_hetMat(nSimplex,
						  nPointsPerSimplex,
						  SHAPE(f)[ND(f)-1],
						  d,
						  IDATA(materialTypes),
						  DDATA(omega),
						  DDATA(alpha_L),
						  DDATA(alpha_T),
						  DDATA(v),
						  DDATA(u),
						  DDATA(m),
						  DDATA(dm),
						  DDATA(f),
						  DDATA(df),
						  DDATA(a));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVariablySaturatedGroundwaterTransportCoefficientsEvaluate_hetMat(PyObject* self, 
													PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double d;
  PyObject *materialTypes,*theta,*alpha_L,*alpha_T,*v,*u,*m,*dm,*f,*df,*a;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOO",
		       &d,
		       &materialTypes,
                       &theta,
                       &alpha_L,
                       &alpha_T,
                       &v,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a))
    return NULL;
  assert(ND(u) == 2 || ND(u) == 3);
  if (ND(u) == 2)
    {
      nSimplex=SHAPE(u)[0]; nPointsPerSimplex=SHAPE(u)[1];
    }
  else
    {
      nSimplex=SHAPE(u)[0]*SHAPE(u)[1]; nPointsPerSimplex=SHAPE(u)[2];
    }
  variablySaturatedGroundwaterTransportCoefficientsEvaluate_hetMat(nSimplex,
								   nPointsPerSimplex,
								   SHAPE(f)[ND(f)-1],
								   d,
								   IDATA(materialTypes),
								   DDATA(theta),
								   DDATA(alpha_L),
								   DDATA(alpha_T),
								   DDATA(v),
								   DDATA(u),
								   DDATA(m),
								   DDATA(dm),
								   DDATA(f),
								   DDATA(df),
								   DDATA(a));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVariablySaturatedGroundwaterEnergyTransportCoefficientsEvaluate_hetMat(PyObject* self, 
													      PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double density_w,density_n,specificHeat_w,specificHeat_n;
  PyObject *materialTypes,*theta,*omega,*alpha_L,*alpha_T,*density_s,*specificHeat_s,*lambda_sat,*lambda_dry,*lambda_ani,*v,*u,*m,*dm,*f,*df,*a;
  if(!PyArg_ParseTuple(args,"ddddOOOOOOOOOOOOOOOOO",
		       &density_w,
		       &density_n,
		       &specificHeat_w,
		       &specificHeat_n,
		       &materialTypes,
                       &theta,
		       &omega,
                       &alpha_L,
                       &alpha_T,
		       &density_s,
		       &specificHeat_s,
		       &lambda_sat,
		       &lambda_dry,
		       &lambda_ani,
                       &v,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df,
                       &a))
    return NULL;
  assert(ND(u) == 2 || ND(u) == 3);
  if (ND(u) == 2)
    {
      nSimplex=SHAPE(u)[0]; nPointsPerSimplex=SHAPE(u)[1];
    }
  else
    {
      nSimplex=SHAPE(u)[0]*SHAPE(u)[1]; nPointsPerSimplex=SHAPE(u)[2];
    }
  variablySaturatedGroundwaterEnergyTransportCoefficientsEvaluate_hetMat(nSimplex,
									 nPointsPerSimplex,
									 SHAPE(f)[ND(f)-1],
									 density_w,
									 density_n,
									 specificHeat_w,
									 specificHeat_n,
									 IDATA(materialTypes),
									 DDATA(theta),
									 DDATA(omega),
									 DDATA(alpha_L),
									 DDATA(alpha_T),
									 DDATA(density_s),
									 DDATA(specificHeat_s),
									 DDATA(lambda_sat),
									 DDATA(lambda_dry),
									 DDATA(lambda_ani),
									 DDATA(v),
									 DDATA(u),
									 DDATA(m),
									 DDATA(dm),
									 DDATA(f),
									 DDATA(df),
									 DDATA(a));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsNonlinearADR_pqrstEvaluate(PyObject* self, 
                                                                 PyObject* args)
{
  int i,nPoints=1;
  double M,C,t,p_pow,q_pow,r_pow,s_pow,t_pow;
  PyObject *B,*A,*x,*u,*m,*dm,*f,*df,*a,*da,*phi,*dphi,*r,*dr;
  if(!PyArg_ParseTuple(args,"dOOdddddddOOOOOOOOOOOO",
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
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  nonlinearADR_pqrstEvaluate(nPoints,
                             SHAPE(f)[ND(f)-1],
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

static PyObject* ctransportCoefficientsNonlinearADR_pqrstDualEvaluate(PyObject* self, 
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

static PyObject* ctransportCoefficientsUnitSquareRotationEvaluate(PyObject* self, 
                                                                 PyObject* args)
{
  int i,nPoints=1;
  PyObject *x,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  unitSquareRotationEvaluate(nPoints,
                             SHAPE(f)[ND(f)-1],
                             DDATA(x),
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsUnitCubeRotationEvaluate(PyObject* self, 
								PyObject* args)
{
  int i,nPoints=1;
  PyObject *x,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  unitCubeRotationEvaluate(nPoints,
			   SHAPE(f)[ND(f)-1],
			   DDATA(x),
			   DDATA(u),
			   DDATA(m),
			   DDATA(dm),
			   DDATA(f),
			   DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsRotatingPulseVelEvaluate(PyObject* self,
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
static PyObject* ctransportCoefficientsUnitSquareVortexEvaluate(PyObject* self, 
								PyObject* args)
{
  int i,nPoints=1;
  double t;
  PyObject *x,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"dOOOOOO",
		       &t,
                       &x,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  unitSquareVortexEvaluate(nPoints,
                             SHAPE(f)[ND(f)-1],
                             t,
			     DDATA(x),
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsDisRotatingPulseVelEvaluate(PyObject* self,
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

static PyObject* ctransportCoefficientsDisVelEvaluate(PyObject* self,
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

static PyObject* ctransportCoefficientsBurgersDiagonalVelEvaluate(PyObject* self,
                                                                 PyObject* args)
{
  int nPoints;
  int nSpace,i;
  double self_a;
  PyObject *self_v,*u,*m,*dm,*f,*df,*a,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"dOOOOOOOOO",
                        &self_a,
			&self_v,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a,
                        &phi,
                        &dphi))
    return NULL;
  nPoints = 1;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  nSpace = SHAPE(f)[ND(f)-1];

  burgersDiagonalVelEvaluate(nPoints,
                             nSpace,
                             self_a,
			     DDATA(self_v),
                             DDATA(u),
                             DDATA(m),
                             DDATA(dm),
                             DDATA(f),
                             DDATA(df),
                             DDATA(a),
                             DDATA(phi),
                             DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsBurgersDiagonalVelHJEvaluate(PyObject* self,
								    PyObject* args)
{
  int nPoints;
  int nSpace,i;
  double self_a;
  PyObject *self_v,*u,*gradu,*m,*dm,*H,*dH,*a,*phi,*dphi;
  if (!PyArg_ParseTuple(args,"dOOOOOOOOOO",
                        &self_a,
			&self_v,
                        &u,
			&gradu,
                        &m,
                        &dm,
                        &H,
                        &dH,
                        &a,
                        &phi,
                        &dphi))
    return NULL;
  nPoints = 1;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  nSpace = SHAPE(dH)[ND(dH)-1];

  burgersDiagonalVelHJEvaluate(nPoints,
			       nSpace,
			       self_a,
			       DDATA(self_v),
			       DDATA(u),
			       DDATA(gradu),
			       DDATA(m),
			       DDATA(dm),
			       DDATA(H),
			       DDATA(dH),
			       DDATA(a),
			       DDATA(phi),
			       DDATA(dphi));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsEvaluateBuckleyLeverettLiuExample(PyObject* self,
									 PyObject* args)
{
  int nPoints;
  int nSpace,i;
  PyObject *x,*u,*m,*dm,*f,*df,*a;
  if (!PyArg_ParseTuple(args,"OOOOOOO",
			&x,
                        &u,
                        &m,
                        &dm,
                        &f,
                        &df,
                        &a))
    return NULL;
  nPoints = 1;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  nSpace = SHAPE(f)[ND(f)-1];

  evaluateBuckleyLeverettLiuExample(nPoints,
				    nSpace,
				    DDATA(x),
				    DDATA(u),
				    DDATA(m),
				    DDATA(dm),
				    DDATA(f),
				    DDATA(df),
				    DDATA(a));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophasePotentialFlowEvaluate(PyObject* self, 
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

static PyObject* ctransportCoefficientsTwophasePotentialFlowUpdateFreeSurface(PyObject* self, 
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

static PyObject* ctransportCoefficientsTwophaseLevelSetCoefficientsUpdateVelocity(PyObject* self, 
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

static PyObject* ctransportCoefficientsTwophaseLevelSetCoefficientsEvaluate(PyObject* self, 
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

static PyObject* ctransportCoefficients_ncLevelSetCoefficientsEvaluate(PyObject* self, 
                                                                       PyObject* args)
{
  int nPoints=1,i;
  PyObject *v,*u,*grad_u,*m,*dm,*H,*dH;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &v,
                       &u,
                       &grad_u,
                       &m,
                       &dm,
                       &H,
                       &dH))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  ncLevelSetCoefficientsEvaluate(nPoints,
                                 SHAPE(dH)[ND(dH)-1],
                                 DDATA(v),
                                 DDATA(u),
                                 DDATA(grad_u),
                                 DDATA(m),
                                 DDATA(dm),
                                 DDATA(H),
                                 DDATA(dH));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_cLevelSetCoefficientsEvaluate(PyObject* self, 
                                                                       PyObject* args)
{
  int nPoints=1,i;
  PyObject *v,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &v,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  cLevelSetCoefficientsEvaluate(nPoints,
                                SHAPE(f)[ND(f)-1],
                                DDATA(v),
                                DDATA(u),
                                DDATA(m),
                                DDATA(dm),
                                DDATA(f),
                                DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficients_VOFCoefficientsEvaluate(PyObject* self, 
                                                                PyObject* args)
{
  int nPoints=1,i;
  double eps;
  PyObject *v,*phi,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"dOOOOOOO",
                       &eps,
                       &v,
                       &phi,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  VOFCoefficientsEvaluate(nPoints,
                          SHAPE(f)[ND(f)-1],
                          eps,
                          DDATA(v),
                          DDATA(phi),
                          DDATA(u),
                          DDATA(m),
                          DDATA(dm),
                          DDATA(f),
                          DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficients_levelSetCurvatureCoefficientsEvaluate(PyObject* self, 
                                                                              PyObject* args)
{
  int nPoints=1,i;
  PyObject *grad_phi,*u,*r,*dr,*f;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &grad_phi,
                       &u,
                       &f,
                       &r,
                       &dr))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  levelSetCurvatureCoefficientsEvaluate(nPoints,
                                        SHAPE(f)[ND(f)-1],
                                        DDATA(grad_phi),
                                        DDATA(u),
                                        DDATA(f),
                                        DDATA(r),
                                        DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_redistanceLevelSetCoefficientsEvaluate(PyObject* self, 
									       PyObject* args)
{
  int nPoints=1,i;
  double eps;
  PyObject *u0,*u,*grad_u,*m,*dm,*H,*dH,*r;
  if(!PyArg_ParseTuple(args,"dOOOOOOOO",
		       &eps,
                       &u0,
                       &u,
                       &grad_u,
                       &m,
                       &dm,
                       &H,
                       &dH,
		       &r))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  redistanceLevelSetCoefficientsEvaluate(nPoints,
					 SHAPE(dH)[ND(dH)-1],
					 eps,
					 DDATA(u0),
					 DDATA(u),
					 DDATA(grad_u),
					 DDATA(m),
					 DDATA(dm),
					 DDATA(H),
					 DDATA(dH),
					 DDATA(r));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_redistanceLevelSetCoefficientsWithWeakPenaltyEvaluate(PyObject* self, 
											      PyObject* args)
{
  int nPoints=1,i;
  double eps,lambda_penalty;
  PyObject *u0,*u,*grad_u,*m,*dm,*H,*dH,*r,*dr;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOO",
		       &eps,
		       &lambda_penalty,
                       &u0,
                       &u,
                       &grad_u,
                       &m,
                       &dm,
                       &H,
                       &dH,
		       &r,
		       &dr))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  redistanceLevelSetCoefficientsWithWeakPenaltyEvaluate(nPoints,
							SHAPE(dH)[ND(dH)-1],
							eps,
							lambda_penalty,
							DDATA(u0),
							DDATA(u),
							DDATA(grad_u),
							DDATA(m),
							DDATA(dm),
							DDATA(H),
							DDATA(dH),
							DDATA(r),
							DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficients_redistanceLevelSetSandFCoefficientsEvaluate(PyObject* self, 
										    PyObject* args)
{
  int nSimplex=1,nPointsPerSimplex=1;
  double eps;
  PyObject *u0,*dV,*u,*grad_u,*m,*dm,*H,*dH,*r;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOO",
		       &eps,
                       &u0,
		       &dV,
                       &u,
                       &grad_u,
                       &m,
                       &dm,
                       &H,
                       &dH,
		       &r))
    return NULL;
  nSimplex = 1;
  assert(ND(u) == 2 || ND(u) == 3);
  if (ND(u) == 2)
    {
      nSimplex=SHAPE(u)[0]; nPointsPerSimplex=SHAPE(u)[1];
    }
  else
    {
      nSimplex=SHAPE(u)[0]*SHAPE(u)[1]; nPointsPerSimplex=SHAPE(u)[2];
    }

  redistanceLevelSetSandFCoefficientsEvaluate(nSimplex,
					      nPointsPerSimplex,
					      SHAPE(dH)[ND(dH)-1],
					      eps,
					      DDATA(u0),
					      DDATA(dV),
					      DDATA(u),
					      DDATA(grad_u),
					      DDATA(m),
					      DDATA(dm),
					      DDATA(H),
					      DDATA(dH),
					      DDATA(r));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_setWeakDirichletConditionsForLevelSet(PyObject* self,
							     PyObject* args)
{
  int nElements_global,nDOF_trial_element;
  double epsilon_freeze_factor;
  PyObject * elementDiameter,
    *u_l2g,
    *u_dof,
    *freeze_nodes_tmp,
    *weakDirichletConditionFlags;

  if (!PyArg_ParseTuple(args,"iidOOOOO",
			&nElements_global,
			&nDOF_trial_element,
			&epsilon_freeze_factor,
			&elementDiameter,
			&u_l2g,
			&u_dof,
			&freeze_nodes_tmp,
			&weakDirichletConditionFlags))
    return NULL;

  setWeakDirichletConditionsForLevelSet(nElements_global,
					nDOF_trial_element,
					epsilon_freeze_factor,
					DDATA(elementDiameter),
					IDATA(u_l2g),
					DDATA(u_dof),
					IDATA(freeze_nodes_tmp),
					IDATA(weakDirichletConditionFlags));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficients_setSimpleWeakDirichletConditionsForLevelSet(PyObject* self,
                                                                                    PyObject* args)
{
  int nElements_global,nDOF_trial_element;
  double epsilon_freeze_factor;
  PyObject * elementDiameter,
    *u_l2g,
    *u_dof,
    *freeze_nodes_tmp,
    *weakDirichletConditionFlags;

  if (!PyArg_ParseTuple(args,"iidOOOOO",
			&nElements_global,
			&nDOF_trial_element,
			&epsilon_freeze_factor,
			&elementDiameter,
			&u_l2g,
			&u_dof,
			&freeze_nodes_tmp,
			&weakDirichletConditionFlags))
    return NULL;

  setSimpleWeakDirichletConditionsForLevelSet(nElements_global,
                                              nDOF_trial_element,
                                              epsilon_freeze_factor,
                                              DDATA(elementDiameter),
                                              IDATA(u_l2g),
                                              DDATA(u_dof),
                                              IDATA(freeze_nodes_tmp),
                                              IDATA(weakDirichletConditionFlags));

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* ctransportCoefficientsTwophaseLevelSetCoefficientsEvaluateCI(PyObject* self, 
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

static PyObject* ctransportCoefficientsTwophaseSignedDistanceCoefficientsUpdateSignFunction(PyObject* self, 
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

static PyObject* ctransportCoefficientsTwophaseSignedDistanceCoefficientsEvaluate(PyObject* self, 
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

/***********************************************************************
  begin mwf unnecessary duplication for two phase interface flow
 ***********************************************************************/
static PyObject* ctransportCoefficients_darcySharpInterfaceFlowEvaluate(PyObject* self, 
									PyObject* args)
{
  int nPoints,i;
  double Km,rhoM,Kp,rhoP,eps;
  PyObject *g_u,*u,*grad_u,*u_ls,*phi,*a,*f,*r,*m,
    *dphi,*da,*df,*dr,*dm;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOO",
		       &Km,
		       &rhoM,
		       &Kp,
		       &rhoP,
		       &eps,
		       &g_u,
		       &u,
		       &grad_u,
		       &u_ls,
		       &phi,
		       &a,
		       &f,
		       &r,
		       &m,
		       &dphi,
		       &da,
		       &df,
		       &dr,
		       &dm))
    return NULL;
  nPoints = 1;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];

  darcySharpInterfaceFlowEvaluate(nPoints,
				  SHAPE(f)[ND(f)-1],
				  Km,rhoM,
				  Kp,rhoP,
				  eps,
				  DDATA(g_u),
				  DDATA(u),
				  DDATA(grad_u),
				  DDATA(u_ls),
				  DDATA(phi),
				  DDATA(a),
				  DDATA(f),
				  DDATA(r),
				  DDATA(m),
				  DDATA(dphi),
				  DDATA(da),
				  DDATA(df),
				  DDATA(dr),
				  DDATA(dm));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_darcySharpInterfaceFlowImEvaluate(PyObject* self, 
									  PyObject* args)
{
  int nPoints,i;
  double Km,rhoM,Kp,rhoP,eps;
  PyObject *g_u,*u,*grad_u,*u_ls,*phi,*a,*f,*r,*m,
    *dphi,*da,*df,*dr,*dm;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOO",
		       &Km,
		       &rhoM,
		       &Kp,
		       &rhoP,
		       &eps,
		       &g_u,
		       &u,
		       &grad_u,
		       &u_ls,
		       &phi,
		       &a,
		       &f,
		       &r,
		       &m,
		       &dphi,
		       &da,
		       &df,
		       &dr,
		       &dm))
    return NULL;
  nPoints = 1;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];

  darcySharpInterfaceFlowImEvaluate(nPoints,
				    SHAPE(f)[ND(f)-1],
				    Km,rhoM,
				    Kp,rhoP,
				    eps,
				    DDATA(g_u),
				    DDATA(u),
				    DDATA(grad_u),
				    DDATA(u_ls),
				    DDATA(phi),
				    DDATA(a),
				    DDATA(f),
				    DDATA(r),
				    DDATA(m),
				    DDATA(dphi),
				    DDATA(da),
				    DDATA(df),
				    DDATA(dr),
				    DDATA(dm));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsConservativeTotalHeadRichardsMualemVanGenuchtenHomEvaluate(PyObject* self, 
												  PyObject* args)
{
  int i,nPoints=1;
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*x,*u,*mass,*dmass,*a,*da,*f,*df,*phi,*dphi;
  if(!PyArg_ParseTuple(args,"dOOddddddOOOOOOOOO",
                       &rho,
                       &gravity,
                       &x,
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeTotalHeadRichardsMualemVanGenuchtenHomEvaluate(nPoints,
							     SHAPE(f)[ND(f)-1],
							     rho,
							     DDATA(gravity),
							     DDATA(x),
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

/***********************************************************************
  end mwf unnecessary duplication for two phase interface flow
 ***********************************************************************/
static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHomEvaluate(PyObject* self, 
                                                                                         PyObject* args)
{
  int i,nPoints=1;
  double rho,beta,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*x,*u,*mass,*dmass,*a,*da,*f,*df,*phi,*dphi;
  if(!PyArg_ParseTuple(args,"ddOOddddddOOOOOOOOO",
                       &rho,
                       &beta,
                       &gravity,
                       &x,
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeHeadRichardsMualemVanGenuchtenHomEvaluate(nPoints,
                                                        SHAPE(f)[ND(f)-1],
                                                        rho,
                                                        beta,
                                                        DDATA(gravity),
                                                        DDATA(x),
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
/***********************************************************************
  mwf more additions for RE
 ***********************************************************************/
static PyObject* ctransportCoefficientsConservativeHeadRichardsL2projMualemVanGenuchtenHomEvaluate(PyObject* self, 
												   PyObject* args)
{
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*w_dV;
  if(!PyArg_ParseTuple(args,"dOddddddOOOOOOOO",
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
                       &m,
                       &thetaR,
                       &thetaSR,
                       &KWs,
		       &w_dV,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
/*   for(i=0;i<ND(f)-2;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/*   if (ND(f) > 2) */
/*     nPointsPerSimplex = SHAPE(f)[ND(f)-2]; */
  assert(ND(f) == 3);
  conservativeHeadRichardsL2projMualemVanGenuchtenHomEvaluate(SHAPE(f)[0],
							      SHAPE(f)[1],
							      SHAPE(f)[2],
							      rho,
							      DDATA(gravity),
							      alpha,
							      n,
							      m,
							      thetaR,
							      thetaSR,
							      KWs,
							      DDATA(w_dV),
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
static PyObject* ctransportCoefficientsConservativeHeadRichardsL2projBndMualemVanGenuchtenHomEvaluate(PyObject* self, 
												      PyObject* args)
{
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*w_dV;
  if(!PyArg_ParseTuple(args,"dOddddddOOOOOOOO",
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
                       &m,
                       &thetaR,
                       &thetaSR,
                       &KWs,
		       &w_dV,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
/*   for(i=0;i<ND(f)-2;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/*   if (ND(f) > 2) */
/*     nPointsPerSimplex = SHAPE(f)[ND(f)-2]; */
  assert(ND(f) == 4);
  conservativeHeadRichardsL2projBndMualemVanGenuchtenHomEvaluate(SHAPE(f)[0],
								 SHAPE(f)[1],
								 SHAPE(f)[2],
								 SHAPE(f)[3],
								 rho,
								 DDATA(gravity),
								 alpha,
								 n,
								 m,
								 thetaR,
								 thetaSR,
								 KWs,
								 DDATA(w_dV),
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
static PyObject* ctransportCoefficientsConservativeHeadRichardsL2projMualemVanGenuchtenHetEvaluate(PyObject* self, 
												   PyObject* args)
{
  int i,nPoints=1,nPointsPerSimplex=1;
  double rho;
  PyObject *alpha,*n,*thetaR,*thetaSR,*KWs,
    *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*dV;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOOOO",
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
                       &thetaR,
                       &thetaSR,
                       &KWs,
		       &dV,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0;i<ND(f)-2;i++)
      nPoints *= SHAPE(f)[i];
  if (ND(f) > 2)
    nPointsPerSimplex = SHAPE(f)[ND(f)-2];
  conservativeHeadRichardsL2projMualemVanGenuchtenHetEvaluate(nPoints,
							      nPointsPerSimplex,
							      SHAPE(f)[ND(f)-1],
							      rho,
							      DDATA(gravity),
							      DDATA(alpha),
							      DDATA(n),
							      DDATA(thetaR),
							      DDATA(thetaSR),
							      DDATA(KWs),
							      DDATA(dV),
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
static PyObject* ctransportCoefficientsL2projectEvaluate(PyObject* self,PyObject* args)
{
  int i,nPoints=1,nPointsPerSimplex=1,rank,offset,dim=1;
  PyObject *r,*dV;
  if(!PyArg_ParseTuple(args,"iOO",
		       &rank,
                       &dV,
		       &r))
    return NULL;
  offset = 1+rank;
  for(i=0;i<ND(r)-offset;i++)
      nPoints *= SHAPE(r)[i];
  if (ND(r) > offset)
    nPointsPerSimplex = SHAPE(r)[ND(r)-offset];
  if (rank > 0)
    dim = SHAPE(r)[ND(r)-1]; 
  if (rank == 0)
    {
      /*scalar*/
      l2projectScalar(nPoints,
		      nPointsPerSimplex,
		      DDATA(dV),
		      DDATA(r));
    }
  else if (rank == 1)
    {
      /*vector*/
      l2projectVector(nPoints,
		      nPointsPerSimplex,
		      dim,
		      DDATA(dV),
		      DDATA(r));
      
    }
  else
    {
      assert(rank == 2);
      /*2 tensor*/
      l2project2Tensor(nPoints,
		       nPointsPerSimplex,
		       dim,
		       DDATA(dV),
		       DDATA(r));
      
    }

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchten_sd_het(PyObject* self, 
											 PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double rho,beta;
  int linearize_at_zero=0;
  double pc_eps=1.0e-4;
  PyObject *rowptr,*colind,*materialTypes,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n,*thetaR,*thetaSR,*vol_frac;
  if(!PyArg_ParseTuple(args,"OOOddOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &alpha,
                       &n,
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
		       &vol_frac,
		       &linearize_at_zero,
		       &pc_eps))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  if (linearize_at_zero == 1)
    {
      conservativeHeadRichardsMualemVanGenuchten_sd_het_linearized_at_saturation(nSimplex,
										 nPointsPerSimplex,
										 SHAPE(f)[ND(f)-1],
										 pc_eps,
										 IDATA(rowptr),
										 IDATA(colind),
										 IDATA(materialTypes),
										 rho,
										 beta,
										 DDATA(gravity),
										 DDATA(alpha),
										 DDATA(n),
										 DDATA(thetaR),
										 DDATA(thetaSR),
										 DDATA(KWs),
										 DDATA(u),
										 DDATA(mass),
										 DDATA(dmass),
										 DDATA(f),
										 DDATA(df),
										 DDATA(a),
										 DDATA(da),
										 DDATA(vol_frac));

    }
  else
    {
      conservativeHeadRichardsMualemVanGenuchten_sd_het(nSimplex,
							nPointsPerSimplex,
							SHAPE(f)[ND(f)-1],
							pc_eps,
							IDATA(rowptr),
							IDATA(colind),
							IDATA(materialTypes),
							rho,
							beta,
							DDATA(gravity),
							DDATA(alpha),
							DDATA(n),
							DDATA(thetaR),
							DDATA(thetaSR),
							DDATA(KWs),
							DDATA(u),
							DDATA(mass),
							DDATA(dmass),
							DDATA(f),
							DDATA(df),
							DDATA(a),
							DDATA(da),
							DDATA(vol_frac));
    }
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(PyObject* self, 
											       PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double rho,beta;
  PyObject *materialTypes,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n,*thetaR,*thetaSR;
  if(!PyArg_ParseTuple(args,"OddOOOOOOOOOOOOO",
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &alpha,
                       &n,
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
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(nSimplex,
                                                          nPointsPerSimplex,
                                                          SHAPE(f)[ND(f)-1],
                                                          IDATA(materialTypes),
                                                          rho,
                                                          beta,
                                                          DDATA(gravity),
                                                          DDATA(alpha),
                                                          DDATA(n),
                                                          DDATA(thetaR),
                                                          DDATA(thetaSR),
                                                          DDATA(KWs),
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

static PyObject* ctransportCoefficientsSeepageBrezis(PyObject* self, 
						     PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double rho,beta,epsFact;
  PyObject *materialTypes,*elementDiameter,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n,*thetaR,*thetaSR;
  if(!PyArg_ParseTuple(args,"OdddOOOOOOOOOOOOOO",
		       &materialTypes,
		       &epsFact,
                       &rho,
                       &beta,
		       &elementDiameter,
                       &gravity,
                       &alpha,
                       &n,
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
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  seepageBrezis(nSimplex,
		nPointsPerSimplex,
		SHAPE(f)[ND(f)-1],
		IDATA(materialTypes),
		epsFact,
		rho,
		beta,
		DDATA(elementDiameter),
		DDATA(gravity),
		DDATA(alpha),
		DDATA(n),
		DDATA(thetaR),
		DDATA(thetaSR),
		DDATA(KWs),
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

static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwind(PyObject* self, 
													 PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1,upwindFlag=0,computeAverages=0;
  double rho,beta;
  PyObject *materialTypes,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n_vg,*thetaR,*thetaSR;
  PyObject *elementBoundaryElementsArray,*quadraturePointToElementBoundary,
    *gradu,*n_global,*dV,*f_avg,*df_avg,*a_avg,*da_avg;
  if(!PyArg_ParseTuple(args,"iiOOOddOOOOOOOOOOOOOOOOOOOO",
		       &upwindFlag,
		       &computeAverages,
		       &elementBoundaryElementsArray,
		       &quadraturePointToElementBoundary,
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &alpha,
                       &n_vg,
                       &thetaR,
                       &thetaSR,
                       &KWs,
                       &u,
		       &gradu,
		       &n_global,
		       &dV,
                       &mass,
                       &dmass,
                       &f_avg,
                       &df_avg,
                       &a_avg,
                       &da_avg,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwind(upwindFlag,
								    computeAverages,
								    nSimplex,
								    nPointsPerSimplex,
								    SHAPE(f)[ND(f)-1],
								    SHAPE(n_global)[1],
								    IDATA(elementBoundaryElementsArray),
								    IDATA(quadraturePointToElementBoundary),
								    IDATA(materialTypes),
								    rho,
								    beta,
								    DDATA(gravity),
								    DDATA(alpha),
								    DDATA(n_vg),
								    DDATA(thetaR),
								    DDATA(thetaSR),
								    DDATA(KWs),
								    DDATA(u),
								    DDATA(gradu),
								    DDATA(n_global),
								    DDATA(dV),
								    DDATA(mass),
								    DDATA(dmass),
								    DDATA(f_avg),
								    DDATA(df_avg),
								    DDATA(a_avg),
								    DDATA(da_avg),
								    DDATA(f),
								    DDATA(df),
								    DDATA(a),
								    DDATA(da));
/*mwf hack try using harmonic average here for conductivity, something wrong with call
   in this example!*/
/*   conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm(upwindFlag, */
/* 									    computeAverages, */
/* 									    nSimplex, */
/* 									    nPointsPerSimplex, */
/* 									    SHAPE(f)[ND(f)-1], */
/* 									    SHAPE(n_global)[1], */
/* 									    IDATA(elementBoundaryElementsArray), */
/* 									    IDATA(quadraturePointToElementBoundary), */
/* 									    IDATA(materialTypes), */
/* 									    rho, */
/* 									    beta, */
/* 									    DDATA(gravity), */
/* 									    DDATA(alpha), */
/* 									    DDATA(n_vg), */
/* 									    DDATA(thetaR), */
/* 									    DDATA(thetaSR), */
/* 									    DDATA(KWs), */
/* 									    DDATA(u), */
/* 									    DDATA(gradu), */
/* 									    DDATA(n_global), */
/* 									    DDATA(dV), */
/* 									    DDATA(mass), */
/* 									    DDATA(dmass), */
/* 									    DDATA(f_avg), */
/* 									    DDATA(df_avg), */
/* 									    DDATA(a_avg), */
/* 									    DDATA(da_avg), */
/* 									    DDATA(f), */
/* 									    DDATA(df), */
/* 									    DDATA(a), */
/* 									    DDATA(da)); */

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm(PyObject* self, 
														PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1,upwindFlag=0,computeAverages=0;
  double rho,beta;
  PyObject *materialTypes,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n_vg,*thetaR,*thetaSR;
  PyObject *elementBoundaryElementsArray,*quadraturePointToElementBoundary,
    *gradu,*n_global,*dV,*f_avg,*df_avg,*a_avg,*da_avg;
  if(!PyArg_ParseTuple(args,"iiOOOddOOOOOOOOOOOOOOOOOOOO",
		       &upwindFlag,
		       &computeAverages,
		       &elementBoundaryElementsArray,
		       &quadraturePointToElementBoundary,
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &alpha,
                       &n_vg,
                       &thetaR,
                       &thetaSR,
                       &KWs,
                       &u,
		       &gradu,
		       &n_global,
		       &dV,
                       &mass,
                       &dmass,
                       &f_avg,
                       &df_avg,
                       &a_avg,
                       &da_avg,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
/*mwf hack try using harmonic average here for conductivity, something wrong with call
   in this example!*/
  conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm(upwindFlag,
									    computeAverages,
									    nSimplex,
									    nPointsPerSimplex,
									    SHAPE(f)[ND(f)-1],
									    SHAPE(n_global)[1],
									    IDATA(elementBoundaryElementsArray),
									    IDATA(quadraturePointToElementBoundary),
									    IDATA(materialTypes),
									    rho,
									    beta,
									    DDATA(gravity),
									    DDATA(alpha),
									    DDATA(n_vg),
									    DDATA(thetaR),
									    DDATA(thetaSR),
									    DDATA(KWs),
									    DDATA(u),
									    DDATA(gradu),
									    DDATA(n_global),
									    DDATA(dV),
									    DDATA(mass),
									    DDATA(dmass),
									    DDATA(f_avg),
									    DDATA(df_avg),
									    DDATA(a_avg),
									    DDATA(da_avg),
									    DDATA(f),
									    DDATA(df),
									    DDATA(a),
									    DDATA(da));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm_sd(PyObject* self, 
														   PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1,upwindFlag=0,computeAverages=0;
  double rho,beta;
  PyObject *materialTypes,*gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n_vg,*thetaR,*thetaSR;
  PyObject *elementBoundaryElementsArray,*quadraturePointToElementBoundary,
    *gradu,*n_global,*dV,*f_avg,*df_avg,*a_avg,*da_avg,*rowptr,*colind;
  
  if(!PyArg_ParseTuple(args,"iiOOOOOddOOOOOOOOOOOOOOOOOOOO",
		       &upwindFlag,
		       &computeAverages,
		       &rowptr,
		       &colind,
		       &elementBoundaryElementsArray,
		       &quadraturePointToElementBoundary,
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &alpha,
                       &n_vg,
                       &thetaR,
                       &thetaSR,
                       &KWs,
                       &u,
		       &gradu,
		       &n_global,
		       &dV,
                       &mass,
                       &dmass,
                       &f_avg,
                       &df_avg,
                       &a_avg,
                       &da_avg,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
/*mwf hack try using harmonic average here for conductivity, something wrong with call
   in this example!*/
  conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm_sd(upwindFlag,
									      computeAverages,
									      nSimplex,
									      nPointsPerSimplex,
									      SHAPE(f)[ND(f)-1],
									      SHAPE(n_global)[1],
									      IDATA(rowptr),
									      IDATA(colind),
									      IDATA(elementBoundaryElementsArray),
									      IDATA(quadraturePointToElementBoundary),
									      IDATA(materialTypes),
									      rho,
									      beta,
									      DDATA(gravity),
									      DDATA(alpha),
									      DDATA(n_vg),
									      DDATA(thetaR),
									      DDATA(thetaSR),
									      DDATA(KWs),
									      DDATA(u),
									      DDATA(gradu),
									      DDATA(n_global),
									      DDATA(dV),
									      DDATA(mass),
									      DDATA(dmass),
									      DDATA(f_avg),
									      DDATA(df_avg),
									      DDATA(a_avg),
									      DDATA(da_avg),
									      DDATA(f),
									      DDATA(df),
									      DDATA(a),
									      DDATA(da));
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsConservativeHeadRichardsJLeverettEvaluate(PyObject* self, 
                                                                         PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double rho,beta;
  PyObject *materialTypes,*gravity,
    *u,*mass,*dmass,*a,*da,*f,*df,
    *phi,*psiD,*ns,*nk,*S_wirr,*S_nwr,*kr0;
  if(!PyArg_ParseTuple(args,"OddOOOOOOOOOOOOOOO",
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &phi,
                       &psiD,
                       &ns,
                       &nk,
                       &S_wirr,
                       &S_nwr,
                       &kr0,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  conservativeHeadRichardsJLeverett(nSimplex,
                                    nPointsPerSimplex,
                                    SHAPE(f)[ND(f)-1],
                                    IDATA(materialTypes),
                                    rho,
                                    beta,
                                    DDATA(gravity),
                                    DDATA(phi),
                                    DDATA(psiD),
                                    DDATA(ns),
                                    DDATA(nk),
                                    DDATA(S_wirr),
                                    DDATA(S_nwr),
                                    DDATA(kr0),
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
static PyObject* ctransportCoefficientsConservativeHeadRichardsJLeverettAniEvaluate(PyObject* self, 
                                                                         PyObject* args)
{
  int i,nSimplex=1,nPointsPerSimplex=1;
  double rho,beta;
  PyObject *materialTypes,*gravity,
    *u,*mass,*dmass,*a,*da,*f,*df,
    *phi,*psiD,*ns,*nk,*S_wirr,*S_nwr,*kr0x,*kr0y,*kr0z;
  if(!PyArg_ParseTuple(args,"OddOOOOOOOOOOOOOOOOO",
		       &materialTypes,
                       &rho,
                       &beta,
                       &gravity,
                       &phi,
                       &psiD,
                       &ns,
                       &nk,
                       &S_wirr,
                       &S_nwr,
                       &kr0x,
                       &kr0y,
                       &kr0z,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0; i < ND(f)-2; i++)
    nSimplex *= SHAPE(f)[i];
  for(i=ND(f)-2;i<ND(f)-1;i++)
      nPointsPerSimplex *= SHAPE(f)[i];
  conservativeHeadRichardsJLeverettAni(nSimplex,
                                    nPointsPerSimplex,
                                    SHAPE(f)[ND(f)-1],
                                    IDATA(materialTypes),
                                    rho,
                                    beta,
                                    DDATA(gravity),
                                    DDATA(phi),
                                    DDATA(psiD),
                                    DDATA(ns),
                                    DDATA(nk),
                                    DDATA(S_wirr),
                                    DDATA(S_nwr),
                                    DDATA(kr0x),
                                    DDATA(kr0y),
                                    DDATA(kr0z),
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


/*mwf end more additions */
static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluate(PyObject* self, 
                                                                                         PyObject* args)
{
  int i,nPoints=1;
  double rho;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*KWs,
    *alpha,*n,*thetaR,*thetaSR;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOOO",
                       &rho,
                       &gravity,
                       &alpha,
                       &n,
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeHeadRichardsMualemVanGenuchtenHetEvaluate(nPoints,
                                                        SHAPE(f)[ND(f)-1],
                                                        rho,
                                                        DDATA(gravity),
                                                        DDATA(alpha),
                                                        DDATA(n),
                                                        DDATA(thetaR),
                                                        DDATA(thetaSR),
                                                        DDATA(KWs),
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

/* static PyObject* ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(PyObject* self,  */
/* 											       PyObject* args) */
/* { */
/*   int i,nPoints=1,nElements=1,nFacePerElement=0; */
/*   double rho; */
/*   PyObject *gravity,*materialTypes,*u,*mass,*dmass,*a,*da,*f,*df,*KWs, */
/*     *alpha,*n,*thetaR,*thetaSR; */
/*   if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOOOO", */
/*                        &rho, */
/*                        &gravity, */
/* 		       &materialTypes, */
/*                        &alpha, */
/*                        &n, */
/*                        &thetaR, */
/*                        &thetaSR, */
/*                        &KWs, */
/*                        &u, */
/*                        &mass, */
/*                        &dmass, */
/*                        &f, */
/*                        &df, */
/*                        &a, */
/*                        &da)) */
/*     return NULL; */
/*   for(i=0;i<ND(f)-1;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/*   nElements = SHAPE(u)[0]; */
/*   if (ND(u) > 2) */
/*     nFacePerElement = SHAPE(u)[1]; */
/*   conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(nElements, */
/* 							  nFacePerElement, */
/* 							  nPoints, */
/* 							  SHAPE(f)[ND(f)-1], */
/* 							  rho, */
/* 							  DDATA(gravity), */
/* 							  IDATA(materialTypes), */
/* 							  DDATA(alpha), */
/* 							  DDATA(n), */
/* 							  DDATA(thetaR), */
/* 							  DDATA(thetaSR), */
/* 							  DDATA(KWs), */
/* 							  DDATA(u), */
/* 							  DDATA(mass), */
/* 							  DDATA(dmass), */
/* 							  DDATA(f), */
/* 							  DDATA(df), */
/* 							  DDATA(a), */
/* 							  DDATA(da)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

static PyObject* ctransportCoefficientsConservativeSatRichardsMualemVanGenuchtenHomEvaluate(PyObject* self, 
                                                                                           PyObject* args)
{
  int i,nPoints=1;
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*x,*u,*mass,*dmass,*a,*da,*f,*df,*phi,*dphi;
  if(!PyArg_ParseTuple(args,"dOOddddddOOOOOOOOO",
                       &rho,
                       &gravity,
                       &x,
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeSatRichardsMualemVanGenuchtenHomEvaluate(nPoints,
                                                       SHAPE(f)[ND(f)-1],
                                                       rho,
                                                       DDATA(gravity),
                                                       DDATA(x),
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

static PyObject* ctransportCoefficientsConservativeHeadRichardsBrooksCoreyBurdineHomEvaluate(PyObject* self, 
											     PyObject* args)
{
  int i,nPoints=1;
  double rho,beta,lambda,pd,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df;
  if(!PyArg_ParseTuple(args,"ddOdddddOOOOOOO",
                       &rho,
		       &beta,
                       &gravity,
                       &lambda,
                       &pd,
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeHeadRichardsBrooksCoreyBurdineHomEvaluate(nPoints,
                                                        SHAPE(f)[ND(f)-1],
                                                        rho,
							beta,
                                                        DDATA(gravity),
                                                        lambda,
                                                        pd,
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

static PyObject* ctransportCoefficientsConservativeHeadRichardsBrooksCoreyBurdineHetEvaluate(PyObject* self, 
											     PyObject* args)
{
  int i,nPoints=1;
  double rho;
  PyObject *lambda, *pd, *thetaR, *thetaS, *KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOOO",
                       &rho,
                       &gravity,
                       &lambda,
                       &pd,
                       &thetaR,
                       &thetaS,
                       &KWs,
                       &u,
                       &mass,
                       &dmass,
                       &f,
                       &df,
                       &a,
                       &da))
    return NULL;
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeHeadRichardsBrooksCoreyBurdineHetEvaluate(nPoints,
                                                        SHAPE(f)[ND(f)-1],
                                                        rho,
                                                        DDATA(gravity),
                                                        DDATA(lambda),
                                                        DDATA(pd),
                                                        DDATA(thetaR),
                                                        DDATA(thetaS),
                                                        DDATA(KWs),
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

static PyObject* ctransportCoefficientsConservativeSatRichardsBrooksCoreyBurdineHomEvaluate(PyObject* self, 
											    PyObject* args)
{
  int nPoints,nSpace;
  double rho,lambda,pd,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df,*phi,*dphi;
  if(!PyArg_ParseTuple(args,"iidOdddddOOOOOOOOO",
                       &nPoints,
                       &nSpace,
                       &rho,
                       &gravity,
                       &lambda,
                       &pd,
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
  conservativeSatRichardsBrooksCoreyBurdineHomEvaluate(nPoints,
                                                       nSpace,
                                                       rho,
                                                       DDATA(gravity),
                                                       lambda,
                                                       pd,
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

static PyObject* ctransportCoefficientsConservativeHeadRichardsBCBfromMVGHomEvaluate(PyObject* self, 
										     PyObject* args)
{
  int i,nPoints=1;
  double rho,alpha,n,m,thetaR,thetaSR,KWs;
  PyObject *gravity,*u,*mass,*dmass,*a,*da,*f,*df;
  if(!PyArg_ParseTuple(args,"dOddddddOOOOOOO",
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
  for(i=0;i<ND(f)-1;i++)
      nPoints *= SHAPE(f)[i];
  conservativeHeadRichardsBCBfromMVGHomEvaluate(nPoints,
						SHAPE(f)[ND(f)-1],
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

static PyObject* ctransportCoefficientsMass_2D_Evaluate(PyObject* self,
							PyObject* args)
{
  int i,nPoints=1;
  double rho;
  PyObject *p,
    *u,
    *v,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v;

  if(!PyArg_ParseTuple(args,"dOOOOOOOOO",
		       &rho,
		       &p,
		       &u,
		       &v,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v))
    
  return NULL;
  for(i=0;i<ND(p);i++)
    nPoints*=SHAPE(p)[i];

  Mass_2D_Evaluate(nPoints,
		   rho,
		   DDATA(p),
		   DDATA(u),
		   DDATA(v),
		   DDATA(mom_p_acc),
		   DDATA(mom_u_acc),
		   DDATA(mom_v_acc),
		   DDATA(dmom_p_acc_p),
		   DDATA(dmom_u_acc_u),
		   DDATA(dmom_v_acc_v));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsMass_3D_Evaluate(PyObject* self,
							PyObject* args)
{
  int i,nPoints=1;
  double rho;
  PyObject *p,
    *u,
    *v,
    *w,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *mom_w_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v,
    *dmom_w_acc_w;

  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOO",
		       &rho,
		       &p,
		       &u,
		       &v,
		       &w,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &mom_w_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v,
		       &dmom_w_acc_w))
    
  return NULL;
  for(i=0;i<ND(p);i++)
    nPoints*=SHAPE(p)[i];

  Mass_3D_Evaluate(nPoints,
		   rho,
		   DDATA(p),
		   DDATA(u),
		   DDATA(v),
		   DDATA(w),
		   DDATA(mom_p_acc),
		   DDATA(mom_u_acc),
		   DDATA(mom_v_acc),
		   DDATA(mom_w_acc),
		   DDATA(dmom_p_acc_p),
		   DDATA(dmom_u_acc_u),
		   DDATA(dmom_v_acc_v),
		   DDATA(dmom_w_acc_w));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseMass_2D_Evaluate(PyObject* self,
								PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseMass_2D_Evaluate(nPoints,
			   eps,
			   rho_0,
			   nu_0,
			   rho_1,
			   nu_1,
			   DDATA(phi),
			   DDATA(p),
			   DDATA(u),
			   DDATA(v),
			   DDATA(mom_p_acc),
			   DDATA(mom_u_acc),
			   DDATA(mom_v_acc),
			   DDATA(dmom_p_acc_p),
			   DDATA(dmom_u_acc_u),
			   DDATA(dmom_v_acc_v));
			   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseMass_3D_Evaluate(PyObject* self,
								PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *w,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *mom_w_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v,
    *dmom_w_acc_w;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &w,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &mom_w_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v,
		       &dmom_w_acc_w))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseMass_3D_Evaluate(nPoints,
			   eps,
			   rho_0,
			   nu_0,
			   rho_1,
			   nu_1,
			   DDATA(phi),
			   DDATA(p),
			   DDATA(u),
			   DDATA(v),
			   DDATA(w),
			   DDATA(mom_p_acc),
			   DDATA(mom_u_acc),
			   DDATA(mom_v_acc),
			   DDATA(mom_w_acc),
			   DDATA(dmom_p_acc_p),
			   DDATA(dmom_u_acc_u),
			   DDATA(dmom_v_acc_v),
			   DDATA(dmom_w_acc_w));
			   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseInvScaledMass_2D_Evaluate(PyObject* self,
									 PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseInvScaledMass_2D_Evaluate(nPoints,
				    eps,
				    rho_0,
				    nu_0,
				    rho_1,
				    nu_1,
				    DDATA(phi),
				    DDATA(p),
				    DDATA(u),
				    DDATA(v),
				    DDATA(mom_p_acc),
				    DDATA(mom_u_acc),
				    DDATA(mom_v_acc),
				    DDATA(dmom_p_acc_p),
				    DDATA(dmom_u_acc_u),
				    DDATA(dmom_v_acc_v));
			   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseMass_mu_2D_Evaluate(PyObject* self,
								   PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseMass_mu_2D_Evaluate(nPoints,
			      eps,
			      rho_0,
			      nu_0,
			      rho_1,
			      nu_1,
			      DDATA(phi),
			      DDATA(p),
			      DDATA(u),
			      DDATA(v),
			      DDATA(mom_p_acc),
			      DDATA(mom_u_acc),
			      DDATA(mom_v_acc),
			      DDATA(dmom_p_acc_p),
			      DDATA(dmom_u_acc_u),
			      DDATA(dmom_v_acc_v));
			   
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* ctransportCoefficientsTwoPhaseInvScaledMass_3D_Evaluate(PyObject* self,
									 PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *w,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *mom_w_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v,
    *dmom_w_acc_w;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &w,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &mom_w_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v,
		       &dmom_w_acc_w))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseInvScaledMass_3D_Evaluate(nPoints,
				    eps,
				    rho_0,
				    nu_0,
				    rho_1,
				    nu_1,
				    DDATA(phi),
				    DDATA(p),
				    DDATA(u),
				    DDATA(v),
				    DDATA(w),
				    DDATA(mom_p_acc),
				    DDATA(mom_u_acc),
				    DDATA(mom_v_acc),
				    DDATA(mom_w_acc),
				    DDATA(dmom_p_acc_p),
				    DDATA(dmom_u_acc_u),
				    DDATA(dmom_v_acc_v),
				    DDATA(dmom_w_acc_w));
			   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseMass_mu_3D_Evaluate(PyObject* self,
								   PyObject* args)
{
  int i, nPoints = 1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *w,
    *mom_p_acc,
    *mom_u_acc,
    *mom_v_acc,
    *mom_w_acc,
    *dmom_p_acc_p,
    *dmom_u_acc_u,
    *dmom_v_acc_v,
    *dmom_w_acc_w;

  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOO",
		       &eps,
		       &rho_0,
		       &nu_0,
		       &rho_1,
		       &nu_1,
		       &phi,
		       &p,
		       &u,
		       &v,
		       &w,
		       &mom_p_acc,
		       &mom_u_acc,
		       &mom_v_acc,
		       &mom_w_acc,
		       &dmom_p_acc_p,
		       &dmom_u_acc_u,
		       &dmom_v_acc_v,
		       &dmom_w_acc_w))

    return NULL;

  for (i=0 ; i<ND(p) ; i++)
    nPoints*=SHAPE(p)[i];

  TwoPhaseMass_mu_3D_Evaluate(nPoints,
			      eps,
			      rho_0,
			      nu_0,
			      rho_1,
			      nu_1,
			      DDATA(phi),
			      DDATA(p),
			      DDATA(u),
			      DDATA(v),
			      DDATA(w),
			      DDATA(mom_p_acc),
			      DDATA(mom_u_acc),
			      DDATA(mom_v_acc),
			      DDATA(mom_w_acc),
			      DDATA(dmom_p_acc_p),
			      DDATA(dmom_u_acc_u),
			      DDATA(dmom_v_acc_v),
			      DDATA(dmom_w_acc_w));
			   
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* ctransportCoefficientsLaplace_2D_Evaluate(PyObject* self,
							PyObject* args)
{
  // Need to double check I'm not adding too many attributes here.
  int i,nPoints=1;
  PyObject *u,*v,*p,*mom_u_diff_ten,*mom_v_diff_ten,*mom_p_diff_ten ;

  if (!PyArg_ParseTuple(args,"OOOOOO",
			&p,
			&u,
			&v,
			&mom_p_diff_ten,
			&mom_u_diff_ten,
			&mom_v_diff_ten))
  return NULL;

  for (i=0;i<ND(p);i++){
    nPoints *= SHAPE(p)[i];
  }

  Laplace_2D_Evaluate(nPoints,
		      DDATA(mom_p_diff_ten),
		      DDATA(mom_u_diff_ten),
		      DDATA(mom_v_diff_ten));

  Py_INCREF(Py_None);
  return Py_None;	
}

static PyObject* ctransportCoefficientsLaplace_3D_Evaluate(PyObject* self,
							   PyObject* args)
{
  int i,nPoints=1;
  PyObject *u,*v,*w,*p,*mom_u_diff_ten,*mom_v_diff_ten,
           *mom_w_diff_ten, *mom_p_diff_ten ;

  if (!PyArg_ParseTuple(args,"OOOOOOOO",
			&p,
			&u,
			&v,
			&w,
			&mom_p_diff_ten,
			&mom_u_diff_ten,
			&mom_v_diff_ten,
			&mom_w_diff_ten))
  return NULL;

  for (i=0;i<ND(p);i++){
    nPoints *= SHAPE(p)[i];
  }

  Laplace_3D_Evaluate(nPoints,
		      DDATA(mom_p_diff_ten),
		      DDATA(mom_u_diff_ten),
		      DDATA(mom_v_diff_ten),
		      DDATA(mom_w_diff_ten));

  Py_INCREF(Py_None);
  return Py_None;	
}

static PyObject* ctransportCoefficientsTwoPhaseInvScaledLaplace_2D_Evaluate(PyObject* self,
									    PyObject* args)
{
  int i,nPoints=1;
  double eps, rho_0, nu_0, rho_1, nu_1 ;
  PyObject *phi, *u,*v,*p,*mom_u_diff_ten,*mom_v_diff_ten,*mom_p_diff_ten ;

  if (!PyArg_ParseTuple(args,"dddddOOOOOOO",
			&eps,
			&rho_0,
			&nu_0,
			&rho_1,
			&nu_1,
			&phi,
			&p,
			&u,
			&v,
			&mom_p_diff_ten,
			&mom_u_diff_ten,
			&mom_v_diff_ten))
  return NULL;

  for (i=0;i<ND(p);i++){
    nPoints *= SHAPE(p)[i];
  }

  TwoPhaseInvScaledLaplace_2D_Evaluate(nPoints,
				       eps,
				       rho_0,
				       nu_0,
				       rho_1,
				       nu_1,
				       DDATA(phi),
				       DDATA(mom_p_diff_ten),
				       DDATA(mom_u_diff_ten),
				       DDATA(mom_v_diff_ten));

  Py_INCREF(Py_None);
  return Py_None;	
}

static PyObject* ctransportCoefficientsTwoPhaseInvScaledLaplace_3D_Evaluate(PyObject* self,
									    PyObject* args)
{
  int i,nPoints=1;
  double eps, rho_0, nu_0, rho_1, nu_1 ;
  PyObject *phi,*u,*v,*w,*p,*mom_u_diff_ten,*mom_v_diff_ten,
           *mom_w_diff_ten, *mom_p_diff_ten ;

  if (!PyArg_ParseTuple(args,"dddddOOOOOOOO",
			&eps,
			&rho_0,
			&nu_0,
			&rho_1,
			&nu_1,
			&phi,
			&p,
			&u,
			&v,
			&w,
			&mom_p_diff_ten,
			&mom_u_diff_ten,
			&mom_v_diff_ten,
			&mom_w_diff_ten))
  return NULL;

  for (i=0;i<ND(p);i++){
    nPoints *= SHAPE(p)[i];
  }

  TwoPhaseInvScaledLaplace_3D_Evaluate(nPoints,
				       eps,
				       rho_0,
				       nu_0,
				       rho_1,
				       nu_1,
				       DDATA(phi),
				       DDATA(mom_p_diff_ten),
				       DDATA(mom_u_diff_ten),
				       DDATA(mom_v_diff_ten),
				       DDATA(mom_w_diff_ten));

  Py_INCREF(Py_None);
  return Py_None;	
}

static PyObject* ctransportCoefficientsAdvection_2D_Evaluate(PyObject* self,
							     PyObject* args)
{
  int i,nPoints=1;
  PyObject *p,
    *u,
    *v,
    *mass_adv,
    *dmass_adv_p,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v;
  
  if (!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
			&p,
			&u,
			&v,
			&mass_adv,
			&mom_u_adv,
			&mom_v_adv,
			&dmass_adv_p,
			&dmass_adv_u,
			&dmass_adv_v,
			&dmom_u_adv_u,
			&dmom_u_adv_v,
			&dmom_v_adv_u,
			&dmom_v_adv_v))
  return NULL;

  for (i=0;i<ND(p);i++)
    nPoints *= SHAPE(p)[i];


  Advection_2D_Evaluate(nPoints,
			DDATA(p),
			DDATA(u),
			DDATA(v),
			DDATA(mass_adv),
			DDATA(dmass_adv_p),
			DDATA(dmass_adv_u),
			DDATA(dmass_adv_v),
			DDATA(mom_u_adv),
			DDATA(dmom_u_adv_u),
			DDATA(dmom_u_adv_v),
			DDATA(mom_v_adv),
			DDATA(dmom_v_adv_u),
			DDATA(dmom_v_adv_v));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsAdvection_3D_Evaluate(PyObject* self,
							     PyObject* args)
{
  int i,nPoints=1;
  PyObject *p,
    *u,
    *v,
    *w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w;
  if (!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOO",
			&mass_adv,
			&mom_u_adv,
			&mom_v_adv,
			&mom_w_adv,
			&dmass_adv_u,
			&dmass_adv_v,
			&dmass_adv_w,
			&dmom_u_adv_u,
			&dmom_u_adv_v,
			&dmom_u_adv_w,
			&dmom_v_adv_u,
			&dmom_v_adv_v,
			&dmom_v_adv_w,
			&dmom_w_adv_u,
			&dmom_w_adv_v,
			&dmom_w_adv_w))
  return NULL;
  for (i=0;i<ND(p);i++)
    nPoints *= SHAPE(p)[i];

  Advection_3D_Evaluate(nPoints,
			DDATA(p),
			DDATA(u),
			DDATA(v),
			DDATA(w),
			DDATA(mass_adv),
			DDATA(dmass_adv_u),
			DDATA(dmass_adv_v),
			DDATA(dmass_adv_w),
			DDATA(mom_u_adv),
			DDATA(dmom_u_adv_u),
			DDATA(dmom_u_adv_v),
			DDATA(dmom_u_adv_w),
			DDATA(mom_v_adv),
			DDATA(dmom_v_adv_u),
			DDATA(dmom_v_adv_v),
			DDATA(dmom_v_adv_w),
			DDATA(mom_w_adv),
			DDATA(dmom_w_adv_u),
			DDATA(dmom_w_adv_v),
			DDATA(dmom_w_adv_w));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseAdvection_2D_Evaluate(PyObject* self,
								     PyObject* args)
{
  int i,nPoints=1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *mass_adv,
    *dmass_adv_p,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v;
  
  if (!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOO",
			&eps,
			&rho_0,
			&nu_0,
			&rho_1,
			&nu_1,
			&phi,
			&p,
			&u,
			&v,
			&mass_adv,
			&mom_u_adv,
			&mom_v_adv,
			&dmass_adv_p,
			&dmass_adv_u,
			&dmass_adv_v,
			&dmom_u_adv_u,
			&dmom_u_adv_v,
			&dmom_v_adv_u,
			&dmom_v_adv_v))
  return NULL;

  for (i=0;i<ND(p);i++)
    nPoints *= SHAPE(p)[i];


  TwoPhaseAdvection_2D_Evaluate(nPoints,
				eps,
				rho_0,
				nu_0,
				rho_1,
				nu_1,
				DDATA(phi),
				DDATA(p),
				DDATA(u),
				DDATA(v),
				DDATA(mass_adv),
				DDATA(dmass_adv_p),
				DDATA(dmass_adv_u),
				DDATA(dmass_adv_v),
				DDATA(mom_u_adv),
				DDATA(dmom_u_adv_u),
				DDATA(dmom_u_adv_v),
				DDATA(mom_v_adv),
				DDATA(dmom_v_adv_u),
				DDATA(dmom_v_adv_v));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* ctransportCoefficientsTwoPhaseAdvection_3D_Evaluate(PyObject* self,
								     PyObject* args)
{
  int i,nPoints=1;
  double eps, rho_0, nu_0, rho_1, nu_1;
  PyObject *phi,
    *p,
    *u,
    *v,
    *w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w;
  
  if (!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOO",
			&eps,
			&rho_0,
			&nu_0,
			&rho_1,
			&nu_1,
			&phi,
			&mass_adv,
			&mom_u_adv,
			&mom_v_adv,
			&mom_w_adv,
			&dmass_adv_u,
			&dmass_adv_v,
			&dmass_adv_w,
			&dmom_u_adv_u,
			&dmom_u_adv_v,
			&dmom_u_adv_w,
			&dmom_v_adv_u,
			&dmom_v_adv_v,
			&dmom_v_adv_w,
			&dmom_w_adv_u,
			&dmom_w_adv_v,
			&dmom_w_adv_w))
  return NULL;
  for (i=0;i<ND(p);i++)
    nPoints *= SHAPE(p)[i];

  TwoPhaseAdvection_3D_Evaluate(nPoints,
				eps,
				rho_0,
				nu_0,
				rho_1,
				nu_1,
				DDATA(phi),
				DDATA(p),
				DDATA(u),
				DDATA(v),
				DDATA(w),
				DDATA(mass_adv),
				DDATA(dmass_adv_u),
				DDATA(dmass_adv_v),
				DDATA(dmass_adv_w),
				DDATA(mom_u_adv),
				DDATA(dmom_u_adv_u),
				DDATA(dmom_u_adv_v),
				DDATA(dmom_u_adv_w),
				DDATA(mom_v_adv),
				DDATA(dmom_v_adv_u),
				DDATA(dmom_v_adv_v),
				DDATA(dmom_v_adv_w),
				DDATA(mom_w_adv),
				DDATA(dmom_w_adv_u),
				DDATA(dmom_w_adv_v),
				DDATA(dmom_w_adv_w));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* ctransportCoefficientsB_2D_Evaluate(PyObject* self,
						     PyObject* args)
{
  int i,nPoints = 1;
  PyObject *p,
    *grad_p,
    *u,
    *v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_ham,
    *mom_v_ham,
    *dmom_u_ham_grad_p,
    *dmom_v_ham_grad_p;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
		       &p,
		       &grad_p,
		       &u,
		       &v,
		       &mass_adv,
		       &dmass_adv_u,
		       &dmass_adv_v,
		       &mom_u_ham,
		       &mom_v_ham,
		       &dmom_u_ham_grad_p,
		       &dmom_v_ham_grad_p))
    return NULL;
  
  for(i=0; i<ND(p); i++)
    nPoints *= SHAPE(p)[i];

  B_2D_Evaluate(nPoints,
  		DDATA(grad_p),
  		DDATA(u),
  		DDATA(v),
  		DDATA(mass_adv),
  		DDATA(dmass_adv_u),
  		DDATA(dmass_adv_v),
  		DDATA(mom_u_ham),
  		DDATA(mom_v_ham),
  		DDATA(dmom_u_ham_grad_p),
  		DDATA(dmom_v_ham_grad_p));

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* ctransportCoefficientsB_3D_Evaluate(PyObject* self,
						     PyObject* args)
{
  int i,nPoints = 1;
  PyObject *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_ham,
    *mom_v_ham,
    *mom_w_ham,
    *dmom_u_ham_grad_p,
    *dmom_v_ham_grad_p,
    *dmom_w_ham_grad_p;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
		       &p,
		       &grad_p,
		       &u,
		       &v,
		       &w,
		       &mass_adv,
		       &dmass_adv_u,
		       &dmass_adv_v,
		       &dmass_adv_w,
		       &mom_u_ham,
		       &mom_v_ham,
		       &mom_w_ham,
		       &dmom_u_ham_grad_p,
		       &dmom_v_ham_grad_p,
		       &dmom_w_ham_grad_p))
    return NULL;
  
  for(i=0; i<ND(p); i++)
    nPoints *= SHAPE(p)[i];
  B_3D_Evaluate(nPoints,
		DDATA(grad_p),
		DDATA(u),
		DDATA(v),
		DDATA(w),
		DDATA(mass_adv),
		DDATA(dmass_adv_u),
		DDATA(dmass_adv_v),
		DDATA(dmass_adv_w),
		DDATA(mom_u_ham),
		DDATA(mom_v_ham),
		DDATA(mom_w_ham),
		DDATA(dmom_u_ham_grad_p),
		DDATA(dmom_v_ham_grad_p),
		DDATA(dmom_w_ham_grad_p));
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* ctransportCoefficientsNavierStokes_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  NavierStokes_2D_Evaluate(nPoints,
                           rho,
                           nu,
                           DDATA(g),
                           DDATA(p),
                           DDATA(grad_p),
                           DDATA(u),
                           DDATA(v),
                           DDATA(mom_u_acc),
                           DDATA(dmom_u_acc_u),
                           DDATA(mom_v_acc),
                           DDATA(dmom_v_acc_v),
                           DDATA(mass_adv),
                           DDATA(dmass_adv_u),
                           DDATA(dmass_adv_v),
                           DDATA(mom_u_adv),
                           DDATA(dmom_u_adv_u),
                           DDATA(dmom_u_adv_v),
                           DDATA(mom_v_adv),
                           DDATA(dmom_v_adv_u),
                           DDATA(dmom_v_adv_v),
                           DDATA(mom_u_diff_ten),
                           DDATA(mom_v_diff_ten),
                           DDATA(mom_u_source),
                           DDATA(mom_v_source),
                           DDATA(mom_u_ham),
                           DDATA(dmom_u_ham_grad_p),
                           DDATA(mom_v_ham),
                           DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* ctransportCoefficientsNavierStokes_3D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  NavierStokes_3D_Evaluate(nPoints,
                           rho,
                           nu,
                           DDATA(g),
                           DDATA(p),
                           DDATA(grad_p),
                           DDATA(u),
                           DDATA(v),
                           DDATA(w),
                           DDATA(mom_u_acc),
                           DDATA(dmom_u_acc_u),
                           DDATA(mom_v_acc),
                           DDATA(dmom_v_acc_v),
                           DDATA(mom_w_acc),
                           DDATA(dmom_w_acc_w),
                           DDATA(mass_adv),
                           DDATA(dmass_adv_u),
                           DDATA(dmass_adv_v),
                           DDATA(dmass_adv_w),
                           DDATA(mom_u_adv),
                           DDATA(dmom_u_adv_u),
                           DDATA(dmom_u_adv_v),
                           DDATA(dmom_u_adv_w),
                           DDATA(mom_v_adv),
                           DDATA(dmom_v_adv_u),
                           DDATA(dmom_v_adv_v),
                           DDATA(dmom_v_adv_w),
                           DDATA(mom_w_adv),
                           DDATA(dmom_w_adv_u),
                           DDATA(dmom_w_adv_v),
                           DDATA(dmom_w_adv_w),
                           DDATA(mom_u_diff_ten),
                           DDATA(mom_v_diff_ten),
                           DDATA(mom_w_diff_ten),
                           DDATA(mom_u_source),
                           DDATA(mom_v_source),
                           DDATA(mom_w_source),
                           DDATA(mom_u_ham),
                           DDATA(dmom_u_ham_grad_p),
                           DDATA(mom_v_ham),
                           DDATA(dmom_v_ham_grad_p),
                           DDATA(mom_w_ham),
                           DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsStokes_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  Stokes_2D_Evaluate(nPoints,
                     rho,
                     nu,
                     DDATA(g),
                     DDATA(p),
                     DDATA(grad_p),
                     DDATA(u),
                     DDATA(v),
                     DDATA(mom_u_acc),
                     DDATA(dmom_u_acc_u),
                     DDATA(mom_v_acc),
                     DDATA(dmom_v_acc_v),
                     DDATA(mass_adv),
                     DDATA(dmass_adv_u),
                     DDATA(dmass_adv_v),
                     DDATA(mom_u_diff_ten),
                     DDATA(mom_v_diff_ten),
                     DDATA(mom_u_source),
                     DDATA(mom_v_source),
                     DDATA(mom_u_ham),
                     DDATA(dmom_u_ham_grad_p),
                     DDATA(mom_v_ham),
                     DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsStokesP_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_p,
    *mom_v_adv,
    *dmom_v_adv_p,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_p,
                       &mom_v_adv,
                       &dmom_v_adv_p,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  StokesP_2D_Evaluate(nPoints,
                     rho,
                     nu,
                     DDATA(g),
                     DDATA(p),
                     DDATA(u),
                     DDATA(v),
                     DDATA(mom_u_acc),
                     DDATA(dmom_u_acc_u),
                     DDATA(mom_v_acc),
                     DDATA(dmom_v_acc_v),
                     DDATA(mass_adv),
                     DDATA(dmass_adv_u),
                     DDATA(dmass_adv_v),
                     DDATA(mom_u_adv),
                     DDATA(dmom_u_adv_p),
                     DDATA(mom_v_adv),
                     DDATA(dmom_v_adv_p),
                     DDATA(mom_u_diff_ten),
                     DDATA(mom_v_diff_ten),
                     DDATA(mom_u_source),
                     DDATA(mom_v_source));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsStokesP_3D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_p,
    *mom_v_adv,
    *dmom_v_adv_p,
    *mom_w_adv,
    *dmom_w_adv_p,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &u,
                       &v,
		       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_p,
                       &mom_v_adv,
                       &dmom_v_adv_p,
                       &mom_w_adv,
                       &dmom_w_adv_p,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
		       &mom_w_source))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
/*   StokesP_3D_Evaluate(nPoints, */
/* 		      rho, */
/* 		      nu, */
/* 		      DDATA(g), */
/* 		      DDATA(p), */
/* 		      DDATA(u), */
/* 		      DDATA(v), */
/* 		      DDATA(w), */
/* 		      DDATA(mom_u_acc), */
/* 		      DDATA(dmom_u_acc_u), */
/* 		      DDATA(mom_v_acc), */
/* 		      DDATA(dmom_v_acc_v), */
/* 		      DDATA(mom_w_acc), */
/* 		      DDATA(dmom_w_acc_w), */
/* 		      DDATA(mass_adv), */
/* 		      DDATA(dmass_adv_u), */
/* 		      DDATA(dmass_adv_v), */
/* 		      DDATA(dmass_adv_w), */
/* 		      DDATA(mom_u_adv), */
/* 		      DDATA(dmom_u_adv_p), */
/* 		      DDATA(mom_v_adv), */
/* 		      DDATA(dmom_v_adv_p), */
/* 		      DDATA(mom_w_adv), */
/* 		      DDATA(dmom_w_adv_p), */
/* 		      DDATA(mom_u_diff_ten), */
/* 		      DDATA(mom_v_diff_ten), */
/* 		      DDATA(mom_w_diff_ten), */
/* 		      DDATA(mom_u_source), */
/* 		      DDATA(mom_v_source), */
/* 		      DDATA(mom_w_source)); */
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsStokes_3D_Evaluate(PyObject* self, 
                                                          PyObject* args)
{
  int i,nPoints=1;
  double rho,nu;
  PyObject *g,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &nu,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  Stokes_3D_Evaluate(nPoints,
                     rho,
                     nu,
                     DDATA(g),
                     DDATA(p),
                     DDATA(grad_p),
                     DDATA(u),
                     DDATA(v),
                     DDATA(w),
                     DDATA(mom_u_acc),
                     DDATA(dmom_u_acc_u),
                     DDATA(mom_v_acc),
                     DDATA(dmom_v_acc_v),
                     DDATA(mom_w_acc),
                     DDATA(dmom_w_acc_w),
                     DDATA(mass_adv),
                     DDATA(dmass_adv_u),
                     DDATA(dmass_adv_v),
                     DDATA(dmass_adv_w),
                     DDATA(mom_u_diff_ten),
                     DDATA(mom_v_diff_ten),
                     DDATA(mom_w_diff_ten),
                     DDATA(mom_u_source),
                     DDATA(mom_v_source),
                     DDATA(mom_w_source),
                     DDATA(mom_u_ham),
                     DDATA(dmom_u_ham_grad_p),
                     DDATA(mom_v_ham),
                     DDATA(dmom_v_ham_grad_p),
                     DDATA(mom_w_ham),
                     DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_LS_SO_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_LS_SO_2D_Evaluate(nPoints,
                                 eps,
                                 rho0,
                                 nu0,
                                 rho1,
                                 nu1,
                                 DDATA(g),
                                 DDATA(phi),
                                 DDATA(p),
                                 DDATA(grad_p),
                                 DDATA(u),
                                 DDATA(v),
                                 DDATA(mom_u_acc),
                                 DDATA(dmom_u_acc_u),
                                 DDATA(mom_v_acc),
                                 DDATA(dmom_v_acc_v),
                                 DDATA(mass_adv),
                                 DDATA(dmass_adv_u),
                                 DDATA(dmass_adv_v),
                                 DDATA(mom_u_adv),
                                 DDATA(dmom_u_adv_u),
                                 DDATA(dmom_u_adv_v),
                                 DDATA(mom_v_adv),
                                 DDATA(dmom_v_adv_u),
                                 DDATA(dmom_v_adv_v),
                                 DDATA(mom_u_diff_ten),
                                 DDATA(mom_v_diff_ten),
                                 DDATA(mom_u_source),
                                 DDATA(mom_v_source),
                                 DDATA(mom_u_ham),
                                 DDATA(dmom_u_ham_grad_p),
                                 DDATA(mom_v_ham),
                                 DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_2D_Evaluate(PyObject* self, 
                                                                                 PyObject* args)
{
  int i,nPoints=1;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_ST_LS_SO_2D_Evaluate(nPoints,
                                            eps_density,
                                            eps_viscosity,
                                            sigma,
                                            rho0,
                                            nu0,
                                            rho1,
                                            nu1,
                                            DDATA(g),
                                            DDATA(phi),
                                            DDATA(n),
                                            DDATA(kappa),
                                            DDATA(p),
                                            DDATA(grad_p),
                                            DDATA(u),
                                            DDATA(v),
                                            DDATA(mom_u_acc),
                                            DDATA(dmom_u_acc_u),
                                            DDATA(mom_v_acc),
                                            DDATA(dmom_v_acc_v),
                                            DDATA(mass_adv),
                                            DDATA(dmass_adv_u),
                                            DDATA(dmass_adv_v),
                                            DDATA(mom_u_adv),
                                            DDATA(dmom_u_adv_u),
                                            DDATA(dmom_u_adv_v),
                                            DDATA(mom_v_adv),
                                            DDATA(dmom_v_adv_u),
                                            DDATA(dmom_v_adv_v),
                                            DDATA(mom_u_diff_ten),
                                            DDATA(mom_v_diff_ten),
                                            DDATA(mom_uv_diff_ten),
                                            DDATA(mom_vu_diff_ten),
                                            DDATA(mom_u_source),
                                            DDATA(mom_v_source),
                                            DDATA(mom_u_ham),
                                            DDATA(dmom_u_ham_grad_p),
                                            DDATA(mom_v_ham),
                                            DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(PyObject* self, 
                                                                                 PyObject* args)
{
  int i,nPoints=1;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(nPoints,
                                            eps_density,
                                            eps_viscosity,
                                            sigma,
                                            rho0,
                                            nu0,
                                            rho1,
                                            nu1,
                                            DDATA(g),
                                            DDATA(phi),
                                            DDATA(n),
                                            DDATA(kappa),
                                            DDATA(p),
                                            DDATA(grad_p),
                                            DDATA(u),
                                            DDATA(v),
                                            DDATA(mom_u_acc),
                                            DDATA(dmom_u_acc_u),
                                            DDATA(mom_v_acc),
                                            DDATA(dmom_v_acc_v),
                                            DDATA(mass_adv),
                                            DDATA(dmass_adv_u),
                                            DDATA(dmass_adv_v),
                                            DDATA(mom_u_adv),
                                            DDATA(dmom_u_adv_u),
                                            DDATA(dmom_u_adv_v),
                                            DDATA(mom_v_adv),
                                            DDATA(dmom_v_adv_u),
                                            DDATA(dmom_v_adv_v),
                                            DDATA(mom_u_diff_ten),
                                            DDATA(mom_v_diff_ten),
                                            DDATA(mom_uv_diff_ten),
                                            DDATA(mom_vu_diff_ten),
                                            DDATA(mom_u_source),
                                            DDATA(mom_v_source),
                                            DDATA(mom_u_ham),
                                            DDATA(dmom_u_ham_grad_p),
                                            DDATA(mom_v_ham),
                                            DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsThreephaseNavierStokes_ST_LS_SO_2D_Evaluate(PyObject* self, 
										   PyObject* args)
{
  int i,nPoints=1;
  double boundaryPenaltyCoef,volumePenaltyCoef,
    eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1,rho_s,nu_s;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *phi_s,
    *n_s,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *mom_v_source,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &boundaryPenaltyCoef,
		       &volumePenaltyCoef,
		       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &rho_s,
		       &nu_s,
                       &g,
                       &phi,
                       &n,
                       &kappa,
		       &phi_s,
		       &n_s,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &dmom_u_source_u,
                       &dmom_u_source_v,
                       &mom_v_source,
                       &dmom_v_source_u,
                       &dmom_v_source_v,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  ThreephaseNavierStokes_ST_LS_SO_2D_Evaluate(nPoints,
					      boundaryPenaltyCoef,
					      volumePenaltyCoef,
					      eps_density,
					      eps_viscosity,
					      sigma,
					      rho0,
					      nu0,
					      rho1,
					      nu1,
					      rho_s,
					      nu_s,
					      DDATA(g),
					      DDATA(phi),
					      DDATA(n),
					      DDATA(kappa),
					      DDATA(phi_s),
					      DDATA(n_s),
					      DDATA(p),
					      DDATA(grad_p),
					      DDATA(u),
					      DDATA(v),
					      DDATA(mom_u_acc),
					      DDATA(dmom_u_acc_u),
					      DDATA(mom_v_acc),
					      DDATA(dmom_v_acc_v),
					      DDATA(mass_adv),
					      DDATA(dmass_adv_u),
					      DDATA(dmass_adv_v),
					      DDATA(mom_u_adv),
					      DDATA(dmom_u_adv_u),
					      DDATA(dmom_u_adv_v),
					      DDATA(mom_v_adv),
					      DDATA(dmom_v_adv_u),
					      DDATA(dmom_v_adv_v),
					      DDATA(mom_u_diff_ten),
					      DDATA(mom_v_diff_ten),
					      DDATA(mom_uv_diff_ten),
					      DDATA(mom_vu_diff_ten),
					      DDATA(mom_u_source),
					      DDATA(dmom_u_source_u),
					      DDATA(dmom_u_source_v),
					      DDATA(mom_v_source),
					      DDATA(dmom_v_source_u),
					      DDATA(dmom_v_source_v),
					      DDATA(mom_u_ham),
					      DDATA(dmom_u_ham_grad_p),
					      DDATA(mom_v_ham),
					      DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_3D_Evaluate(PyObject* self, 
                                                                                 PyObject* args)
{
  int i,nPoints=1;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_ST_LS_SO_3D_Evaluate(nPoints,
                                            eps_density,
                                            eps_viscosity,
                                            sigma,
                                            rho0,
                                            nu0,
                                            rho1,
                                            nu1,
                                            DDATA(g),
                                            DDATA(phi),
                                            DDATA(n),
                                            DDATA(kappa),
                                            DDATA(p),
                                            DDATA(grad_p),
                                            DDATA(u),
                                            DDATA(v),
                                            DDATA(w),
                                            DDATA(mom_u_acc),
                                            DDATA(dmom_u_acc_u),
                                            DDATA(mom_v_acc),
                                            DDATA(dmom_v_acc_v),
                                            DDATA(mom_w_acc),
                                            DDATA(dmom_w_acc_w),
                                            DDATA(mass_adv),
                                            DDATA(dmass_adv_u),
                                            DDATA(dmass_adv_v),
                                            DDATA(dmass_adv_w),
                                            DDATA(mom_u_adv),
                                            DDATA(dmom_u_adv_u),
                                            DDATA(dmom_u_adv_v),
                                            DDATA(dmom_u_adv_w),
                                            DDATA(mom_v_adv),
                                            DDATA(dmom_v_adv_u),
                                            DDATA(dmom_v_adv_v),
                                            DDATA(dmom_v_adv_w),
                                            DDATA(mom_w_adv),
                                            DDATA(dmom_w_adv_u),
                                            DDATA(dmom_w_adv_v),
                                            DDATA(dmom_w_adv_w),
                                            DDATA(mom_u_diff_ten),
                                            DDATA(mom_v_diff_ten),
                                            DDATA(mom_w_diff_ten),
                                            DDATA(mom_uv_diff_ten),
                                            DDATA(mom_uw_diff_ten),
                                            DDATA(mom_vu_diff_ten),
                                            DDATA(mom_vw_diff_ten),
                                            DDATA(mom_wu_diff_ten),
                                            DDATA(mom_wv_diff_ten),
                                            DDATA(mom_u_source),
                                            DDATA(mom_v_source),
                                            DDATA(mom_w_source),
                                            DDATA(mom_u_ham),
                                            DDATA(dmom_u_ham_grad_p),
                                            DDATA(mom_v_ham),
                                            DDATA(dmom_v_ham_grad_p),
                                            DDATA(mom_w_ham),
                                            DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(PyObject* self, 
                                                                                 PyObject* args)
{
  int i,nPoints=1;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(nPoints,
                                            eps_density,
                                            eps_viscosity,
                                            sigma,
                                            rho0,
                                            nu0,
                                            rho1,
                                            nu1,
                                            DDATA(g),
                                            DDATA(phi),
                                            DDATA(n),
                                            DDATA(kappa),
                                            DDATA(p),
                                            DDATA(grad_p),
                                            DDATA(u),
                                            DDATA(v),
                                            DDATA(w),
                                            DDATA(mom_u_acc),
                                            DDATA(dmom_u_acc_u),
                                            DDATA(mom_v_acc),
                                            DDATA(dmom_v_acc_v),
                                            DDATA(mom_w_acc),
                                            DDATA(dmom_w_acc_w),
                                            DDATA(mass_adv),
                                            DDATA(dmass_adv_u),
                                            DDATA(dmass_adv_v),
                                            DDATA(dmass_adv_w),
                                            DDATA(mom_u_adv),
                                            DDATA(dmom_u_adv_u),
                                            DDATA(dmom_u_adv_v),
                                            DDATA(dmom_u_adv_w),
                                            DDATA(mom_v_adv),
                                            DDATA(dmom_v_adv_u),
                                            DDATA(dmom_v_adv_v),
                                            DDATA(dmom_v_adv_w),
                                            DDATA(mom_w_adv),
                                            DDATA(dmom_w_adv_u),
                                            DDATA(dmom_w_adv_v),
                                            DDATA(dmom_w_adv_w),
                                            DDATA(mom_u_diff_ten),
                                            DDATA(mom_v_diff_ten),
                                            DDATA(mom_w_diff_ten),
                                            DDATA(mom_uv_diff_ten),
                                            DDATA(mom_uw_diff_ten),
                                            DDATA(mom_vu_diff_ten),
                                            DDATA(mom_vw_diff_ten),
                                            DDATA(mom_wu_diff_ten),
                                            DDATA(mom_wv_diff_ten),
                                            DDATA(mom_u_source),
                                            DDATA(mom_v_source),
                                            DDATA(mom_w_source),
                                            DDATA(mom_u_ham),
                                            DDATA(dmom_u_ham_grad_p),
                                            DDATA(mom_v_ham),
                                            DDATA(dmom_v_ham_grad_p),
                                            DDATA(mom_w_ham),
                                            DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsThreephaseNavierStokes_ST_LS_SO_3D_Evaluate(PyObject* self, 
										   PyObject* args)
{
  int i,nPoints=1;
  double boundaryPenaltyCoef,volumePenaltyCoef,eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1,rho_s,nu_s;
  PyObject *g,
    *phi,
    *n,
    *kappa,
    *phi_s,
    *n_s,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_u_source_w,
    *mom_v_source,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *dmom_v_source_w,
    *mom_w_source,
    *dmom_w_source_u,
    *dmom_w_source_v,
    *dmom_w_source_w,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &boundaryPenaltyCoef,
		       &volumePenaltyCoef,
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &rho_s,
                       &nu_s,
		       &g,
                       &phi,
                       &n,
                       &kappa,
                       &phi_s,
                       &n_s,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
		       &dmom_u_source_w,
                       &mom_v_source,
		       &dmom_v_source_u,
		       &dmom_v_source_v,
		       &dmom_v_source_w,
                       &mom_w_source,
		       &dmom_w_source_u,
		       &dmom_w_source_v,
		       &dmom_w_source_w,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  ThreephaseNavierStokes_ST_LS_SO_3D_Evaluate(nPoints,
					      boundaryPenaltyCoef,
					      volumePenaltyCoef,
					      eps_density,
					      eps_viscosity,
					      sigma,
					      rho0,
					      nu0,
					      rho1,
					      nu1,
					      rho_s,
					      nu_s,
					      DDATA(g),
					      DDATA(phi),
					      DDATA(n),
					      DDATA(kappa),
					      DDATA(phi_s),
					      DDATA(n_s),
					      DDATA(p),
					      DDATA(grad_p),
					      DDATA(u),
					      DDATA(v),
					      DDATA(w),
					      DDATA(mom_u_acc),
					      DDATA(dmom_u_acc_u),
					      DDATA(mom_v_acc),
					      DDATA(dmom_v_acc_v),
					      DDATA(mom_w_acc),
					      DDATA(dmom_w_acc_w),
					      DDATA(mass_adv),
					      DDATA(dmass_adv_u),
					      DDATA(dmass_adv_v),
					      DDATA(dmass_adv_w),
					      DDATA(mom_u_adv),
					      DDATA(dmom_u_adv_u),
					      DDATA(dmom_u_adv_v),
					      DDATA(dmom_u_adv_w),
					      DDATA(mom_v_adv),
					      DDATA(dmom_v_adv_u),
					      DDATA(dmom_v_adv_v),
					      DDATA(dmom_v_adv_w),
					      DDATA(mom_w_adv),
					      DDATA(dmom_w_adv_u),
					      DDATA(dmom_w_adv_v),
					      DDATA(dmom_w_adv_w),
					      DDATA(mom_u_diff_ten),
					      DDATA(mom_v_diff_ten),
					      DDATA(mom_w_diff_ten),
					      DDATA(mom_uv_diff_ten),
					      DDATA(mom_uw_diff_ten),
					      DDATA(mom_vu_diff_ten),
					      DDATA(mom_vw_diff_ten),
					      DDATA(mom_wu_diff_ten),
					      DDATA(mom_wv_diff_ten),
					      DDATA(mom_u_source),
					      DDATA(dmom_u_source_u),
					      DDATA(dmom_u_source_v),
					      DDATA(dmom_u_source_w),
					      DDATA(mom_v_source),
					      DDATA(dmom_v_source_u),
					      DDATA(dmom_v_source_v),
					      DDATA(dmom_v_source_w),
					      DDATA(mom_w_source),
					      DDATA(dmom_w_source_u),
					      DDATA(dmom_w_source_v),
					      DDATA(dmom_w_source_w),
					      DDATA(mom_u_ham),
					      DDATA(dmom_u_ham_grad_p),
					      DDATA(mom_v_ham),
					      DDATA(dmom_v_ham_grad_p),
					      DDATA(mom_w_ham),
					      DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_LS_SO_3D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &g,
                       &phi,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_LS_SO_3D_Evaluate(nPoints,
					 eps,
					 rho0,
					 nu0,
					 rho1,
					 nu1,
					 DDATA(g),
					 DDATA(phi),
					 DDATA(p),
					 DDATA(grad_p),
					 DDATA(u),
					 DDATA(v),
					 DDATA(w),
					 DDATA(mom_u_acc),
					 DDATA(dmom_u_acc_u),
					 DDATA(mom_v_acc),
					 DDATA(dmom_v_acc_v),
					 DDATA(mom_w_acc),
					 DDATA(dmom_w_acc_w),
					 DDATA(mass_adv),
					 DDATA(dmass_adv_u),
					 DDATA(dmass_adv_v),
					 DDATA(dmass_adv_w),
					 DDATA(mom_u_adv),
					 DDATA(dmom_u_adv_u),
					 DDATA(dmom_u_adv_v),
					 DDATA(dmom_u_adv_w),
					 DDATA(mom_v_adv),
					 DDATA(dmom_v_adv_u),
					 DDATA(dmom_v_adv_v),
					 DDATA(dmom_v_adv_w),
					 DDATA(mom_w_adv),
					 DDATA(dmom_w_adv_u),
					 DDATA(dmom_w_adv_v),
					 DDATA(dmom_w_adv_w),
					 DDATA(mom_u_diff_ten),
					 DDATA(mom_v_diff_ten),
					 DDATA(mom_w_diff_ten),
					 DDATA(mom_u_source),
					 DDATA(mom_v_source),
					 DDATA(mom_w_source),
                                         DDATA(mom_u_ham),
					 DDATA(dmom_u_ham_grad_p),
					 DDATA(mom_v_ham),
					 DDATA(dmom_v_ham_grad_p),
					 DDATA(mom_w_ham),
					 DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseStokes_LS_SO_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseStokes_LS_SO_2D_Evaluate(nPoints,
                                   eps,
                                   rho0,
                                   nu0,
                                   rho1,
                                   nu1,
                                   DDATA(g),
                                   DDATA(phi),
                                   DDATA(p),
                                   DDATA(grad_p),
                                   DDATA(u),
                                   DDATA(v),
                                   DDATA(mom_u_acc),
                                   DDATA(dmom_u_acc_u),
                                   DDATA(mom_v_acc),
                                   DDATA(dmom_v_acc_v),
                                   DDATA(mass_adv),
                                   DDATA(dmass_adv_u),
                                   DDATA(dmass_adv_v),
                                   DDATA(mom_u_diff_ten),
                                   DDATA(mom_v_diff_ten),
                                   DDATA(mom_u_source),
                                   DDATA(mom_v_source),
                                   DDATA(mom_u_ham),
                                   DDATA(dmom_u_ham_grad_p),
                                   DDATA(mom_v_ham),
                                   DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseStokes_LS_SO_3D_Evaluate(PyObject* self, 
                                                          PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *phi,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &phi,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseStokes_LS_SO_3D_Evaluate(nPoints,
                           eps,
                           rho0,
                           nu0,
                           rho1,
                           nu1,
                           DDATA(g),
                           DDATA(phi),
                           DDATA(p),
                                   DDATA(grad_p),
                           DDATA(u),
                           DDATA(v),
                           DDATA(w),
                           DDATA(mom_u_acc),
                           DDATA(dmom_u_acc_u),
                           DDATA(mom_v_acc),
                           DDATA(dmom_v_acc_v),
                           DDATA(mom_w_acc),
                           DDATA(dmom_w_acc_w),
                           DDATA(mass_adv),
                           DDATA(dmass_adv_u),
                           DDATA(dmass_adv_v),
                           DDATA(dmass_adv_w),
                           DDATA(mom_u_diff_ten),
                           DDATA(mom_v_diff_ten),
                           DDATA(mom_w_diff_ten),
                           DDATA(mom_u_source),
                           DDATA(mom_v_source),
                           DDATA(mom_w_source),
                                   DDATA(mom_u_ham),
                           DDATA(dmom_u_ham_grad_p),
                           DDATA(mom_v_ham),
                           DDATA(dmom_v_ham_grad_p),
                           DDATA(mom_w_ham),
                           DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_VOF_SO_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *vof,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &vof,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_VOF_SO_2D_Evaluate(nPoints,
                                 eps,
                                 rho0,
                                 nu0,
                                 rho1,
                                 nu1,
                                 DDATA(g),
                                 DDATA(vof),
                                 DDATA(p),
                                 DDATA(grad_p),
                                 DDATA(u),
                                 DDATA(v),
                                 DDATA(mom_u_acc),
                                 DDATA(dmom_u_acc_u),
                                 DDATA(mom_v_acc),
                                 DDATA(dmom_v_acc_v),
                                 DDATA(mass_adv),
                                 DDATA(dmass_adv_u),
                                 DDATA(dmass_adv_v),
                                 DDATA(mom_u_adv),
                                 DDATA(dmom_u_adv_u),
                                 DDATA(dmom_u_adv_v),
                                 DDATA(mom_v_adv),
                                 DDATA(dmom_v_adv_u),
                                 DDATA(dmom_v_adv_v),
                                 DDATA(mom_u_diff_ten),
                                 DDATA(mom_v_diff_ten),
                                 DDATA(mom_u_source),
                                 DDATA(mom_v_source),
                                 DDATA(mom_u_ham),
                                 DDATA(dmom_u_ham_grad_p),
                                 DDATA(mom_v_ham),
                                 DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseNavierStokes_VOF_SO_3D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *vof,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &vof,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseNavierStokes_VOF_SO_3D_Evaluate(nPoints,
                                 eps,
                                 rho0,
                                 nu0,
                                 rho1,
                                 nu1,
                                 DDATA(g),
                                 DDATA(vof),
                                 DDATA(p),
                                 DDATA(grad_p),
                                 DDATA(u),
                                 DDATA(v),
                                 DDATA(w),
                                 DDATA(mom_u_acc),
                                 DDATA(dmom_u_acc_u),
                                 DDATA(mom_v_acc),
                                 DDATA(dmom_v_acc_v),
                                 DDATA(mom_w_acc),
                                 DDATA(dmom_w_acc_w),
                                 DDATA(mass_adv),
                                 DDATA(dmass_adv_u),
                                 DDATA(dmass_adv_v),
                                 DDATA(dmass_adv_w),
                                 DDATA(mom_u_adv),
                                 DDATA(dmom_u_adv_u),
                                 DDATA(dmom_u_adv_v),
                                 DDATA(dmom_u_adv_w),
                                 DDATA(mom_v_adv),
                                 DDATA(dmom_v_adv_u),
                                 DDATA(dmom_v_adv_v),
                                 DDATA(dmom_v_adv_w),
                                 DDATA(mom_w_adv),
                                 DDATA(dmom_w_adv_u),
                                 DDATA(dmom_w_adv_v),
                                 DDATA(dmom_w_adv_w),
                                 DDATA(mom_u_diff_ten),
                                 DDATA(mom_v_diff_ten),
                                 DDATA(mom_w_diff_ten),
                                 DDATA(mom_u_source),
                                 DDATA(mom_v_source),
                                 DDATA(mom_w_source),
                                          DDATA(mom_u_ham),
                                 DDATA(dmom_u_ham_grad_p),
                                 DDATA(mom_v_ham),
                                 DDATA(dmom_v_ham_grad_p),
                                 DDATA(mom_w_ham),
                                 DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseStokes_VOF_SO_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *vof,
    *p,
    *grad_p,
    *u,
    *v,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &vof,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseStokes_VOF_SO_2D_Evaluate(nPoints,
                           eps,
                           rho0,
                           nu0,
                           rho1,
                           nu1,
                           DDATA(g),
                           DDATA(vof),
                           DDATA(p),
                           DDATA(grad_p),
                           DDATA(u),
                           DDATA(v),
                           DDATA(mom_u_acc),
                           DDATA(dmom_u_acc_u),
                           DDATA(mom_v_acc),
                           DDATA(dmom_v_acc_v),
                           DDATA(mass_adv),
                           DDATA(dmass_adv_u),
                           DDATA(dmass_adv_v),
                           DDATA(mom_u_diff_ten),
                           DDATA(mom_v_diff_ten),
                           DDATA(mom_u_source),
                           DDATA(mom_v_source),
                           DDATA(mom_u_ham),
                           DDATA(dmom_u_ham_grad_p),
                           DDATA(mom_v_ham),
                           DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsTwophaseStokes_VOF_SO_3D_Evaluate(PyObject* self, 
                                                          PyObject* args)
{
  int i,nPoints=1;
  double eps,rho0,nu0,rho1,nu1;
  PyObject *g,
    *vof,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &eps,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
                       &g,
                       &vof,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  TwophaseStokes_VOF_SO_3D_Evaluate(nPoints,
                           eps,
                           rho0,
                           nu0,
                           rho1,
                           nu1,
                           DDATA(g),
                           DDATA(vof),
                           DDATA(p),
                                    DDATA(grad_p),
                           DDATA(u),
                           DDATA(v),
                           DDATA(w),
                           DDATA(mom_u_acc),
                           DDATA(dmom_u_acc_u),
                           DDATA(mom_v_acc),
                           DDATA(dmom_v_acc_v),
                           DDATA(mom_w_acc),
                           DDATA(dmom_w_acc_w),
                           DDATA(mass_adv),
                           DDATA(dmass_adv_u),
                           DDATA(dmass_adv_v),
                           DDATA(dmass_adv_w),
                           DDATA(mom_u_diff_ten),
                           DDATA(mom_v_diff_ten),
                           DDATA(mom_w_diff_ten),
                           DDATA(mom_u_source),
                           DDATA(mom_v_source),
                           DDATA(mom_w_source),
                                    DDATA(mom_u_ham),
                           DDATA(dmom_u_ham_grad_p),
                           DDATA(mom_v_ham),
                           DDATA(dmom_v_ham_grad_p),
                           DDATA(mom_w_ham),
                           DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsConstantVelocityLevelSetEvaluate(PyObject* self, 
									PyObject* args)
{
  int i,nPoints=1;
  PyObject *b,*x,*u,*gradu,*m,*dm,*f,*df,*H,*dH;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &b,
                       &x,
                       &u,
		       &gradu,
                       &m,
                       &dm,
                       &f,
                       &df,
		       &H,
		       &dH))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  constantVelocityLevelSetEvaluate(nPoints,
				   SHAPE(f)[ND(f)-1],
				   DDATA(b),
				   DDATA(x),
				   DDATA(u),
				   DDATA(gradu),
				   DDATA(m),
				   DDATA(dm),
				   DDATA(f),
				   DDATA(df),
				   DDATA(H),
				   DDATA(dH));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsConstantNormalVelocityLevelSetEvaluate(PyObject* self, 
									      PyObject* args)
{
  int i,nPoints=1;
  PyObject *x,*u,*gradu,*m,*dm,*f,*df,*H,*dH;
  double ns;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOO",
		       &ns,
                       &x,
                       &u,
		       &gradu,
                       &m,
                       &dm,
                       &f,
                       &df,
		       &H,
		       &dH))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  constantNormalVelocityLevelSetEvaluate(nPoints,
					 SHAPE(f)[ND(f)-1],
					 ns,
					 DDATA(x),
					 DDATA(u),
					 DDATA(gradu),
					 DDATA(m),
					 DDATA(dm),
					 DDATA(f),
					 DDATA(df),
					 DDATA(H),
					 DDATA(dH));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsUnitSquareVortexLevelSetEvaluate(PyObject* self, 
									PyObject* args)
{
  int i,nPoints=1;
  double t;
  PyObject *x,*u,*gradu,*m,*dm,*f,*df,*H,*dH;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOO",
		       &t,
                       &x,
                       &u,
		       &gradu,
                       &m,
                       &dm,
                       &f,
                       &df,
		       &H,
		       &dH))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  unitSquareVortexLevelSetEvaluate(nPoints,
				   SHAPE(f)[ND(f)-1],
				   t,
				   DDATA(x),
				   DDATA(u),
				   DDATA(gradu),
				   DDATA(m),
				   DDATA(dm),
				   DDATA(f),
				   DDATA(df),
				   DDATA(H),
				   DDATA(dH));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsUnitSquareRotationLevelSetEvaluate(PyObject* self, 
									PyObject* args)
{
  int i,nPoints=1;
  double t;
  PyObject *x,*u,*gradu,*m,*dm,*f,*df,*H,*dH;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOO",
		       &t,
                       &x,
                       &u,
		       &gradu,
                       &m,
                       &dm,
                       &f,
                       &df,
		       &H,
		       &dH))
    return NULL;
  for(i=0;i<ND(x)-1;i++)
      nPoints *= SHAPE(x)[i];
  unitSquareRotationLevelSetEvaluate(nPoints,
				   SHAPE(f)[ND(f)-1],
				   t,
				   DDATA(x),
				   DDATA(u),
				   DDATA(gradu),
				   DDATA(m),
				   DDATA(dm),
				   DDATA(f),
				   DDATA(df),
				   DDATA(H),
				   DDATA(dH));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsHJBurgersEvaluate(PyObject* self, 
							 PyObject* args)
{
  int i,nPoints=1;
  double offset;
  PyObject *u,*gradu,*m,*dm,*H,*dH;
  if(!PyArg_ParseTuple(args,"dOOOOOO",
		       &offset,
                       &u,
		       &gradu,
                       &m,
                       &dm,
		       &H,
		       &dH))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  HJBurgersEvaluate(nPoints,
		    SHAPE(dH)[ND(dH)-1],
		    offset,
		    DDATA(u),
		    DDATA(gradu),
		    DDATA(m),
		    DDATA(dm),
		    DDATA(H),
		    DDATA(dH));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsEikonalEquationEvaluate(PyObject* self, 
							       PyObject* args)
{
  int i,nPoints=1;
  double rhsval;
  PyObject *u,*gradu,*m,*dm,*H,*dH,*r;
  if(!PyArg_ParseTuple(args,"dOOOOOOO",
		       &rhsval,
                       &u,
		       &gradu,
                       &m,
                       &dm,
		       &H,
		       &dH,
		       &r))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  eikonalEquationEvaluate(nPoints,
			  SHAPE(dH)[ND(dH)-1],
			  rhsval,
			  DDATA(u),
			  DDATA(gradu),
			  DDATA(m),
			  DDATA(dm),
			  DDATA(H),
			  DDATA(dH),
			  DDATA(r));

  Py_INCREF(Py_None); 
  return Py_None;
}


/* /\* jcc Begin add for two phase fractional flow *\/  */
/* static PyObject* ctransportCoefficientsFractionalFlowPhaseForm_saturationEvaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,nc,pskModelFlag,nPoints=1; */
/*   double  Kbar, rhon, rhow, alpha,gMag, sw_min,sw_max,M,R,Temp,p_o, b,bc_lambda, bc_pd,mvg_n,mvg_m,omega,mun,muw; */
/*   PyObject *g,*u,*m,*dm,*phi,*dphi,*f,*df,*a,*da,*q_t,*psiw; */
/*   if(!PyArg_ParseTuple(args,"iidddOddddddddddddddddOOOOOOOOOOO",		   */
/* 		       &nc, */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &gMag, */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &omega,  */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &sw_min, */
/* 		       &sw_max, */
/* 		       &M, */
/* 		       &R, */
/* 		       &Temp, */
/* 		       &p_o, */
/* 		       &b,    */
/*                        &u, */
/*                        &m, */
/*                        &dm, */
/* 		       &phi, */
/* 		       &dphi, */
/*                        &f, */
/*                        &df, */
/*                        &a, */
/*                        &da, */
/* 		       &q_t, */
/* 		       &psiw)) */
/*     return NULL; */
/*   for(i=0;i<ND(f)-1;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/* /\*   FractionalFlowPhaseForm_saturationEvaluate(nPoints, *\/ */
/* /\*                                      SHAPE(f)[ND(f)-1], *\/ */
/* /\* 				     nc, *\/ */
/* /\* 				     pskModelFlag, *\/ */
/* /\* 		     		     Kbar, *\/ */
/* /\* 			     	     rhon, *\/ */
/* /\* 		     		     rhow, *\/ */
/* /\* 		     		     DDATA(g), *\/ */
/* /\* 				     gMag, *\/ */
/* /\* 		     		     alpha, *\/ */
/* /\* 		     		     bc_lambda, *\/ */
/* /\* 		   		     bc_pd, *\/ */
/* /\* 		    		     mvg_n, *\/ */
/* /\* 		     		     mvg_m, *\/ */
/* /\* 		     		     omega, *\/ */
/* /\* 		     		     mun, *\/ */
/* /\* 				     muw, *\/ */
/* /\* 				     sw_min, *\/ */
/* /\* 				     sw_max, *\/ */
/* /\* 				     M, *\/ */
/* /\* 				     R, *\/ */
/* /\* 				     Temp, *\/ */
/* /\* 				     p_o, *\/ */
/* /\* 				     b, *\/ */
/* /\* 		   		     DDATA(u), *\/ */
/* /\* 		    		     DDATA(m), *\/ */
/* /\* 		     		     DDATA(dm), *\/ */
/* /\* 		     		     DDATA(phi), *\/ */
/* /\* 		     		     DDATA(dphi), *\/ */
/* /\* 		 		     DDATA(f), *\/ */
/* /\* 		     		     DDATA(df), *\/ */
/* /\* 		     		     DDATA(a), *\/ */
/* /\* 		     		     DDATA(da), *\/ */
/* /\* 		     		     DDATA(q_t), *\/ */
/* /\* 				     DDATA(psiw)); *\/ */
/*   FractionalFlowPhaseForm_saturationEvaluateV3(nPoints, */
/*                                      SHAPE(f)[ND(f)-1], */
/* 				     nc, */
/* 				     pskModelFlag, */
/* 		     		     Kbar, */
/* 			     	     rhon, */
/* 		     		     rhow, */
/* 		     		     DDATA(g), */
/* 				     gMag, */
/* 		     		     alpha, */
/* 		     		     bc_lambda, */
/* 		   		     bc_pd, */
/* 		    		     mvg_n, */
/* 		     		     mvg_m, */
/* 		     		     omega, */
/* 		     		     mun, */
/* 				     muw, */
/* 				     sw_min, */
/* 				     sw_max, */
/* 				     M, */
/* 				     R, */
/* 				     Temp, */
/* 				     p_o, */
/* 				     b, */
/* 		   		     DDATA(u), */
/* 		    		     DDATA(m), */
/* 		     		     DDATA(dm), */
/* 		     		     DDATA(phi), */
/* 		     		     DDATA(dphi), */
/* 		 		     DDATA(f), */
/* 		     		     DDATA(df), */
/* 		     		     DDATA(a), */
/* 		     		     DDATA(da), */
/* 		     		     DDATA(q_t), */
/* 				     DDATA(psiw)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsFractionalFlowPhaseForm_potentialEvaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nc,nPoints=1; */
/*   double  Kbar, rhon, rhow,gMag, alpha, sw_min,sw_max, M,R,Temp,p_o,b, bc_lambda, bc_pd,mvg_n,mvg_m,omega,mun,muw; */
/*   PyObject *g,*u,*m,*dm,*phi,*dphi,*f,*df,*a,*da,*s_w,*grad_psic; */
/*   if(!PyArg_ParseTuple(args,"iidddOddddddddddddddddOOOOOOOOOOO",		   */
/*                        &nc, */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &gMag, */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &omega,  */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &sw_min, */
/* 		       &sw_max, */
/* 		       &M, */
/* 		       &R, */
/* 		       &Temp, */
/* 		       &p_o, */
/* 		       &b,    */
/*                        &u, */
/*                        &m, */
/*                        &dm, */
/* 		       &phi, */
/* 		       &dphi, */
/*                        &f, */
/*                        &df, */
/*                        &a, */
/*                        &da, */
/* 		       &s_w, */
/* 		       &grad_psic)) */
/*     return NULL; */
/*   for(i=0;i<ND(f)-1;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/* /\*   FractionalFlowPhaseForm_potentialEvaluate(nPoints, *\/ */
/* /\*                                      SHAPE(f)[ND(f)-1], *\/ */
/* /\* 				     nc, *\/ */
/* /\* 				     pskModelFlag, *\/ */
/* /\* 		     		     Kbar, *\/ */
/* /\* 			     	     rhon, *\/ */
/* /\* 		     		     rhow, *\/ */
/* /\* 		     		     DDATA(g), *\/ */
/* /\* 				     gMag, *\/ */
/* /\* 		     		     alpha, *\/ */
/* /\* 		     		     bc_lambda, *\/ */
/* /\* 		   		     bc_pd, *\/ */
/* /\* 		    		     mvg_n, *\/ */
/* /\* 		     		     mvg_m, *\/ */
/* /\* 		     		     omega, *\/ */
/* /\* 		     		     mun, *\/ */
/* /\* 				     muw, *\/ */
/* /\* 				     sw_min, *\/ */
/* /\* 				     sw_max, *\/ */
/* /\* 				     M, *\/ */
/* /\* 				     R, *\/ */
/* /\* 				     Temp, *\/ */
/* /\* 				     p_o, *\/ */
/* /\* 				     b, *\/ */
/* /\* 		   		     DDATA(u), *\/ */
/* /\* 		    		     DDATA(m), *\/ */
/* /\* 		     		     DDATA(dm), *\/ */
/* /\* 		     		     DDATA(phi), *\/ */
/* /\* 		     		     DDATA(dphi), *\/ */
/* /\* 		 		     DDATA(f), *\/ */
/* /\* 		     		     DDATA(df), *\/ */
/* /\* 		     		     DDATA(a), *\/ */
/* /\* 		     		     DDATA(da), *\/ */
/* /\* 		     		     DDATA(s_w), *\/ */
/* /\* 				     DDATA(grad_psic)); *\/ */
/*   FractionalFlowPhaseForm_pressureEvaluateV3(nPoints, */
/*                                      SHAPE(f)[ND(f)-1], */
/* 				     nc, */
/* 				     pskModelFlag, */
/* 		     		     Kbar, */
/* 			     	     rhon, */
/* 		     		     rhow, */
/* 		     		     DDATA(g), */
/* 				     gMag, */
/* 		     		     alpha, */
/* 		     		     bc_lambda, */
/* 		   		     bc_pd, */
/* 		    		     mvg_n, */
/* 		     		     mvg_m, */
/* 		     		     omega, */
/* 		     		     mun, */
/* 				     muw, */
/* 				     sw_min, */
/* 				     sw_max, */
/* 				     M, */
/* 				     R, */
/* 				     Temp, */
/* 				     p_o, */
/* 				     b, */
/* 		   		     DDATA(u), */
/* 		    		     DDATA(m), */
/* 		     		     DDATA(dm), */
/* 		     		     DDATA(phi), */
/* 		     		     DDATA(dphi), */
/* 		 		     DDATA(f), */
/* 		     		     DDATA(df), */
/* 		     		     DDATA(a), */
/* 		     		     DDATA(da), */
/* 		     		     DDATA(s_w), */
/* 				     DDATA(grad_psic)); */
/*   printf("done");  */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */
/* /\* For the heterogenious case *\/  */
/* static PyObject* ctransportCoefficientsFractionalFlowPhaseForm_saturationHetEvaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,nc,pskModelFlag,nPoints=1; */
/*   double  rhon, rhow, b, mun,muw; */
/*   PyObject *Kbar,*alpha,*bc_lambda, *bc_pd,*mvg_m,*thetaS, *thetaR,*g,*u,*m,*dm,*phi,*dphi,*f,*df,*a,*da,*q_t; */
/*   if(!PyArg_ParseTuple(args,"iiOddOOOOOOOdddOOOOOOOOOO",		   */
/* 		       &nc, */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_m, */
/* 		       &thetaS, */
/* 		       &thetaR, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b, */
/*                        &u, */
/*                        &m, */
/*                        &dm, */
/* 		       &phi, */
/* 		       &dphi, */
/*                        &f, */
/*                        &df, */
/*                        &a, */
/*                        &da, */
/* 		       &q_t)) */
/*     return NULL; */
/*   for(i=0;i<ND(f)-1;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/*   FractionalFlowPhaseForm_saturationHetEvaluate(nPoints, */
/*                                      SHAPE(f)[ND(f)-1], */
/* 				     nc, */
/* 				     pskModelFlag, */
/* 		     		     DDATA(Kbar), */
/* 			     	     rhon, */
/* 		     		     rhow, */
/* 		     		     DDATA(g),  */
/* 		     		     DDATA(alpha), */
/* 		     		     DDATA(bc_lambda), */
/* 		   		     DDATA(bc_pd),  */
/* 		     		     DDATA(mvg_m), */
/* 				     DDATA(thetaS), */
/* 		     		     DDATA(thetaR),  */
/* 		     		     mun, */
/* 				     muw, */
/* 				     b,    */
/* 		   		     DDATA(u), */
/* 		    		     DDATA(m), */
/* 		     		     DDATA(dm), */
/* 		     		     DDATA(phi), */
/* 		     		     DDATA(dphi), */
/* 		 		     DDATA(f), */
/* 		     		     DDATA(df), */
/* 		     		     DDATA(a), */
/* 		     		     DDATA(da), */
/* 		     		     DDATA(q_t));			                                                  */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsFractionalFlowPhaseForm_potentialHetEvaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nc,nPoints=1; */
/*   double  rhon, rhow, b, mun, muw; */
/*   PyObject *Kbar,*g,*alpha,*bc_lambda, *bc_pd,*mvg_m,*thetaS, *thetaR, *u,*m,*dm,*phi,*dphi,*f,*df,*a,*da,*s_w,*grad_psic; */
/*   if(!PyArg_ParseTuple(args,"iiOddOOOOOOOdddOOOOOOOOOOO",		   */
/*                        &nc, */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_m, */
/* 		       &thetaS, */
/* 		       &thetaR, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b,    */
/*                        &u, */
/*                        &m, */
/*                        &dm, */
/* 		       &phi, */
/* 		       &dphi, */
/*                        &f, */
/*                        &df, */
/*                        &a, */
/*                        &da, */
/* 		       &s_w, */
/* 		       &grad_psic)) */
/*     return NULL; */
/*   for(i=0;i<ND(f)-1;i++) */
/*       nPoints *= SHAPE(f)[i]; */
/*   FractionalFlowPhaseForm_potentialHetEvaluate(nPoints, */
/*                                      SHAPE(f)[ND(f)-1], */
/* 				     nc, */
/* 				     pskModelFlag, */
/* 		     		     DDATA(Kbar), */
/* 			     	     rhon, */
/* 		     		     rhow, */
/* 		     		     DDATA(g),  */
/* 		     		     DDATA(alpha), */
/* 		     		     DDATA(bc_lambda), */
/* 		   		     DDATA(bc_pd),  */
/* 		     		     DDATA(mvg_m), */
/* 		     		     DDATA(thetaS), */
/* 				     DDATA(thetaR),  */
/* 		     		     mun, */
/* 				     muw, */
/* 				     b,    */
/* 		   		     DDATA(u), */
/* 		    		     DDATA(m), */
/* 		     		     DDATA(dm), */
/* 		     		     DDATA(phi), */
/* 		     		     DDATA(dphi), */
/* 		 		     DDATA(f), */
/* 		     		     DDATA(df), */
/* 		     		     DDATA(a), */
/* 		     		     DDATA(da), */
/* 		     		     DDATA(s_w), */
/* 				     DDATA(grad_psic));			                                                  */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */
/* /\* jcc end add for two phase fractional flow *\/ */

/* static PyObject* ctransportCoefficientsTwophaseDarcyFC_Evaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nPoints=1; */
/*   double  Kbar, rhon, rhow, alpha,b, bc_lambda, bc_pd,mvg_n,mvg_m,omega,omega_r,mun,muw; */
/*   PyObject *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw, */
/*     *mn,*dmn, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *phi_psin,*dphi_psin_dpsiw,*dphi_psin_dsw, */
/*     *fw,*dfw, */
/*     *fn,*dfn, */
/*     *aw,*daw, */
/*     *an,*dan; */
/*   if(!PyArg_ParseTuple(args,"idddOOddddddddddOOOOOOOOOOOOOOOOOOO", */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g, */
/* 		       &x, */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd, */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &omega, */
/* 		       &omega_r, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b, */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw, */
/*                        &mn, */
/*                        &dmn, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/* 		       &phi_psin, */
/* 		       &dphi_psin_dpsiw, */
/* 		       &dphi_psin_dsw, */
/*                        &fw, */
/*                        &dfw, */
/*                        &fn, */
/*                        &dfn, */
/*                        &aw, */
/*                        &daw, */
/*                        &an, */
/*                        &dan)) */
/*     return NULL; */
/*   for(i=0;i<ND(fw)-1;i++) */
/*       nPoints *= SHAPE(fw)[i]; */
/*   TwophaseDarcyFC_Evaluate(nPoints, */
/*                            SHAPE(fw)[ND(fw)-1], */
/*                            pskModelFlag, */
/*                            Kbar, */
/*                            rhon, */
/*                            rhow, */
/*                            DDATA(g), */
/*                            DDATA(x), */
/*                            alpha, */
/*                            bc_lambda, */
/*                            bc_pd, */
/*                            mvg_n, */
/*                            mvg_m, */
/*                            omega, */
/* 			   omega_r, */
/*                            mun, */
/*                            muw, */
/*                            b, */
/*                            DDATA(sw), */
/*                            DDATA(psiw), */
/*                            DDATA(mw), */
/*                            DDATA(dmw), */
/*                            DDATA(mn), */
/*                            DDATA(dmn), */
/*                            DDATA(phi_psiw), */
/*                            DDATA(dphi_psiw_dpsiw), */
/*                            DDATA(phi_psin), */
/*                            DDATA(dphi_psin_dpsiw), */
/*                            DDATA(dphi_psin_dsw), */
/*                            DDATA(fw), */
/*                            DDATA(dfw), */
/*                            DDATA(fn), */
/*                            DDATA(dfn), */
/*                            DDATA(aw), */
/*                            DDATA(daw), */
/*                            DDATA(an), */
/*                            DDATA(dan)); */
/*   Py_INCREF(Py_None); */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsTwophaseDarcyFCHet_EvaluateV2(PyObject* self,PyObject* args) */
/* { */
/*   int i,nSimplex=1,nPointsPerSimplex=1,pskModelFlag,nPoints=1,nTypes=1; */
/*   double  rhon, rhow, mun, muw, b; */
/*   PyObject *materialTypes, *Ksw,  */
/*     *mvg_alpha, *mvg_m, *mvg_n, */
/*     *bc_pd, *bc_lambda, */
/*     *thetaS, *thetaR,  */
/*     *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw, */
/*     *mn,*dmn, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *phi_psin,*dphi_psin_dpsiw,*dphi_psin_dsw, */
/*     *fw,*dfw, */
/*     *fn,*dfn, */
/*     *aw,*daw, */
/*     *an,*dan; */
/*   if(!PyArg_ParseTuple(args,"iOOdddOOOOOOOOOddOOOOOOOOOOOOOOOOOOO", */
/* 		       &pskModelFlag, */
/* 		       &materialTypes, */
/*                        &Ksw, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &b, */
/* 		       &g,  */
/* 		       &x,  */
/* 		       &mvg_alpha, */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &bc_pd,  */
/* 		       &bc_lambda, */
/* 		       &thetaS, */
/* 		       &thetaR, */
/* 		       &mun, */
/* 		       &muw, */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw, */
/*                        &mn, */
/*                        &dmn, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/* 		       &phi_psin, */
/* 		       &dphi_psin_dpsiw, */
/* 		       &dphi_psin_dsw, */
/*                        &fw, */
/*                        &dfw, */
/*                        &fn, */
/*                        &dfn, */
/*                        &aw, */
/*                        &daw, */
/*                        &an, */
/*                        &dan)) */
/*     return NULL; */
/*   for(i=0; i < ND(sw)-1; i++) */
/*     nSimplex *= SHAPE(sw)[i]; */
/*   for(i=ND(sw)-1;i<ND(sw);i++) */
/*       nPointsPerSimplex *= SHAPE(sw)[i]; */
/*   nTypes = SHAPE(Ksw)[0]; */
/*   TwophaseDarcyFCHet_EvaluateV2(nSimplex, */
/* 				nPointsPerSimplex, */
/* 				SHAPE(fw)[ND(fw)-1], */
/* 				nTypes, */
/* 				pskModelFlag, */
/* 				IDATA(materialTypes), */
/* 				DDATA(Ksw), */
/* 				rhon, */
/* 				rhow, */
/* 				b, */
/* 				DDATA(g),  */
/* 				DDATA(x),  */
/* 				DDATA(mvg_alpha), */
/* 				DDATA(mvg_n), */
/* 				DDATA(mvg_m), */
/* 				DDATA(bc_pd), */
/* 				DDATA(bc_lambda), */
/* 				DDATA(thetaS), */
/* 				DDATA(thetaR), */
/* 				mun, */
/* 				muw, */
/* 				DDATA(sw), */
/* 				DDATA(psiw), */
/* 				DDATA(mw), */
/* 				DDATA(dmw), */
/* 				DDATA(mn), */
/* 				DDATA(dmn), */
/* 				DDATA(phi_psiw), */
/* 				DDATA(dphi_psiw_dpsiw), */
/* 				DDATA(phi_psin), */
/* 				DDATA(dphi_psin_dpsiw), */
/* 				DDATA(dphi_psin_dsw), */
/* 				DDATA(fw), */
/* 				DDATA(dfw), */
/* 				DDATA(fn), */
/* 				DDATA(dfn), */
/* 				DDATA(aw), */
/* 				DDATA(daw), */
/* 				DDATA(an), */
/* 				DDATA(dan)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsTwophaseFFDarcyFCHet_EvaluateV2(PyObject* self,PyObject* args) */
/* { */
/*   int i,nSimplex=1,nPointsPerSimplex=1,pskModelFlag,nPoints=1,nTypes=1; */
/*   double  rhon, rhow, mun, muw, b; */
/*   PyObject *materialTypes, *Ksw,  */
/*     *mvg_alpha, *mvg_m, *mvg_n, */
/*     *bc_pd, *bc_lambda, */
/*     *thetaS, *thetaR,  */
/*     *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw_dsw, */
/*     *mm,*dmm_dsw, */
/*     *phi_psic,*dphi_psic_dsw, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *fw,*dfw_dsw, */
/*     *fm,*dfm_dsw, */
/*     *aw_psiw,*daw_psiw_dsw, */
/*     *am_psiw,*dam_psiw_dsw, */
/*     *am_psic,*dam_psic_dsw; */
/*   if(!PyArg_ParseTuple(args,"iOOdddOOOOOOOOOddOOOOOOOOOOOOOOOOOOOO", */
/* 		       &pskModelFlag, */
/* 		       &materialTypes, */
/*                        &Ksw, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &b, */
/* 		       &g,  */
/* 		       &x,  */
/* 		       &mvg_alpha, */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &bc_pd,  */
/* 		       &bc_lambda, */
/* 		       &thetaS, */
/* 		       &thetaR, */
/* 		       &mun, */
/* 		       &muw, */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw_dsw, */
/*                        &mm, */
/*                        &dmm_dsw, */
/* 		       &phi_psic, */
/* 		       &dphi_psic_dsw, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/*                        &fw, */
/*                        &dfw_dsw, */
/*                        &fm, */
/*                        &dfm_dsw, */
/*                        &aw_psiw, */
/*                        &daw_psiw_dsw, */
/*                        &am_psiw, */
/*                        &dam_psiw_dsw, */
/*                        &am_psic, */
/*                        &dam_psic_dsw)) */
/*     return NULL; */
/*   for(i=0; i < ND(sw)-1; i++) */
/*     nSimplex *= SHAPE(sw)[i]; */
/*   for(i=ND(sw)-1;i<ND(sw);i++) */
/*       nPointsPerSimplex *= SHAPE(sw)[i]; */
/*   nTypes = SHAPE(Ksw)[0]; */
/*   TwophaseFFDarcyFCHet_EvaluateV2(nSimplex, */
/* 				  nPointsPerSimplex, */
/* 				  SHAPE(fw)[ND(fw)-1], */
/* 				  nTypes, */
/* 				  pskModelFlag, */
/* 				  IDATA(materialTypes), */
/* 				  DDATA(Ksw), */
/* 				  rhon, */
/* 				  rhow, */
/* 				  b, */
/* 				  DDATA(g),  */
/* 				  DDATA(x),  */
/* 				  DDATA(mvg_alpha), */
/* 				  DDATA(mvg_n), */
/* 				  DDATA(mvg_m), */
/* 				  DDATA(bc_pd), */
/* 				  DDATA(bc_lambda), */
/* 				  DDATA(thetaS), */
/* 				  DDATA(thetaR), */
/* 				  mun, */
/* 				  muw, */
/* 				  DDATA(sw), */
/* 				  DDATA(psiw), */
/* 				  DDATA(mw), */
/* 				  DDATA(dmw_dsw), */
/* 				  DDATA(mm), */
/* 				  DDATA(dmm_dsw), */
/* 				  DDATA(phi_psic), */
/* 				  DDATA(dphi_psic_dsw), */
/* 				  DDATA(phi_psiw), */
/* 				  DDATA(dphi_psiw_dpsiw), */
/* 				  DDATA(fw), */
/* 				  DDATA(dfw_dsw), */
/* 				  DDATA(fm), */
/* 				  DDATA(dfm_dsw), */
/* 				  DDATA(aw_psiw), */
/* 				  DDATA(daw_psiw_dsw), */
/* 				  DDATA(am_psiw), */
/* 				  DDATA(dam_psiw_dsw), */
/* 				  DDATA(am_psic), */
/* 				  DDATA(dam_psic_dsw)); */

/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */


/* static PyObject* ctransportCoefficientsTwophaseFFDarcyFC_Evaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nPoints=1; */
/*   double  Kbar, rhon, rhow, alpha,b, bc_lambda, bc_pd,mvg_n,mvg_m,omega,mun,muw,omega_r; */
/*   PyObject *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw_dsw, */
/*     *mm,*dmm_dsw, */
/*     *phi_psic,*dphi_psic_dsw, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *fm,*dfm_dsw, */
/*     *aw_psiw,*daw_psiw_dsw, */
/*     *am_psiw,*dam_psiw_dsw, */
/*     *am_psic,*dam_psic_dsw; */
/*   if(!PyArg_ParseTuple(args,"idddOOddddddddddOOOOOOOOOOOOOOOOOO", */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g, */
/* 		       &x, */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd, */
/* 		       &mvg_n, */
/* 		       &mvg_m, */
/* 		       &omega, */
/* 		       &omega_r, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b, */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw_dsw, */
/*                        &mm, */
/*                        &dmm_dsw, */
/* 		       &phi_psic, */
/* 		       &dphi_psic_dsw, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/*                        &fm, */
/*                        &dfm_dsw, */
/*                        &aw_psiw, */
/*                        &daw_psiw_dsw, */
/*                        &am_psiw, */
/*                        &dam_psiw_dsw, */
/*                        &am_psic, */
/*                        &dam_psic_dsw)) */
/*     return NULL; */
/*   for(i=0;i<ND(fm)-1;i++) */
/*       nPoints *= SHAPE(fm)[i]; */
/*   TwophaseFFDarcyFC_Evaluate(nPoints, */
/*                              SHAPE(fm)[ND(fm)-1], */
/*                              pskModelFlag, */
/*                              Kbar, */
/*                              rhon, */
/*                              rhow, */
/*                              DDATA(g), */
/*                              DDATA(x), */
/*                              alpha, */
/*                              bc_lambda, */
/*                              bc_pd, */
/*                              mvg_n, */
/*                              mvg_m, */
/*                              omega, */
/* 			     omega_r, */
/*                              mun, */
/*                              muw, */
/*                              b, */
/*                              DDATA(sw), */
/*                              DDATA(psiw), */
/*                              DDATA(mw), */
/*                              DDATA(dmw_dsw), */
/*                              DDATA(mm), */
/*                              DDATA(dmm_dsw), */
/*                              DDATA(phi_psic), */
/*                              DDATA(dphi_psic_dsw), */
/*                              DDATA(phi_psiw), */
/*                              DDATA(dphi_psiw_dpsiw), */
/*                              DDATA(fm), */
/*                              DDATA(dfm_dsw), */
/*                              DDATA(aw_psiw), */
/*                              DDATA(daw_psiw_dsw), */
/*                              DDATA(am_psiw), */
/*                              DDATA(dam_psiw_dsw), */
/*                              DDATA(am_psic), */
/*                              DDATA(dam_psic_dsw)); */
/*   Py_INCREF(Py_None); */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsTwophaseDarcyFCHet_Evaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nPoints=1; */
/*   double   rhon, rhow, b, mun,muw; */
/*   PyObject *Kbar, *alpha, *bc_lambda, *bc_pd, *mvg_m, *omega, *omega_r, */
/*     *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw, */
/*     *mn,*dmn, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *phi_psin,*dphi_psin_dpsiw,*dphi_psin_dsw, */
/*     *fw,*dfw, */
/*     *fn,*dfn, */
/*     *aw,*daw, */
/*     *an,*dan; */
/*   if(!PyArg_ParseTuple(args,"iOddOOOOOOOOdddOOOOOOOOOOOOOOOOOOO",		   */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &x,  */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_m, */
/* 		       &omega,  */
/* 		       &omega_r, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b,    */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw, */
/*                        &mn, */
/*                        &dmn, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/* 		       &phi_psin, */
/* 		       &dphi_psin_dpsiw, */
/* 		       &dphi_psin_dsw, */
/*                        &fw, */
/*                        &dfw, */
/*                        &fn, */
/*                        &dfn, */
/*                        &aw, */
/*                        &daw, */
/*                        &an, */
/*                        &dan)) */
/*     return NULL; */
/*   for(i=0;i<ND(fw)-1;i++) */
/*       nPoints *= SHAPE(fw)[i]; */
/*   TwophaseDarcyFCHet_Evaluate(nPoints, */
/*                            SHAPE(fw)[ND(fw)-1], */
/*                            pskModelFlag, */
/*                            DDATA(Kbar), */
/*                            rhon, */
/*                            rhow, */
/*                            DDATA(g),  */
/*                            DDATA(x),  */
/*                            DDATA(alpha), */
/*                            DDATA(bc_lambda), */
/*                            DDATA(bc_pd),  */
/*                            DDATA(mvg_m), */
/*                            DDATA(omega),  */
/*                            DDATA(omega_r),  */
/*                            mun, */
/*                            muw, */
/*                            b,    */
/*                            DDATA(sw), */
/*                            DDATA(psiw), */
/*                            DDATA(mw), */
/*                            DDATA(dmw), */
/*                            DDATA(mn), */
/*                            DDATA(dmn), */
/*                            DDATA(phi_psiw), */
/*                            DDATA(dphi_psiw_dpsiw), */
/*                            DDATA(phi_psin), */
/*                            DDATA(dphi_psin_dpsiw), */
/*                            DDATA(dphi_psin_dsw), */
/*                            DDATA(fw), */
/*                            DDATA(dfw), */
/*                            DDATA(fn), */
/*                            DDATA(dfn), */
/*                            DDATA(aw), */
/*                            DDATA(daw), */
/*                            DDATA(an), */
/*                            DDATA(dan)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

/* static PyObject* ctransportCoefficientsTwophaseFFDarcyFCHet_Evaluate(PyObject* self,PyObject* args) */
/* { */
/*   int i,pskModelFlag,nPoints=1; */
/*   double  rhon, rhow, b,mun,muw; */
/*   PyObject *Kbar, *alpha, *bc_lambda, *bc_pd, *mvg_m, *omega, *omega_r, */
/*     *g,*x, */
/*     *sw,*psiw, */
/*     *mw,*dmw_dsw, */
/*     *mm,*dmm_dsw, */
/*     *phi_psic,*dphi_psic_dsw, */
/*     *phi_psiw,*dphi_psiw_dpsiw, */
/*     *fm,*dfm_dsw, */
/*     *aw_psiw,*daw_psiw_dsw, */
/*     *am_psiw,*dam_psiw_dsw, */
/*     *am_psic,*dam_psic_dsw; */
/*   if(!PyArg_ParseTuple(args,"iOddOOOOOOOOdddOOOOOOOOOOOOOOOOOO",		   */
/* 		       &pskModelFlag, */
/*                        &Kbar, */
/* 		       &rhon, */
/* 		       &rhow, */
/* 		       &g,  */
/* 		       &x,  */
/* 		       &alpha, */
/* 		       &bc_lambda, */
/* 		       &bc_pd,  */
/* 		       &mvg_m, */
/* 		       &omega, */
/* 		       &omega_r, */
/* 		       &mun, */
/* 		       &muw, */
/* 		       &b,  */
/*                        &sw, */
/*                        &psiw, */
/*                        &mw, */
/*                        &dmw_dsw, */
/*                        &mm, */
/*                        &dmm_dsw, */
/* 		       &phi_psic, */
/* 		       &dphi_psic_dsw, */
/* 		       &phi_psiw, */
/* 		       &dphi_psiw_dpsiw, */
/*                        &fm, */
/*                        &dfm_dsw, */
/*                        &aw_psiw, */
/*                        &daw_psiw_dsw, */
/*                        &am_psiw, */
/*                        &dam_psiw_dsw, */
/*                        &am_psic, */
/*                        &dam_psic_dsw)) */
/*     return NULL; */
/*   for(i=0;i<ND(fm)-1;i++) */
/*       nPoints *= SHAPE(fm)[i]; */
/*   TwophaseFFDarcyFCHet_Evaluate(nPoints, */
/*                              SHAPE(fm)[ND(fm)-1], */
/*                              pskModelFlag, */
/*                              DDATA(Kbar), */
/*                              rhon, */
/*                              rhow, */
/*                              DDATA(g),  */
/*                              DDATA(x),  */
/*                              DDATA(alpha), */
/*                              DDATA(bc_lambda), */
/*                              DDATA(bc_pd),  */
/*                              DDATA(mvg_m), */
/*                              DDATA(omega),  */
/*                              DDATA(omega_r),  */
/*                              mun, */
/*                              muw, */
/*                              b,    */
/*                              DDATA(sw), */
/*                              DDATA(psiw), */
/*                              DDATA(mw), */
/*                              DDATA(dmw_dsw), */
/*                              DDATA(mm), */
/*                              DDATA(dmm_dsw), */
/*                              DDATA(phi_psic), */
/*                              DDATA(dphi_psic_dsw), */
/*                              DDATA(phi_psiw), */
/*                              DDATA(dphi_psiw_dpsiw), */
/*                              DDATA(fm), */
/*                              DDATA(dfm_dsw), */
/*                              DDATA(aw_psiw), */
/*                              DDATA(daw_psiw_dsw), */
/*                              DDATA(am_psiw), */
/*                              DDATA(dam_psiw_dsw), */
/*                              DDATA(am_psic), */
/*                              DDATA(dam_psic_dsw)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */

static PyObject* ctransportCoefficientsLinearElasticity_1D_Evaluate(PyObject* self, 
                                                                    PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *u,
    *uu_diff_ten,
    *u_force;
  if(!PyArg_ParseTuple(args,"ddOOOO",
                       &E,
                       &nu,
                       &g,
                       &u,
                       &uu_diff_ten,
                       &u_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  LinearElasticity_1D_Evaluate(nPoints,
                               E,
                               nu,
                               DDATA(g),
                               DDATA(u),
                               DDATA(uu_diff_ten),
                               DDATA(u_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsLinearElasticity_2D_Evaluate(PyObject* self, 
                                                                    PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *u,
    *v,
    *uu_diff_ten,*uv_diff_ten,
    *vu_diff_ten,*vv_diff_ten,
    *u_force,
    *v_force;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOO",
                       &E,
                       &nu,
                       &g,
                       &u,
                       &v,
                       &uu_diff_ten,&uv_diff_ten,
                       &vu_diff_ten,&vv_diff_ten,
                       &u_force,
                       &v_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  LinearElasticity_2D_Evaluate(nPoints,
                               E,
                               nu,
                               DDATA(g),
                               DDATA(u),
                               DDATA(v),
                               DDATA(uu_diff_ten),DDATA(uv_diff_ten),
                               DDATA(vu_diff_ten),DDATA(vv_diff_ten),
                               DDATA(u_force),
                               DDATA(v_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsLinearElasticity_3D_Evaluate(PyObject* self, 
                                                                    PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *u,
    *v,
    *w,
    *uu_diff_ten,*uv_diff_ten,*uw_diff_ten,
    *vu_diff_ten,*vv_diff_ten,*vw_diff_ten,
    *wu_diff_ten,*wv_diff_ten,*ww_diff_ten,
    *u_force,
    *v_force,
    *w_force;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOO",
                       &E,
                       &nu,
                       &g,
                       &u,
                       &v,
                       &w,
                       &uu_diff_ten,&uv_diff_ten,&uw_diff_ten,
                       &vu_diff_ten,&vv_diff_ten,&vw_diff_ten,
                       &wu_diff_ten,&wv_diff_ten,&ww_diff_ten,
                       &u_force,
                       &v_force,
                       &w_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  LinearElasticity_3D_Evaluate(nPoints,
                               E,
                               nu,
                               DDATA(g),
                               DDATA(u),
                               DDATA(v),
                               DDATA(w),
                               DDATA(uu_diff_ten),DDATA(uv_diff_ten),DDATA(uw_diff_ten),
                               DDATA(vu_diff_ten),DDATA(vv_diff_ten),DDATA(vw_diff_ten),
                               DDATA(wu_diff_ten),DDATA(wv_diff_ten),DDATA(ww_diff_ten),
                               DDATA(u_force),
                               DDATA(v_force),
                               DDATA(w_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsMovingMesh_1D_Evaluate(PyObject* self, 
							      PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *det_J,
    *u,
    *uu_diff_ten,
    *u_force;
  if(!PyArg_ParseTuple(args,"ddOOOOO",
                       &E,
                       &nu,
                       &g,
                       &det_J,
		       &u,
                       &uu_diff_ten,
                       &u_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  MovingMesh_1D_Evaluate(nPoints,
                               E,
                               nu,
                               DDATA(g),
			       DDATA(det_J),
                               DDATA(u),
                               DDATA(uu_diff_ten),
                               DDATA(u_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsMovingMesh_2D_Evaluate(PyObject* self, 
                                                                    PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *det_J,
    *u,
    *v,
    *uu_diff_ten,*uv_diff_ten,
    *vu_diff_ten,*vv_diff_ten,
    *u_force,
    *v_force;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOO",
                       &E,
                       &nu,
                       &g,
		       &det_J,
                       &u,
                       &v,
                       &uu_diff_ten,&uv_diff_ten,
                       &vu_diff_ten,&vv_diff_ten,
                       &u_force,
                       &v_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  MovingMesh_2D_Evaluate(nPoints,
                               E,
                               nu,
                               DDATA(g),
			       DDATA(det_J),
                               DDATA(u),
                               DDATA(v),
                               DDATA(uu_diff_ten),DDATA(uv_diff_ten),
                               DDATA(vu_diff_ten),DDATA(vv_diff_ten),
                               DDATA(u_force),
                               DDATA(v_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsMovingMesh_3D_Evaluate(PyObject* self, 
							      PyObject* args)
{
  int i,nPoints=1;
  double E,nu;
  PyObject *g,
    *det_J,
    *u,
    *v,
    *w,
    *uu_diff_ten,*uv_diff_ten,*uw_diff_ten,
    *vu_diff_ten,*vv_diff_ten,*vw_diff_ten,
    *wu_diff_ten,*wv_diff_ten,*ww_diff_ten,
    *u_force,
    *v_force,
    *w_force;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOO",
                       &E,
                       &nu,
                       &g,
		       &det_J,
                       &u,
                       &v,
                       &w,
                       &uu_diff_ten,&uv_diff_ten,&uw_diff_ten,
                       &vu_diff_ten,&vv_diff_ten,&vw_diff_ten,
                       &wu_diff_ten,&wv_diff_ten,&ww_diff_ten,
                       &u_force,
                       &v_force,
                       &w_force))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  MovingMesh_3D_Evaluate(nPoints,
                               E,
                               nu,
			 DDATA(det_J),
                               DDATA(g),
                               DDATA(u),
                               DDATA(v),
                               DDATA(w),
                               DDATA(uu_diff_ten),DDATA(uv_diff_ten),DDATA(uw_diff_ten),
                               DDATA(vu_diff_ten),DDATA(vv_diff_ten),DDATA(vw_diff_ten),
                               DDATA(wu_diff_ten),DDATA(wv_diff_ten),DDATA(ww_diff_ten),
                               DDATA(u_force),
                               DDATA(v_force),
                               DDATA(w_force));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsLevelSetConservationCoefficientsEvaluate(PyObject* self, 
                                                                                PyObject* args)
{
  int nPoints=1,i;
  double epsHeaviside,epsDirac,epsDiffusion;
  PyObject *u_ls,*H_vof,*u,*r,*dr,*a;
  if(!PyArg_ParseTuple(args,"dddOOOOOO",
                       &epsHeaviside,
                       &epsDirac,
                       &epsDiffusion,
                       &u_ls,
                       &H_vof,
                       &u,
                       &r,
                       &dr,
                       &a))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  levelSetConservationCoefficientsEvaluate(nPoints,
                                           SHAPE(a)[ND(a)-1],
                                           epsHeaviside,
                                           epsDirac,
                                           epsDiffusion,
                                           DDATA(u_ls),
                                           DDATA(H_vof),
                                           DDATA(u),
                                           DDATA(r),
                                           DDATA(dr),
                                           DDATA(a));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsLevelSetConservationCoefficientsEvaluate_sd(PyObject* self, 
                                                                                PyObject* args)
{
  int nPoints=1,i;
  double epsHeaviside,epsDirac;
  PyObject *u_ls,*H_vof,*u,*r,*dr;
  if(!PyArg_ParseTuple(args,"ddOOOOO",
                       &epsHeaviside,
                       &epsDirac,
                       &u_ls,
                       &H_vof,
                       &u,
                       &r,
                       &dr))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  levelSetConservationCoefficientsEvaluate_sd(nPoints,
					      epsHeaviside,
					      epsDirac,
					      DDATA(u_ls),
					      DDATA(H_vof),
					      DDATA(u),
					      DDATA(r),
					      DDATA(dr));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsVolumeAveragedNavierStokesFullDevStress_2D_Evaluate(PyObject* self, 
											   PyObject* args)
{
  int i,nPoints=1;
  double rho,mu;
  PyObject *meanGrainSize,
    *g,
    *p,
    *grad_p,
    *u,
    *v,
    *porosity,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &mu,
		       &meanGrainSize,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
                       &dmom_v_source_u,
                       &dmom_v_source_v,
		       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedNavierStokesFullDevStress_2D_Evaluate(nPoints,
						      rho,
						      mu,
						      DDATA(meanGrainSize),
						      DDATA(g),
						      DDATA(p),
						      DDATA(grad_p),
						      DDATA(u),
						      DDATA(v),
						      DDATA(porosity),
						      DDATA(mom_u_acc),
						      DDATA(dmom_u_acc_u),
						      DDATA(mom_v_acc),
						      DDATA(dmom_v_acc_v),
						      DDATA(mass_adv),
						      DDATA(dmass_adv_u),
						      DDATA(dmass_adv_v),
						      DDATA(mom_u_adv),
						      DDATA(dmom_u_adv_u),
						      DDATA(dmom_u_adv_v),
						      DDATA(mom_v_adv),
						      DDATA(dmom_v_adv_u),
						      DDATA(dmom_v_adv_v),
						      DDATA(mom_u_diff_ten),
						      DDATA(mom_v_diff_ten),
						      DDATA(mom_uv_diff_ten),
						      DDATA(mom_vu_diff_ten),
						      DDATA(mom_u_source),
						      DDATA(mom_v_source),
						      DDATA(dmom_u_source_u),
						      DDATA(dmom_u_source_v),
						      DDATA(dmom_v_source_u),
						      DDATA(dmom_v_source_v),
						      DDATA(mom_u_ham),
						      DDATA(dmom_u_ham_grad_p),
						      DDATA(mom_v_ham),
						      DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVolumeAveragedNavierStokesFullDevStress_3D_Evaluate(PyObject* self, 
											   PyObject* args)
{
  int i,nPoints=1;
  double rho,mu;
  PyObject *meanGrainSize,
    *g,
    *p,
    *grad_p,
    *u,
    *v,
    *w,
    *porosity,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_u_source_w,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *dmom_v_source_w,
    *dmom_w_source_u,
    *dmom_w_source_v,
    *dmom_w_source_w,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &rho,
                       &mu,
                       &meanGrainSize,
                       &g,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
		       &dmom_u_source_w,
                       &dmom_v_source_u,
                       &dmom_v_source_v,
                       &dmom_v_source_w,
                       &dmom_w_source_u,
                       &dmom_w_source_v,
                       &dmom_w_source_w,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedNavierStokesFullDevStress_3D_Evaluate(nPoints,
						      rho,
						      mu,
						      DDATA(meanGrainSize),
						      DDATA(g),
						      DDATA(p),
						      DDATA(grad_p),
						      DDATA(u),
						      DDATA(v),
						      DDATA(w),
						      DDATA(porosity),
						      DDATA(mom_u_acc),
						      DDATA(dmom_u_acc_u),
						      DDATA(mom_v_acc),
						      DDATA(dmom_v_acc_v),
						      DDATA(mom_w_acc),
						      DDATA(dmom_w_acc_w),
						      DDATA(mass_adv),
						      DDATA(dmass_adv_u),
						      DDATA(dmass_adv_v),
						      DDATA(dmass_adv_w),
						      DDATA(mom_u_adv),
						      DDATA(dmom_u_adv_u),
						      DDATA(dmom_u_adv_v),
						      DDATA(dmom_u_adv_w),
						      DDATA(mom_v_adv),
						      DDATA(dmom_v_adv_u),
						      DDATA(dmom_v_adv_v),
						      DDATA(dmom_v_adv_w),
						      DDATA(mom_w_adv),
						      DDATA(dmom_w_adv_u),
						      DDATA(dmom_w_adv_v),
						      DDATA(dmom_w_adv_w),
						      DDATA(mom_u_diff_ten),
						      DDATA(mom_v_diff_ten),
						      DDATA(mom_w_diff_ten),
						      DDATA(mom_uv_diff_ten),
						      DDATA(mom_uw_diff_ten),
						      DDATA(mom_vu_diff_ten),
						      DDATA(mom_vw_diff_ten),
						      DDATA(mom_wu_diff_ten),
						      DDATA(mom_wv_diff_ten),
						      DDATA(mom_u_source),
						      DDATA(mom_v_source),
						      DDATA(mom_w_source),
						      DDATA(dmom_u_source_u),
						      DDATA(dmom_u_source_v),
						      DDATA(dmom_u_source_w),
						      DDATA(dmom_v_source_u),
						      DDATA(dmom_v_source_v),
						      DDATA(dmom_v_source_w),
						      DDATA(dmom_w_source_u),
						      DDATA(dmom_w_source_v),
						      DDATA(dmom_w_source_w),
						      DDATA(mom_u_ham),
						      DDATA(dmom_u_ham_grad_p),
						      DDATA(mom_v_ham),
						      DDATA(dmom_v_ham_grad_p),
						      DDATA(mom_w_ham),
						      DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate(PyObject* self, 
											       PyObject* args)
{
  int i,nPoints=1,killNonlinearDrag;
  double eps_rho,eps_mu,sigma,rho0,nu0,rho1,nu1;
  PyObject *meanGrainSize,
    *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *porosity,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"idddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &killNonlinearDrag,
                       &eps_rho,
		       &eps_mu,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &meanGrainSize,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
                       &dmom_v_source_u,
                       &dmom_v_source_v,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate(nPoints,
							  killNonlinearDrag,
							  eps_rho,
							  eps_mu,
							  sigma,
							  rho0,
							  nu0,
							  rho1,
							  nu1,
							  DDATA(meanGrainSize),
							  DDATA(g),
							  DDATA(phi),
							  DDATA(n),
							  DDATA(kappa),
							  DDATA(p),
							  DDATA(grad_p),
							  DDATA(u),
							  DDATA(v),
							  DDATA(porosity),
							  DDATA(mom_u_acc),
							  DDATA(dmom_u_acc_u),
							  DDATA(mom_v_acc),
							  DDATA(dmom_v_acc_v),
							  DDATA(mass_adv),
							  DDATA(dmass_adv_u),
							  DDATA(dmass_adv_v),
							  DDATA(mom_u_adv),
							  DDATA(dmom_u_adv_u),
							  DDATA(dmom_u_adv_v),
							  DDATA(mom_v_adv),
							  DDATA(dmom_v_adv_u),
							  DDATA(dmom_v_adv_v),
							  DDATA(mom_u_diff_ten),
							  DDATA(mom_v_diff_ten),
							  DDATA(mom_uv_diff_ten),
							  DDATA(mom_vu_diff_ten),
							  DDATA(mom_u_source),
							  DDATA(mom_v_source),
							  DDATA(dmom_u_source_u),
							  DDATA(dmom_u_source_v),
							  DDATA(dmom_v_source_u),
							  DDATA(dmom_v_source_v),
							  DDATA(mom_u_ham),
							  DDATA(dmom_u_ham_grad_p),
							  DDATA(mom_v_ham),
							  DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(PyObject* self, 
												  PyObject* args)
{
  int i,nPoints=1,killNonlinearDrag;
  double eps_rho,eps_mu,sigma,rho0,nu0,rho1,nu1;
  PyObject *meanGrainSize,
    *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *porosity,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p;
  if(!PyArg_ParseTuple(args,"idddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &killNonlinearDrag,
                       &eps_rho,
		       &eps_mu,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &meanGrainSize,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
                       &dmom_v_source_u,
                       &dmom_v_source_v,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd(nPoints,
							     killNonlinearDrag,
							     eps_rho,
							     eps_mu,
							     sigma,
							     rho0,
							     nu0,
							     rho1,
							     nu1,
							     DDATA(meanGrainSize),
							     DDATA(g),
							     DDATA(phi),
							     DDATA(n),
							     DDATA(kappa),
							     DDATA(p),
							     DDATA(grad_p),
							     DDATA(u),
							     DDATA(v),
							     DDATA(porosity),
							     DDATA(mom_u_acc),
							     DDATA(dmom_u_acc_u),
							     DDATA(mom_v_acc),
							     DDATA(dmom_v_acc_v),
							     DDATA(mass_adv),
							     DDATA(dmass_adv_u),
							     DDATA(dmass_adv_v),
							     DDATA(mom_u_adv),
							     DDATA(dmom_u_adv_u),
							     DDATA(dmom_u_adv_v),
							     DDATA(mom_v_adv),
							     DDATA(dmom_v_adv_u),
							     DDATA(dmom_v_adv_v),
							     DDATA(mom_u_diff_ten),
							     DDATA(mom_v_diff_ten),
							     DDATA(mom_uv_diff_ten),
							     DDATA(mom_vu_diff_ten),
							     DDATA(mom_u_source),
							     DDATA(mom_v_source),
							     DDATA(dmom_u_source_u),
							     DDATA(dmom_u_source_v),
							     DDATA(dmom_v_source_u),
							     DDATA(dmom_v_source_v),
							     DDATA(mom_u_ham),
							     DDATA(dmom_u_ham_grad_p),
							     DDATA(mom_v_ham),
							     DDATA(dmom_v_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate(PyObject* self, 
											       PyObject* args)
{
  int i,nPoints=1,killNonlinearDrag;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *meanGrainSize,
    *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *porosity,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_u_source_w,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *dmom_v_source_w,
    *dmom_w_source_u,
    *dmom_w_source_v,
    *dmom_w_source_w,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"idddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &killNonlinearDrag,
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &meanGrainSize,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
		       &dmom_u_source_w,
		       &dmom_v_source_u,
		       &dmom_v_source_v,
		       &dmom_v_source_w,
		       &dmom_w_source_u,
		       &dmom_w_source_v,
		       &dmom_w_source_w,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate(nPoints,
							  killNonlinearDrag,
							  eps_density,
							  eps_viscosity,
							  sigma,
							  rho0,
							  nu0,
							  rho1,
							  nu1,
							  DDATA(meanGrainSize),
							  DDATA(g),
							  DDATA(phi),
							  DDATA(n),
							  DDATA(kappa),
							  DDATA(p),
							  DDATA(grad_p),
							  DDATA(u),
							  DDATA(v),
							  DDATA(w),
							  DDATA(porosity),
							  DDATA(mom_u_acc),
							  DDATA(dmom_u_acc_u),
							  DDATA(mom_v_acc),
							  DDATA(dmom_v_acc_v),
							  DDATA(mom_w_acc),
							  DDATA(dmom_w_acc_w),
							  DDATA(mass_adv),
							  DDATA(dmass_adv_u),
							  DDATA(dmass_adv_v),
							  DDATA(dmass_adv_w),
							  DDATA(mom_u_adv),
							  DDATA(dmom_u_adv_u),
							  DDATA(dmom_u_adv_v),
							  DDATA(dmom_u_adv_w),
							  DDATA(mom_v_adv),
							  DDATA(dmom_v_adv_u),
							  DDATA(dmom_v_adv_v),
							  DDATA(dmom_v_adv_w),
							  DDATA(mom_w_adv),
							  DDATA(dmom_w_adv_u),
							  DDATA(dmom_w_adv_v),
							  DDATA(dmom_w_adv_w),
							  DDATA(mom_u_diff_ten),
							  DDATA(mom_v_diff_ten),
							  DDATA(mom_w_diff_ten),
							  DDATA(mom_uv_diff_ten),
							  DDATA(mom_uw_diff_ten),
							  DDATA(mom_vu_diff_ten),
							  DDATA(mom_vw_diff_ten),
							  DDATA(mom_wu_diff_ten),
							  DDATA(mom_wv_diff_ten),
							  DDATA(mom_u_source),
							  DDATA(mom_v_source),
							  DDATA(mom_w_source),
							  DDATA(dmom_u_source_u),
							  DDATA(dmom_u_source_v),
							  DDATA(dmom_u_source_w),
							  DDATA(dmom_v_source_u),
							  DDATA(dmom_v_source_v),
							  DDATA(dmom_v_source_w),
							  DDATA(dmom_w_source_u),
							  DDATA(dmom_w_source_v),
							  DDATA(dmom_w_source_w),
							  DDATA(mom_u_ham),
							  DDATA(dmom_u_ham_grad_p),
							  DDATA(mom_v_ham),
							  DDATA(dmom_v_ham_grad_p),
							  DDATA(mom_w_ham),
							  DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(PyObject* self, 
												  PyObject* args)
{
  int i,nPoints=1,killNonlinearDrag;
  double eps_density,eps_viscosity,sigma,rho0,nu0,rho1,nu1;
  PyObject *meanGrainSize,
    *g,
    *phi,
    *n,
    *kappa,
    *p,
    *grad_p,
    *u,
    *v,
    *porosity,
    *w,
    *mom_u_acc,
    *dmom_u_acc_u,
    *mom_v_acc,
    *dmom_v_acc_v,
    *mom_w_acc,
    *dmom_w_acc_w,
    *mass_adv,
    *dmass_adv_u,
    *dmass_adv_v,
    *dmass_adv_w,
    *mom_u_adv,
    *dmom_u_adv_u,
    *dmom_u_adv_v,
    *dmom_u_adv_w,
    *mom_v_adv,
    *dmom_v_adv_u,
    *dmom_v_adv_v,
    *dmom_v_adv_w,
    *mom_w_adv,
    *dmom_w_adv_u,
    *dmom_w_adv_v,
    *dmom_w_adv_w,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source,
    *dmom_u_source_u,
    *dmom_u_source_v,
    *dmom_u_source_w,
    *dmom_v_source_u,
    *dmom_v_source_v,
    *dmom_v_source_w,
    *dmom_w_source_u,
    *dmom_w_source_v,
    *dmom_w_source_w,
    *mom_u_ham,
    *dmom_u_ham_grad_p,
    *mom_v_ham,
    *dmom_v_ham_grad_p,
    *mom_w_ham,
    *dmom_w_ham_grad_p;
  if(!PyArg_ParseTuple(args,"idddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &killNonlinearDrag,
                       &eps_density,
                       &eps_viscosity,
                       &sigma,
                       &rho0,
                       &nu0,
                       &rho1,
                       &nu1,
		       &meanGrainSize,
                       &g,
                       &phi,
                       &n,
                       &kappa,
                       &p,
                       &grad_p,
                       &u,
                       &v,
                       &w,
		       &porosity,
                       &mom_u_acc,
                       &dmom_u_acc_u,
                       &mom_v_acc,
                       &dmom_v_acc_v,
                       &mom_w_acc,
                       &dmom_w_acc_w,
                       &mass_adv,
                       &dmass_adv_u,
                       &dmass_adv_v,
                       &dmass_adv_w,
                       &mom_u_adv,
                       &dmom_u_adv_u,
                       &dmom_u_adv_v,
                       &dmom_u_adv_w,
                       &mom_v_adv,
                       &dmom_v_adv_u,
                       &dmom_v_adv_v,
                       &dmom_v_adv_w,
                       &mom_w_adv,
                       &dmom_w_adv_u,
                       &dmom_w_adv_v,
                       &dmom_w_adv_w,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source,
		       &dmom_u_source_u,
		       &dmom_u_source_v,
		       &dmom_u_source_w,
		       &dmom_v_source_u,
		       &dmom_v_source_v,
		       &dmom_v_source_w,
		       &dmom_w_source_u,
		       &dmom_w_source_v,
		       &dmom_w_source_w,
                       &mom_u_ham,
                       &dmom_u_ham_grad_p,
                       &mom_v_ham,
                       &dmom_v_ham_grad_p,
                       &mom_w_ham,
                       &dmom_w_ham_grad_p))
    return NULL;
  for(i=0;i<ND(p);i++)
      nPoints *= SHAPE(p)[i];
  VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd(nPoints,
							     killNonlinearDrag,
							     eps_density,
							     eps_viscosity,
							     sigma,
							     rho0,
							     nu0,
							     rho1,
							     nu1,
							     DDATA(meanGrainSize),
							     DDATA(g),
							     DDATA(phi),
							     DDATA(n),
							     DDATA(kappa),
							     DDATA(p),
							     DDATA(grad_p),
							     DDATA(u),
							     DDATA(v),
							     DDATA(w),
							     DDATA(porosity),
							     DDATA(mom_u_acc),
							     DDATA(dmom_u_acc_u),
							     DDATA(mom_v_acc),
							     DDATA(dmom_v_acc_v),
							     DDATA(mom_w_acc),
							     DDATA(dmom_w_acc_w),
							     DDATA(mass_adv),
							     DDATA(dmass_adv_u),
							     DDATA(dmass_adv_v),
							     DDATA(dmass_adv_w),
							     DDATA(mom_u_adv),
							     DDATA(dmom_u_adv_u),
							     DDATA(dmom_u_adv_v),
							     DDATA(dmom_u_adv_w),
							     DDATA(mom_v_adv),
							     DDATA(dmom_v_adv_u),
							     DDATA(dmom_v_adv_v),
							     DDATA(dmom_v_adv_w),
							     DDATA(mom_w_adv),
							     DDATA(dmom_w_adv_u),
							     DDATA(dmom_w_adv_v),
							     DDATA(dmom_w_adv_w),
							     DDATA(mom_u_diff_ten),
							     DDATA(mom_v_diff_ten),
							     DDATA(mom_w_diff_ten),
							     DDATA(mom_uv_diff_ten),
							     DDATA(mom_uw_diff_ten),
							     DDATA(mom_vu_diff_ten),
							     DDATA(mom_vw_diff_ten),
							     DDATA(mom_wu_diff_ten),
							     DDATA(mom_wv_diff_ten),
							     DDATA(mom_u_source),
							     DDATA(mom_v_source),
							     DDATA(mom_w_source),
							     DDATA(dmom_u_source_u),
							     DDATA(dmom_u_source_v),
							     DDATA(dmom_u_source_w),
							     DDATA(dmom_v_source_u),
							     DDATA(dmom_v_source_v),
							     DDATA(dmom_v_source_w),
							     DDATA(dmom_w_source_u),
							     DDATA(dmom_w_source_v),
							     DDATA(dmom_w_source_w),
							     DDATA(mom_u_ham),
							     DDATA(dmom_u_ham_grad_p),
							     DDATA(mom_v_ham),
							     DDATA(dmom_v_ham_grad_p),
							     DDATA(mom_w_ham),
							     DDATA(dmom_w_ham_grad_p));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficients_VolumeAveragedVOFCoefficientsEvaluate(PyObject* self, 
									      PyObject* args)
{
  int nPoints=1,i;
  double eps;
  PyObject *v,*phi,*porosity,*u,*m,*dm,*f,*df;
  if(!PyArg_ParseTuple(args,"dOOOOOOOO",
                       &eps,
                       &v,
                       &phi,
		       &porosity,
                       &u,
                       &m,
                       &dm,
                       &f,
                       &df))
    return NULL;
  for(i=0;i<ND(u);i++)
    nPoints *= SHAPE(u)[i];
  VolumeAveragedVOFCoefficientsEvaluate(nPoints,
					SHAPE(f)[ND(f)-1],
					eps,
					DDATA(v),
					DDATA(phi),
					DDATA(porosity),
					DDATA(u),
					DDATA(m),
					DDATA(dm),
					DDATA(f),
					DDATA(df));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficients_kEpsilon_2D_Evaluate(PyObject* self, 
							     PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *m_e,
    *dm_e,
    *phi_k,
    *dphi_k,
    *phi_e,
    *dphi_e,
    *f_k,
    *df_k,
    *f_e,
    *df_e,
    *a_k,
    *da_k_dk,
    *da_k_de,
    *a_e,
    *da_e_dk,
    *da_e_de,
    *r_k,
    *dr_k_dk,
    *dr_k_de,
    *r_e,
    *dr_e_dk,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &m_e,
		       &dm_e,
		       &phi_k,
		       &dphi_k,
		       &phi_e,
		       &dphi_e,
		       &f_k,
		       &df_k,
		       &f_e,
		       &df_e,
		       &a_k,
		       &da_k_dk,
		       &da_k_de,
		       &a_e,
		       &da_e_dk,
		       &da_e_de,
		       &r_k,
		       &dr_k_dk,
		       &dr_k_de,
		       &r_e,
		       &dr_e_dk,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_2D_Evaluate(nPoints,
		       SHAPE(f_k)[ND(f_k)-1],
		       sigma_k,
		       sigma_e,
		       c_1,
		       c_2,
		       c_mu,
		       c_e,
		       nu,
		       DDATA(velocity),
		       DDATA(gradu),
		       DDATA(gradv),
		       DDATA(k),
		       DDATA(epsilon),
		       DDATA(m_k),
		       DDATA(dm_k),
		       DDATA(m_e),
		       DDATA(dm_e),
		       DDATA(phi_k),
		       DDATA(dphi_k),
		       DDATA(phi_e),
		       DDATA(dphi_e),
		       DDATA(f_k),
		       DDATA(df_k),
		       DDATA(f_e),
		       DDATA(df_e),
		       DDATA(a_k),
		       DDATA(da_k_dk),
		       DDATA(da_k_de),
		       DDATA(a_e),
		       DDATA(da_e_dk),
		       DDATA(da_e_de),
		       DDATA(r_k),
		       DDATA(dr_k_dk),
		       DDATA(dr_k_de),
		       DDATA(r_e),
		       DDATA(dr_e_dk),
		       DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}
static PyObject* ctransportCoefficients_kEpsilon_2D_Evaluate_sd(PyObject* self, 
								PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *m_e,
    *dm_e,
    *phi_k,
    *dphi_k,
    *phi_e,
    *dphi_e,
    *f_k,
    *df_k,
    *f_e,
    *df_e,
    *a_k,
    *da_k_dk,
    *da_k_de,
    *a_e,
    *da_e_dk,
    *da_e_de,
    *r_k,
    *dr_k_dk,
    *dr_k_de,
    *r_e,
    *dr_e_dk,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &m_e,
		       &dm_e,
		       &phi_k,
		       &dphi_k,
		       &phi_e,
		       &dphi_e,
		       &f_k,
		       &df_k,
		       &f_e,
		       &df_e,
		       &a_k,
		       &da_k_dk,
		       &da_k_de,
		       &a_e,
		       &da_e_dk,
		       &da_e_de,
		       &r_k,
		       &dr_k_dk,
		       &dr_k_de,
		       &r_e,
		       &dr_e_dk,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_2D_Evaluate_sd(nPoints,
			  SHAPE(f_k)[ND(f_k)-1],
			  sigma_k,
			  sigma_e,
			  c_1,
			  c_2,
			  c_mu,
			  c_e,
			  nu,
			  DDATA(velocity),
			  DDATA(gradu),
			  DDATA(gradv),
			  DDATA(k),
			  DDATA(epsilon),
			  DDATA(m_k),
			  DDATA(dm_k),
			  DDATA(m_e),
			  DDATA(dm_e),
			  DDATA(phi_k),
			  DDATA(dphi_k),
			  DDATA(phi_e),
			  DDATA(dphi_e),
			  DDATA(f_k),
			  DDATA(df_k),
			  DDATA(f_e),
			  DDATA(df_e),
			  DDATA(a_k),
			  DDATA(da_k_dk),
			  DDATA(da_k_de),
			  DDATA(a_e),
			  DDATA(da_e_dk),
			  DDATA(da_e_de),
			  DDATA(r_k),
			  DDATA(dr_k_dk),
			  DDATA(dr_k_de),
			  DDATA(r_e),
			  DDATA(dr_e_dk),
			  DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject* ctransportCoefficients_kEpsilon_k_2D_Evaluate_sd(PyObject* self, 
								  PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,c_mu,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *phi_k,
    *dphi_k,
    *f_k,
    *df_k,
    *a_k,
    *da_k_dk,
    *r_k,
    *dr_k_dk;

  if(!PyArg_ParseTuple(args,"dddOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &c_mu,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &phi_k,
		       &dphi_k,
		       &f_k,
		       &df_k,
		       &a_k,
		       &da_k_dk,
		       &r_k,
		       &dr_k_dk))

    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_k_2D_Evaluate_sd(nPoints,
			    SHAPE(f_k)[ND(f_k)-1],
			    sigma_k,
			    c_mu,
			    nu,
			    DDATA(velocity),
			    DDATA(gradu),
			    DDATA(gradv),
			    DDATA(k),
			    DDATA(epsilon),
			    DDATA(m_k),
			    DDATA(dm_k),
			    DDATA(phi_k),
			    DDATA(dphi_k),
			    DDATA(f_k),
			    DDATA(df_k),
			    DDATA(a_k),
			    DDATA(da_k_dk),
			    DDATA(r_k),
			    DDATA(dr_k_dk));
			  

  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject* ctransportCoefficients_kEpsilon_epsilon_2D_Evaluate_sd(PyObject* self, 
									PyObject* args)
{
  int nPoints=1,i;
  double sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *k,
    *epsilon,
    *m_e,
    *dm_e,
    *phi_e,
    *dphi_e,
    *f_e,
    *df_e,
    *a_e,
    *da_e_de,
    *r_e,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"ddddddOOOOOOOOOOOOOOO",
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &k,
		       &epsilon,
		       &m_e,
		       &dm_e,
		       &phi_e,
		       &dphi_e,
		       &f_e,
		       &df_e,
		       &a_e,
		       &da_e_de,
		       &r_e,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(epsilon)[i];

  kEpsilon_epsilon_2D_Evaluate_sd(nPoints,
				  SHAPE(f_e)[ND(f_e)-1],
				  sigma_e,
				  c_1,
				  c_2,
				  c_mu,
				  c_e,
				  nu,
				  DDATA(velocity),
				  DDATA(gradu),
				  DDATA(gradv),
				  DDATA(k),
				  DDATA(epsilon),
				  DDATA(m_e),
				  DDATA(dm_e),
				  DDATA(phi_e),
				  DDATA(dphi_e),
				  DDATA(f_e),
				  DDATA(df_e),
				  DDATA(a_e),
				  DDATA(da_e_de),
				  DDATA(r_e),
				  DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}
static PyObject* ctransportCoefficients_kEpsilon_k_3D_Evaluate_sd(PyObject* self, 
								  PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,c_mu,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *gradw,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *phi_k,
    *dphi_k,
    *f_k,
    *df_k,
    *a_k,
    *da_k_dk,
    *r_k,
    *dr_k_dk;

  if(!PyArg_ParseTuple(args,"dddOOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &c_mu,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &gradw,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &phi_k,
		       &dphi_k,
		       &f_k,
		       &df_k,
		       &a_k,
		       &da_k_dk,
		       &r_k,
		       &dr_k_dk))

    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_k_3D_Evaluate_sd(nPoints,
			    SHAPE(f_k)[ND(f_k)-1],
			    sigma_k,
			    c_mu,
			    nu,
			    DDATA(velocity),
			    DDATA(gradu),
			    DDATA(gradv),
			    DDATA(gradw),
			    DDATA(k),
			    DDATA(epsilon),
			    DDATA(m_k),
			    DDATA(dm_k),
			    DDATA(phi_k),
			    DDATA(dphi_k),
			    DDATA(f_k),
			    DDATA(df_k),
			    DDATA(a_k),
			    DDATA(da_k_dk),
			    DDATA(r_k),
			    DDATA(dr_k_dk));
			  

  Py_INCREF(Py_None); 
  return Py_None;

}
static PyObject* ctransportCoefficients_kEpsilon_epsilon_3D_Evaluate_sd(PyObject* self, 
									PyObject* args)
{
  int nPoints=1,i;
  double sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *gradw,
    *k,
    *epsilon,
    *m_e,
    *dm_e,
    *phi_e,
    *dphi_e,
    *f_e,
    *df_e,
    *a_e,
    *da_e_de,
    *r_e,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"ddddddOOOOOOOOOOOOOOOO",
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &gradw,
		       &k,
		       &epsilon,
		       &m_e,
		       &dm_e,
		       &phi_e,
		       &dphi_e,
		       &f_e,
		       &df_e,
		       &a_e,
		       &da_e_de,
		       &r_e,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(epsilon)[i];

  kEpsilon_epsilon_3D_Evaluate_sd(nPoints,
				  SHAPE(f_e)[ND(f_e)-1],
				  sigma_e,
				  c_1,
				  c_2,
				  c_mu,
				  c_e,
				  nu,
				  DDATA(velocity),
				  DDATA(gradu),
				  DDATA(gradv),
				  DDATA(gradw),
				  DDATA(k),
				  DDATA(epsilon),
				  DDATA(m_e),
				  DDATA(dm_e),
				  DDATA(phi_e),
				  DDATA(dphi_e),
				  DDATA(f_e),
				  DDATA(df_e),
				  DDATA(a_e),
				  DDATA(da_e_de),
				  DDATA(r_e),
				  DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject* ctransportCoefficients_kEpsilon_3D_Evaluate_sd(PyObject* self, 
								PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *gradw,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *m_e,
    *dm_e,
    *phi_k,
    *dphi_k,
    *phi_e,
    *dphi_e,
    *f_k,
    *df_k,
    *f_e,
    *df_e,
    *a_k,
    *da_k_dk,
    *da_k_de,
    *a_e,
    *da_e_dk,
    *da_e_de,
    *r_k,
    *dr_k_dk,
    *dr_k_de,
    *r_e,
    *dr_e_dk,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &gradw,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &m_e,
		       &dm_e,
		       &phi_k,
		       &dphi_k,
		       &phi_e,
		       &dphi_e,
		       &f_k,
		       &df_k,
		       &f_e,
		       &df_e,
		       &a_k,
		       &da_k_dk,
		       &da_k_de,
		       &a_e,
		       &da_e_dk,
		       &da_e_de,
		       &r_k,
		       &dr_k_dk,
		       &dr_k_de,
		       &r_e,
		       &dr_e_dk,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_3D_Evaluate_sd(nPoints,
			  SHAPE(f_k)[ND(f_k)-1],
			  sigma_k,
			  sigma_e,
			  c_1,
			  c_2,
			  c_mu,
			  c_e,
			  nu,
			  DDATA(velocity),
			  DDATA(gradu),
			  DDATA(gradv),
			  DDATA(gradw),
			  DDATA(epsilon),
			  DDATA(m_k),
			  DDATA(dm_k),
			  DDATA(m_e),
			  DDATA(dm_e),
			  DDATA(phi_k),
			  DDATA(dphi_k),
			  DDATA(phi_e),
			  DDATA(dphi_e),
			  DDATA(f_k),
			  DDATA(df_k),
			  DDATA(f_e),
			  DDATA(df_e),
			  DDATA(a_k),
			  DDATA(da_k_dk),
			  DDATA(da_k_de),
			  DDATA(a_e),
			  DDATA(da_e_dk),
			  DDATA(da_e_de),
			  DDATA(r_k),
			  DDATA(dr_k_dk),
			  DDATA(dr_k_de),
			  DDATA(r_e),
			  DDATA(dr_e_dk),
			  DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject* ctransportCoefficients_kEpsilon_3D_Evaluate(PyObject* self, 
							     PyObject* args)
{
  int nPoints=1,i;
  double sigma_k,sigma_e,c_1,c_2,c_mu,c_e,nu;
  PyObject *velocity,
    *gradu,
    *gradv,
    *gradw,
    *k,
    *epsilon,
    *m_k,
    *dm_k,
    *m_e,
    *dm_e,
    *phi_k,
    *dphi_k,
    *phi_e,
    *dphi_e,
    *f_k,
    *df_k,
    *f_e,
    *df_e,
    *a_k,
    *da_k_dk,
    *da_k_de,
    *a_e,
    *da_e_dk,
    *da_e_de,
    *r_k,
    *dr_k_dk,
    *dr_k_de,
    *r_e,
    *dr_e_dk,
    *dr_e_de;

  if(!PyArg_ParseTuple(args,"dddddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &sigma_k,
		       &sigma_e,
		       &c_1,
		       &c_2,
		       &c_mu,
		       &c_e,
		       &nu,
		       &velocity,
		       &gradu,
		       &gradv,
		       &gradw,
		       &k,
		       &epsilon,
		       &m_k,
		       &dm_k,
		       &m_e,
		       &dm_e,
		       &phi_k,
		       &dphi_k,
		       &phi_e,
		       &dphi_e,
		       &f_k,
		       &df_k,
		       &f_e,
		       &df_e,
		       &a_k,
		       &da_k_dk,
		       &da_k_de,
		       &a_e,
		       &da_e_dk,
		       &da_e_de,
		       &r_k,
		       &dr_k_dk,
		       &dr_k_de,
		       &r_e,
		       &dr_e_dk,
		       &dr_e_de))
    
    return NULL;
  for(i=0;i<ND(k);i++)
    nPoints *= SHAPE(k)[i];

  kEpsilon_3D_Evaluate(nPoints,
		       SHAPE(f_k)[ND(f_k)-1],
		       sigma_k,
		       sigma_e,
		       c_1,
		       c_2,
		       c_mu,
		       c_e,
		       nu,
		       DDATA(velocity),
		       DDATA(gradu),
		       DDATA(gradv),
		       DDATA(gradw),
		       DDATA(epsilon),
		       DDATA(m_k),
		       DDATA(dm_k),
		       DDATA(m_e),
		       DDATA(dm_e),
		       DDATA(phi_k),
		       DDATA(dphi_k),
		       DDATA(phi_e),
		       DDATA(dphi_e),
		       DDATA(f_k),
		       DDATA(df_k),
		       DDATA(f_e),
		       DDATA(df_e),
		       DDATA(a_k),
		       DDATA(da_k_dk),
		       DDATA(da_k_de),
		       DDATA(a_e),
		       DDATA(da_e_dk),
		       DDATA(da_e_de),
		       DDATA(r_k),
		       DDATA(dr_k_dk),
		       DDATA(dr_k_de),
		       DDATA(r_e),
		       DDATA(dr_e_dk),
		       DDATA(dr_e_de));

  Py_INCREF(Py_None); 
  return Py_None;

}


static PyObject* ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_2D_Update(PyObject* self, 
										       PyObject* args)
{
  int i,nPoints=1;
  double rho,nu,c_mu;
  PyObject *k,
    *grad_k,
    *epsilon,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOOO",
		       &rho,
                       &nu,
                       &c_mu,
		       &k,
		       &grad_k,
		       &epsilon,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source))
    return NULL;
  for(i=0;i<ND(k);i++)
      nPoints *= SHAPE(k)[i];
  ReynoldsAveragedNavierStokes_kEpsilon_2D_Update(nPoints,
						  rho,
						  nu,
						  c_mu,
						  DDATA(k),
						  DDATA(grad_k),
						  DDATA(epsilon),
						  DDATA(mom_u_diff_ten),
						  DDATA(mom_v_diff_ten),
						  DDATA(mom_uv_diff_ten),
						  DDATA(mom_vu_diff_ten),
						  DDATA(mom_u_source),
						  DDATA(mom_v_source));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_2D_Update_sd(PyObject* self, 
										       PyObject* args)
{
  int i,nPoints=1;
  double rho,nu,c_mu;
  PyObject *k,
    *grad_k,
    *epsilon,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten,
    *mom_u_source,
    *mom_v_source;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOOO",
		       &rho,
                       &nu,
                       &c_mu,
		       &k,
		       &grad_k,
		       &epsilon,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_u_source,
                       &mom_v_source))
    return NULL;
  for(i=0;i<ND(k);i++)
      nPoints *= SHAPE(k)[i];
  ReynoldsAveragedNavierStokes_kEpsilon_2D_Update_sd(nPoints,
						     rho,
						     nu,
						     c_mu,
						     DDATA(k),
						     DDATA(grad_k),
						     DDATA(epsilon),
						     DDATA(mom_u_diff_ten),
						     DDATA(mom_v_diff_ten),
						     DDATA(mom_uv_diff_ten),
						     DDATA(mom_vu_diff_ten),
						     DDATA(mom_u_source),
						     DDATA(mom_v_source));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_3D_Update(PyObject* self, 
										       PyObject* args)
{
  int i,nPoints=1;
  double rho,nu,c_mu;
  PyObject *k,
    *grad_k,
    *epsilon,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOOOOOOOOO",
		       &rho,
                       &nu,
                       &c_mu,
		       &k,
		       &grad_k,
		       &epsilon,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source))
    return NULL;
  for(i=0;i<ND(k);i++)
      nPoints *= SHAPE(k)[i];
  ReynoldsAveragedNavierStokes_kEpsilon_3D_Update(nPoints,
						  rho,
						  nu,
						  c_mu,
						  DDATA(k),
						  DDATA(grad_k),
						  DDATA(epsilon),
						  DDATA(mom_u_diff_ten),
						  DDATA(mom_v_diff_ten),
						  DDATA(mom_w_diff_ten),
						  DDATA(mom_uv_diff_ten),
						  DDATA(mom_uw_diff_ten),
						  DDATA(mom_vu_diff_ten),
						  DDATA(mom_vw_diff_ten),
						  DDATA(mom_wu_diff_ten),
						  DDATA(mom_wv_diff_ten),
						  DDATA(mom_u_source),
						  DDATA(mom_v_source),
						  DDATA(mom_w_source));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_3D_Update_sd(PyObject* self, 
											  PyObject* args)
{
  int i,nPoints=1;
  double rho,nu,c_mu;
  PyObject *k,
    *grad_k,
    *epsilon,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten,
    *mom_u_source,
    *mom_v_source,
    *mom_w_source;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOOOOOOOOO",
		       &rho,
                       &nu,
                       &c_mu,
		       &k,
		       &grad_k,
		       &epsilon,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_w_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_uw_diff_ten,
                       &mom_vu_diff_ten,
                       &mom_vw_diff_ten,
                       &mom_wu_diff_ten,
                       &mom_wv_diff_ten,
                       &mom_u_source,
                       &mom_v_source,
                       &mom_w_source))
    return NULL;
  for(i=0;i<ND(k);i++)
      nPoints *= SHAPE(k)[i];
  ReynoldsAveragedNavierStokes_kEpsilon_3D_Update_sd(nPoints,
						     rho,
						     nu,
						     c_mu,
						     DDATA(k),
						     DDATA(grad_k),
						     DDATA(epsilon),
						     DDATA(mom_u_diff_ten),
						     DDATA(mom_v_diff_ten),
						     DDATA(mom_w_diff_ten),
						     DDATA(mom_uv_diff_ten),
						     DDATA(mom_uw_diff_ten),
						     DDATA(mom_vu_diff_ten),
						     DDATA(mom_vw_diff_ten),
						     DDATA(mom_wu_diff_ten),
						     DDATA(mom_wv_diff_ten),
						     DDATA(mom_u_source),
						     DDATA(mom_v_source),
						     DDATA(mom_w_source));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsEddyViscosity_2D_Update_sd(PyObject* self, 
								  PyObject* args)
{
  int i,nPoints=1;
  PyObject *nu_t,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &nu_t,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten))
    return NULL;
  for(i=0;i<ND(nu_t);i++)
      nPoints *= SHAPE(nu_t)[i];
  eddyViscosity_2D_Update_sd(nPoints,
			     DDATA(nu_t),
			     DDATA(mom_u_diff_ten),
			     DDATA(mom_v_diff_ten),
			     DDATA(mom_uv_diff_ten),
			     DDATA(mom_vu_diff_ten));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsEddyViscosity_2D_Update(PyObject* self, 
							       PyObject* args)
{
  int i,nPoints=1;
  PyObject *nu_t,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_uv_diff_ten,
    *mom_vu_diff_ten;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &nu_t,
                       &mom_u_diff_ten,
                       &mom_v_diff_ten,
                       &mom_uv_diff_ten,
                       &mom_vu_diff_ten))
    return NULL;
  for(i=0;i<ND(nu_t);i++)
      nPoints *= SHAPE(nu_t)[i];
  eddyViscosity_2D_Update(nPoints,
			  DDATA(nu_t),
			  DDATA(mom_u_diff_ten),
			  DDATA(mom_v_diff_ten),
			  DDATA(mom_uv_diff_ten),
			  DDATA(mom_vu_diff_ten));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsEddyViscosity_3D_Update_sd(PyObject* self, 
								  PyObject* args)
{
  int i,nPoints=1;
  PyObject *nu_t,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &nu_t,
		       &mom_u_diff_ten,
		       &mom_v_diff_ten,
		       &mom_w_diff_ten,
		       &mom_uv_diff_ten,
		       &mom_uw_diff_ten,
		       &mom_vu_diff_ten,
		       &mom_vw_diff_ten,
		       &mom_wu_diff_ten,
		       &mom_wv_diff_ten))
    return NULL;
  for(i=0;i<ND(nu_t);i++)
      nPoints *= SHAPE(nu_t)[i];
  eddyViscosity_3D_Update_sd(nPoints,
			     DDATA(nu_t),
			     DDATA(mom_u_diff_ten),
			     DDATA(mom_v_diff_ten),
			     DDATA(mom_w_diff_ten),
			     DDATA(mom_uv_diff_ten),
			     DDATA(mom_uw_diff_ten),
			     DDATA(mom_vu_diff_ten),
			     DDATA(mom_vw_diff_ten),
			     DDATA(mom_wu_diff_ten),
			     DDATA(mom_wv_diff_ten));
  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctransportCoefficientsEddyViscosity_3D_Update(PyObject* self, 
								  PyObject* args)
{
  int i,nPoints=1;
  PyObject *nu_t,
    *mom_u_diff_ten,
    *mom_v_diff_ten,
    *mom_w_diff_ten,
    *mom_uv_diff_ten,
    *mom_uw_diff_ten,
    *mom_vu_diff_ten,
    *mom_vw_diff_ten,
    *mom_wu_diff_ten,
    *mom_wv_diff_ten;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &nu_t,
		       mom_u_diff_ten,
		       mom_v_diff_ten,
		       mom_w_diff_ten,
		       mom_uv_diff_ten,
		       mom_uw_diff_ten,
		       mom_vu_diff_ten,
		       mom_vw_diff_ten,
		       mom_wu_diff_ten,
		       mom_wv_diff_ten))
    return NULL;
  for(i=0;i<ND(nu_t);i++)
      nPoints *= SHAPE(nu_t)[i];
  eddyViscosity_3D_Update(nPoints,
			  DDATA(nu_t),
			  DDATA(mom_u_diff_ten),
			  DDATA(mom_v_diff_ten),
			  DDATA(mom_w_diff_ten),
			  DDATA(mom_uv_diff_ten),
			  DDATA(mom_uw_diff_ten),
			  DDATA(mom_vu_diff_ten),
			  DDATA(mom_vw_diff_ten),
			  DDATA(mom_wu_diff_ten),
			  DDATA(mom_wv_diff_ten));
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsCalculateEddyViscosity_Smagorinsky_2D(PyObject* self, 
									     PyObject* args)
{
  int i;
  double smagorinskyCoefficient;
  PyObject *h_e,
    *grad_u,
    *grad_v,
    *nu_t;

  if(!PyArg_ParseTuple(args,"dOOOO",
		       &smagorinskyCoefficient,
		       &h_e,
		       &grad_u,
		       &grad_v,
		       &nu_t))
    return NULL;
  calculateEddyViscosity_Smagorinsky_2D(SHAPE(grad_u)[0],
					SHAPE(grad_u)[1],
					smagorinskyCoefficient,
					DDATA(h_e),
					DDATA(grad_u),
					DDATA(grad_v),
					DDATA(nu_t));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* ctransportCoefficientsCalculateEddyViscosity_Smagorinsky_3D(PyObject* self, 
									     PyObject* args)
{
  int i;
  double smagorinskyCoefficient;
  PyObject *h_e,
    *grad_u,
    *grad_v,
    *grad_w,
    *nu_t;

  if(!PyArg_ParseTuple(args,"dOOOOO",
		       &smagorinskyCoefficient,
		       &h_e,
		       &grad_u,
		       &grad_v,
		       &grad_w,
		       &nu_t))
    return NULL;
  calculateEddyViscosity_Smagorinsky_3D(SHAPE(grad_u)[0],
					SHAPE(grad_u)[1],
					smagorinskyCoefficient,
					DDATA(h_e),
					DDATA(grad_u),
					DDATA(grad_v),
					DDATA(grad_w),
					DDATA(nu_t));
						  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsCalculateEddyViscosity_Smagorinsky2P_2D(PyObject* self, 
									       PyObject* args)
{
  int i;
  double smagorinskyCoefficient_0,smagorinskyCoefficient_1,eps;
  PyObject *phi_ls,
    *h_e,
    *grad_u,
    *grad_v,
    *nu_t;

  if(!PyArg_ParseTuple(args,"dddOOOOO",
		       &smagorinskyCoefficient_0,
		       &smagorinskyCoefficient_1,
		       &eps,
		       &phi_ls,
		       &h_e,
		       &grad_u,
		       &grad_v,
		       &nu_t))
    return NULL;
  calculateEddyViscosity_Smagorinsky2P_2D(SHAPE(grad_u)[0],
					  SHAPE(grad_u)[1],
					  smagorinskyCoefficient_0,
					  smagorinskyCoefficient_1,
					  eps,
					  DDATA(phi_ls),
					  DDATA(h_e),
					  DDATA(grad_u),
					  DDATA(grad_v),
					  DDATA(nu_t));
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsCalculateEddyViscosity_Smagorinsky2P_3D(PyObject* self, 
									       PyObject* args)
{
  int i;
  double smagorinskyCoefficient_0,smagorinskyCoefficient_1,eps;
  PyObject *phi_ls,
    *h_e,
    *grad_u,
    *grad_v,
    *grad_w,
    *nu_t;

  if(!PyArg_ParseTuple(args,"dddOOOOOO",
		       &smagorinskyCoefficient_0,
		       &smagorinskyCoefficient_1,
		       &eps,
		       &phi_ls,
		       &h_e,
		       &grad_u,
		       &grad_v,
		       &grad_w,
		       &nu_t))
    return NULL;
  calculateEddyViscosity_Smagorinsky2P_3D(SHAPE(grad_u)[0],
					  SHAPE(grad_u)[1],
					  smagorinskyCoefficient_0,
					  smagorinskyCoefficient_1,
					  eps,
					  DDATA(phi_ls),
					  DDATA(h_e),
					  DDATA(grad_u),
					  DDATA(grad_v),
					  DDATA(grad_w),
					  DDATA(nu_t));
  
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* ctransportCoefficientsScriptedSphereMotionSignedDistance(PyObject* self, 
									  PyObject* args)
{
  int i,nSpheres,nSpace,nPoints=1;
  double t;
  PyObject *radii,
    *centers,
    *phi,
    *n;

  if(!PyArg_ParseTuple(args,"diiOOOO",
 		       &t,
		       &nSpace,
                       &nSpheres,
                       &radii,
                       &centers,
                       &phi,
		       &n))
    return NULL;
  for(i=0;i<ND(phi);i++)
      nPoints *= SHAPE(phi)[i];
  scriptedSphereMotionSignedDistance(nPoints,
				     t,
				     SHAPE(n)[ND(n)-1],
				     SHAPE(radii)[0],
				     DDATA(radii),
				     DDATA(centers),
				     DDATA(phi),
				     DDATA(n));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsShallowWater_1D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double h_eps,g,bedFrictionCoefficient,bedFrictionPower,eddyViscosity;
  PyObject *h,
    *hu,
    *H,
    *mass_acc,
    *dmass_acc_dh,
    *mom_acc,
    *dmom_acc_dhu,
    *mass_adv,
    *dmass_adv_dhu,
    *mom_adv,
    *dmom_adv_dh,
    *dmom_adv_dhu,
    *x,
    *db_dx,
    *mom_source,
    *dmom_source_dh,
    *dmom_source_dhu,
    *mom_diff;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOO",
                       &h_eps,
                       &g,
		       &bedFrictionCoefficient,
		       &bedFrictionPower,
		       &eddyViscosity,
		       &x,
		       &db_dx,
		       &h,
		       &hu,
		       &H,
		       &mass_acc,
		       &dmass_acc_dh,
		       &mom_acc,
		       &dmom_acc_dhu,
		       &mass_adv,
		       &dmass_adv_dhu,
		       &mom_adv,
		       &dmom_adv_dh,
		       &dmom_adv_dhu,
		       &mom_diff,
		       &mom_source,
		       &dmom_source_dh,
		       &dmom_source_dhu))
    return NULL;
  for(i=0;i<ND(h);i++)
    nPoints *= SHAPE(h)[i];
  shallowWater_1D_Evaluate(nPoints,
                           h_eps,
			   g,
			   bedFrictionCoefficient,
			   bedFrictionPower,
			   eddyViscosity,
			   DDATA(x),
			   DDATA(db_dx),
			   DDATA(h),
			   DDATA(hu),
			   DDATA(H),
			   DDATA(mass_acc),
			   DDATA(dmass_acc_dh),
			   DDATA(mom_acc),
			   DDATA(dmom_acc_dhu),
			   DDATA(mass_adv),
			   DDATA(dmass_adv_dhu),
			   DDATA(mom_adv),
			   DDATA(dmom_adv_dh),
			   DDATA(dmom_adv_dhu),
			   DDATA(mom_diff),
			   DDATA(mom_source),
			   DDATA(dmom_source_dh),
			   DDATA(dmom_source_dhu));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsShallowWater_2D_Evaluate(PyObject* self, 
                                                                PyObject* args)
{
  int i,nPoints=1;
  double h_eps,g,bedFrictionCoefficient,bedFrictionPower,eddyViscosity;
  PyObject *h,
    *hu,
    *hv,
    *H,
    *mass_acc,
    *dmass_acc_dh,
    *mom_u_acc,
    *dmom_u_acc_dhu,
    *mom_v_acc,
    *dmom_v_acc_dhv,
    *mass_adv,
    *dmass_adv_dhu,
    *dmass_adv_dhv,
    *mom_u_adv,
    *dmom_u_adv_dh,
    *dmom_u_adv_dhu,
    *dmom_u_adv_dhv,
    *mom_v_adv,
    *dmom_v_adv_dh,
    *dmom_v_adv_dhu,
    *dmom_v_adv_dhv,
    *mom_u_diff,
    *mom_v_diff,
    *x,
    *grad_b,
    *mom_u_source,
    *dmom_u_source_dh,
    *dmom_u_source_dhu,
    *dmom_u_source_dhv,
    *mom_v_source,
    *dmom_v_source_dh,
    *dmom_v_source_dhu,
    *dmom_v_source_dhv;
  if(!PyArg_ParseTuple(args,"dddddOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &h_eps,
                       &g,
		       &bedFrictionCoefficient,
		       &bedFrictionPower,
		       &eddyViscosity,
		       &x,
		       &grad_b,
		       &h,
		       &hu,
		       &hv,
		       &H,
		       &mass_acc,
		       &dmass_acc_dh,
		       &mom_u_acc,
		       &dmom_u_acc_dhu,
		       &mom_v_acc,
		       &dmom_v_acc_dhv,
		       &mass_adv,
		       &dmass_adv_dhu,
		       &dmass_adv_dhv,
		       &mom_u_adv,
		       &dmom_u_adv_dh,
		       &dmom_u_adv_dhu,
		       &dmom_u_adv_dhv,
		       &mom_v_adv,
		       &dmom_v_adv_dh,
		       &dmom_v_adv_dhu,
		       &dmom_v_adv_dhv,
		       &mom_u_diff,
		       &mom_v_diff,
		       &mom_u_source,
		       &dmom_u_source_dh,
		       &dmom_u_source_dhu,
		       &dmom_u_source_dhv,
		       &mom_v_source,
		       &dmom_v_source_dh,
		       &dmom_v_source_dhu,
		       &dmom_v_source_dhv))
    return NULL;
  for(i=0;i<ND(h);i++)
    nPoints *= SHAPE(h)[i];
  shallowWater_2D_Evaluate(nPoints,
                           h_eps,
			   g,
			   bedFrictionCoefficient,
			   bedFrictionPower,
			   eddyViscosity,
			   DDATA(x),
			   DDATA(grad_b),
			   DDATA(h),
			   DDATA(hu),
			   DDATA(hv),
			   DDATA(H),
			   DDATA(mass_acc),
			   DDATA(dmass_acc_dh),
			   DDATA(mom_u_acc),
			   DDATA(dmom_u_acc_dhu),
			   DDATA(mom_v_acc),
			   DDATA(dmom_v_acc_dhv),
			   DDATA(mass_adv),
			   DDATA(dmass_adv_dhu),
			   DDATA(dmass_adv_dhv),
			   DDATA(mom_u_adv),
			   DDATA(dmom_u_adv_dh),
			   DDATA(dmom_u_adv_dhu),
			   DDATA(dmom_u_adv_dhv),
			   DDATA(mom_v_adv),
			   DDATA(dmom_v_adv_dh),
			   DDATA(dmom_v_adv_dhu),
			   DDATA(dmom_v_adv_dhv),
			   DDATA(mom_u_diff),
			   DDATA(mom_v_diff),
			   DDATA(mom_u_source),
			   DDATA(dmom_u_source_dh),
			   DDATA(dmom_u_source_dhu),
			   DDATA(dmom_u_source_dhv),
			   DDATA(mom_v_source),
			   DDATA(dmom_v_source_dh),
			   DDATA(dmom_v_source_dhu),
			   DDATA(dmom_v_source_dhv));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctransportCoefficientsSmoothedHeaviside(PyObject* self, 
							 PyObject* args)
{
  double eps,phi,H;
  if(!PyArg_ParseTuple(args,"dd",
                       &eps,
                       &phi))
    return NULL;

  H = smoothedHeaviside(eps,phi);
  
  return Py_BuildValue("d",H);
}

static PyObject* ctransportCoefficientsSmoothedHeaviside_integral(PyObject* self, 
								  PyObject* args)
{
  double eps=0.0,phi=0.0,HI=0.0;
  if(!PyArg_ParseTuple(args,"dd",
                       &eps,
                       &phi))
    return NULL;
  HI = smoothedHeaviside_integral(eps,phi);
  return Py_BuildValue("d",HI);
}

static PyObject* ctransportCoefficientsSmoothedDirac(PyObject* self, 
						     PyObject* args)
{
  double eps,phi,dH;
  if(!PyArg_ParseTuple(args,"dd",
                       &eps,
                       &phi))
    return NULL;

  dH = smoothedDirac(eps,phi);
  
  return Py_BuildValue("d",dH);
}

static PyObject* ctransportCoefficientsDiffusiveWave1DCoefficientsEvaluate(PyObject* self, PyObject* args)

{
  int nPoints=1,i;
  double alpha,gamma,epsilon;
  PyObject *x,*u,*grad_u,*m,*dm,*a,*da;
  if(!PyArg_ParseTuple(args,"dddOOOOOOO",
		       &alpha,
		       &gamma,
		       &epsilon,
		       &x,
		       &u,
		       &grad_u,
		       &m,
		       &dm,
		       &a,
		       &da))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  diffusiveWave1DEvaluate(nPoints,
			  alpha,
			  gamma,
			  epsilon,
			  DDATA(x),
			  DDATA(u),
			  DDATA(grad_u),
			  DDATA(m),
			  DDATA(dm),
			  DDATA(a),
			  DDATA(da));
  
  Py_INCREF(Py_None);
  return Py_None; 
}


static PyObject* ctransportCoefficientsDiffusiveWave2DCoefficientsEvaluate(PyObject* self, PyObject* args)

{
  int nPoints=1,i,nd;
  double alpha,gamma,epsilon;
  PyObject *x,*u,*grad_u,*m,*dm,*a,*da;
  if(!PyArg_ParseTuple(args,"idddOOOOOOO",
		       &nd,
		       &alpha,
		       &gamma,
		       &epsilon,
		       &x,
		       &u,
		       &grad_u,
		       &m,
		       &dm,
		       &a,
		       &da))
    return NULL;
  for(i=0;i<ND(u);i++)
      nPoints *= SHAPE(u)[i];
  diffusiveWave2DEvaluate(nd,
			  nPoints,
			  alpha,
			  gamma,
			  epsilon,
			  DDATA(x),
			  DDATA(u),
			  DDATA(grad_u),
			  DDATA(m),
			  DDATA(dm),
			  DDATA(a),
			  DDATA(da));
  
  Py_INCREF(Py_None);
  return Py_None; 
}
static PyObject*
ctransportCoefficientsCalculateWaveFunction3d_ref(PyObject* self, 
						  PyObject* args)
{
  double t,epsFact,waveHeight,waveCelerity,waveFrequency,waveNumber,waterDepth;
  int waveFlag;
  PyObject *mesh_trial_ref,*mesh_dof,*mesh_l2g,*elementDiametersArray,
    *omega_s_x,*omega_s_y,*omega_s_z,*source;
  if(!PyArg_ParseTuple(args,"OOOOOOOdiddddddO",
		       &mesh_trial_ref,
		       &mesh_dof,
		       &mesh_l2g,
		       &elementDiametersArray,
		       &omega_s_x,
		       &omega_s_y,
		       &omega_s_z,
                       &t,
                       &waveFlag,
                       &epsFact,
		       &waveHeight,
		       &waveCelerity,
		       &waveFrequency,
		       &waveNumber,
		       &waterDepth,
                       &source))
    return NULL;
  calculateWaveFunction3d_ref(SHAPE(mesh_l2g)[0],
			      SHAPE(mesh_l2g)[1],
			      SHAPE(mesh_trial_ref)[0],
			      DDATA(mesh_trial_ref),
			      DDATA(mesh_dof),
			      IDATA(mesh_l2g),
			      DDATA(elementDiametersArray),
			      DDATA(omega_s_x),
			      DDATA(omega_s_y),
			      DDATA(omega_s_z),
			      t,
			      waveFlag,
			      epsFact,
			      waveHeight,
			      waveCelerity,
			      waveFrequency,
			      waveNumber,
			      waterDepth,
			      DDATA(source));
  Py_INCREF(Py_None); 
  return Py_None;
}

/* static PyObject* ctransportCoefficientsPiecewiseLinearTableLookup(PyObject* self, PyObject* args) */

/* { */
/*   int start=0; */
/*   double x,y,dy; */
/*   PyObject *xv,*yv; */
/*   if(!PyArg_ParseTuple(args,"dOO|i", */
/* 		       &x, */
/* 		       &xv, */
/* 		       &yv, */
/* 		       &start)) */
/*     return NULL; */
/*   piecewiseLinearTableLookup(x, */
/* 			     SHAPE(xv)[0], */
/* 			     &start, */
/* 			     &y, */
/* 			     &dy, */
/* 			     DDATA(xv), */
/* 			     DDATA(yv)); */
/*   return Py_BuildValue("ddi",y,dy,start); */
/* } */

static PyMethodDef ctransportCoefficientsMethods[] = {
  { "linearADR_ConstantCoefficientsEvaluate", 
    ctransportCoefficientsLinearADR_ConstatCoefficientsEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar vtransport equation for the linear, constant coefficients case"}, 
  { "groundwaterTransportCoefficientsEvaluate", 
      ctransportCoefficientsGroundwaterTransportCoefficientsEvaluate,
      METH_VARARGS, 
      "Evaluate  the coefficients of linear advective-diffusion transport in porous media"}, 
  { "groundwaterBiodegradation01EvaluateFC", 
    ctransportCoefficientsGroundwaterBiodegradation01EvaluateFC,
    METH_VARARGS, 
    "Evaluate  the coefficients of nonlinear advective-diffusion transport in porous media with toy biodegradation system"}, 
  { "groundwaterBryantDawsonIonExEvaluateFC", 
    ctransportCoefficientsGroundwaterBryantDawsonIonExEvaluateFC,
    METH_VARARGS, 
    "Evaluate  the coefficients of nonlinear advective-diffusion transport in porous media with model ion-exchange system"}, 
  { "groundwaterTransportCoefficientsEvaluate_hetMat", 
      ctransportCoefficientsGroundwaterTransportCoefficientsEvaluate_hetMat,
      METH_VARARGS, 
      "Evaluate  the coefficients of linear advective-diffusion transport in porous media"}, 
  { "variablySaturatedGroundwaterTransportCoefficientsEvaluate_hetMat", 
      ctransportCoefficientsVariablySaturatedGroundwaterTransportCoefficientsEvaluate_hetMat,
      METH_VARARGS, 
      "Evaluate  the coefficients of linear advective-diffusion transport in variably saturated porous media"}, 
  { "variablySaturatedGroundwaterEnergyTransportCoefficientsEvaluate_hetMat", 
      ctransportCoefficientsVariablySaturatedGroundwaterEnergyTransportCoefficientsEvaluate_hetMat,
      METH_VARARGS, 
      "Evaluate  the coefficients of linear advective-diffusion transport in variably saturated porous media"}, 
  { "nonlinearADR_pqrstEvaluate", 
    ctransportCoefficientsNonlinearADR_pqrstEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar vtransport equation for simple monomial coefficients"}, 
  { "nonlinearADR_pqrstDualEvaluate", 
    ctransportCoefficientsNonlinearADR_pqrstDualEvaluate, 
    METH_VARARGS, 
    "Coefficients for doubly degenerate pqrst equation"},
  { "unitSquareRotationEvaluate", 
    ctransportCoefficientsUnitSquareRotationEvaluate, 
    METH_VARARGS, 
    "Coefficients for rotating velocity field on the unit square"}, 
  { "unitCubeRotationEvaluate", 
    ctransportCoefficientsUnitCubeRotationEvaluate, 
    METH_VARARGS, 
    "Coefficients for rotating velocity field on the unit cube"}, 
  { "unitSquareVortexEvaluate", 
    ctransportCoefficientsUnitSquareVortexEvaluate, 
    METH_VARARGS, 
    "Coefficients for vortex field on the unit square"}, 
  { "rotatingPulseVelEvaluate", 
    ctransportCoefficientsRotatingPulseVelEvaluate, 
    METH_VARARGS, 
    "coefficients advection-diffusion with a  rotating velocity field"}, 
  { "disRotatingPulseVelEvaluate", 
    ctransportCoefficientsDisRotatingPulseVelEvaluate, 
    METH_VARARGS, 
    "vtransport coefficients for a discontinuous rotating pulse"}, 
  { "disVelEvaluate", 
    ctransportCoefficientsDisVelEvaluate, 
    METH_VARARGS, 
    "vtransport coefficients for a discontinuous velocity field"}, 
  { "burgersDiagonalVelEvaluate", 
    ctransportCoefficientsBurgersDiagonalVelEvaluate, 
    METH_VARARGS, 
    "Burgers equations coefficients for a diagonal velocity field"}, 
  { "burgersDiagonalVelHJEvaluate", 
    ctransportCoefficientsBurgersDiagonalVelHJEvaluate, 
    METH_VARARGS, 
    "Burgers equations coefficients for a diagonal velocity field"}, 
  { "evaluateBuckleyLeverettLiuExample", 
    ctransportCoefficientsEvaluateBuckleyLeverettLiuExample, 
    METH_VARARGS, 
    "Simplified Buckley Leverett example with potential for velocity field"}, 
  { "twophasePotentialFlowEvaluate", 
    ctransportCoefficientsTwophasePotentialFlowEvaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients of a scalar vtransport equation for spatially variable coefficients"}, 
  { "twophasePotentialFlowUpdateFreeSurface", 
    ctransportCoefficientsTwophasePotentialFlowUpdateFreeSurface, 
    METH_VARARGS, 
    "Update  the coefficients of a scalar vtransport equation based on the location of a free surface"}, 
  { "twophaseLevelSetCoefficientsUpdateVelocity", 
    ctransportCoefficientsTwophaseLevelSetCoefficientsUpdateVelocity, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseLevelSetCoefficientsEvaluate", 
    ctransportCoefficientsTwophaseLevelSetCoefficientsEvaluate, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseLevelSetCoefficientsEvaluateCI", 
    ctransportCoefficientsTwophaseLevelSetCoefficientsEvaluateCI, 
    METH_VARARGS, 
    "Update  the level set velocity"}, 
  { "twophaseSignedDistanceUpdateSignFunction", 
    ctransportCoefficientsTwophaseSignedDistanceCoefficientsUpdateSignFunction, 
    METH_VARARGS, 
    "Update  the level set sign function"}, 
  { "twophaseSignedDistanceCoefficientsEvaluate", 
    ctransportCoefficientsTwophaseSignedDistanceCoefficientsEvaluate,
    METH_VARARGS, 
    "Update  the eikonal equation coefficients"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHomEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsL2projMualemVanGenuchtenHomEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsL2projMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of Richards' equation w/ L_2 proj"}, 
  { "conservativeHeadRichardsL2projBndMualemVanGenuchtenHomEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsL2projBndMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of Richards' equation w/ L_2 proj on element boundaries"}, 
  { "conservativeHeadRichardsL2projMualemVanGenuchtenHetEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsL2projMualemVanGenuchtenHetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of Richards' equation w/ L_2 proj"}, 
  { "conservativeTotalHeadRichardsMualemVanGenuchtenHomEvaluate", 
    ctransportCoefficientsConservativeTotalHeadRichardsMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation in terms of h=psi+z"}, 
  { "L2projectEvaluate", 
    ctransportCoefficientsL2projectEvaluate,
    METH_VARARGS, 
    "just project scalar quantity onto constants for a mesh entity"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHetEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsMualemVanGenuchten_sd_het", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchten_sd_het,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation for block heterogeneity with sparse diffusion rep for hydraulic conductivity"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "seepageBrezis", 
    ctransportCoefficientsSeepageBrezis,
    METH_VARARGS, 
    "evaluate the coefficients  of the seepage free boundary problem"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwind", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwind,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation with upwinding and harmonic average for Ks"}, 
  { "conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm_sd", 
    ctransportCoefficientsConservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2withUpwindAndHarm_sd,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation with upwinding and harmonic average for Ks"}, 
  { "conservativeSatRichardsMualemVanGenuchtenHomEvaluate", 
    ctransportCoefficientsConservativeSatRichardsMualemVanGenuchtenHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsBrooksCoreyBurdineHomEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsBrooksCoreyBurdineHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsBrooksCoreyBurdineHetEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsBrooksCoreyBurdineHetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
   { "conservativeSatRichardsBrooksCoreyBurdineHomEvaluate", 
    ctransportCoefficientsConservativeSatRichardsBrooksCoreyBurdineHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsBCBfromMVGHomEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsBCBfromMVGHomEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsJLeverettEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsJLeverettEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "conservativeHeadRichardsJLeverettAniEvaluate", 
    ctransportCoefficientsConservativeHeadRichardsJLeverettAniEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of Richards' equation"}, 
  { "Mass_2D_Evaluate", 
    ctransportCoefficientsMass_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the a 2D mass matrix"}, 
  { "Mass_3D_Evaluate", 
    ctransportCoefficientsMass_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of a 3D mass matrix"},
  { "TwoPhaseMass_2D_Evaluate",
    ctransportCoefficientsTwoPhaseMass_2D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 2D mass matrix"},
  { "TwoPhaseMass_3D_Evaluate",
    ctransportCoefficientsTwoPhaseMass_3D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 3D mass matrix"},
  { "TwoPhaseInvScaledMass_2D_Evaluate",
    ctransportCoefficientsTwoPhaseInvScaledMass_2D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 2D mass matrix"},
  { "TwoPhaseMass_mu_2D_Evaluate",
    ctransportCoefficientsTwoPhaseMass_mu_2D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 2D mass matrix"},
  { "TwoPhaseInvScaledMass_3D_Evaluate",
    ctransportCoefficientsTwoPhaseInvScaledMass_3D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 3D mass matrix"},
   { "TwoPhaseMass_mu_3D_Evaluate",
    ctransportCoefficientsTwoPhaseMass_mu_3D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of a two-phase 3D mass matrix"},
  { "Advection_2D_Evaluate", 
    ctransportCoefficientsAdvection_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Advection equations"}, 
  { "Advection_3D_Evaluate", 
    ctransportCoefficientsAdvection_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Advection equations"},
  { "TwoPhaseAdvection_2D_Evaluate", 
    ctransportCoefficientsTwoPhaseAdvection_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Advection equations"},
  { "TwoPhaseAdvection_3D_Evaluate", 
    ctransportCoefficientsTwoPhaseAdvection_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Advection equations"},
  { "NavierStokes_2D_Evaluate", 
    ctransportCoefficientsNavierStokes_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations"}, 
  { "NavierStokes_3D_Evaluate", 
    ctransportCoefficientsNavierStokes_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Navier-Stokes equations"}, 
  { "Stokes_2D_Evaluate", 
    ctransportCoefficientsStokes_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Stokes equations"},
  { "Laplace_2D_Evaluate",
    ctransportCoefficientsLaplace_2D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of the distcrete 2D Laplace operator"},
  { "Laplace_3D_Evaluate",
    ctransportCoefficientsLaplace_3D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients of the discrete 3D Laplace operator"},
  { "TwoPhaseInvScaledLaplace_2D_Evaluate",
    ctransportCoefficientsTwoPhaseInvScaledLaplace_2D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients for a two-phase inverse scaled 2D Laplace operator"},
   { "TwoPhaseInvScaledLaplace_3D_Evaluate",
    ctransportCoefficientsTwoPhaseInvScaledLaplace_3D_Evaluate,
    METH_VARARGS,
    "evaluate the coefficients for a two-phase inverse scaled 3D Laplace operator"},
  { "B_2D_Evaluate", 
    ctransportCoefficientsB_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the 2D B-operator"}, 
  { "B_3D_Evaluate", 
    ctransportCoefficientsB_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the 2D B-operator"}, 
  { "StokesP_2D_Evaluate", 
    ctransportCoefficientsStokesP_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Stokes equations"}, 
  { "StokesP_3D_Evaluate", 
    ctransportCoefficientsStokesP_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Stokes equations"}, 
  { "Stokes_3D_Evaluate", 
    ctransportCoefficientsStokes_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Stokes equations"}, 
  { "TwophaseNavierStokes_LS_SO_2D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_LS_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations"}, 
  { "TwophaseNavierStokes_ST_LS_SO_2D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations with surface  tension"}, 
  { "TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd", 
    ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations with surface  tension"}, 
  { "ThreephaseNavierStokes_ST_LS_SO_2D_Evaluate", 
    ctransportCoefficientsThreephaseNavierStokes_ST_LS_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations with surface  tension and a possible mobile rigid solid phase"}, 
  { "TwophaseNavierStokes_ST_LS_SO_3D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Navier-Stokes equations with surface  tension"}, 
  { "TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd", 
    ctransportCoefficientsTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Navier-Stokes equations with surface  tension"}, 
  { "ThreephaseNavierStokes_ST_LS_SO_3D_Evaluate", 
    ctransportCoefficientsThreephaseNavierStokes_ST_LS_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations with surface  tension and a possible mobile rigid solid phase"}, 
  { "TwophaseNavierStokes_LS_SO_3D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_LS_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Navier-Stokes equations"}, 
  { "TwophaseStokes_LS_SO_2D_Evaluate", 
    ctransportCoefficientsTwophaseStokes_LS_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Stokes equations"}, 
  { "TwophaseStokes_LS_SO_3D_Evaluate", 
    ctransportCoefficientsTwophaseStokes_LS_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Stokes equations"}, 
  { "TwophaseNavierStokes_VOF_SO_2D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_VOF_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Navier-Stokes equations"}, 
  { "TwophaseNavierStokes_VOF_SO_3D_Evaluate", 
    ctransportCoefficientsTwophaseNavierStokes_VOF_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Navier-Stokes equations"}, 
  { "TwophaseStokes_VOF_SO_2D_Evaluate", 
    ctransportCoefficientsTwophaseStokes_VOF_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 2D Stokes equations"}, 
  { "TwophaseStokes_VOF_SO_3D_Evaluate", 
    ctransportCoefficientsTwophaseStokes_VOF_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients  of the 3D Stokes equations"}, 
  { "constantVelocityLevelSetEvaluate", 
    ctransportCoefficientsConstantVelocityLevelSetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set equation with constant velocity (1,2,3)D"}, 
  { "constantNormalVelocityLevelSetEvaluate", 
    ctransportCoefficientsConstantNormalVelocityLevelSetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set equation with constant normal velocity (1,2,3)D"}, 
  { "unitSquareVortexLevelSetEvaluate", 
    ctransportCoefficientsUnitSquareVortexLevelSetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set equation with oscillating vortex 2D"}, 
  { "unitSquareRotationLevelSetEvaluate", 
    ctransportCoefficientsUnitSquareRotationLevelSetEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set equation with rotating velocity 2D"}, 
  { "HJBurgersEvaluate", 
    ctransportCoefficientsHJBurgersEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients Hamilton Jacobi version of Burgers equation (not same as noncons Burgers)"}, 
  { "eikonalEquationEvaluate", 
    ctransportCoefficientsEikonalEquationEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for Eikonal Equation"}, 
  { "ncLevelSetCoefficientsEvaluate", 
    ctransportCoefficients_ncLevelSetCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the non-conservative level set advection equation"}, 
  { "cLevelSetCoefficientsEvaluate", 
    ctransportCoefficients_cLevelSetCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the conservative level set advection equation"}, 
  { "VOFCoefficientsEvaluate", 
    ctransportCoefficients_VOFCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the conservative level set advection equation"}, 
  { "levelSetCurvatureCoefficientsEvaluate", 
    ctransportCoefficients_levelSetCurvatureCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the level set curvature equation"}, 
  { "redistanceLevelSetCoefficientsEvaluate", 
    ctransportCoefficients_redistanceLevelSetCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set redistancing equation"}, 
  { "redistanceLevelSetCoefficientsWithWeakPenaltyEvaluate", 
    ctransportCoefficients_redistanceLevelSetCoefficientsWithWeakPenaltyEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set redistancing equation including a weak penalty"}, 
  { "redistanceLevelSetSandFCoefficientsEvaluate", 
    ctransportCoefficients_redistanceLevelSetSandFCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for level set redistancing equation using Sussman and Fatemi volume correction"}, 
  { "setWeakDirichletConditionsForLevelSet",
    ctransportCoefficients_setWeakDirichletConditionsForLevelSet,
    METH_VARARGS, 
   "Determine which unknowns to freeze for the zero level set contour"},
  { "setSimpleWeakDirichletConditionsForLevelSet",
    ctransportCoefficients_setSimpleWeakDirichletConditionsForLevelSet,
    METH_VARARGS, 
   "Determine which unknowns to freeze for the zero level set contour"},
  { "darcySharpInterfaceFlowEvaluate", 
    ctransportCoefficients_darcySharpInterfaceFlowEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for two phase interface flow with Darcy's law for mom."}, 
  { "darcySharpInterfaceFlowImEvaluate", 
    ctransportCoefficients_darcySharpInterfaceFlowImEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients for two phase interface flow with Darcy's law for mom. and immersed interface type approx"}, 
/*   { "FractionalFlowPhaseForm_saturationEvaluate",  */
/*     ctransportCoefficientsFractionalFlowPhaseForm_saturationEvaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase fractional flow in the saturation phase"}, */
/*   { "FractionalFlowPhaseForm_potentialEvaluate",  */
/*     ctransportCoefficientsFractionalFlowPhaseForm_potentialEvaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase fractional flow in the potential phase"}, */
/*   { "FractionalFlowPhaseForm_saturationHetEvaluate",  */
/*     ctransportCoefficientsFractionalFlowPhaseForm_saturationHetEvaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase fractional flow in the saturation phase heterogeneous case"}, */
/*   { "FractionalFlowPhaseForm_potentialHetEvaluate",  */
/*     ctransportCoefficientsFractionalFlowPhaseForm_potentialHetEvaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase fractional flow in the potential phase heterogeneous case"},    */
/*   { "TwophaseDarcyFC_Evaluate", */
/*     ctransportCoefficientsTwophaseDarcyFC_Evaluate, */
/*     METH_VARARGS, */
/*     "evaluate the coefficients for two phase flow in fully coupled primitive form"}, */
/*   { "TwophaseFFDarcyFC_Evaluate", */
/*     ctransportCoefficientsTwophaseFFDarcyFC_Evaluate, */
/*     METH_VARARGS, */
/*     "evaluate the coefficients for two phase flow in fully coupled fractional flow form"}, */
/*   { "TwophaseDarcyFCHet_Evaluate",  */
/*     ctransportCoefficientsTwophaseDarcyFCHet_Evaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase flow in fully coupled primitive form heterogeneous case"},    */
/*   { "TwophaseDarcyFCHet_EvaluateV2",  */
/*     ctransportCoefficientsTwophaseDarcyFCHet_EvaluateV2, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase flow in fully coupled primitive form heterogeneous case"}, */
/*   { "TwophaseFFDarcyFCHet_EvaluateV2",  */
/*     ctransportCoefficientsTwophaseFFDarcyFCHet_EvaluateV2, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase flow in fully coupled, fractional flow version heterogeneous case"}, */
/*   { "TwophaseFFDarcyFCHet_Evaluate",  */
/*     ctransportCoefficientsTwophaseFFDarcyFCHet_Evaluate, */
/*     METH_VARARGS,  */
/*     "evaluate the coefficients for two phase flow in fully coupled fractional flow form heterogeneous case"},     */
   { "LinearElasticity_1D_Evaluate", 
    ctransportCoefficientsLinearElasticity_1D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 1D"},
   { "LinearElasticity_2D_Evaluate", 
    ctransportCoefficientsLinearElasticity_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 2D"},
   { "LinearElasticity_3D_Evaluate", 
    ctransportCoefficientsLinearElasticity_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 3D"},
   { "MovingMesh_1D_Evaluate", 
    ctransportCoefficientsMovingMesh_1D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 1D"},
   { "MovingMesh_2D_Evaluate", 
    ctransportCoefficientsMovingMesh_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 2D"},
   { "MovingMesh_3D_Evaluate", 
    ctransportCoefficientsMovingMesh_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the linear elasticity equation in 3D"},
   { "levelSetConservationCoefficientsEvaluate", 
    ctransportCoefficientsLevelSetConservationCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the volume conservation correction for level set methods"},
   { "levelSetConservationCoefficientsEvaluate_sd", 
    ctransportCoefficientsLevelSetConservationCoefficientsEvaluate_sd,
    METH_VARARGS, 
    "evaluate the coefficients of the volume conservation correction for level set methods"},
  { "VolumeAveragedNavierStokesFullDevStress_2D_Evaluate", 
    ctransportCoefficientsVolumeAveragedNavierStokesFullDevStress_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the volume-averaged 2D Navier-Stokes equations"}, 
  { "VolumeAveragedNavierStokesFullDevStress_3D_Evaluate", 
    ctransportCoefficientsVolumeAveragedNavierStokesFullDevStress_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the volume-averaged 3D Navier-Stokes equations"}, 
  { "VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate", 
    ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the Volume-Averaged 2D Navier-Stokes equations with surface  tension"}, 
  { "VolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd", 
    ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd,
    METH_VARARGS, 
    "evaluate the coefficients of the Volume-Averaged 2D Navier-Stokes equations with surface  tension"}, 
  { "VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate", 
    ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the Volume-Averaged 3D Navier-Stokes equations with surface  tension"}, 
  { "VolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd", 
    ctransportCoefficientsVolumeAveragedTwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd,
    METH_VARARGS, 
    "evaluate the coefficients of the Volume-Averaged 3D Navier-Stokes equations with surface  tension"}, 
  { "VolumeAveragedVOFCoefficientsEvaluate", 
    ctransportCoefficients_VolumeAveragedVOFCoefficientsEvaluate,
    METH_VARARGS, 
    "evaluate the coefficients of the conservative level set advection equation with variable porosity"}, 
  { "kEpsilon_2D_Evaluate", 
    ctransportCoefficients_kEpsilon_2D_Evaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_2D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_2D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_k_2D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_k_2D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for k assuming lagged epsilon for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_epsilon_2D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_epsilon_2D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for epsilon assuming lagged epsilon for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_3D_Evaluate", 
    ctransportCoefficients_kEpsilon_3D_Evaluate, 
    METH_VARARGS, 
    "Evaluate  the coefficients for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_3D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_3D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for standard incompressible flow k-epsilon model 3D"}, 
  { "kEpsilon_k_3D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_k_3D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for k assuming lagged epsilon for standard incompressible flow k-epsilon model"}, 
  { "kEpsilon_epsilon_3D_Evaluate_sd", 
    ctransportCoefficients_kEpsilon_epsilon_3D_Evaluate_sd, 
    METH_VARARGS, 
    "Evaluate  the coefficients for epsilon assuming lagged epsilon for standard incompressible flow k-epsilon model"}, 
  { "ReynoldsAveragedNavierStokes_kEpsilon_2D_Update", 
    ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_2D_Update, 
    METH_VARARGS, 
    "Update the NS coefficients for RANS terms in standard incompressible flow k-epsilon model"}, 
  { "ReynoldsAveragedNavierStokes_kEpsilon_2D_Update_sd", 
    ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_2D_Update_sd, 
    METH_VARARGS, 
    "Update the NS coefficients for RANS terms in standard incompressible flow k-epsilon model"}, 
  { "ReynoldsAveragedNavierStokes_kEpsilon_3D_Update", 
    ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_3D_Update, 
    METH_VARARGS, 
    "Update the NS coefficients for RANS terms in standard incompressible flow k-epsilon model"}, 
  { "ReynoldsAveragedNavierStokes_kEpsilon_3D_Update_sd", 
    ctransportCoefficientsReynoldsAveragedNavierStokes_kEpsilon_3D_Update_sd, 
    METH_VARARGS, 
    "Update the NS coefficients for RANS terms in standard incompressible flow k-epsilon model"}, 
  { "eddyViscosity_2D_Update", 
    ctransportCoefficientsEddyViscosity_2D_Update, 
    METH_VARARGS, 
    "Update the NS coefficients with eddy viscosity from some closure model"}, 
  { "eddyViscosity_2D_Update_sd", 
    ctransportCoefficientsEddyViscosity_2D_Update_sd, 
    METH_VARARGS, 
    "Update the NS coefficients with eddy viscosity from some closure model"}, 
   { "eddyViscosity_3D_Update", 
    ctransportCoefficientsEddyViscosity_3D_Update, 
    METH_VARARGS, 
    "Update the NS coefficients with eddy viscosity from some closure model"}, 
  { "eddyViscosity_3D_Update_sd", 
    ctransportCoefficientsEddyViscosity_3D_Update_sd, 
    METH_VARARGS, 
    "Update the NS coefficients with eddy viscosity from some closure model"}, 
 { "calculateEddyViscosity_Smagorinsky_2D", 
    ctransportCoefficientsCalculateEddyViscosity_Smagorinsky_2D, 
    METH_VARARGS, 
    "Calculate simple Smagorinsky eddy viscosity"}, 
 { "calculateEddyViscosity_Smagorinsky_3D", 
    ctransportCoefficientsCalculateEddyViscosity_Smagorinsky_3D, 
    METH_VARARGS, 
    "Calculate simple Smagorinsky eddy viscosity"}, 
 { "calculateEddyViscosity_Smagorinsky2P_2D", 
    ctransportCoefficientsCalculateEddyViscosity_Smagorinsky2P_2D, 
    METH_VARARGS, 
    "Calculate simple Smagorinsky eddy viscosity with smoothed Heaviside for 2p calculation"}, 
 { "calculateEddyViscosity_Smagorinsky2P_3D", 
    ctransportCoefficientsCalculateEddyViscosity_Smagorinsky2P_3D, 
    METH_VARARGS, 
    "Calculate simple Smagorinsky eddy viscosity with smoothed Heaviside for 2p calculation"}, 
  { "scriptedSphereMotionSignedDistance", 
    ctransportCoefficientsScriptedSphereMotionSignedDistance, 
    METH_VARARGS, 
    "evaluate distance to a set of spheres with varying radii and distributed in space according to centers"}, 
  { "shallowWater_1D_Evaluate", 
    ctransportCoefficientsShallowWater_1D_Evaluate, 
    METH_VARARGS, 
    "Evaluate the coefficients of the shallow water equations in 1D"}, 
  { "shallowWater_2D_Evaluate", 
    ctransportCoefficientsShallowWater_2D_Evaluate, 
    METH_VARARGS, 
    "Evaluate the coefficients of the shallow water equations in 2D"}, 
  { "smoothedHeaviside", 
    ctransportCoefficientsSmoothedHeaviside, 
    METH_VARARGS, 
    "Evaluate smoothed Heaviside function with smearing eps"}, 
  { "smoothedHeaviside_integral", 
    ctransportCoefficientsSmoothedHeaviside_integral, 
    METH_VARARGS, 
    "Evaluate the definite integral of the smoothed Heaviside function w.r.t. phi, with smearing eps, from -eps to phi"}, 
  { "smoothedDirac", 
    ctransportCoefficientsSmoothedDirac, 
    METH_VARARGS, 
    "Evaluate smoothed Dirac delta function with smearing eps"}, 
  { "applyContactLineSlip", 
    ctransportCoefficientsApplyContactLineSlip, 
    METH_VARARGS, 
    "Apply contact line slip condition"}, 
  { "applyContactLineSlipJacobian", 
    ctransportCoefficientsApplyContactLineSlipJacobian, 
    METH_VARARGS, 
    "Apply contact line slip condition to Jacobian"}, 
  { "diffusiveWave1DCoefficientsEvaluate",
    ctransportCoefficientsDiffusiveWave1DCoefficientsEvaluate,
    METH_VARARGS, 
    "Evaluate the coefficients of the 1D diffusive wave equation"}, 
  { "diffusiveWave2DCoefficientsEvaluate",
    ctransportCoefficientsDiffusiveWave2DCoefficientsEvaluate,
    METH_VARARGS, 
    "Evaluate the coefficients of the 2D diffusive wave equation"},
  { "calculateWaveFunction3d_ref",
    ctransportCoefficientsCalculateWaveFunction3d_ref,
    METH_VARARGS,
    "generate continuity source term for driving waves following Liu etal"},
  { NULL,NULL,0,NULL}
};


PyMODINIT_FUNC initctransportCoefficients(void)
{
  PyObject *m,*d;
  m = Py_InitModule("ctransportCoefficients", ctransportCoefficientsMethods);
  d = PyModule_GetDict(m);
  import_array();
}















/** @} */


