#ifndef MCorrV2_H
#define MCorrV2_H
#include "CompKernel.h"
#include <cmath>
#include <iostream>

#define MCorrV2_SPACE_DIM 3
#if (MCorrV2_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 5
 #define nDOF_mesh_trial_element 4
 #define nDOF_trial_element 4
 #define nDOF_test_element 4
 #define nDOF_test_X_trial_element 16
 #define nQuadraturePoints_elementBoundary 4
 #define MCORRV2_NAME MCorrV23D
 #define MCORRV2_RES calculateResidual_MCorrV23D
 #define MCORRV2_JAC calculateJacobian_MCorrV23D
#else
 #if (MCorrV2_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 6
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 4
  #define MCORRV2_NAME MCorrV22D
  #define MCORRV2_RES calculateResidual_MCorrV22D
  #define MCORRV2_JAC calculateJacobian_MCorrV22D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 2
  #define nDOF_test_element 2
  #define nDOF_test_X_trial_element 4
  #define nQuadraturePoints_elementBoundary 1
  #define MCORRV2_NAME MCorrV21D
  #define MCORRV2_RES calculateResidual_MCorrV21D
  #define MCORRV2_JAC calculateJacobian_MCorrV21D
 #endif
#endif

namespace MCORRV2_NAME // MCorrV2
{
inline double smoothedHeaviside(double eps, double phi)
{
  double H;
  if (phi > eps)
    H=1.0;
  else if (phi < -eps)
    H=0.0;
  else if (phi==0.0)
    H=0.5;
  else
    H = 0.5*(1.0 + phi/eps + sin(M_PI*phi/eps)/M_PI);
  return H;
}

inline double smoothedHeaviside_integral(double eps, double phi)
{
  double HI;
  if (phi > eps)
    {
      HI= phi - eps + 	0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  else if (phi < -eps)
    {
      HI=0.0;
    }
  else
    {
      HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  return HI;
}

inline double smoothedDirac(double eps, double phi)
{
  double d;
  if (phi > eps)
    d=0.0;
  else if (phi < -eps)
    d=0.0;
  else
    d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
  return d;
}
inline
void evaluateCoefficients(const double& epsHeaviside,
                          const double& epsDirac,
                          const double& phi,
                          const double& H,
                          const double& u,
                          double& r,
                          double& dr)
{
  r = smoothedHeaviside(epsHeaviside,u + phi) - H;
  dr = smoothedDirac(epsDirac,u + phi);
}

// inline
// double NumericalDiffusion_c(const double& numDiff,
// 			    const double grad_u[nSpace],
// 			    const double grad_w_dV[nSpace])
// {
//   double tmp=0.0;
//   for (int I=0;I<nSpace;I++)
//     tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
//   return tmp;
// }

// inline
// double NumericalDiffusionJacobian_c(const double& numDiff,
// 					const double grad_v[nSpace],
// 					const double grad_w_dV[nSpace])
// {
//   double tmp=0.0;
//   for (int I=0;I<nSpace;I++)
//     tmp += numDiff*grad_v[I]*grad_w_dV[I];
//   return tmp;
// }

}//MCorrV2 namespace
extern "C"
{
  //void calculateResidual_MCorrV2(int nElements_global,
  void MCORRV2_RES (//element
                            double* mesh_trial_ref,
                            double* mesh_grad_trial_ref,
                            double* mesh_dof,
                            int* mesh_l2g,
                            double* dV_ref,
                            double* u_trial_ref,
                            double* u_grad_trial_ref,
                            double* u_test_ref,
                            double* u_grad_test_ref,
                            //element boundary
                            double* mesh_trial_trace_ref,
                            double* mesh_grad_trial_trace_ref,
                            double* dS_ref,
                            double* u_trial_trace_ref,
                            double* u_grad_trial_trace_ref,
                            double* u_test_trace_ref,
                            double* u_grad_test_trace_ref,
                            double* normal_ref,
                            double* boundaryJac_ref,
                            //physics
                            int nElements_global,
			 double epsHeaviside,
			 double epsDirac,
			 double epsDiffusion,
			 int* u_l2g, 
			 double* elementDiameter,
			 double* u_dof,
			 double* u_trial, 
			 double* u_grad_trial, 
			 double* u_test_dV, 
			 double* u_grad_test_dV, 
			 double* q_phi,
			 double* q_H,
			 double* q_u,
			 double* q_r,
			 int offset_u, int stride_u, 
			 double* globalResidual);
  //void calculateJacobian_MCorrV2(int nElements_global,
  void MCORRV2_JAC(//element
                          double* mesh_trial_ref,
                          double* mesh_grad_trial_ref,
                          double* mesh_dof,
                          int* mesh_l2g,
                          double* dV_ref,
                          double* u_trial_ref,
                          double* u_grad_trial_ref,
                          double* u_test_ref,
                          double* u_grad_test_ref,
                          //element boundary
                          double* mesh_trial_trace_ref,
                          double* mesh_grad_trial_trace_ref,
                          double* dS_ref,
                          double* u_trial_trace_ref,
                          double* u_grad_trial_trace_ref,
                          double* u_test_trace_ref,
                          double* u_grad_test_trace_ref,
                          double* normal_ref,
                          double* boundaryJac_ref,
                          //physics
                          int nElements_global,
			 double epsHeaviside,
			 double epsDirac,
			 double epsDiffusion,
			 int* u_l2g,
			 double* elementDiameter,
			 double* u_dof, 
			 double* u_trial, 
			 double* u_grad_trial, 
			 double* u_test_dV, 
			 double* u_grad_test_dV, 
			 double* q_phi,
			 double* q_H,
			 int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			 double* globalJacobian);
}//extern "C"
#endif
