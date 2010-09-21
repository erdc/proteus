#ifndef RDLSV2_H
#define RDLSV2_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"

#define RDLSV2_SPACE_DIM 3
#if (RDLSV2_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 5
 #define nDOF_mesh_trial_element 4
 #define nDOF_trial_element 4
 #define nDOF_test_element 4
 #define nDOF_test_X_trial_element 16
 #define nQuadraturePoints_elementBoundary 4
 #define RDLSV2_NAME RDLSV23D
 #define RDLSV2_RES calculateResidual_RDLSV23D
 #define RDLSV2_JAC calculateJacobian_RDLSV23D
#else
 #if (RDLSV2_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 6
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 4
  #define RDLSV2_NAME RDLSV22D
  #define RDLSV2_RES calculateResidual_RDLSV22D
  #define RDLSV2_JAC calculateJacobian_RDLSV22D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 2
  #define nDOF_test_element 2
  #define nDOF_test_X_trial_element 4
  #define nQuadraturePoints_elementBoundary 1
  #define RDLSV2_NAME RDLSV21D
  #define RDLSV2_RES calculateResidual_RDLSV21D
  #define RDLSV2_JAC calculateJacobian_RDLSV21D
 #endif
#endif
namespace RDLSV2_NAME
{
  inline
  double smoothedHeaviside(double eps, double phi)
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

inline
void evaluateCoefficients(const double& eps,
			    const double& u_levelSet,
			    const double& u,
			    const double grad_u[nSpace],
			    double& m,
			    double& dm,
			    double& H,
			    double dH[nSpace],
			    double& r)
{
  int I;
  double normGradU=0.0,Si=0.0;
  m = u;
  dm=1.0;
  H = 0.0;
  Si= -1.0+2.0*smoothedHeaviside(eps,u_levelSet);
  r = -Si;
  for (I=0; I < nSpace; I++)
    {
      normGradU += grad_u[I]*grad_u[I];
    }
  normGradU = sqrt(normGradU);
  H = Si*normGradU;
  for (I=0; I < nSpace; I++)
    {
      dH[I] = Si*grad_u[I]/(normGradU+1.0e-12);
    }
}

inline
void calculateSubgridError_tau(const double& elementDiameter,
                               const double& dmt,
                               const double dH[nSpace],
                               double& cfl,
                               double& tau)
{
  double h,nrm_v,oneByAbsdt;
  h = elementDiameter;
  nrm_v=0.0;
  for(int I=0;I<nSpace;I++)
    nrm_v+=dH[I]*dH[I];
  nrm_v = sqrt(nrm_v);
  cfl = nrm_v/h;
  oneByAbsdt =  dmt;
  //'1'
  //tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
  //'2'
  tau = 1.0/sqrt(4.0*nrm_v*nrm_v/(h*h) + oneByAbsdt*oneByAbsdt + 1.0e-8);
  //mwf debug
  //std::cout<<"tau calc h= "<<h<<" dmt= "<<dmt<<" dH[0]= "<<dH[0]<<" cfl= "<<cfl<<" tau= "<<tau<<std::endl;
  
}
}//RDLSV2 namespace

extern "C"
{
  //void calculateResidual_RDLSV2(int nElements_global,
  void RDLSV2_RES(//element
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
                  double dt,
                  double epsilon_redist,
                  int freezeLevelSet,
                  int useTimeIntegration,
                  int lag_shockCapturing, 
                  int lag_subgridError, //0 nothing lagged
                  //1 dH lagged in tau
                  //2 dH lagged in tau and Residual, adjoint calculations
                  double shockCapturingDiffusion,
                  int* u_l2g, 
                  double* elementDiameter,
                  double* u_dof,
                  double* u_trial, 
                  double* u_grad_trial, 
                  double* u_test_dV, 
                  double* u_grad_test_dV, 
                  double* phi_ls,
                  double* q_m,
                  double* q_u,
                  double* q_dH,
                  double* u_weak_internal_bc_dofs,
                  double* q_m_last,
                  double* q_dH_last,
                  double* q_cfl,
                  double* q_numDiff_u, 
                  double* q_numDiff_u_last, 
                  double* q_elementResidual_u,
                  //mwf for debugging
                  double* q_m_t,
                  double* q_r,
                  double* q_subgridError,
                  double* q_Lstar,
                  double* q_tau_last,
                  //mwf end debugging
                  int * weakDirichletConditionFlags,
                  int offset_u, int stride_u, 
                  double* globalResidual,
                  int nExteriorElementBoundaries_global,
                  int* exteriorElementBoundariesArray,
                  int* elementBoundaryElementsArray,
                  int* elementBoundaryLocalElementBoundariesArray,
                  double* u_trial_ext,
                  double* u_grad_trial_ext,
                  double* ebqe_phi_ls_ext,
                  double* ebqe_n_ext,
                  int* isDOFBoundary_u,
                  double* ebqe_bc_u_ext,
                  double* u_test_dS_ext,
                  double* ebqe_u);
  
  //void calculateJacobian_RDLSV2(int nElements_global,
  void RDLSV2_JAC(//element
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
                  double dt,
                  double epsilon_redist,
                  int freezeLevelSet,
                  int useTimeIntegration,
                  int lag_shockCapturing, 
                  int lag_subgridError, 
                  double shockCapturingDiffusion,
                  int* u_l2g,
                  double* elementDiameter,
                  double* u_dof, 
                  double* u_trial, 
                  double* u_grad_trial, 
                  double* u_test_dV, 
                  double* u_grad_test_dV, 
                  double* phi_ls,
                  double* q_m_last,
                  double* q_dH_last,
                  double* q_cfl,
                  double* q_numDiff_u,
                  double* q_numDiff_u_last,
                  int* weakDirichletConditionFlags,
                  int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                  double* globalJacobian,
                  int nExteriorElementBoundaries_global,
                  int* exteriorElementBoundariesArray,
                  int* elementBoundaryElementsArray,
                  int* elementBoundaryLocalElementBoundariesArray,
                  double* u_trial_ext,
                  double* u_grad_trial_ext,
                  double* ebqe_phi_ls_ext,
                  double* ebqe_n,
                  int* isDOFBoundary_u,
                  double* ebqe_bc_u_ext,
                  double* u_test_dS_ext,
                  int* csrColumnOffsets_eb_u_u);

}//extern "C"
#endif
