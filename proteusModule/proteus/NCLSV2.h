#ifndef NCLSV2_H
#define NCLSV2_H
#include <cmath>
#include <iostream>

#define NCLSV2_SPACE_DIM 3
#if (NCLSV2_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 5
 #define nDOF_mesh_trial_element 4
 #define nDOF_trial_element 4
 #define nDOF_test_element 4
 #define nDOF_test_X_trial_element 16
 #define nQuadraturePoints_elementBoundary 4
 #define NCLSV2_NAME NCLSV23D
 #define NCLSV2_RES calculateResidual_NCLSV23D
 #define NCLSV2_JAC calculateJacobian_NCLSV23D
#else
 #if (NCLSV2_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 6
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 4
  #define NCLSV2_NAME NCLSV22D
  #define NCLSV2_RES calculateResidual_NCLSV22D
  #define NCLSV2_JAC calculateJacobian_NCLSV22D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 2
  #define nDOF_test_element 2
  #define nDOF_test_X_trial_element 4
  #define nQuadraturePoints_elementBoundary 1
  #define NCLSV2_NAME NCLSV21D
  #define NCLSV2_RES calculateResidual_NCLSV21D
  #define NCLSV2_JAC calculateJacobian_NCLSV21D
 #endif
#endif

namespace NCLSV2_NAME
{
  inline void evaluateCoefficients(const double v[nSpace],
                                   const double& u,
                                   const double grad_u[nSpace],
                                   double& m,
                                   double& dm,
                                   double& H,
                                   double dH[nSpace])
  {
    int I;
    m = u;
    dm=1.0;
    H = 0.0;
    for (I=0; I < nSpace; I++)
      {
        H += v[I]*grad_u[I];
        dH[I] = v[I];
      }
  }
  
  inline void calculateSubgridError_tau(const double& elementDiameter,
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
    oneByAbsdt =  fabs(dmt);
    tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
  }  
}//NCLSV2

extern "C"
{
  //void calculateResidual_NCLSV2(int nElements_global,
  void NCLSV2_RES (//element
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
                   int lag_shockCapturing, /*mwf not used yet*/
                   double shockCapturingDiffusion,
                   int* u_l2g, 
                   double* elementDiameter,
                   double* u_dof,
                   double* u_trial, 
                   double* u_grad_trial, 
                   double* u_test_dV, 
                   double* u_grad_test_dV, 
                   double* velocity,
                   double* q_m,
                   double* q_u,
                   double* q_dH,
                   double* q_m_last,
                   double* q_cfl,
                   double* q_numDiff_u, 
                   double* q_numDiff_u_last, 
                   double* q_elementResidual_u, 
                   int offset_u, int stride_u, 
                   double* globalResidual,
                   int nExteriorElementBoundaries_global,
                   int* exteriorElementBoundariesArray,
                   int* elementBoundaryElementsArray,
                   int* elementBoundaryLocalElementBoundariesArray,
                   double* u_trial_ext,
                   double* u_grad_trial_ext,
                   double* ebqe_velocity_ext,
                   double* ebqe_n_ext,
                   int* isDOFBoundary_u,
                   double* ebqe_bc_u_ext,
                   double* u_test_dS_ext,
                   double* ebqe_u);
  
  //  void calculateJacobian_NCLSV2(int nElements_global,
  void NCLSV2_JAC (//element
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
                 int lag_shockCapturing,/*mwf not used yet*/
                 double shockCapturingDiffusion,
                 int* u_l2g,
                 double* elementDiameter,
                 double* u_dof, 
                 double* u_trial, 
                 double* u_grad_trial, 
                 double* u_test_dV, 
                 double* u_grad_test_dV, 
                 double* velocity,
                 double* q_m_last, 
                 double* q_cfl,
                 double* q_numDiff_u_last, 
                 int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                 double* globalJacobian,
                 int nExteriorElementBoundaries_global,
                 int* exteriorElementBoundariesArray,
                 int* elementBoundaryElementsArray,
                 int* elementBoundaryLocalElementBoundariesArray,
                 double* u_trial_ext,
                 double* u_grad_trial_ext,
                 double* ebqe_velocity_ext,
                 double* ebqe_n,
                 int* isDOFBoundary_u,
                 double* ebqe_bc_u_ext,
                 double* u_test_dS_ext,
                 int* csrColumnOffsets_eb_u_u);
}//extern "C"
#endif
