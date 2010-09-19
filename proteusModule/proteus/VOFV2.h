#ifndef VOFV2_H
#define VOFV2_H
#include <cmath>
#include <iostream>

#define VOFV2_SPACE_DIM 3
#if (VOFV2_SPACE_DIM == 3)
#define nSpace 3
#define nQuadraturePoints_element 5
#define nDOF_mesh_trial_element 4
#define nDOF_trial_element 4
#define nDOF_test_element 4
#define nDOF_test_X_trial_element 16
#define nQuadraturePoints_elementBoundary 4
#define VOFV2_NAME VOFV23D
#define VOFV2_RES calculateResidual_VOFV23D
#define VOFV2_JAC calculateJacobian_VOFV23D
#else
#if (VOFV2_SPACE_DIM == 2)
#define nSpace 2
#define nQuadraturePoints_element 6
#define nDOF_trial_element 3
#define nDOF_test_element 3
#define nDOF_test_X_trial_element 9
#define nQuadraturePoints_elementBoundary 4
#define VOFV2_NAME VOFV22D
#define VOFV2_RES calculateResidual_VOFV22D
#define VOFV2_JAC calculateJacobian_VOFV22D
#else 
#define nSpace 1
#define nQuadraturePoints_element 3
#define nDOF_trial_element 2
#define nDOF_test_element 2
#define nDOF_test_X_trial_element 4
#define nQuadraturePoints_elementBoundary 1
#define VOFV2_NAME VOFV21D
#define VOFV2_RES calculateResidual_VOFV21D
#define VOFV2_JAC calculateJacobian_VOFV21D
#endif
#endif

namespace VOFV2_NAME
{
  inline
    void evaluateCoefficients(const double v[nSpace],
				const double& u,
				double& m,
				double& dm,
				double f[nSpace],
				double df[nSpace])
  {
    m = u;
    dm = 1.0;
    for (int I=0; I < nSpace; I++)
      {
	f[I] = v[I]*u;
	df[I] = v[I];
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
    oneByAbsdt =  fabs(dmt);
    tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
  }

  inline 
  void calculateNumericalDiffusion(const double& shockCapturingDiffusion,
				       const double& elementDiameter,
				       const double& strong_residual,
				       const double grad_u[nSpace],
				       double& numDiff)
  {
    double h,
      num,
      den,
      n_grad_u;
    h = elementDiameter;
    n_grad_u = 0.0;
    for (int I=0;I<nSpace;I++)
      n_grad_u += grad_u[I]*grad_u[I];
    num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
    den = sqrt(n_grad_u) + 1.0e-8;
    numDiff = num/den;
  }

  inline
    void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_u,
					  const int& isFluxBoundary_u,
					  const double n[nSpace],
					  const double& bc_u,
					  const double& bc_flux_u,
					  const double& u,
					  const double velocity[nSpace],
					  double& flux)
  {

    double flow=0.0;
    for (int I=0; I < nSpace; I++)
      flow += n[I]*velocity[I];
    //std::cout<<" isDOFBoundary_u= "<<isDOFBoundary_u<<" flow= "<<flow<<std::endl;
    if (isDOFBoundary_u == 1)
      {
	//std::cout<<"Dirichlet boundary u and bc_u "<<u<<'\t'<<bc_u<<std::endl;
	if (flow >= 0.0)
	  {
	    flux = u*flow;
	    //flux = flow;
	  }
	else
	  {
	    flux = bc_u*flow;
	    //flux = flow;
	  }
      }
    else if (isFluxBoundary_u == 1)
      {
	flux = bc_flux_u;
	//std::cout<<"Flux boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
      }
    else
      {
	//std::cout<<"No BC boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
	if (flow >= 0.0)
	  {
	    flux = u*flow;
	  }
	else
	  {
	    std::cout<<"warning: open boundary with no external trace, setting to zero for inflow"<<std::endl;
	    flux = 0.0;
	  }

      }
    //flux = flow;
    //std::cout<<"flux error "<<flux-flow<<std::endl;
    //std::cout<<"flux in computationa"<<flux<<std::endl;
  }

  inline
  void exteriorNumericalAdvectiveFluxDerivative(const int& isDOFBoundary_u,
                                                const int& isFluxBoundary_u,
                                                const double n[nSpace],
                                                const double velocity[nSpace],
                                                double& dflux)
  {
    double flow=0.0;
    for (int I=0; I < nSpace; I++)
      flow += n[I]*velocity[I];
    //double flow=n[0]*velocity[0]+n[1]*velocity[1]+n[2]*velocity[2];
    dflux=0.0;//default to no flux
    if (isDOFBoundary_u == 1)
      {
	if (flow >= 0.0)
	  {
	    dflux = flow;
	  }
	else
	  {
	    dflux = 0.0;
	  }
      }
    if (isFluxBoundary_u == 1)
      {
	dflux = 0.0;
      }
  }

}//VOFV2 namespace
extern "C"
{
  //void calculateResidual_VOFV2(int nElements_global,
  void VOFV2_RES (//element
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
                          int nElements_global,
		double alphaBDF,
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
		double* q_m_betaBDF,
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
		int* isFluxBoundary_u,
		double* ebqe_bc_flux_u_ext,
		double* u_test_dS_ext,
                double* phi,double eps,
		double* ebqe_u,
		double* ebqe_flux);
  //void calculateJacobian_VOFV2(int nElements_global,
  void VOFV2_JAC(//element
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
	       double alphaBDF,
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
	       double* q_m_betaBDF, 
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
	       int* isFluxBoundary_u,
	       double* ebqe_bc_flux_u_ext,
	       double* u_test_dS_ext,
	       int* csrColumnOffsets_eb_u_u);
}//extern "C"
#endif
