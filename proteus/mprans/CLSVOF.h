#ifndef CLSVOF_H
#define CLSVOF_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

// True characteristic functions
#define heaviside(z) (z>0 ? 1. : (z<0 ? 0. : 0.5))
#define sign(z) (z>0 ? 1. : (z<0 ? -1. : 0.))

// These functions are used just on experimental code
#define betaNormGrad(gradu2,beta2) std::sqrt(gradu2+beta2)
#define rabs(z) std::sqrt(z*z+1E-15)
#define ISOTROPIC_REGULARIZATION 1 // to use with ..._experimental 

namespace proteus
{
  class CLSVOF_base
  {
    //The base class defining the interface
  public:
    virtual ~CLSVOF_base(){}
    virtual void calculateResidual_Monolithic(//element
					      double dt,
					      double* mesh_trial_ref,
					      double* mesh_grad_trial_ref,
					      double* mesh_dof,
					      double* mesh_velocity_dof,
					      double MOVING_DOMAIN,
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
					      double useMetrics,
					      double alphaBDF,
					      //VRANS
					      const double* q_porosity,
					      const double* porosity_dof, /////
					      //
					      int* u_l2g, 
					      double* elementDiameter,
					      double* nodeDiametersArray,
					      int degree_polynomial,
					      double* u_dof,
					      double* u_dof_old,
					      double* velocity,
					      double* velocity_old,
					      double* q_m,
					      double* q_u,
					      double* q_m_betaBDF,
					      double* q_dV,
					      double* q_dV_last,
					      double* cfl,
					      int offset_u, int stride_u, 
					      double* globalResidual,
					      int nExteriorElementBoundaries_global,
					      int* exteriorElementBoundariesArray,
					      int* elementBoundaryElementsArray,
					      int* elementBoundaryLocalElementBoundariesArray,
					      double* ebqe_velocity_ext,
					      //VRANS
					      const double* ebqe_porosity_ext,
					      //
					      int* isDOFBoundary_u,
					      double* ebqe_bc_u_ext,
					      int* isFluxBoundary_u,
					      double* ebqe_bc_flux_u_ext,
					      double* ebqe_phi,double epsFact,
					      double* ebqe_u,
					      double* ebqe_flux,
					      // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					      int timeOrder,
					      int timeStage,
					      int useFullNewton,      
					      double epsFactHeaviside,
					      double epsFactDirac,
					      double epsFactDiffusion,
					      double* phin_dof,
					      double* phiHat_dof,
					      // interface locator
					      double* norm_factor,
					      double norm_factor_lagged,
					      double* interface_locator,
					      double* interface_locator_lagged,
					      // normal reconstruction
					      double* lumped_wx,
					      double* lumped_wy,
					      double* lumped_wz,
					      double* lumped_wx_tStar,
					      double* lumped_wy_tStar,
					      double* lumped_wz_tStar,
					      // AUX QUANTITIES OF INTEREST
					      double* quantDOFs)=0;
    virtual void calculateResidual_non_Monolithic(//element
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
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
						  double useMetrics,
						  double alphaBDF,
						  //VRANS
						  const double* q_porosity,
						  const double* porosity_dof,
						  //
						  int* u_l2g, 
						  double* elementDiameter,
						  double* nodeDiametersArray,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* velocity,
						  double* velocity_old,
						  double* q_m,
						  double* q_u,
						  double* q_m_betaBDF,
						  double* q_dV,
						  double* q_dV_last,
						  double* cfl,
						  int offset_u, int stride_u, 
						  double* globalResidual,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  double* ebqe_phi,double epsFact,
						  double* ebqe_u,
						  double* ebqe_flux,
						  // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
						  int timeOrder,
						  int timeStage,
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  // interface locator
						  double* norm_factor,
						  double norm_factor_lagged,
						  double* interface_locator,
						  double* interface_locator_lagged,
						  // normal reconstruction
						  double* lumped_wx,
						  double* lumped_wy,
						  double* lumped_wz,
						  double* lumped_wx_tStar,
						  double* lumped_wy_tStar,
						  double* lumped_wz_tStar,
						  // AUX QUANTITIES OF INTEREST
						  double* quantDOFs)=0;
    virtual void calculateResidual_experimental(//element 
						double dt,
						double* mesh_trial_ref,
						double* mesh_grad_trial_ref,
						double* mesh_dof,
						double* mesh_velocity_dof,
						double MOVING_DOMAIN,
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
						double useMetrics,
						double alphaBDF,
						//VRANS
						const double* q_porosity,
						const double* porosity_dof,
						//
						int* u_l2g, 
						double* elementDiameter,
						double* nodeDiametersArray,
						int degree_polynomial,
						double* u_dof,
						double* u_dof_old,
						double* velocity,
						double* velocity_old,
						double* q_m,
						double* q_u,
						double* q_m_betaBDF,
						double* q_dV,
						double* q_dV_last,
						double* cfl,
						int offset_u, int stride_u, 
						double* globalResidual,
						int nExteriorElementBoundaries_global,
						int* exteriorElementBoundariesArray,
						int* elementBoundaryElementsArray,
						int* elementBoundaryLocalElementBoundariesArray,
						double* ebqe_velocity_ext,
						//VRANS
						const double* ebqe_porosity_ext,
						//
						int* isDOFBoundary_u,
						double* ebqe_bc_u_ext,
						int* isFluxBoundary_u,
						double* ebqe_bc_flux_u_ext,
						double* ebqe_phi,double epsFact,
						double* ebqe_u,
						double* ebqe_flux,
						// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
						int timeOrder,
						int timeStage,
						int useFullNewton,
						double epsFactHeaviside,
						double epsFactDirac,
						double epsFactDiffusion,
						double* phin_dof,
						double* phiHat_dof,
						// interface locator
						double* norm_factor,
						double norm_factor_lagged,
						double* interface_locator,
						double* interface_locator_lagged,
						// normal reconstruction
						double* lumped_wx,
						double* lumped_wy,
						double* lumped_wz,
						double* lumped_wx_tStar,
						double* lumped_wy_tStar,
						double* lumped_wz_tStar,
						// AUX QUANTITIES OF INTEREST
						double* quantDOFs)=0;
    virtual void calculateJacobian_Monolithic(//element
					      double dt,
					      double* mesh_trial_ref,
					      double* mesh_grad_trial_ref,
					      double* mesh_dof,
					      double* mesh_velocity_dof,
					      double MOVING_DOMAIN,
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
					      double useMetrics,
					      double alphaBDF,
					      //VRANS
					      const double* q_porosity,
					      //
					      int* u_l2g,
					      double* elementDiameter,
					      double* nodeDiametersArray,
					      int degree_polynomial,
					      double* u_dof,
					      double* u_dof_old,
					      double* velocity,
					      double* q_m_betaBDF, 
					      double* cfl,
					      int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					      double* globalJacobian,
					      int nExteriorElementBoundaries_global,
					      int* exteriorElementBoundariesArray,
					      int* elementBoundaryElementsArray,
					      int* elementBoundaryLocalElementBoundariesArray,
					      double* ebqe_velocity_ext,
					      //VRANS
					      const double* ebqe_porosity_ext,
					      //
					      int* isDOFBoundary_u,
					      double* ebqe_bc_u_ext,
					      int* isFluxBoundary_u,
					      double* ebqe_bc_flux_u_ext,
					      int* csrColumnOffsets_eb_u_u,
					      // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					      int timeOrder,
					      int timeStage,
					      int useFullNewton,
					      double epsFactHeaviside,
					      double epsFactDirac,
					      double epsFactDiffusion,
					      double* phin_dof,
					      double* phiHat_dof,
					      // interface locator
					      double norm_factor_lagged,
					      double* interface_locator_lagged,
					      // normal reconstruction
					      double* lumped_wx,
					      double* lumped_wy,
					      double* lumped_wz,
					      double* lumped_wx_tStar,
					      double* lumped_wy_tStar,
					      double* lumped_wz_tStar
					      )=0;
    virtual void calculateJacobian_non_Monolithic(//element
						  double dt,
						  double* mesh_trial_ref,
						  double* mesh_grad_trial_ref,
						  double* mesh_dof,
						  double* mesh_velocity_dof,
						  double MOVING_DOMAIN,
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
						  double useMetrics,
						  double alphaBDF,
						  //VRANS
						  const double* q_porosity,
						  //
						  int* u_l2g,
						  double* elementDiameter,
						  double* nodeDiametersArray,
						  int degree_polynomial,
						  double* u_dof,
						  double* u_dof_old,
						  double* velocity,
						  double* q_m_betaBDF, 
						  double* cfl,
						  int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
						  double* globalJacobian,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  double* ebqe_velocity_ext,
						  //VRANS
						  const double* ebqe_porosity_ext,
						  //
						  int* isDOFBoundary_u,
						  double* ebqe_bc_u_ext,
						  int* isFluxBoundary_u,
						  double* ebqe_bc_flux_u_ext,
						  int* csrColumnOffsets_eb_u_u,
						  // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
						  int timeOrder,
						  int timeStage,
						  int useFullNewton,
						  double epsFactHeaviside,
						  double epsFactDirac,
						  double epsFactDiffusion,
						  double* phin_dof,
						  double* phiHat_dof,
						  // interface locator
						  double norm_factor_lagged,
						  double* interface_locator_lagged,
						  // normal reconstruction
						  double* lumped_wx,
						  double* lumped_wy,
						  double* lumped_wz,
						  double* lumped_wx_tStar,
						  double* lumped_wy_tStar,
						  double* lumped_wz_tStar)=0;
    virtual void calculateJacobian_experimental(//element
						double dt,
						double* mesh_trial_ref,
						double* mesh_grad_trial_ref,
						double* mesh_dof,
						double* mesh_velocity_dof,
						double MOVING_DOMAIN,
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
						double useMetrics,
						double alphaBDF,
						//VRANS
						const double* q_porosity,
						//
						int* u_l2g,
						double* elementDiameter,
						double* nodeDiametersArray,
						int degree_polynomial,
						double* u_dof,
						double* u_dof_old,
						double* velocity,
						double* q_m_betaBDF, 
						double* cfl,
						int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
						double* globalJacobian,
						int nExteriorElementBoundaries_global,
						int* exteriorElementBoundariesArray,
						int* elementBoundaryElementsArray,
						int* elementBoundaryLocalElementBoundariesArray,
						double* ebqe_velocity_ext,
						//VRANS
						const double* ebqe_porosity_ext,
						//
						int* isDOFBoundary_u,
						double* ebqe_bc_u_ext,
						int* isFluxBoundary_u,
						double* ebqe_bc_flux_u_ext,
						int* csrColumnOffsets_eb_u_u,
						// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
						int timeOrder,
						int timeStage,
						int useFullNewton,
						double epsFactHeaviside,
						double epsFactDirac,
						double epsFactDiffusion,
						double* phin_dof,
						double* phiHat_dof,
						// interface locator
						double norm_factor_lagged,
						double* interface_locator_lagged,
						// normal reconstruction
						double* lumped_wx,
						double* lumped_wy,
						double* lumped_wz,
						double* lumped_wx_tStar,
						double* lumped_wy_tStar,
						double* lumped_wz_tStar)=0;
    virtual void calculateMetricsAtEOS( //EOS=End Of Simulation
				       double* mesh_trial_ref,
				       double* mesh_grad_trial_ref,
				       double* mesh_dof,
				       int* mesh_l2g,
				       double* dV_ref,
				       double* u_trial_ref,
				       double* u_grad_trial_ref,
				       double* u_test_ref,
				       //physics
				       int nElements_global,
				       int* u_l2g, 
				       double* elementDiameter,
				       //double* nodeDiametersArray,
				       double degree_polynomial,
				       double epsFactHeaviside,
				       double* u_dof,
				       double* u0_dof,
				       double* u_exact,
				       int offset_u, int stride_u,
				       double* global_I_err,
				       double* global_Ieps_err,
				       double* global_V_err,
				       double* global_Veps_err,
				       double* global_D_err)=0;
    virtual void calculateMetricsAtETS( //ETS=Every Time Step
				       double dt,
				       double* mesh_trial_ref,
				       double* mesh_grad_trial_ref,
				       double* mesh_dof,
				       int* mesh_l2g,
				       double* dV_ref,
				       double* u_trial_ref,
				       double* u_grad_trial_ref,
				       double* u_test_ref,
				       //physics
				       int nElements_global,
				       int* u_l2g, 
				       double* elementDiameter,
				       //double* nodeDiametersArray,
				       double degree_polynomial,
				       double epsFactHeaviside,
				       double* u_dof,
				       double* u_dof_old,
				       double* u0_dof,
				       double* velocity,
				       int offset_u, int stride_u,
				       int numDOFs,
				       double* global_R,
				       double* global_Reps,
				       double* global_V_err,
				       double* global_Veps_err,
				       double* global_D_err)=0;    
    virtual void normalReconstruction(double* mesh_trial_ref,
				      double* mesh_grad_trial_ref,
				      double* mesh_dof,
				      int* mesh_l2g,
				      double* dV_ref,
				      double* u_trial_ref,
				      double* u_grad_trial_ref,
				      double* u_test_ref,
				      int nElements_global,
				      int* u_l2g, 
				      double* elementDiameter,
				      double* phi_dof,
				      int offset_u, int stride_u, 
				      int numDOFs,
				      double* lumped_wx,
				      double* lumped_wy,
				      double* lumped_wz)=0;
  };

  template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class CLSVOF : public CLSVOF_base
    {
    public:
      const int nDOF_test_X_trial_element;
      CompKernelType ck;
    CLSVOF():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
	ck()
	  {}
      
      inline
	void calculateCFL(const double& elementDiameter,
			  const double df[nSpace],
			  double& cfl)
      {
	double h,nrm_v;
	h = elementDiameter;
	nrm_v=0.0;
	for(int I=0;I<nSpace;I++)
	  nrm_v+=df[I]*df[I];
	nrm_v = sqrt(nrm_v);
	cfl = nrm_v/h;
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
		std::cout<<"warning: CLSVOF open boundary with no external trace, setting to zero for inflow"<<std::endl;
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
	else if (isFluxBoundary_u == 1)
	  {
	    dflux = 0.0;
	  }
	else
	  {
	    if (flow >= 0.0)
	      {
		dflux = flow;
	      }
	  }
      }

      inline void Mult(const double mat[nSpace][nSpace],
		       const double vec[nSpace],
		       double mat_times_vector[nSpace])
      // Matrix*vector. This is used to project grad_u onto normal direction
      {
	for (int I=0; I<nSpace; I++)
	  {
	    mat_times_vector[I] = 0.;
	    for (int J=0; J<nSpace; J++)
	      mat_times_vector[I] += mat[I][J] * vec[J];
	  }
      }

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

      inline double smoothedSign(double eps, double phi)
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
	return 2*H-1;
      }

      inline double smoothedDerSign(double eps, double phi)
      {
	double d;
	if (phi > eps)
	  d=0.0;
	else if (phi < -eps)
	  d=0.0;
	else
	  d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
	return 2*d; 
      }      

      void calculateResidual_Monolithic(//element
					double dt,
					double* mesh_trial_ref,
					double* mesh_grad_trial_ref,
					double* mesh_dof,
					double* mesh_velocity_dof,
					double MOVING_DOMAIN,
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
					double useMetrics,
					double alphaBDF,
					//VRANS
					const double* q_porosity,
					const double* porosity_dof,
					//
					int* u_l2g, 
					double* elementDiameter,
					double* nodeDiametersArray,
					int degree_polynomial,
					double* u_dof,
					double* u_dof_old,
					double* velocity,
					double* velocity_old,
					double* q_m,
					double* q_u,
					double* q_m_betaBDF,
					double* q_dV,
					double* q_dV_last,
					double* cfl,
					int offset_u, int stride_u, 
					double* globalResidual,
					int nExteriorElementBoundaries_global,
					int* exteriorElementBoundariesArray,
					int* elementBoundaryElementsArray,
					int* elementBoundaryLocalElementBoundariesArray,
					double* ebqe_velocity_ext,
					//VRANS
					const double* ebqe_porosity_ext,
					//
					int* isDOFBoundary_u,
					double* ebqe_bc_u_ext,
					int* isFluxBoundary_u,
					double* ebqe_bc_flux_u_ext,
					double* ebqe_phi,double epsFact,
					double* ebqe_u,
					double* ebqe_flux,
					// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					int timeOrder,
					int timeStage,
					int useFullNewton,
					double epsFactHeaviside,
					double epsFactDirac,
					double epsFactDiffusion,
					double* phin_dof,
					double* phiHat_dof,
					// interface locator
					double* norm_factor,
					double norm_factor_lagged,
					double* interface_locator,
					double* interface_locator_lagged,
					// normal reconstruction
					double* lumped_wx,
					double* lumped_wy,
					double* lumped_wz,
					double* lumped_wx_tStar,
					double* lumped_wy_tStar,
					double* lumped_wz_tStar,
					// AUX QUANTITIES OF INTEREST 
					double* quantDOFs)
      {
	double min_distance = 1E10;
	double max_distance = -1E10;
	double mean_distance = 0.;
	
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    interface_locator[eN] = 1.0;
	    //declare local storage for local contributions and initialize
	    register double elementResidual_u[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementResidual_u[i]=0.0;
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double 
		  //for mass matrix contributions
		  un, grad_un[nSpace], grad_unHalf[nSpace], 
		  u, grad_u[nSpace],
		  normalReconstruction[nSpace],
		  qxn, qyn, qzn, qxnStar, qynStar, qznStar,
		  relative_velocity[nSpace], relative_velocity_old[nSpace],
		  fnp1[nSpace], fnHalf[nSpace], //f=velocity*H(phi)
		  u_test_dV[nDOF_trial_element], 
		  u_grad_trial[nDOF_trial_element*nSpace], 
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z,xt,yt,zt,h_phi;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
                ck.calculateH_element(eN,
                                      k,
                                      nodeDiametersArray,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      h_phi);		
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		dV = fabs(jacDet)*dV_ref[k];
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		// get the components of the normal reconstruction
		ck.valFromDOF(lumped_wx,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxn);
		ck.valFromDOF(lumped_wy,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qyn);
		ck.valFromDOF(lumped_wz,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qzn);
		ck.valFromDOF(lumped_wx_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxnStar);
		ck.valFromDOF(lumped_wy_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qynStar);
		ck.valFromDOF(lumped_wz_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qznStar);
		// get the solution (of Newton's solver)
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		// get old solution
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      un);
		//get the solution gradients at quad points	      
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_u);
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_un);
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		//calculate time derivative at quadrature points
		if (q_dV_last[eN_k] <= -100)
		  q_dV_last[eN_k] = dV;
		q_dV[eN_k] = dV;

		/////////////////
		// MOVING MESH //
		/////////////////
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double lambda = epsFactDiffusion;	      
		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		double Hn = smoothedSign(epsHeaviside,un);
		double Hnp1 = smoothedSign(epsHeaviside,u);
		
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    relative_velocity_old[I] = (velocity_old[eN_k_nSpace+I]
						-MOVING_DOMAIN*mesh_velocity[I]);
		    fnp1[I] = relative_velocity[I]*Hnp1; //implicit advection via BDF
		    fnHalf[I] = 0.5*(relative_velocity[I]*Hnp1
				     +relative_velocity_old[I]*Hn); //implicit advection via CN
		    grad_unHalf[I] = 0.5*(grad_u[I] + grad_un[I]);
		  }
		//////////////////////////////
		// CALCULATE CELL BASED CFL //
		//////////////////////////////
		calculateCFL(elementDiameter[eN]/degree_polynomial,relative_velocity,cfl[eN_k]);

		/////////////////////
		// TIME DERIVATIVE //
		/////////////////////
		double time_derivative_residual = (Hnp1-Hn)/dt;

		///////////////////////
		// INTERFACE LOCATOR //
		///////////////////////
		interface_locator[eN] *= u;

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		//double gamma = interface_locator_lagged[eN] > 0 ? 1. : 0.;
		double gamma = (useMetrics*h_phi + (1.0-useMetrics)
				*elementDiameter[eN])/degree_polynomial/norm_factor_lagged;

		// CALCULATE min, max and mean distance
		min_distance = fmin(min_distance,u);
		max_distance = fmax(max_distance,u);
		mean_distance += u*dV;

		///////////////////////////
		// NORMAL RECONSTRUCTION //
		///////////////////////////
		if (timeOrder == 2 && timeStage == 2)
		  {
		    normalReconstruction[0] = 0.5*(qxnStar+qxn);
		    normalReconstruction[1] = 0.5*(qynStar+qyn);
#if nSpace==3
		    normalReconstruction[2] = 0.5*(qznStar+qzn);
#endif
		  }
		else //timeOrder == 1 or timeStage==1
		  {
		    normalReconstruction[0] = qxn;
		    normalReconstruction[1] = qyn;
#if nSpace==3
		    normalReconstruction[2] = qzn;
#endif
		  }
		//////////////////
		// LOOP ON DOFs //
		//////////////////
		for(int i=0;i<nDOF_test_element;i++) 
		  { 
		    register int i_nSpace=i*nSpace;
		    if (timeOrder==1)
		      {
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnp1,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + gamma*lambda*ck.NumericalDiffusion(1.0,
							       grad_u,
							       &u_grad_test_dV[i_nSpace])
			  // TARGET for PENALIZATION. This is EXPLICIT
			  - gamma*lambda*ck.NumericalDiffusion(1.0,
							       normalReconstruction,
							       &u_grad_test_dV[i_nSpace]);
		      }
		    else // timeOrder=2
		      elementResidual_u[i] +=
			// TIME DERIVATIVE
			time_derivative_residual*u_test_dV[i]
			// ADVECTION TERM. This is IMPLICIT
			+ ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			// REGULARIZATION TERM. This is IMPLICIT
			+ gamma*lambda*ck.NumericalDiffusion(1.0,
							     grad_unHalf,
							     &u_grad_test_dV[i_nSpace])
			// TARGET for PENALIZATION. This is EXPLICIT 
			- gamma*lambda*ck.NumericalDiffusion(1.0,
							     normalReconstruction,
							     &u_grad_test_dV[i_nSpace]);
		  }//i
		//save solution for other models 
		q_u[eN_k] = u;
		q_m[eN_k] = u;//porosity*u;
	      }
	    // NORMALIZE INTERFACE LOCATOR //
	    interface_locator[eN] /= fabs(interface_locator[eN]);
	    /////////////////
	    // DISTRIBUTE // load cell based element into global residual
	    ////////////////
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		// distribute global residual for (lumped) mass matrix
		globalResidual[gi] += elementResidual_u[i];
	      }//i
	  }//elements

	norm_factor[0] = fmax(fabs(max_distance - mean_distance),
			      fabs(mean_distance - min_distance));
	  
	//////////////
	// BOUNDARY //
	//////////////
	// END OF BOUNDARY //
      }

      // piece of code to incorporate: |Seps(phi)| with (\nabla\phi-q)
		///////////////////////////
		// NORMAL RECONSTRUCTION //
		///////////////////////////
		/*
		if (timeOrder == 2 && timeStage == 2)
		  {
		    normalReconstruction[0] = 0.5*(rabs(Hnp1)*qxnStar+rabs(Hn)*qxn);
		    normalReconstruction[1] = 0.5*(rabs(Hnp1)*qynStar+rabs(Hn)*qyn);
#if nSpace==3
		    normalReconstruction[2] = 0.5*(rabs(Hnp1)*qznStar+rabs(Hn)*qzn);
#endif
		  }
		else //timeOrder == 1 or timeStage==1
		  {
		    normalReconstruction[0] = rabs(Hn)*qxn;
		    normalReconstruction[1] = rabs(Hn)*qyn;
#if nSpace==3
		    normalReconstruction[2] = rabs(Hn)*qzn;
#endif
		  }
		//////////////////
		// LOOP ON DOFs //
		//////////////////
		for(int i=0;i<nDOF_test_element;i++) 
		  { 
		    register int i_nSpace=i*nSpace;
		    if (timeOrder==1)
		      {
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnp1,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + rabs(Hnp1)*gamma*lambda*ck.NumericalDiffusion(1.0,
									grad_u,
									&u_grad_test_dV[i_nSpace])
			  // TARGET for PENALIZATION. This is EXPLICIT
			  - rabs(Hn)*gamma*lambda*ck.NumericalDiffusion(1.0,
									normalReconstruction,
									&u_grad_test_dV[i_nSpace]);
		      }
		    else // timeOrder=2
		      elementResidual_u[i] +=
			// TIME DERIVATIVE
			time_derivative_residual*u_test_dV[i]
			// ADVECTION TERM. This is IMPLICIT
			+ ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			// REGULARIZATION TERM. This is IMPLICIT
			+ gamma*lambda*ck.NumericalDiffusion(1.0,
							     grad_unHalf,
							     &u_grad_test_dV[i_nSpace])
			// TARGET for PENALIZATION. This is EXPLICIT 
			- gamma*lambda*ck.NumericalDiffusion(1.0,
							     normalReconstruction,
							     &u_grad_test_dV[i_nSpace]);
		  }//i
		*/
      
      void calculateResidual_non_Monolithic(//element
					    double dt,
					    double* mesh_trial_ref,
					    double* mesh_grad_trial_ref,
					    double* mesh_dof,
					    double* mesh_velocity_dof,
					    double MOVING_DOMAIN,
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
					    double useMetrics,
					    double alphaBDF,
					    //VRANS
					    const double* q_porosity,
					    const double* porosity_dof,
					    //
					    int* u_l2g, 
					    double* elementDiameter,
					    double* nodeDiametersArray,
					    int degree_polynomial,
					    double* u_dof,
					    double* u_dof_old,
					    double* velocity,
					    double* velocity_old,
					    double* q_m,
					    double* q_u,
					    double* q_m_betaBDF,
					    double* q_dV,
					    double* q_dV_last,
					    double* cfl,
					    int offset_u, int stride_u, 
					    double* globalResidual,
					    int nExteriorElementBoundaries_global,
					    int* exteriorElementBoundariesArray,
					    int* elementBoundaryElementsArray,
					    int* elementBoundaryLocalElementBoundariesArray,
					    double* ebqe_velocity_ext,
					    //VRANS
					    const double* ebqe_porosity_ext,
					    //
					    int* isDOFBoundary_u,
					    double* ebqe_bc_u_ext,
					    int* isFluxBoundary_u,
					    double* ebqe_bc_flux_u_ext,
					    double* ebqe_phi,double epsFact,
					    double* ebqe_u,
					    double* ebqe_flux,
					    // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					    int timeOrder,
					    int timeStage,
					    int useFullNewton,
					    double epsFactHeaviside,
					    double epsFactDirac,
					    double epsFactDiffusion,
					    double* phin_dof,
					    double* phiHat_dof,
					    // interface locator
					    double* norm_factor,
					    double norm_factor_lagged,
					    double* interface_locator,
					    double* interface_locator_lagged,
					    // normal reconstruction
					    double* lumped_wx,
					    double* lumped_wy,
					    double* lumped_wz,
					    double* lumped_wx_tStar,
					    double* lumped_wy_tStar,
					    double* lumped_wz_tStar,
					    // AUX QUANTITIES OF INTEREST 
					    double* quantDOFs)
      {
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double elementResidual_u[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementResidual_u[i]=0.0;
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double 
		  //for mass matrix contributions
		  u, grad_u[nSpace], relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		  phiHatnp1, phin,
		  u_test_dV[nDOF_trial_element], 
		  u_grad_trial[nDOF_trial_element*nSpace], 
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z,xt,yt,zt;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		dV = fabs(jacDet)*dV_ref[k];
		//get the solution (of Newton's solver)
		ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
		//get the solution gradients at quad points
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
		// get phin and phiHatnp1 at quad points
		ck.valFromDOF(phin_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phin);
		ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  }
		//calculate time derivative at quadrature points
		if (q_dV_last[eN_k] <= -100)
		  q_dV_last[eN_k] = dV;
		q_dV[eN_k] = dV;
		//
		//moving mesh
		//
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		double epsDiffusion = epsFactDiffusion*elementDiameter[eN]/degree_polynomial; //kappa = const*h
		double Hn = smoothedSign(epsHeaviside,phin);
		double Hnp1 = smoothedSign(epsHeaviside,phiHatnp1+u);
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    f[I] = relative_velocity[I]*Hn; // Explicit
		    //f[I] = relative_velocity[I]*Hnp1; // Implicit
		  }
		//////////////////////////////
		// CALCULATE CELL BASED CFL //
		//////////////////////////////
		calculateCFL(elementDiameter[eN]/degree_polynomial,relative_velocity,cfl[eN_k]);
		
		double time_derivative_residual=(Hnp1-Hn)/dt;
		//////////////
		// ith-LOOP //
		//////////////	      
		for(int i=0;i<nDOF_test_element;i++) 
		  { 
		    register int i_nSpace=i*nSpace;
		    elementResidual_u[i] += 
		      time_derivative_residual*u_test_dV[i]
		      + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
		      + ck.NumericalDiffusion(epsDiffusion/dt,grad_u,&u_grad_test_dV[i_nSpace]);
		  }//i
		//save solution for other models 
		q_u[eN_k] = u;
		q_m[eN_k] = u;//porosity*u;
	      }
	    /////////////////
	    // DISTRIBUTE // load cell based element into global residual
	    ////////////////
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		// distribute global residual for (lumped) mass matrix
		globalResidual[gi] += elementResidual_u[i];
	      }//i
	  }//elements
	//////////////
	// BOUNDARY //
	//////////////
	// END OF BOUNDARY //
      }

      void calculateResidual_experimental(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
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
					  double useMetrics,
					  double alphaBDF,
					  //VRANS
					  const double* q_porosity,
					  const double* porosity_dof,
					  //
					  int* u_l2g, 
					  double* elementDiameter,
					  double* nodeDiametersArray,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* velocity,
					  double* velocity_old,
					  double* q_m,
					  double* q_u,
					  double* q_m_betaBDF,
					  double* q_dV,
					  double* q_dV_last,
					  double* cfl,
					  int offset_u, int stride_u, 
					  double* globalResidual,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  double* ebqe_phi,double epsFact,
					  double* ebqe_u,
					  double* ebqe_flux,
					  // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					  int timeOrder,
					  int timeStage,
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  // interface locator
					  double* norm_factor,
					  double norm_factor_lagged,
					  double* interface_locator,
					  double* interface_locator_lagged,
					  // normal reconstruction
					  double* lumped_wx,
					  double* lumped_wy,
					  double* lumped_wz,
					  double* lumped_wx_tStar,
					  double* lumped_wy_tStar,
					  double* lumped_wz_tStar,
					  // AUX QUANTITIES OF INTEREST 
					  double* quantDOFs)
      {
	double min_distance = 1E10;
	double max_distance = -1E10;
	double mean_distance = 0.;
	
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double elementResidual_u[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementResidual_u[i]=0.0;
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double 
		  //for mass matrix contributions
		  u, grad_u[nSpace], proj_times_grad_u[nSpace],
		  normalReconstruction[nSpace], wx, wy,
		  un, grad_un[nSpace], 
		  relative_velocity[nSpace], f[nSpace], //f=velocity*H(phi)
		  phiHatnp1, phin,
		  u_test_dV[nDOF_trial_element], 
		  u_grad_trial[nDOF_trial_element*nSpace], 
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z,xt,yt,zt;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		dV = fabs(jacDet)*dV_ref[k];
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		// get the components of the normal reconstruction
		ck.valFromDOF(lumped_wx,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      wx);
		ck.valFromDOF(lumped_wy,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      wy);		
		// get the solution (of Newton's solver)
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		// get old solutions
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      un);
		//get the solution gradients at quad points	      
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_u);
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_un);
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  }
		//calculate time derivative at quadrature points
		if (q_dV_last[eN_k] <= -100)
		  q_dV_last[eN_k] = dV;
		q_dV[eN_k] = dV;
		//
		//moving mesh
		//
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double gradu2 = 0;
		for(int I=0;I<nSpace;I++)
		  gradu2 += grad_u[I]*grad_u[I];
		double beta_norm_grad_u = betaNormGrad(gradu2,1.E-15);

		double gradun2 = 0;
		for(int I=0;I<nSpace;I++)
		  gradun2 += grad_un[I]*grad_un[I];
		double beta_norm_grad_un = betaNormGrad(gradun2,1.E-15);

		double lambda = epsFactDiffusion;	      
		//double lambda = epsFactDiffusion*elementDiameter[eN]/degree_polynomial/dt;

		//double coeffFullNewton = -1./beta_norm_grad_u; // single potential
		double coeffFullNewton = 2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u; //d.pot
       
		//double coeffLinNewton = -1./beta_norm_grad_un; // single potential
		double coeffLinNewton = 2*std::pow(beta_norm_grad_un,2)-3*beta_norm_grad_un; //d.pot

		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		double Hn = smoothedSign(epsHeaviside,un);
		double Hnp1 = smoothedSign(epsHeaviside,u);
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    f[I] = relative_velocity[I]*Hn;
		    //f[I] = relative_velocity[I]*Hnp1;
		    //f[I] = 0.5*relative_velocity[I]*(Hnp1+Hn); //implicit advection via CN
		  }
		//////////////////////////////
		// CALCULATE CELL BASED CFL //
		//////////////////////////////
		calculateCFL(elementDiameter[eN]/degree_polynomial,relative_velocity,cfl[eN_k]);
		
		double time_derivative_residual = (Hnp1-Hn)/dt;

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		double gamma = 1.0;  //elementDiameter[eN]/degree_polynomial/norm_factor_lagged;

		// CALCULATE min, max and mean distance
		min_distance = fmin(min_distance,u);
		max_distance = fmax(max_distance,u);
		mean_distance += u*dV;
		
		//////////////
		// ith-LOOP //
		//////////////
		if (useFullNewton==false)
		  for(int i=0;i<nDOF_test_element;i++) 
		    { 
		      register int i_nSpace=i*nSpace;
		      elementResidual_u[i] += 
			time_derivative_residual*u_test_dV[i]
			+ ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
			+ gamma*lambda*ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
			+ gamma*lambda*ck.NumericalDiffusion(coeffLinNewton,
							     grad_un,&u_grad_test_dV[i_nSpace]);
		    }//i
		else
		  {
		    for(int i=0;i<nDOF_test_element;i++) 
		      { 
			register int i_nSpace=i*nSpace;
			elementResidual_u[i] += 
			  time_derivative_residual*u_test_dV[i]
			  + ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])
			  + gamma*lambda*ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
			  + gamma*lambda*ck.NumericalDiffusion(coeffFullNewton,
							       grad_u,&u_grad_test_dV[i_nSpace]);
		      }//i
		  }
		//save solution for other models 
		q_u[eN_k] = u;
		q_m[eN_k] = u;//porosity*u;
	      }
	    /////////////////
	    // DISTRIBUTE // load cell based element into global residual
	    ////////////////
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		// distribute global residual for (lumped) mass matrix
		globalResidual[gi] += elementResidual_u[i];
	      }//i
	  }//elements

	norm_factor[0] = fmax(fabs(max_distance - mean_distance),
			      fabs(mean_distance - min_distance));
	//////////////
	// BOUNDARY //
	//////////////
	// END OF BOUNDARY //
      }

      /*
//      CODE TO CONSIDER ANISOTROPIC REGULARIZATION/PENALIZATION
	void calculateResidual_experimental(//element
					double dt,
					double* mesh_trial_ref,
					double* mesh_grad_trial_ref,
					double* mesh_dof,
					double* mesh_velocity_dof,
					double MOVING_DOMAIN,
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
					double useMetrics,
					double alphaBDF,
					//VRANS
					const double* q_porosity,
					const double* porosity_dof,
					//
					int* u_l2g, 
					double* elementDiameter,
					double* nodeDiametersArray,
					int degree_polynomial,
					double* u_dof,
					double* u_dof_old,
					double* velocity,
					double* velocity_old,
					double* q_m,
					double* q_u,
					double* q_m_betaBDF,
					double* q_dV,
					double* q_dV_last,
					double* cfl,
					int offset_u, int stride_u, 
					double* globalResidual,
					int nExteriorElementBoundaries_global,
					int* exteriorElementBoundariesArray,
					int* elementBoundaryElementsArray,
					int* elementBoundaryLocalElementBoundariesArray,
					double* ebqe_velocity_ext,
					//VRANS
					const double* ebqe_porosity_ext,
					//
					int* isDOFBoundary_u,
					double* ebqe_bc_u_ext,
					int* isFluxBoundary_u,
					double* ebqe_bc_flux_u_ext,
					double* ebqe_phi,double epsFact,
					double* ebqe_u,
					double* ebqe_flux,
					// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					int timeOrder,
					int timeStage,
					int useFullNewton,
					double epsFactHeaviside,
					double epsFactDirac,
					double epsFactDiffusion,
					double* phin_dof,
					double* phiHat_dof,
					// interface locator
					double* norm_factor,
					double norm_factor_lagged,
					double* interface_locator,
					double* interface_locator_lagged,
					// normal reconstruction
					double* lumped_wx,
					double* lumped_wy,
					double* lumped_wz,
					double* lumped_wx_tStar,
					double* lumped_wy_tStar,
					double* lumped_wz_tStar,
					// AUX QUANTITIES OF INTEREST 
					double* quantDOFs)
      {
	double min_distance = 1E10;
	double max_distance = -1E10;
	double mean_distance = 0.;
	
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    interface_locator[eN] = 1.0;
	    //declare local storage for local contributions and initialize
	    register double elementResidual_u[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementResidual_u[i]=0.0;
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double 
		  //for mass matrix contributions
		  un, grad_un[nSpace], grad_unHalf[nSpace], 
		  u, grad_u[nSpace],		  
		  qn[nSpace], qnStar[nSpace],
		  qxn, qyn, qzn, qxnStar, qynStar, qznStar,
		  projn_times_grad_un[nSpace],
		  projn_times_grad_unHalf[nSpace],
		  projn_times_grad_u[nSpace],
		  projnStar_times_grad_u[nSpace],
		  projn[nSpace][nSpace], projnStar[nSpace][nSpace],		  
		  relative_velocity[nSpace],
		  fnp1[nSpace], fnHalf[nSpace], //f=velocity*H(phi)
		  u_test_dV[nDOF_trial_element], 
		  u_grad_trial[nDOF_trial_element*nSpace], 
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z,xt,yt,zt,h_phi;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
                ck.calculateH_element(eN,
                                      k,
                                      nodeDiametersArray,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      h_phi);		
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		dV = fabs(jacDet)*dV_ref[k];
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		// get the components of the normal reconstruction
		ck.valFromDOF(lumped_wx,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxn);
		ck.valFromDOF(lumped_wy,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qyn);
		ck.valFromDOF(lumped_wz,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qzn);
		ck.valFromDOF(lumped_wx_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxnStar);
		ck.valFromDOF(lumped_wy_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qynStar);
		ck.valFromDOF(lumped_wz_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qznStar);
		// get the solution (of Newton's solver)
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		// get old solution
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      un);
		//get the solution gradients at quad points	      
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_u);
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_un);
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		//calculate time derivative at quadrature points
		if (q_dV_last[eN_k] <= -100)
		  q_dV_last[eN_k] = dV;
		q_dV[eN_k] = dV;

		/////////////////
		// MOVING MESH //
		/////////////////
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double lambda = epsFactDiffusion;	      
		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial; 	    
		double Hn = smoothedSign(epsHeaviside,un);
		double Hnp1 = smoothedSign(epsHeaviside,u);
		
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    fnp1[I] = relative_velocity[I]*Hnp1; //implicit advection via BDF
		    fnHalf[I] = 0.5*relative_velocity[I]*(Hnp1+Hn); //implicit advection via CN
		    grad_unHalf[I] = 0.5*(grad_u[I] + grad_un[I]);
		  }
		//////////////////////////////
		// CALCULATE CELL BASED CFL //
		//////////////////////////////
		calculateCFL(elementDiameter[eN]/degree_polynomial,relative_velocity,cfl[eN_k]);

		/////////////////////
		// TIME DERIVATIVE //
		/////////////////////
		double time_derivative_residual = (Hnp1-Hn)/dt;

		///////////////////////
		// INTERFACE LOCATOR //
		///////////////////////
		interface_locator[eN] *= u;

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		//double gamma = interface_locator_lagged[eN] > 0 ? 1. : 0.;
		double gamma = (useMetrics*h_phi
				+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial/norm_factor_lagged;

		// CALCULATE min, max and mean distance
		min_distance = fmin(min_distance,u);
		max_distance = fmax(max_distance,u);
		mean_distance += u*dV;

		///////////////////////////
		// NORMAL RECONSTRUCTION //
		///////////////////////////
		qn[0] = qxn; 
		qn[1] = qyn;
		qnStar[0] = qxnStar;
		qnStar[1] = qynStar;
#if nSpace==3
		qn[2] = qzn;
		qnStar[2] = qznStar;
#endif

		///////////////////////////////////
		// PROJECTOR TO NORMAL COMPONENT //
		///////////////////////////////////
		double norm2_qn = 0; 
		double norm2_qnStar = 0;
		for (int I=0; I<nSpace; I++)
		  {
		    norm2_qn += qn[I]*qn[I];
		    norm2_qnStar += qnStar[I]*qnStar[I];
		  }
		double norm_qn = std::sqrt(norm2_qn);
		double norm_qnStar = std::sqrt(norm2_qnStar);
		for (int I = 0; I < nSpace; ++I)
		  for (int J = 0; J < nSpace; ++J)
		    {
		      projn[I][J] = qn[I]*qn[J]/(norm2_qn+1E-15);
		      projnStar[I][J] = qnStar[I]*qnStar[J]/(norm2_qnStar+1E-15);
		    }
		Mult(projn,grad_un,projn_times_grad_un);
		Mult(projn,grad_unHalf,projn_times_grad_unHalf);
		Mult(projn,grad_u,projn_times_grad_u);
		Mult(projnStar,grad_u,projnStar_times_grad_u);
		//////////////////
		// LOOP ON DOFs //
		//////////////////
		for(int i=0;i<nDOF_test_element;i++) 
		  { 
		    register int i_nSpace=i*nSpace;
		    if (timeOrder==1)
		      {
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnp1,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + gamma*lambda*(ISOTROPIC_REGULARIZATION == 1 ? 
					  ck.NumericalDiffusion(1.0,
			  					grad_u,
			  					&u_grad_test_dV[i_nSpace])
					  :
					  +ck.NumericalDiffusion(norm_qn,
			  					 projn_times_grad_u,
			  					 &u_grad_test_dV[i_nSpace])
			  		  +ck.NumericalDiffusion(1.0-norm_qn,
			  					 grad_u,
			  					 &u_grad_test_dV[i_nSpace]))
			  // TARGET for PENALIZATION. This is EXPLICIT
			  - gamma*lambda*ck.NumericalDiffusion(1.0,
							       qn,
							       &u_grad_test_dV[i_nSpace]);
		      }
		    else // timeOrder=2
		      if (timeStage==1)
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + gamma*lambda*(ISOTROPIC_REGULARIZATION == 1 ?
					  ck.NumericalDiffusion(1.0,
								grad_unHalf,
								&u_grad_test_dV[i_nSpace])
					  :
					  // at time tn
					  +ck.NumericalDiffusion(norm_qn,
								 projn_times_grad_unHalf,
								 &u_grad_test_dV[i_nSpace])
					  +ck.NumericalDiffusion(1.0-norm_qn,
								 grad_unHalf,
								 &u_grad_test_dV[i_nSpace]))
			  // TARGET for PENALIZATION. This is EXPLICIT 
			  - gamma*lambda*ck.NumericalDiffusion(1.0,
							       qn,
							       &u_grad_test_dV[i_nSpace]);
		      else
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + gamma*lambda*(ISOTROPIC_REGULARIZATION == 1 ?
			  		  ck.NumericalDiffusion(1.0,
			  					grad_unHalf,
			  					&u_grad_test_dV[i_nSpace])
					  // at time tn
			  		  :
			  		  +0.5*(ck.NumericalDiffusion(norm_qn,
			  					      projn_times_grad_un,
			  					      &u_grad_test_dV[i_nSpace])
			  			+ck.NumericalDiffusion(1.0-norm_qn,
			  					       grad_un,
			  					       &u_grad_test_dV[i_nSpace]))
					  // at time tnp1
			  		  +0.5*(ck.NumericalDiffusion(norm_qnStar,
								      projnStar_times_grad_u,
			  					      &u_grad_test_dV[i_nSpace])
			  			+ck.NumericalDiffusion(1.0-norm_qnStar, 
								       grad_u,
			  					       &u_grad_test_dV[i_nSpace])))
			  // TARGET for PENALIZATION. This is EXPLICIT
		    	  - 0.5*gamma*lambda*(ck.NumericalDiffusion(1.0,
		    						    qn,
		    						    &u_grad_test_dV[i_nSpace])
		    			      +ck.NumericalDiffusion(1.0,
		    						     qnStar,
		    						     &u_grad_test_dV[i_nSpace]));
		  }//i
		//save solution for other models 
		q_u[eN_k] = u;
		q_m[eN_k] = u;//porosity*u;
	      }
	    // NORMALIZE INTERFACE LOCATOR //
	    interface_locator[eN] /= fabs(interface_locator[eN]);
	    /////////////////
	    // DISTRIBUTE // load cell based element into global residual
	    ////////////////
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		// distribute global residual for (lumped) mass matrix
		globalResidual[gi] += elementResidual_u[i];
	      }//i
	  }//elements

	norm_factor[0] = fmax(fabs(max_distance - mean_distance),
			      fabs(mean_distance - min_distance));
	  
	//////////////
	// BOUNDARY //
	//////////////
	// END OF BOUNDARY //
      }
       */
      
      void calculateJacobian_Monolithic(//element
					double dt,
					double* mesh_trial_ref,
					double* mesh_grad_trial_ref,
					double* mesh_dof,
					double* mesh_velocity_dof,
					double MOVING_DOMAIN,
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
					double useMetrics,
					double alphaBDF,
					//VRANS
					const double* q_porosity,
					//
					int* u_l2g,
					double* elementDiameter,
					double* nodeDiametersArray,
					int degree_polynomial,
					double* u_dof,
					double* u_dof_old,
					double* velocity,
					double* q_m_betaBDF, 
					double* cfl,
					int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					double* globalJacobian,
					int nExteriorElementBoundaries_global,
					int* exteriorElementBoundariesArray,
					int* elementBoundaryElementsArray,
					int* elementBoundaryLocalElementBoundariesArray,
					double* ebqe_velocity_ext,
					//VRANS
					const double* ebqe_porosity_ext,
					//
					int* isDOFBoundary_u,
					double* ebqe_bc_u_ext,
					int* isFluxBoundary_u,
					double* ebqe_bc_flux_u_ext,
					int* csrColumnOffsets_eb_u_u,
					// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					int timeOrder,
					int timeStage,
					int useFullNewton,
					double epsFactHeaviside,
					double epsFactDirac,
					double epsFactDiffusion,
					double* phin_dof,
					double* phiHat_dof,
					// interface locator
					double norm_factor_lagged,
					double* interface_locator_lagged,
					// normal reconstruction
					double* lumped_wx,
					double* lumped_wy,
					double* lumped_wz,
					double* lumped_wx_tStar,
					double* lumped_wy_tStar,
					double* lumped_wz_tStar)
      {
	double timeCoeff=1.0;
	if (timeOrder==2)
	  timeCoeff=0.5;
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      for (int j=0;j<nDOF_trial_element;j++)
		elementJacobian_u_u[i][j]=0.0;	      
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
		//declare local storage
		register double
		  u, un, u_grad_trial[nDOF_trial_element*nSpace],		  
		  relative_velocity[nSpace], df[nSpace], 		  
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		  u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		  dV, x,y,z,xt,yt,zt,h_phi;
		//get jacobian, etc for mapping reference element
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
                ck.calculateH_element(eN,
                                      k,
                                      nodeDiametersArray,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      h_phi);		
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		//get the physical integration weight
		dV = fabs(jacDet)*dV_ref[k];
		//get the trial function gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		//get the solution 	
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      un);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		  }

		/////////////////
		// MOVING MESH //
		/////////////////
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double lambda = epsFactDiffusion;
		//double lambda = epsFactDiffusion*elementDiameter[eN]/degree_polynomial/dt;
		double epsDirac = epsFactDirac*elementDiameter[eN]/degree_polynomial;
		double dHnp1 = smoothedDerSign(epsDirac,u); //derivative of smoothed sign
		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
				
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    df[I] = relative_velocity[I]*dHnp1;
		  }

		/////////////////////
		// TIME DERIVATIVE //
		/////////////////////
		double time_derivative_jacobian = dHnp1/dt;

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		//double gamma = interface_locator_lagged[eN] > 0 ? 1. : 0.;
		double gamma = (useMetrics*h_phi + (1.0-useMetrics)
				*elementDiameter[eN])/degree_polynomial/norm_factor_lagged;

		//////////////////
		// LOOP ON DOFs //
		//////////////////
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;

			elementJacobian_u_u[i][j] +=
			  // TIME DERIVATIVE 
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  // IMPLICIT TERMS: ADVECTION, DIFFUSION
			  + timeCoeff*
			  (ck.AdvectionJacobian_weak(df,
			   			     u_trial_ref[k*nDOF_trial_element+j],
			   			     &u_grad_test_dV[i_nSpace])
			   + gamma*lambda*ck.NumericalDiffusionJacobian(1.0,
									&u_grad_trial[j_nSpace],
									&u_grad_test_dV[i_nSpace]));
		      }//j
		  }//i
	      }//k
	    //
	    //load into element Jacobian into global Jacobian
	    //
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		      elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements
      }//computeJacobian for MCorr with CLSVOF      

      // piece of code to incorporate: |Seps(phi)| with (\nabla\phi-q)
      /*
                double Hnp1 = smoothedSign(epsHeaviside,u);
		
		//////////////////
		// LOOP ON DOFs //
		//////////////////
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;

			elementJacobian_u_u[i][j] +=
			  // TIME DERIVATIVE 
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  // IMPLICIT TERMS: ADVECTION, DIFFUSION
			  + timeCoeff*
			  (ck.AdvectionJacobian_weak(df,
						     u_trial_ref[k*nDOF_trial_element+j],
						     &u_grad_test_dV[i_nSpace])
			   + rabs(Hnp1)*gamma*lambda*
			   ck.NumericalDiffusionJacobian(1.0,
							 &u_grad_trial[j_nSpace],
							 &u_grad_test_dV[i_nSpace]));
			//+ Hnp1*dHnp1/rabs(Hnp1)*u_test_dV[j]*
			// ck.NumericalDiffusionJacobian(1.0,
			//				 grad_u,
			//				 &u_grad_trial[i_nSpace]));
		      }//j
		  }//i
      */
      
      void calculateJacobian_non_Monolithic(//element
					    double dt,
					    double* mesh_trial_ref,
					    double* mesh_grad_trial_ref,
					    double* mesh_dof,
					    double* mesh_velocity_dof,
					    double MOVING_DOMAIN,
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
					    double useMetrics,
					    double alphaBDF,
					    //VRANS
					    const double* q_porosity,
					    //
					    int* u_l2g,
					    double* elementDiameter,
					    double* nodeDiametersArray,
					    int degree_polynomial,
					    double* u_dof,
					    double* u_dof_old,
					    double* velocity,
					    double* q_m_betaBDF, 
					    double* cfl,
					    int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					    double* globalJacobian,
					    int nExteriorElementBoundaries_global,
					    int* exteriorElementBoundariesArray,
					    int* elementBoundaryElementsArray,
					    int* elementBoundaryLocalElementBoundariesArray,
					    double* ebqe_velocity_ext,
					    //VRANS
					    const double* ebqe_porosity_ext,
					    //
					    int* isDOFBoundary_u,
					    double* ebqe_bc_u_ext,
					    int* isFluxBoundary_u,
					    double* ebqe_bc_flux_u_ext,
					    int* csrColumnOffsets_eb_u_u,
					    // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					    int timeOrder,
					    int timeStage,
					    int useFullNewton,
					    double epsFactHeaviside,
					    double epsFactDirac,
					    double epsFactDiffusion,
					    double* phin_dof,
					    double* phiHat_dof,
					    // interface locator
					    double norm_factor_lagged,
					    double* interface_locator_lagged,
					    // normal reconstruction
					    double* lumped_wx,
					    double* lumped_wy,
					    double* lumped_wz,
					    double* lumped_wx_tStar,
					    double* lumped_wy_tStar,
					    double* lumped_wz_tStar)
      {
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      for (int j=0;j<nDOF_trial_element;j++)
		elementJacobian_u_u[i][j]=0.0;	      
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
		//declare local storage
		register double u, phiHatnp1, u_grad_trial[nDOF_trial_element*nSpace],
		  df[nSpace], relative_velocity[nSpace],
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		  dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		  x,y,z,xt,yt,zt;
		//get jacobian, etc for mapping reference element
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		//get the physical integration weight
		dV = fabs(jacDet)*dV_ref[k];
		//get the trial function gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		//get the solution 	
		ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
		//get phiHat at tnp1
		ck.valFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],phiHatnp1);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		  }
		double epsDiffusion = epsFactDiffusion*elementDiameter[eN]/degree_polynomial; //kappa*h
		double epsDirac = epsFactDirac*elementDiameter[eN]/degree_polynomial; //kappa*h
		double time_derivative_jacobian = smoothedDerSign(epsDirac,phiHatnp1+u)/dt;
		
		///////////////////////////////////
		// FOR IMPLICIT NONLINEAR CLSVOF // (implicit advection term)
		///////////////////////////////////
		//double mesh_velocity[3];
		//mesh_velocity[0] = xt;
		//mesh_velocity[1] = yt;
		//mesh_velocity[2] = zt;
		//double dHnp1 = smoothedDerSign(epsDirac,phiHatnp1+u);
		//for (int I=0;I<nSpace;I++)
		//{
		//relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		//df[I] = relative_velocity[I]*dHnp1;
		//}
		////////////////////////////////
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;		
			elementJacobian_u_u[i][j] +=
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  //+ ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
			  //			    &u_grad_test_dV[i_nSpace])
			  + ck.NumericalDiffusionJacobian(epsDiffusion/dt,
							  &u_grad_trial[j_nSpace],
							  &u_grad_test_dV[i_nSpace]);
		      }//j
		  }//i
	      }//k
	    //
	    //load into element Jacobian into global Jacobian
	    //
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		      elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements
      }//computeJacobian for MCorr with CLSVOF

      void calculateJacobian_experimental(//element
					  double dt,
					  double* mesh_trial_ref,
					  double* mesh_grad_trial_ref,
					  double* mesh_dof,
					  double* mesh_velocity_dof,
					  double MOVING_DOMAIN,
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
					  double useMetrics,
					  double alphaBDF,
					  //VRANS
					  const double* q_porosity,
					  //
					  int* u_l2g,
					  double* elementDiameter,
					  double* nodeDiametersArray,
					  int degree_polynomial,
					  double* u_dof,
					  double* u_dof_old,
					  double* velocity,
					  double* q_m_betaBDF, 
					  double* cfl,
					  int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					  double* globalJacobian,
					  int nExteriorElementBoundaries_global,
					  int* exteriorElementBoundariesArray,
					  int* elementBoundaryElementsArray,
					  int* elementBoundaryLocalElementBoundariesArray,
					  double* ebqe_velocity_ext,
					  //VRANS
					  const double* ebqe_porosity_ext,
					  //
					  int* isDOFBoundary_u,
					  double* ebqe_bc_u_ext,
					  int* isFluxBoundary_u,
					  double* ebqe_bc_flux_u_ext,
					  int* csrColumnOffsets_eb_u_u,
					  // FOR NONLINEAR CLSVOF; i.e., MCorr with CLSVOF
					  int timeOrder,
					  int timeStage,
					  int useFullNewton,
					  double epsFactHeaviside,
					  double epsFactDirac,
					  double epsFactDiffusion,
					  double* phin_dof,
					  double* phiHat_dof,
					  // interface locator
					  double norm_factor_lagged,
					  double* interface_locator_lagged,
					  // normal reconstruction
					  double* lumped_wx,
					  double* lumped_wy,
					  double* lumped_wz,
					  double* lumped_wx_tStar,
					  double* lumped_wy_tStar,
					  double* lumped_wz_tStar)
      {
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      for (int j=0;j<nDOF_trial_element;j++)
		elementJacobian_u_u[i][j]=0.0;	      
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
		//declare local storage
		register double u, grad_u[nSpace], phiHatnp1, u_grad_trial[nDOF_trial_element*nSpace],
		  normalReconstruction[nSpace], wx, wy, proj_times_grad_trial[nSpace],
		  df[nSpace], relative_velocity[nSpace],
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		  dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		  x,y,z,xt,yt,zt;
		//get jacobian, etc for mapping reference element
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		//get the physical integration weight
		dV = fabs(jacDet)*dV_ref[k];
		//get the trial function gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		// get the components of the normal reconstruction
		ck.valFromDOF(lumped_wx,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      wx);
		ck.valFromDOF(lumped_wy,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      wy);
		//get the solution 	
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_u);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		  }

		double gradu2 = 0;
		for(int I=0;I<nSpace;I++)
		  gradu2 += grad_u[I]*grad_u[I];
		double beta_norm_grad_u = betaNormGrad(gradu2,1E-15);
	      
		double lambda = epsFactDiffusion;
		//double lambda = epsFactDiffusion*elementDiameter[eN]/degree_polynomial/dt;

		// single potential
		//double coeff1FullNonlinear = 1.-1./beta_norm_grad_u; 
		//double coeff2FullNonlinear = 1./std::pow(beta_norm_grad_u,3);
	      
		// double potential
		double coeff1FullNonlinear=fmax(1E-10,
						1+2*std::pow(beta_norm_grad_u,2)-3*beta_norm_grad_u);
		double coeff2FullNonlinear=fmax(1E-10,
						4.-3./beta_norm_grad_u);

		double epsDirac = epsFactDirac*elementDiameter[eN]/degree_polynomial;
		double time_derivative_jacobian = smoothedDerSign(epsDirac,u)/dt; //der of sign
		//double time_derivative_jacobian = 1/dt;

		///////////////////////////////////
		// FOR IMPLICIT NONLINEAR CLSVOF // (implicit advection term)
		///////////////////////////////////
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;
		double dHnp1 = smoothedDerSign(epsDirac,u); // derivative of sign
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    df[I] = relative_velocity[I]*dHnp1;
		  }

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		double gamma = 1.0;//elementDiameter[eN]/degree_polynomial/norm_factor_lagged;		

		if (useFullNewton==false)
		  for(int i=0;i<nDOF_test_element;i++)
		    {
		      for(int j=0;j<nDOF_trial_element;j++)
			{
			  int j_nSpace = j*nSpace;
			  int i_nSpace = i*nSpace;		
			  elementJacobian_u_u[i][j] +=
			    time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			    //+ ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
			    //				  &u_grad_test_dV[i_nSpace])
			    + gamma*lambda*ck.NumericalDiffusionJacobian(1.0,
									 &u_grad_trial[j_nSpace],
									 &u_grad_test_dV[i_nSpace]);
			}//j
		    }//i
		else //use full newton method
		  for(int i=0;i<nDOF_test_element;i++)
		    {
		      for(int j=0;j<nDOF_trial_element;j++)
			{
			  int j_nSpace = j*nSpace;
			  int i_nSpace = i*nSpace;		
			  elementJacobian_u_u[i][j] +=
			    time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			    //+ ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],
			    //				  &u_grad_test_dV[i_nSpace])
			    + gamma*lambda*ck.NumericalDiffusionJacobian(coeff1FullNonlinear,
								   &u_grad_trial[j_nSpace],
								   &u_grad_test_dV[i_nSpace])
			    + gamma*lambda*( coeff2FullNonlinear*dV*
				       ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[i_nSpace])*
				       ck.NumericalDiffusion(1.0,grad_u,&u_grad_trial[j_nSpace]) );
			}//j
		    }//i
	      }//k
	    //
	    //load into element Jacobian into global Jacobian
	    //
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		      elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements
      }//computeJacobian for MCorr with CLSVOF

      /*
//      CODE TO CONSIDER ANISOTROPIC REGULARIZATION/PENALIZATION
	void calculateJacobian_experimental(//element
					double dt,
					double* mesh_trial_ref,
					double* mesh_grad_trial_ref,
					double* mesh_dof,
					double* mesh_velocity_dof,
					double MOVING_DOMAIN,
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
					double useMetrics,
					double alphaBDF,
					//VRANS
					const double* q_porosity,
					//
					int* u_l2g,
					double* elementDiameter,
					double* nodeDiametersArray,
					int degree_polynomial,
					double* u_dof,
					double* u_dof_old, 
					double* velocity,
					double* q_m_betaBDF, 
					double* cfl,
					int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					double* globalJacobian,
					int nExteriorElementBoundaries_global,
					int* exteriorElementBoundariesArray,
					int* elementBoundaryElementsArray,
					int* elementBoundaryLocalElementBoundariesArray,
					double* ebqe_velocity_ext,
					//VRANS
					const double* ebqe_porosity_ext,
					//
					int* isDOFBoundary_u,
					double* ebqe_bc_u_ext,
					int* isFluxBoundary_u,
					double* ebqe_bc_flux_u_ext,
					int* csrColumnOffsets_eb_u_u,
					// FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
					int timeOrder,
					int timeStage,
					int useFullNewton,
					double epsFactHeaviside,
					double epsFactDirac,
					double epsFactDiffusion,
					double* phin_dof,
					double* phiHat_dof,
					// interface locator
					double norm_factor_lagged,
					double* interface_locator_lagged,
					// normal reconstruction
					double* lumped_wx,
					double* lumped_wy,
					double* lumped_wz,
					double* lumped_wx_tStar,
					double* lumped_wy_tStar,
					double* lumped_wz_tStar)
      {
	double timeCoeff=1.0;
	if (timeOrder==2)
	  timeCoeff=0.5;
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      for (int j=0;j<nDOF_trial_element;j++)
		elementJacobian_u_u[i][j]=0.0;	      
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
		//declare local storage
		register double
		  u, grad_u[nSpace], u_grad_trial[nDOF_trial_element*nSpace],		  
		  normalReconstruction[nSpace],
		  qxn, qyn, qzn, qxnStar, qynStar, qznStar,		  
		  proj_times_grad_trial[nSpace], projector[nSpace][nSpace],		  
		  relative_velocity[nSpace], df[nSpace], 		  
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],		
		  u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
		  dV, x,y,z,xt,yt,zt,h_phi;
		//get jacobian, etc for mapping reference element
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
                ck.calculateH_element(eN,
                                      k,
                                      nodeDiametersArray,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      h_phi);		
		ck.calculateMappingVelocity_element(eN,
						    k,
						    mesh_velocity_dof,
						    mesh_l2g,
						    mesh_trial_ref,
						    xt,yt,zt);	      
		//get the physical integration weight
		dV = fabs(jacDet)*dV_ref[k];
		// get the components of the normal reconstruction
		ck.valFromDOF(lumped_wx,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxn);
		ck.valFromDOF(lumped_wy,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qyn);
		ck.valFromDOF(lumped_wz,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qzn);		
		ck.valFromDOF(lumped_wx_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qxnStar);
		ck.valFromDOF(lumped_wy_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qynStar);
		ck.valFromDOF(lumped_wz_tStar,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      qznStar);
		//get the trial function gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,u_grad_trial);
		//get the solution 	
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_u);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		  }

		/////////////////
		// MOVING MESH //
		/////////////////
		double mesh_velocity[3];
		mesh_velocity[0] = xt;
		mesh_velocity[1] = yt;
		mesh_velocity[2] = zt;

		double lambda = epsFactDiffusion;
		//double lambda = epsFactDiffusion*elementDiameter[eN]/degree_polynomial/dt;
		double epsDirac = epsFactDirac*elementDiameter[eN]/degree_polynomial;		
		double dHnp1 = smoothedDerSign(epsDirac,u); //derivative of sign
		
		for (int I=0;I<nSpace;I++)
		  {
		    relative_velocity[I] = (velocity[eN_k_nSpace+I]-MOVING_DOMAIN*mesh_velocity[I]);
		    df[I] = relative_velocity[I]*dHnp1;
		  }

		/////////////////////
		// TIME DERIVATIVE //
		/////////////////////
		double time_derivative_jacobian = dHnp1/dt;

		///////////////////
		// COMPUTE GAMMA // = h/norm_factor_lagged
		///////////////////
		//double gamma = interface_locator_lagged[eN] > 0 ? 1. : 0.;
		double gamma = (useMetrics*h_phi
				+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial/norm_factor_lagged;

		if (timeOrder == 2 && timeStage == 2)
		  {
		    normalReconstruction[0] = qxnStar;
		    normalReconstruction[1] = qynStar;
#if nSpace==3
		    normalReconstruction[2] = qznStar;
#endif
		  }
		else //timeOrder == 1 or timeStage==1
		  {
		    normalReconstruction[0] = qxn;
		    normalReconstruction[1] = qyn;
#if nSpace==3
		    normalReconstruction[2] = qzn;
#endif
		  }		  
		///////////////////////////////////
		// PROJECTOR TO NORMAL COMPONENT //
		///////////////////////////////////
		// NOTE: projector is build based on normalReconstruction
		//       which changes depending on timeOrder and timeStage		
		double norm2_q = 0;
		for (int I=0; I<nSpace; I++)
		  norm2_q += normalReconstruction[I]*normalReconstruction[I];
		double norm_q = std::sqrt(norm2_q);

		for (int I = 0; I < nSpace; ++I)
		  for (int J = 0; J < nSpace; ++J)
		    projector[I][J] =
		      (normalReconstruction[I]*normalReconstruction[J])/(norm2_q+1E-15);
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;
			Mult(projector,&u_grad_trial[j_nSpace],proj_times_grad_trial);
			
			elementJacobian_u_u[i][j] +=
			  // TIME DERIVATIVE 
			  time_derivative_jacobian*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i]
			  // IMPLICIT TERMS: ADVECTION, DIFFUSION
			  + timeCoeff*ck.AdvectionJacobian_weak(df,
								u_trial_ref[k*nDOF_trial_element+j],
								&u_grad_test_dV[i_nSpace])
			  + timeCoeff*gamma*lambda*
			  (ISOTROPIC_REGULARIZATION == 1 ? 
			    ck.NumericalDiffusionJacobian(1.0,
							  &u_grad_trial[j_nSpace],
							  &u_grad_test_dV[i_nSpace])
			   :
			   + ck.NumericalDiffusionJacobian(norm_q,
			   				   proj_times_grad_trial,
			   				   &u_grad_test_dV[i_nSpace])
			   + ck.NumericalDiffusionJacobian(1.0-norm_q,
			   				   &u_grad_trial[j_nSpace],
			   				   &u_grad_test_dV[i_nSpace]));
		      }//j
		  }//i
	      }//k
	    //
	    //load into element Jacobian into global Jacobian
	    //
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		      elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements
      }//computeJacobian for MCorr with CLSVOF      
       */

      void calculateMetricsAtEOS( //EOS=End Of Simulation
				 double* mesh_trial_ref, //
				 double* mesh_grad_trial_ref, //
				 double* mesh_dof, //
				 int* mesh_l2g, //
				 double* dV_ref, //
				 double* u_trial_ref, 
				 double* u_grad_trial_ref,
				 double* u_test_ref, //
				 //physics
				 int nElements_global, //
				 int* u_l2g, //
				 double* elementDiameter,
				 //double* nodeDiametersArray,
				 double degree_polynomial,
				  double epsFactHeaviside,
				 double* u_dof,
				 double* u0_dof,
				 double* u_exact,
				 int offset_u, int stride_u, 
				 double* global_I_err,
				 double* global_Ieps_err,
				 double* global_V_err,
				 double* global_Veps_err,
				 double* global_D_err)
      {
	double global_V = 0.;
	double global_V0 = 0.;
	double global_sV = 0.;
	double global_sV0 = 0.;
	*global_I_err = 0.0;
	*global_Ieps_err = 0.0;
	*global_V_err = 0.0;
	*global_Veps_err = 0.0;
	*global_D_err = 0.0;
	//////////////////////
	// ** LOOP IN CELLS //
	//////////////////////
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double 
	      elementResidual_u[nDOF_test_element];
	    double cell_mass_error = 0., cell_mass_exact = 0.,
	      cell_I_err = 0., cell_Ieps_err = 0.,
	      cell_V = 0., cell_V0 = 0., cell_sV = 0., cell_sV0 = 0.,
	      cell_D_err = 0.;

	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  u, u0, uh,
		  u_grad_trial[nDOF_trial_element*nSpace],
		  grad_uh[nSpace], 
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		dV = fabs(jacDet)*dV_ref[k];
		// get functions at quad points
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      uh);
		ck.valFromDOF(u0_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u0);
		u = u_exact[eN_k];
		// get gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_uh);

		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		// compute (smoothed) heaviside functions //		
		double Hu0 = heaviside(u0);
		double Hu = heaviside(u);
		double Huh = heaviside(uh);
		double sHu0 = smoothedHeaviside(epsHeaviside,u0);
		double sHu = smoothedHeaviside(epsHeaviside,u);
		double sHuh = smoothedHeaviside(epsHeaviside,uh);

		// compute cell metrics //
		cell_I_err += fabs(Hu - Huh)*dV;
		cell_Ieps_err += fabs(sHu - sHuh)*dV;

		cell_V   += Huh*dV;
		cell_V0  += Hu0*dV;  
		cell_sV  += sHuh*dV;
		cell_sV0 += sHu0*dV; 

		double norm2_grad_uh = 0.;
		for (int I=0; I<nSpace; I++)
		  norm2_grad_uh += grad_uh[I]*grad_uh[I];
		cell_D_err += std::pow(std::sqrt(norm2_grad_uh) - 1, 2.)*dV;
	      }
	    global_V += cell_V;
	    global_V0 += cell_V0;
	    global_sV += cell_sV;
	    global_sV0 += cell_sV0;
	    // metrics //
	    *global_I_err    += cell_I_err;
	    *global_Ieps_err += cell_Ieps_err;
	    *global_D_err    += cell_D_err;	    
	  }//elements
	*global_V_err = fabs(global_V0 - global_V)/global_V0;
	*global_Veps_err = fabs(global_sV0 - global_sV)/global_sV0;
	*global_D_err *= 0.5;
      }

      void calculateMetricsAtETS( // ETS=Every Time Step
				 double dt,
				 double* mesh_trial_ref, //
				 double* mesh_grad_trial_ref, //
				 double* mesh_dof, //
				 int* mesh_l2g, //
				 double* dV_ref, //
				 double* u_trial_ref, 
				 double* u_grad_trial_ref,
				 double* u_test_ref, //
				 //physics
				 int nElements_global, //
				 int* u_l2g, //
				 double* elementDiameter,
				 //double* nodeDiametersArray,
				 double degree_polynomial,
				 double epsFactHeaviside,
				 double* u_dof,
				 double* u_dof_old,
				 double* u0_dof,
				 double* velocity,
				 int offset_u, int stride_u,
				 int numDOFs,				 
				 double* global_R,
				 double* global_Reps,				 
				 double* global_V_err,
				 double* global_Veps_err,
				 double* global_D_err)
      {
	register double R_vector[numDOFs], Reps_vector[numDOFs];
	for (int i=0; i<numDOFs; i++)
	  {
	    R_vector[i] = 0.;
	    Reps_vector[i] = 0.;
	  }
	
	double global_V = 0.;
	double global_V0 = 0.;
	double global_sV = 0.;
	double global_sV0 = 0.;
	*global_R = 0.0;
	*global_Reps = 0.0;
	*global_V_err = 0.0;
	*global_Veps_err = 0.0;
	*global_D_err = 0.0;
	//////////////////////////////////////////////
	// ** LOOP IN CELLS FOR CELL BASED TERMS ** //
	//////////////////////////////////////////////
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double element_R[nDOF_test_element], element_Reps[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		element_R[i] = 0.;
		element_Reps[i] = 0.;
	      }
	    
	    double cell_mass_error = 0., cell_mass_exact = 0.,
	      cell_R = 0., cell_Reps = 0.,
	      cell_V = 0., cell_V0 = 0., cell_sV = 0., cell_sV0 = 0.,
	      cell_D_err = 0.;

	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  unp1, un, u0,
		  grad_unp1[nSpace], sFlux_np1[nSpace], Flux_np1[nSpace],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  u_test_dV[nDOF_trial_element], 
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		dV = fabs(jacDet)*dV_ref[k];
		// get functions at quad points
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      unp1);
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      un);
		ck.valFromDOF(u0_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u0);     
		// get gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_unp1);

		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		
		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		// compute (smoothed) heaviside functions //		
		double Hu0 = heaviside(u0);
		double Hunp1 = heaviside(unp1);
		double sHu0 = smoothedHeaviside(epsHeaviside,u0);
		double sHunp1 = smoothedHeaviside(epsHeaviside,unp1);

		// compute cell metrics //
		cell_V   += Hunp1*dV;
		cell_V0  += Hu0*dV;  
		cell_sV  += sHunp1*dV;
		cell_sV0 += sHu0*dV; 

		double norm2_grad_unp1 = 0.;
		for (int I=0; I<nSpace; I++)
		  norm2_grad_unp1 += grad_unp1[I]*grad_unp1[I];
		cell_D_err += std::pow(std::sqrt(norm2_grad_unp1) - 1, 2.)*dV;

		double Sunp1 = sign(unp1);
		double Sun = sign(un);
		double sSunp1 = smoothedSign(epsHeaviside,unp1);
		double sSun = smoothedSign(epsHeaviside,un);
		for (int I=0; I<nSpace; I++)
		  {
		    Flux_np1[I] = velocity[eN_k_nSpace+I]*Sunp1;
		    sFlux_np1[I] = velocity[eN_k_nSpace+I]*sSunp1;
		  }
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    element_R[i] += ((Sunp1-Sun)/dt*u_test_dV[i]
				     + ck.Advection_weak(Flux_np1,&u_grad_test_dV[i_nSpace]));
		    element_Reps[i] += ((sSunp1-sSun)/dt*u_test_dV[i]
					+ ck.Advection_weak(sFlux_np1,&u_grad_test_dV[i_nSpace]));
		  }
	      }
	    // DISTRIBUTE //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index	      
		R_vector[gi] += element_R[i];
		Reps_vector[gi] += element_Reps[i];
	      }
	    global_V += cell_V;
	    global_V0 += cell_V0;
	    global_sV += cell_sV;
	    global_sV0 += cell_sV0;	    
	    // metrics //
	    *global_D_err    += cell_D_err;	    
	  }//elements

	for (int i=0; i<numDOFs; i++)
	  {
	    *global_R += R_vector[i]*R_vector[i];
	    *global_Reps += Reps_vector[i]*Reps_vector[i];
	  }
	
	*global_V_err = fabs(global_V0 - global_V)/global_V0;
	*global_Veps_err = fabs(global_sV0 - global_sV)/global_sV0;
	*global_D_err *= 0.5;
      }      

      /*
      void calculateMetrics(
			   double* mesh_trial_ref, //
			   double* mesh_grad_trial_ref, //
			   double* mesh_dof, //
			   int* mesh_l2g, //
			   double* dV_ref, //
			   double* u_trial_ref, 
			   double* u_grad_trial_ref,
			   double* u_test_ref, //
			   //physics
			   int nElements_global, //
			   int* u_l2g, //
			   double* elementDiameter,
			   //double* nodeDiametersArray,
			   double degree_polynomial,
			   double epsFactHeaviside,
			   double* u_dof,
			   double* phiHat_dof,
			   double* phiExact_dof,
			   int offset_u, int stride_u, 
			   double* globalResidual,
			   double* global_mass_error,
			   double* global_L2_interface,
			   double* global_H1_interface,
			   double* global_L2_Hinterface,
			   double* global_H1_Hinterface,
			   double* global_L2_u,
			   double* global_H1_u)			 
      {
	double global_mass_exact = 0.0;
	*global_mass_error = 0.0;
	*global_L2_interface = 0.0;
	*global_H1_interface = 0.0;
	*global_L2_Hinterface = 0.0;
	*global_H1_Hinterface = 0.0;
	*global_L2_u = 0.0;
	*global_H1_u = 0.0;
	//////////////////////////////////////////////
	// ** LOOP IN CELLS FOR CELL BASED TERMS ** //
	//////////////////////////////////////////////
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double 
	      elementResidual_u[nDOF_test_element];
	    double cell_mass_error = 0., cell_mass_exact = 0.,
	      cell_L2_u = 0., cell_L2_int = 0., cell_L2_Hint = 0.,
	      cell_H1_u = 0., cell_H1_int = 0., cell_H1_Hint = 0.;
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementResidual_u[i]=0.0;
	      }
	  
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  u, phiHat, phiExact,
		  u_grad_trial[nDOF_trial_element*nSpace],
		  grad_u[nSpace], grad_phiHat[nSpace], grad_phiExact[nSpace],
		  grad_int[nSpace], grad_Hint[nSpace],
		  u_test_dV[nDOF_trial_element],
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		dV = fabs(jacDet)*dV_ref[k];
		// get functions at quad points
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      u);
		ck.valFromDOF(phiHat_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      phiHat);
		ck.valFromDOF(phiExact_dof,
			      &u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],
			      phiExact);
		// get gradients
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
		ck.gradFromDOF(phiHat_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_phiHat);
		ck.gradFromDOF(phiExact_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_phiExact);
		      
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

		double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;
		double epsDirac = epsFactHeaviside*elementDiameter[eN]/degree_polynomial;

		// gradients of errors in the interfaces
		for (int I=0;I<nSpace;I++)
		  {
		    grad_int[I] = grad_phiHat[I]+grad_u[I] - grad_phiExact[I];
		    grad_Hint[I] = (smoothedDirac(epsDirac,phiHat+u)*(grad_phiHat[I]+grad_u[I])
				    -smoothedDirac(epsDirac,phiExact)*(grad_phiExact[I]));
		  }
		// cell mass error
		cell_mass_error += smoothedHeaviside(epsHeaviside,phiHat+u)*dV;
		cell_mass_exact += smoothedHeaviside(epsHeaviside,phiExact)*dV;
		// L2 component
		double L2_u_tmp = u*u*dV;
		double L2_int_tmp = std::pow(phiHat+u - phiExact,2)*dV; 
		double L2_Hint_tmp = std::pow(smoothedHeaviside(epsHeaviside,phiHat+u) -
					      smoothedHeaviside(epsHeaviside,phiExact),2)*dV;
		// H1 Semi norm component
		double H1Semi_u_tmp = ck.NumericalDiffusion(dV,grad_u,grad_u);
		double H1Semi_int_tmp = ck.NumericalDiffusion(dV,grad_int,grad_int);
		double H1Semi_Hint_tmp = ck.NumericalDiffusion(dV,grad_Hint,grad_Hint);
		// cell L2 norms
		cell_L2_u    += L2_u_tmp;
		cell_L2_int  += L2_int_tmp;
		cell_L2_Hint += L2_Hint_tmp;
		// cell H1 norms
		cell_H1_u    += L2_u_tmp    + H1Semi_u_tmp;
		cell_H1_int  += L2_int_tmp  + H1Semi_int_tmp;
		cell_H1_Hint += L2_Hint_tmp + H1Semi_Hint_tmp;
	      
		// ith-LOOP //	      
		for(int i=0;i<nDOF_test_element;i++)
		  elementResidual_u[i] += smoothedHeaviside(epsHeaviside,phiHat+u)*u_test_dV[i];
	      }
	    *global_mass_error += cell_mass_error;
	    global_mass_exact += cell_mass_exact;
	    *global_L2_interface += cell_L2_int;
	    *global_H1_interface += cell_H1_int;
	    *global_L2_Hinterface += cell_L2_Hint;
	    *global_H1_Hinterface += cell_H1_Hint;
	    *global_L2_u += cell_L2_u;
	    *global_H1_u += cell_H1_u;
	    /////////////////
	    // DISTRIBUTE // load cell based element into global residual
	    ////////////////
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		// distribute global residual
		globalResidual[gi] += elementResidual_u[i];
	      }//i
	  }//elements
	*global_mass_error -= global_mass_exact;
      }
       */
    
      void normalReconstruction(//element
				double* mesh_trial_ref,//
				double* mesh_grad_trial_ref,
				double* mesh_dof, //
				int* mesh_l2g,//
				double* dV_ref,//
				double* u_trial_ref,
				double* u_grad_trial_ref,
				double* u_test_ref,
				//physics
				int nElements_global,//
				int* u_l2g, //
				double* elementDiameter,//
				double* phi_dof,//
				int offset_u, int stride_u, 
				// PARAMETERS FOR EDGE VISCOSITY 
				int numDOFs,
				double* lumped_wx,
				double* lumped_wy,
				double* lumped_wz)
      //double* rhs_mass_correction,
      //double* lumped_L2p,
      //double* lumped_mass_matrix,
      // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
      //double epsFactHeaviside,
      //double epsFactDiffusion,
      //double* phiHat_dof)
      {
	register double
	  weighted_lumped_mass_matrix[numDOFs],
	  rhsx_normal_reconstruction[numDOFs],
	  rhsy_normal_reconstruction[numDOFs];
	for (int i=0; i<numDOFs; i++)
	  {
	    lumped_wx[i]=0.;
	    lumped_wy[i]=0.;
	    weighted_lumped_mass_matrix[i]=0.;
	    rhsx_normal_reconstruction[i]=0.;
	    rhsy_normal_reconstruction[i]=0.;
	  }
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for local contributions and initialize
	    register double
	      element_weighted_lumped_mass_matrix[nDOF_test_element],
	      element_rhsx_normal_reconstruction[nDOF_test_element],
	      element_rhsy_normal_reconstruction[nDOF_test_element];	    
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		element_weighted_lumped_mass_matrix[i]=0.0;
		element_rhsx_normal_reconstruction[i]=0.0;
		element_rhsy_normal_reconstruction[i]=0.0;
	      }
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double 
		  //for mass matrix contributions
		  grad_phi[nSpace],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_test_dV[nDOF_trial_element], 
		  //for general use
		  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
		  dV,x,y,z;
		//get the physical integration weight
		ck.calculateMapping_element(eN,
					    k,
					    mesh_dof,
					    mesh_l2g,
					    mesh_trial_ref,
					    mesh_grad_trial_ref,
					    jac,
					    jacDet,
					    jacInv,
					    x,y,z);
		dV = fabs(jacDet)*dV_ref[k];
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		ck.gradFromDOF(phi_dof,
			       &u_l2g[eN_nDOF_trial_element],u_grad_trial,
			       grad_phi);	      
		//precalculate test function products with integration weights for mass matrix terms
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	      
		double rhsx = grad_phi[0];
		double rhsy = grad_phi[1];
		double grad_phi2 = 0;
		for (int I=0;I<nSpace; I++)
		  grad_phi2 += grad_phi[I]*grad_phi[I];
		double beta_norm_grad_phi = betaNormGrad(grad_phi2,1E-10);
	      
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    element_weighted_lumped_mass_matrix[i] += beta_norm_grad_phi*u_test_dV[i];
		    element_rhsx_normal_reconstruction[i] += rhsx*u_test_dV[i];
		    element_rhsy_normal_reconstruction[i] += rhsy*u_test_dV[i];
		  }
	      } //k
	    // DISTRIBUTE //
	    for(int i=0;i<nDOF_test_element;i++) 
	      { 
		int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index	      
		weighted_lumped_mass_matrix[gi] += element_weighted_lumped_mass_matrix[i];
		rhsx_normal_reconstruction[gi] += element_rhsx_normal_reconstruction[i];
		rhsy_normal_reconstruction[gi] += element_rhsy_normal_reconstruction[i];
	      }//i
	  }//elements
	// COMPUTE LUMPED L2 PROJECTION
	for (int i=0; i<numDOFs; i++)
	  {	 
	    double mi = weighted_lumped_mass_matrix[i];
	    lumped_wx[i] = 1./mi*rhsx_normal_reconstruction[i];
	    lumped_wy[i] = 1./mi*rhsy_normal_reconstruction[i];
	  }
      }
    
    };//CLSVOF

  inline CLSVOF_base* newCLSVOF(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<CLSVOF_base,CLSVOF,CompKernel>(nSpaceIn,
										       nQuadraturePoints_elementIn,
										       nDOF_mesh_trial_elementIn,
										       nDOF_trial_elementIn,
										       nDOF_test_elementIn,
										       nQuadraturePoints_elementBoundaryIn,
										       CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<CLSVOF_base,CLSVOF,CompKernel>(nSpaceIn,
										     nQuadraturePoints_elementIn,
										     nDOF_mesh_trial_elementIn,
										     nDOF_trial_elementIn,
										     nDOF_test_elementIn,
										     nQuadraturePoints_elementBoundaryIn,
										     CompKernelFlag);
  }
}//proteus
#endif
