#ifndef CLSVOF_H
#define CLSVOF_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"

#define USE_SIGN_FUNCTION 1
#define IMPLICIT_BCs 0
#define LAMBDA_SCALING 0

namespace proteus
{
// True characteristic functions
  inline double heaviside(const double& z){
    return (z>0 ? 1. : (z<0 ? 0. : 0.5));
  }
  inline double Sign(const double& z){
    return (z>0 ? 1. : (z<0 ? -1. : 0.));
  }
}

namespace proteus
{
  class CLSVOF_base
  {
    //The base class defining the interface
  public:
    std::valarray<double> Rpos, Rneg, FluxCorrectionMatrix;
    virtual ~CLSVOF_base(){}
    virtual void calculateResidual(//element
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
                                   int nElements_owned,
                                   double useMetrics,
                                   //VRANS
                                   const double* q_vos,
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
				   double* q_n,
				   double* q_H,
				   double* q_mH,
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
                                   const double* ebqe_vos_ext,
                                   //
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isFluxBoundary_u,
                                   double* ebqe_bc_flux_u_ext,
                                   double* ebqe_u,
				   double* ebqe_n,
				   double* ebqe_H,
                                   double* ebqe_flux,
                                   // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
                                   int timeOrder,
                                   int timeStage,
                                   double epsFactHeaviside,
                                   double epsFactDirac,
				   double epsFactRedist,
                                   double lambdaFact,
                                   // normalization factor
                                   double* min_distance,
                                   double* max_distance,
                                   double* mean_distance,
				   double* volume_domain,
                                   double norm_factor_lagged,
				   double VelMax,
                                   // normal reconstruction
                                   double* projected_qx_tn,
                                   double* projected_qy_tn,
                                   double* projected_qz_tn,
                                   double* projected_qx_tStar,
                                   double* projected_qy_tStar,
                                   double* projected_qz_tStar,
				   // TO COMPUTE H
				   int numDOFs,
				   double* lumped_mass_matrix,
				   double* H_dof,
				   // PRE REDISTANCING
				   int preRedistancingStage,
				   double* interface_locator,
				   double alpha)=0;
    virtual void calculateJacobian(//element
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
                                   //VRANS
                                   const double* q_vos,
                                   //
                                   int* u_l2g,
                                   double* elementDiameter,
                                   double* nodeDiametersArray,
                                   int degree_polynomial,
                                   double* u_dof,
                                   double* u_dof_old,
                                   double* velocity,
                                   double* cfl,
                                   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                                   double* globalJacobian,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   //VRANS
                                   const double* ebqe_vos_ext,
                                   //
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isFluxBoundary_u,
                                   double* ebqe_bc_flux_u_ext,
                                   int* csrColumnOffsets_eb_u_u,
                                   // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
                                   int timeOrder,
                                   int timeStage,
                                   double epsFactHeaviside,
                                   double epsFactDirac,
				   double epsFactRedist,
                                   double lambdaFact,
                                   // normalization factor
				   int preRedistancingStage,
                                   double norm_factor_lagged,
				   double alpha)=0;
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
                                       int nElements_owned,
                                       int useMetrics,
                                       int* u_l2g,
                                       double* elementDiameter,
                                       double* nodeDiametersArray,
                                       double degree_polynomial,
                                       double epsFactHeaviside,
                                       double* u_dof,
                                       double* u0_dof,
                                       double* u_exact,
                                       int offset_u, int stride_u,
                                       double* global_I_err,
                                       double* global_sI_err,
                                       double* global_V,
                                       double* global_V0,
                                       double* global_sV,
                                       double* global_sV0,
                                       double* global_D_err,
				       double* global_L2_err,
				       double* global_L2Banded_err,
				       double* global_area_band,
				       double* global_sH_L2_err)=0;
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
                                       int nElements_owned,
                                       int useMetrics,
				       double* q_vos,
                                       int* u_l2g,
                                       double* elementDiameter,
                                       double* nodeDiametersArray,
                                       double degree_polynomial,
                                       double epsFactHeaviside,
                                       double* u_dof,
                                       double* u_dof_old,
                                       double* u0_dof,
                                       double* velocity,
                                       int offset_u, int stride_u,
                                       int numDOFs,
                                       double* R_vector,
                                       double* sR_vector,
                                       double* global_V,
                                       double* global_V0,
                                       double* global_sV,
                                       double* global_sV0,
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
                                      double* u_dof,
                                      int offset_u, int stride_u,
                                      int numDOFs,
				      double* weighted_lumped_mass_matrix,
				      // normal reconstruction via consistent mass matrix
				      double* rhs_qx,
				      double* rhs_qy,
				      double* rhs_qz,
				      int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				      double* weighted_mass_matrix)=0;
    virtual void calculateRhsL2Proj(double* mesh_trial_ref,
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
				    double he_for_disc_ICs,
				    double* u_dof,
				    int offset_u, int stride_u,
				    int numDOFs,
				    double* rhs_qx)=0;
    virtual void calculateLumpedMassMatrix(double* mesh_trial_ref,
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
					   double* lumped_mass_matrix,
					   int offset_u, int stride_u)=0;
    virtual void assembleSpinUpSystem(//element
				      double* mesh_trial_ref,
				      double* mesh_grad_trial_ref,
				      double* mesh_dof,
				      int* mesh_l2g,
				      double* dV_ref,
				      double* u_trial_ref,
				      double* u_test_ref,
				      //physics
				      int nElements_global,
				      int* u_l2g,
				      double* uInitial,
				      int offset_u, int stride_u,
				      double* globalResidual,
				      int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				      double* globalMassMatrix)=0;
    virtual void FCTStep(int NNZ, //number on non-zero entries on sparsity pattern
                         int numDOFs, //number of DOFs
                         double* lumped_mass_matrix, //lumped mass matrix (as vector)
			 double* soln,
                         double* solH, //DOFs of high order solution at tnp1
                         double* solL,
                         double* limited_solution,
                         int* csrRowIndeces_DofLoops, //csr row indeces
                         int* csrColumnOffsets_DofLoops, //csr column offsets
                         double* MassMatrix //mass matrix
                         )=0;
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

      inline void calculateNonlinearCFL(const double& elementDiameter,
					const double df[nSpace],
					const double norm_factor_lagged,
					const double epsFactHeaviside,
					const double lambdaFact,
					double& cfl)
      {
        double h,nrm2_v;
        h = elementDiameter;
        nrm2_v=0.0;
        for(int I=0;I<nSpace;I++)
          nrm2_v+=df[I]*df[I];
        cfl = nrm2_v*norm_factor_lagged/(epsFactHeaviside*lambdaFact*h*h*h);
	cfl = std::sqrt(cfl);
      }

      inline void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_u,
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
        if (isDOFBoundary_u == 1)
          {
            if (flow >= 0.0)
              {
                flux = u*flow;
              }
            else
              {
		flux = bc_u*flow;
              }
          }
        else if (isFluxBoundary_u == 1)
          {
            flux = bc_flux_u;
          }
        else
          {
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
      }

      inline void exteriorNumericalAdvectiveFluxDerivative(const int& isDOFBoundary_u,
							   const int& isFluxBoundary_u,
							   const double n[nSpace],
							   const double& du,
							   const double velocity[nSpace],
							   double& dflux)
      {
	double flow=0.0;
	for (int I=0; I < nSpace; I++)
	  flow += n[I]*velocity[I];
	dflux=0.0;//default to no flux
	if (isDOFBoundary_u == 1)
	  {
	    if (flow >= 0.0)
	      {
		dflux = du*flow;
	      }
	    else
	      {
		dflux = 0.0; //zero since inflow BC is given by data (so independent on soln)
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
	    else
	      {
		dflux = 0.;
	      }
	  }
      }

      inline double smoothedHeaviside(double eps, double u)
      {
        double H;
        if (u > eps)
          H=1.0;
        else if (u < -eps)
          H=0.0;
        else if (u==0.0)
          H=0.5;
        else
          H = 0.5*(1.0 + u/eps + sin(M_PI*u/eps)/M_PI);
        return H;
      }

      inline double smoothedDirac(double eps, double u)
      {
        double d;
        if (u > eps)
          d=0.0;
        else if (u < -eps)
          d=0.0;
        else
          d = 0.5*(1.0 + cos(M_PI*u/eps))/eps;
        return d;
      }

      inline double smoothedNormalizedDirac(double eps, double u)
      {
        double d;
        if (u > eps)
          d=0.0;
        else if (u < -eps)
          d=0.0;
        else
          d = 0.5*(1.0 + cos(M_PI*u/eps));
        return d;
      }

      inline double smoothedSign(double eps, double u)
      {
        double H;
        if (u > eps)
          H=1.0;
        else if (u < -eps)
          H=0.0;
        else if (u==0.0)
          H=0.5;
        else
          H = 0.5*(1.0 + u/eps + sin(M_PI*u/eps)/M_PI);
	if (USE_SIGN_FUNCTION==1)
	  return 2*H-1;
	else
	  return H;
      }

      inline double smoothedDerivativeSign(double eps, double u)
      {
        double d;
        if (u > eps)
          d=0.0;
        else if (u < -eps)
          d=0.0;
        else
          d = 0.5*(1.0 + cos(M_PI*u/eps))/eps;
        if (USE_SIGN_FUNCTION==1)
	  return 2*d;
	else
	  return d;
      }

      void calculateResidual(//element
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
                             int nElements_owned,
                             double useMetrics,
                             //VRANS
                             const double* q_vos,
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
			     double* q_n,
			     double* q_H,
			     double* q_mH,
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
                             const double* ebqe_vos_ext,
                             //
                             int* isDOFBoundary_u,
                             double* ebqe_bc_u_ext,
                             int* isFluxBoundary_u,
                             double* ebqe_bc_flux_u_ext,
                             double* ebqe_u,
			     double* ebqe_n,
			     double* ebqe_H,
                             double* ebqe_flux,
                             // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
                             int timeOrder,
                             int timeStage,
                             double epsFactHeaviside,
                             double epsFactDirac,
			     double epsFactRedist,
                             double lambdaFact,
                             // normalization factor
                             double* min_distance,
                             double* max_distance,
                             double* mean_distance,
			     double* volume_domain,
                             double norm_factor_lagged,
			     double VelMax,
                             // normal reconstruction
                             double* projected_qx_tn,
                             double* projected_qy_tn,
                             double* projected_qz_tn,
                             double* projected_qx_tStar,
                             double* projected_qy_tStar,
                             double* projected_qz_tStar,
			     // TO COMPUTE H
			     int numDOFs,
			     double* lumped_mass_matrix,
			     double* H_dof,
			     // PRE REDISTANCING
			     int preRedistancingStage,
			     double* interface_locator,
			     double alpha)
      {
        min_distance[0] = 1E10;
        max_distance[0] = -1E10;
        mean_distance[0] = 0.;
        volume_domain[0] = 0.;

        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double elementResidual_u[nDOF_test_element], element_rhs_L2proj_H[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
	      {
		elementResidual_u[i]=0.0;
		element_rhs_L2proj_H[i]=0.0;
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
                  u, un, grad_u[nSpace], grad_un[nSpace],
		  qxn, qyn, qzn, //qxnStar, qynStar, qznStar,
                  normalReconstruction[3], // assume 3D always
                  u_test_dV[nDOF_trial_element],
                  u_grad_trial[nDOF_trial_element*nSpace],
                  u_grad_test_dV[nDOF_test_element*nSpace],
                  //for general use
                  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                  dV,x,y,z,xt,yt,zt,h_phi,
		  porosity;
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
				    jacInv,
				    u_grad_trial);
                // get the components of the normal reconstruction //
                ck.valFromDOF(projected_qx_tn,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      qxn);
                ck.valFromDOF(projected_qy_tn,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      qyn);
                ck.valFromDOF(projected_qz_tn,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      qzn);
                //ck.valFromDOF(projected_qx_tStar,
		//	      &u_l2g[eN_nDOF_trial_element],
		//	      &u_trial_ref[k*nDOF_trial_element],
		//	      qxnStar);
                //ck.valFromDOF(projected_qy_tStar,
		//	      &u_l2g[eN_nDOF_trial_element],
		//	      &u_trial_ref[k*nDOF_trial_element],
		//	      qynStar);
                //ck.valFromDOF(projected_qz_tStar,
		//	      &u_l2g[eN_nDOF_trial_element],
		//	      &u_trial_ref[k*nDOF_trial_element],
		//	      qznStar);
                // get the solution (of Newton's solver)
                ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      u);
                // get old solution
                ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un);
                //get the solution gradients at quad points
                ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_u);
                ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_un);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
                  }
		double hK=(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial;
		//VRANS
		porosity=1.0-q_vos[eN_k];

		///////////////////////////
                // NORMAL RECONSTRUCTION //
                ///////////////////////////
                //if (timeOrder == 2 && timeStage == 2) // this is never reached
		//{
		//normalReconstruction[0] = 0.5*(qxnStar+qxn);
		//normalReconstruction[1] = 0.5*(qynStar+qyn);
		//if (nSpace==3)
		//normalReconstruction[2] = 0.5*(qznStar+qzn);
		//else
		//normalReconstruction[2] = 0.;
		//}
		//else //timeOrder == 1 or timeStage==1
		//{
		normalReconstruction[0] = qxn;
		normalReconstruction[1] = qyn;
		if (nSpace==3)
		  normalReconstruction[2] = qzn;
		else
		  normalReconstruction[2] = 0.;
		//}

		/////////////////////////////////////////
		// ADJUSTMENT ON dV DUE TO MESH MOTION //
		/////////////////////////////////////////
		if (q_dV_last[eN_k] <= -100)
		  q_dV_last[eN_k] = dV;
		q_dV[eN_k] = dV;

		double delta, residualEikonal, tau, backgroundDissipation=0.1*hK;
		double time_derivative_residual, fnHalf[nSpace], lambda, Hnp1;
		double sign;
		int same_sign=1;
		if (preRedistancingStage==1)
		  {
		    //////////////////////////////////////////////////////
		    // *************** EIKONAL EQUATION *************** //
		    //////////////////////////////////////////////////////
		    double norm_grad_u=0, norm_grad_un=0;
		    for (int I=0;I<nSpace; I++)
		      {
			norm_grad_u += grad_u[I]*grad_u[I];
			norm_grad_un += grad_un[I]*grad_un[I];
		      }
		    norm_grad_u = std::sqrt(norm_grad_u)+1E-10;
		    norm_grad_un = std::sqrt(norm_grad_un)+1E-10;
		    // residual of Eikonal equation
		    double epsRedist = epsFactRedist*hK;
		    double Si = -1.0+2.0*smoothedHeaviside(epsRedist,un);
		    residualEikonal = Si*(norm_grad_u-1.0);
		    delta = smoothedDirac(epsRedist,un);

		    // compute (lagged) velocity for redistancing //
		    double Un[nSpace];
		    double normUn=0;
		    for (int I=0; I < nSpace; I++)
		      {
			Un[I]  = Si*grad_un[I]/norm_grad_un;
			normUn += Un[I]*Un[I];
		      }
		    normUn = sqrt(normUn)+1E-10;
		    // compute coefficient for stabilization
		    tau = 0.5*hK;///normUn;
		  }
		else //clsvof model
		  {
		    //////////////////////////////////////////////////
		    // *************** CLSVOF MODEL *************** //
		    //////////////////////////////////////////////////
		    /////////////////
		    // MOVING MESH //
		    /////////////////
		    double mesh_velocity[3];
		    mesh_velocity[0] = xt;
		    mesh_velocity[1] = yt;
		    mesh_velocity[2] = zt;

		    ///////////////////
		    // GENERAL STUFF //
		    ///////////////////
		    double epsHeaviside = epsFactHeaviside*hK;
		    double Sn = smoothedSign(epsHeaviside,un);
		    double Snp1 = smoothedSign(epsHeaviside,u);
		    Hnp1 = smoothedHeaviside(epsHeaviside,u);

		    ////////////
		    // LAMBDA //
		    ////////////
		    lambda = lambdaFact*hK/norm_factor_lagged;
		    if (LAMBDA_SCALING==1)
		      {
			double deltaHat = fmax(smoothedNormalizedDirac(2*epsHeaviside,un),1E-6);
			lambda = lambdaFact*deltaHat;
		      }

		    /////////////////////
		    // TIME DERIVATIVE //
		    /////////////////////
		    time_derivative_residual = porosity*(Snp1-Sn)/dt;

		    ////////////////////
		    // ADVECTIVE TERM //
		    ////////////////////
		    double relative_velocity[nSpace], relative_velocity_old[nSpace];
		    for (int I=0;I<nSpace;I++)
		      {
			// compute relative velocity //
			relative_velocity[I] = (velocity[eN_k_nSpace+I]
						-MOVING_DOMAIN*mesh_velocity[I]);
			relative_velocity_old[I] = (velocity_old[eN_k_nSpace+I]
						    -MOVING_DOMAIN*mesh_velocity[I]);
			// compute advection term //
			//fnp1[I] = relative_velocity[I]*Snp1; //implicit advection via BDF1
			//fnHalf[I] = 0.5*(relative_velocity[I]*Snp1
			//		 +relative_velocity_old[I]*Sn); //implicit advection via CN
			fnHalf[I] = 0.5*relative_velocity[I]*porosity*(Snp1+Sn);
		      }

		    //////////////////////////////
		    // CALCULATE CELL BASED CFL //
		    //////////////////////////////
		    calculateCFL(elementDiameter[eN]/degree_polynomial,relative_velocity,cfl[eN_k]);
		    //calculateNonlinearCFL(elementDiameter[eN]/degree_polynomial,
		    //		      relative_velocity,
		    //		      norm_factor_lagged,
		    //		      epsFactHeaviside,
		    //		      lambdaFact,
		    //		      cfl[eN_k]);

		    ///////////////
		    // CALCULATE min, max and mean distance //
		    ///////////////
		    if (eN<nElements_owned) // locally owned?
		      {
			min_distance[0] = fmin(min_distance[0],fabs(u));
			max_distance[0] = fmax(max_distance[0],fabs(u));
			mean_distance[0] += fabs(u)*dV;
			//min_distance[0] = fmin(min_distance[0],u);
			//max_distance[0] = fmax(max_distance[0],u);
			//mean_distance[0] += u*dV;
			volume_domain[0] += dV;
		      }
		    ///////////////////
		    // SAVE SOLUTION // for other models
		    ///////////////////
		    q_u[eN_k] = u;
		    q_m[eN_k] = porosity*Hnp1;
		    q_H[eN_k] = Hnp1;
		    q_mH[eN_k] = porosity*Hnp1; //porosity*H(\phi)=(1-q_vos)*H(\phi)
		    // gradient //
		    for (int I=0;I<nSpace;I++)
		      q_n[eN_k_nSpace+I]  = grad_u[I];
		  }

		//////////////////////////////////////////////////
                // *************** LOOP ON DOFs *************** //
                //////////////////////////////////////////////////
                for(int i=0;i<nDOF_test_element;i++)
                  {
		    int eN_i=eN*nDOF_test_element+i;
		    int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
                    register int i_nSpace=i*nSpace;
		    if (preRedistancingStage==1)
		      {
			//////////////////////
			// EIKONAL RESIDUAL //
			//////////////////////
			elementResidual_u[i] +=
			  alpha*(u_dof[gi]-u_dof_old[gi])*delta*u_test_dV[i] // BCs
			  + residualEikonal*u_test_dV[i] // Eikonal eqn
			  // Dmitri's ~SUPG + background dissipation //
			  + ck.NumericalDiffusion(tau+backgroundDissipation,
						  grad_u,
						  &u_grad_test_dV[i_nSpace])
			  - ck.NumericalDiffusion(tau,
						  normalReconstruction,
						  &u_grad_test_dV[i_nSpace]);
		      }
		    else // clsvof
		      {
			//////////////////////////
			// LOCATE THE INTERFACE //
			//////////////////////////
			if (i==0)
			  sign = Sign(u_dof[gi]);
			else if (same_sign==1)
			  {
			    same_sign = sign == Sign(u_dof[gi]) ? 1 : 0;
			    sign = Sign(u_dof[gi]);
			  }
			////////////////////////////////////
			// for visualization of VOF field //
			////////////////////////////////////
			element_rhs_L2proj_H[i] += Hnp1*u_test_dV[i];
			//if (timeOrder==1)
			//{
			/////////////////////
			// CLSVOF RESIDUAL //
			/////////////////////
			elementResidual_u[i] +=
			  // TIME DERIVATIVE
			  time_derivative_residual*u_test_dV[i]
			  // ADVECTION TERM. This is IMPLICIT
			  + ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			  // REGULARIZATION TERM. This is IMPLICIT
			  + lambda*(ck.NumericalDiffusion(1.0,
							  grad_u,
							  &u_grad_test_dV[i_nSpace])
				    // TARGET for PENALIZATION. This is EXPLICIT
				    - ck.NumericalDiffusion(1.0,
							    normalReconstruction,
							    &u_grad_test_dV[i_nSpace]));
			//}
			//else // timeOrder=2
			//{
			//elementResidual_u[i] +=
			// TIME DERIVATIVE
			//time_derivative_residual*u_test_dV[i]
			// ADVECTION TERM. This is IMPLICIT
			//+ ck.Advection_weak(fnHalf,&u_grad_test_dV[i_nSpace])
			// REGULARIZATION TERM. This is IMPLICIT
			//+ lambda*ck.NumericalDiffusion(1.0,
			//			     grad_unHalf,
			//				     &u_grad_test_dV[i_nSpace])
			// TARGET for PENALIZATION. This is EXPLICIT
			//- lambda*ck.NumericalDiffusion(1.0,
			//			     normalReconstruction,
			//				     &u_grad_test_dV[i_nSpace]);
			//}
		      }
                  }//i
		if (preRedistancingStage==0)
		  {
		    if (same_sign == 0) // This cell contains the interface
		      {
			for(int i=0;i<nDOF_test_element;i++)
			  {
			    int gi = offset_u+stride_u*u_l2g[eN*nDOF_test_element+i];
			    interface_locator[gi] = 1.0;
			  }
		      }
		  }
              } //k
            /////////////////
            // DISTRIBUTE // load cell based element into global residual
            ////////////////
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
                // distribute global residual for (lumped) mass matrix
                globalResidual[gi] += elementResidual_u[i];
		if (preRedistancingStage==0)
		  H_dof[gi] += element_rhs_L2proj_H[i]; // int(H*wi*dx)
              }//i
          }//elements
	if (preRedistancingStage==0)
	  {
	    // COMPUTE LUMPED L2 PROJECTION
	    for (int i=0; i<numDOFs; i++)
	      H_dof[i] /= lumped_mass_matrix[i];
	  }

        //////////////
        // BOUNDARY //
        //////////////
        //ebNE is the Exterior element boundary INdex
        //ebN is the element boundary INdex
        //eN is the element index
	if (preRedistancingStage==0)
        for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
          {
            register int ebN = exteriorElementBoundariesArray[ebNE],
              eN  = elementBoundaryElementsArray[ebN*2+0],
              ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
              eN_nDOF_trial_element = eN*nDOF_trial_element;
            register double elementResidual_u[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                elementResidual_u[i]=0.0;
              }
            for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
              {
                register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                  ebNE_kb_nSpace = ebNE_kb*nSpace,
                  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                register double
                  u_ext=0.0,un_ext=0.0,grad_u_ext[nSpace],
                  df_ext[nSpace],
                  flux_ext=0.0,
                  bc_u_ext=0.0,
                  jac_ext[nSpace*nSpace],
                  jacDet_ext,
                  jacInv_ext[nSpace*nSpace],
                  boundaryJac[nSpace*(nSpace-1)],
                  metricTensor[(nSpace-1)*(nSpace-1)],
                  metricTensorDetSqrt,
                  dS,
                  u_test_dS[nDOF_test_element],
		  u_grad_trial_trace[nDOF_trial_element*nSpace],
                  normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                  //VRANS
                  porosity_ext;
                //
                //calculate the solution and gradients at quadrature points
                //
                //compute information about mapping from reference element to physical element
                ck.calculateMapping_elementBoundary(eN,
                                                    ebN_local,
                                                    kb,
                                                    ebN_local_kb,
                                                    mesh_dof,
                                                    mesh_l2g,
                                                    mesh_trial_trace_ref,
                                                    mesh_grad_trial_trace_ref,
                                                    boundaryJac_ref,
                                                    jac_ext,
                                                    jacDet_ext,
                                                    jacInv_ext,
                                                    boundaryJac,
                                                    metricTensor,
                                                    metricTensorDetSqrt,
                                                    normal_ref,
                                                    normal,
                                                    x_ext,y_ext,z_ext);
                ck.calculateMappingVelocity_elementBoundary(eN,
                                                            ebN_local,
                                                            kb,
                                                            ebN_local_kb,
                                                            mesh_velocity_dof,
                                                            mesh_l2g,
                                                            mesh_trial_trace_ref,
                                                            xt_ext,yt_ext,zt_ext,
                                                            normal,
                                                            boundaryJac,
                                                            metricTensor,
                                                            integralScaling);
                dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt
		      + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
		ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],
				    jacInv_ext,
				    u_grad_trial_trace);
                //solution at quad points
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],
			      u_ext);
                ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],
			      un_ext);
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial_trace,
			       grad_u_ext);
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                //
                //load the boundary values
                //
		double hK = elementDiameter[eN]/degree_polynomial;
		double epsHeaviside = epsFactHeaviside*hK;
		double Su_ext = smoothedSign(epsHeaviside,u_ext); //Sign(u_ext)
		double Sun_ext = smoothedSign(epsHeaviside,un_ext); //Sign(un_ext)
		// NOTE: ebqe_bc_u_ext is provided by the user as BCs for VOF (i.e., 0 or 1)
		// SuDBC = Sign(uBC)
		double SuBC = (USE_SIGN_FUNCTION == 0 ? ebqe_bc_u_ext[ebNE_kb] :
			       2*ebqe_bc_u_ext[ebNE_kb] - 1);
		if (IMPLICIT_BCs==1)
		  bc_u_ext = (isDOFBoundary_u[ebNE_kb]*SuBC
			      +(1-isDOFBoundary_u[ebNE_kb])*Su_ext);
		else
		  bc_u_ext = (isDOFBoundary_u[ebNE_kb]*SuBC
			      +(1-isDOFBoundary_u[ebNE_kb])*Sun_ext);
                //VRANS
                porosity_ext = 1.-ebqe_vos_ext[ebNE_kb];

                //
                //moving mesh
                //
                double mesh_velocity[3];
                mesh_velocity[0] = xt_ext;
                mesh_velocity[1] = yt_ext;
                mesh_velocity[2] = zt_ext;

                for (int I=0;I<nSpace;I++)
                  df_ext[I] = porosity_ext*(ebqe_velocity_ext[ebNE_kb_nSpace+I]
					    - MOVING_DOMAIN*mesh_velocity[I]);
                //
                //calculate the numerical fluxes
                //
                exteriorNumericalAdvectiveFlux(isDOFBoundary_u[ebNE_kb],
                                               isFluxBoundary_u[ebNE_kb],
                                               normal,
                                               bc_u_ext, //{-1,1} or {0,1}
                                               ebqe_bc_flux_u_ext[ebNE_kb],
					       IMPLICIT_BCs == 1 ? Su_ext : Sun_ext,
                                               df_ext, //VRANS includes porosity
                                               flux_ext);
                ebqe_flux[ebNE_kb] = flux_ext;

		///////////////////
		// save solution // for other models? cek need to be consistent with numerical flux
		///////////////////
		ebqe_u[ebNE_kb] = u_ext;
		// TODO: do I need ebqe_m?
		ebqe_H[ebNE_kb] = smoothedHeaviside(epsHeaviside,u_ext);
		// gradient //
		for (int I=0;I<nSpace;I++)
		  ebqe_n[ebNE_kb_nSpace+I]  = grad_u_ext[I];

                //
                //update residuals
                //
                for (int i=0;i<nDOF_test_element;i++)
                  {
                    elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]);
                  }//i
              }//kb
            //
            //update the element and global residual storage
            //
            for (int i=0;i<nDOF_test_element;i++)
              {
                int eN_i = eN*nDOF_test_element+i;
                globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
              }//i
          }//ebNE
        // END OF BOUNDARY //
      }

      void calculateJacobian(//element
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
                             //VRANS
                             const double* q_vos,
                             //
                             int* u_l2g,
                             double* elementDiameter,
                             double* nodeDiametersArray,
                             int degree_polynomial,
                             double* u_dof,
                             double* u_dof_old,
                             double* velocity,
                             double* cfl,
                             int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                             double* globalJacobian,
                             int nExteriorElementBoundaries_global,
                             int* exteriorElementBoundariesArray,
                             int* elementBoundaryElementsArray,
                             int* elementBoundaryLocalElementBoundariesArray,
                             double* ebqe_velocity_ext,
                             //VRANS
                             const double* ebqe_vos_ext,
                             //
                             int* isDOFBoundary_u,
                             double* ebqe_bc_u_ext,
                             int* isFluxBoundary_u,
                             double* ebqe_bc_flux_u_ext,
                             int* csrColumnOffsets_eb_u_u,
                             // FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
                             int timeOrder,
                             int timeStage,
                             double epsFactHeaviside,
                             double epsFactDirac,
			     double epsFactRedist,
                             double lambdaFact,
                             // normalization factor
			     int preRedistancingStage,
                             double norm_factor_lagged,
			     double alpha)
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
                register double
                  u, un, u_grad_trial[nDOF_trial_element*nSpace],
		  grad_u[nSpace], grad_un[nSpace],
                  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                  u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element*nSpace],
                  dV, x,y,z,xt,yt,zt,h_phi,
		  porosity;
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
				    jacInv,
				    u_grad_trial);
                //get the solution
                ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      u);
                ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un);
		//get the solution gradients
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_u);
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_un);
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
                  }
		double hK=(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial;
		//VRANS
		porosity=1.0-q_vos[eN_k];

		double delta, dH[nSpace], tau, backgroundDissipation=0.1*hK;
		double lambda, time_derivative_jacobian, df[nSpace];
		if (preRedistancingStage==1)
		  {
		    // ************************************************ //
		    // *************** EIKONAL EQUATION *************** //
		    // ************************************************ //
		    double norm_grad_u=0, norm_grad_un=0;
		    for (int I=0;I<nSpace; I++)
		      {
			norm_grad_u += grad_u[I]*grad_u[I];
			norm_grad_un += grad_un[I]*grad_un[I];
		      }
		    norm_grad_u = std::sqrt(norm_grad_u)+1E-10;
		    norm_grad_un = std::sqrt(norm_grad_un)+1E-10;
		    // derivative of residual of Eikonal equation //
		    double epsRedist = epsFactRedist*hK;
		    double Si = -1.0+2.0*smoothedHeaviside(epsRedist,un);
		    for (int I=0; I<nSpace;I++)
		      dH[I] = Si*grad_u[I]/norm_grad_u;
		    delta = smoothedDirac(epsRedist,un);

		    // compute lagged velocity of redistancing
		    double Un[nSpace];
		    double normUn = 0.;
		    for (int I=0; I < nSpace; I++)
		      {
			Un[I] = Si*grad_un[I]/norm_grad_un;
			normUn += Un[I]*Un[I];
		      }
		    normUn = sqrt(normUn)+1E-10;
		    // compute tau coefficient
		    tau = 0.5*hK;///normUn;
		  }
		else
		  {
		    // ******************************************** //
		    // *************** CLSVOF MODEL *************** //
		    // ******************************************** //
		    /////////////////
		    // MOVING MESH //
		    /////////////////
		    double mesh_velocity[3];
		    mesh_velocity[0] = xt;
		    mesh_velocity[1] = yt;
		    mesh_velocity[2] = zt;

		    ///////////////////
		    // GENERAL STUFF //
		    ///////////////////
		    double epsDirac = epsFactDirac*hK;
		    double dSnp1 = smoothedDerivativeSign(epsDirac,u); //derivative of smoothed sign

		    ////////////
		    // LAMBDA //
		    ////////////
		    lambda = lambdaFact*hK/norm_factor_lagged;
		    if (LAMBDA_SCALING==1)
		      {
			double deltaHat = fmax(smoothedNormalizedDirac(2*epsDirac,un),1E-6);
			lambda = lambdaFact*deltaHat;
		      }

		    /////////////////////
		    // TIME DERIVATIVE //
		    /////////////////////
		    time_derivative_jacobian = porosity*dSnp1/dt;

		    ////////////////////
		    // ADVECTIVE TERM //
		    ////////////////////
		    double relative_velocity[nSpace];
		    for (int I=0;I<nSpace;I++)
		      {
			relative_velocity[I] = (velocity[eN_k_nSpace+I]
						-MOVING_DOMAIN*mesh_velocity[I]);
			df[I] = relative_velocity[I]*porosity*dSnp1;
		      }
		  }

                //////////////////
                // LOOP ON DOFs //
                //////////////////
                for(int i=0;i<nDOF_test_element;i++)
                  {
                    for(int j=0;j<nDOF_trial_element;j++)
                      {
                        int j_nSpace = j*nSpace;
                        int i_nSpace = i*nSpace;

			if (preRedistancingStage==1)
			  {
			    //////////////////////
			    // EIKONAL JACOBIAN //
			    //////////////////////
			    elementJacobian_u_u[i][j] +=
			      (i == j ? alpha*delta*u_test_dV[i] : 0.) // BCs
			      + ck.HamiltonianJacobian_weak(dH, // Eikonal equation
							    &u_grad_trial[j_nSpace],
							    u_test_dV[i])
			      // Dmitri's ~SUPG + background dissipation //
			      +ck.NumericalDiffusionJacobian(tau+backgroundDissipation,
							     &u_grad_trial[j_nSpace],
							     &u_grad_test_dV[i_nSpace]);
			  }
			else // clsvof model
			  {
			    /////////////////////
			    // CLSVOF JACOBIAN //
			    /////////////////////
			    elementJacobian_u_u[i][j] +=
			      // TIME DERIVATIVE
			      time_derivative_jacobian*(u_trial_ref[k*nDOF_trial_element+j]
							*u_test_dV[i])
			      // IMPLICIT TERMS: ADVECTION, DIFFUSION
			      + 0.5*ck.AdvectionJacobian_weak(df,
							      u_trial_ref[k*nDOF_trial_element+j],
							      &u_grad_test_dV[i_nSpace])
			      + lambda*ck.NumericalDiffusionJacobian(1.0,
								     &u_grad_trial[j_nSpace],
								     &u_grad_test_dV[i_nSpace]);
			  }
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

        ///////////////////
        // BOUNDARY LOOP //
        ///////////////////
	if (IMPLICIT_BCs==1 && preRedistancingStage==0)
	for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
	  {
	    register int ebN = exteriorElementBoundariesArray[ebNE];
	    register int eN  = elementBoundaryElementsArray[ebN*2+0],
	      ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	      eN_nDOF_trial_element = eN*nDOF_trial_element;
	    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
	      {
		register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		  ebNE_kb_nSpace = ebNE_kb*nSpace,
		  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
		register double u_ext=0.0,
		  grad_u_ext[nSpace],
		  m_ext=0.0,
		  dm_ext=0.0,
		  f_ext[nSpace],
		  df_ext[nSpace],
		  dflux_u_u_ext=0.0,
		  //bc_grad_u_ext[nSpace],
		  bc_m_ext=0.0,
		  bc_dm_ext=0.0,
		  bc_f_ext[nSpace],
		  bc_df_ext[nSpace],
		  fluxJacobian_u_u[nDOF_trial_element],
		  jac_ext[nSpace*nSpace],
		  jacDet_ext,
		  jacInv_ext[nSpace*nSpace],
		  boundaryJac[nSpace*(nSpace-1)],
		  metricTensor[(nSpace-1)*(nSpace-1)],
		  metricTensorDetSqrt,
		  dS,
		  u_test_dS[nDOF_test_element],
		  u_grad_trial_trace[nDOF_trial_element*nSpace],
		  normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		  //VRANS
		  porosity_ext;
		ck.calculateMapping_elementBoundary(eN,
						    ebN_local,
						    kb,
						    ebN_local_kb,
						    mesh_dof,
						    mesh_l2g,
						    mesh_trial_trace_ref,
						    mesh_grad_trial_trace_ref,
						    boundaryJac_ref,
						    jac_ext,
						    jacDet_ext,
						    jacInv_ext,
						    boundaryJac,
						    metricTensor,
						    metricTensorDetSqrt,
						    normal_ref,
						    normal,
						    x_ext,y_ext,z_ext);
		ck.calculateMappingVelocity_elementBoundary(eN,
							    ebN_local,
							    kb,
							    ebN_local_kb,
							    mesh_velocity_dof,
							    mesh_l2g,
							    mesh_trial_trace_ref,
							    xt_ext,yt_ext,zt_ext,
							    normal,
							    boundaryJac,
							    metricTensor,
							    integralScaling);
		dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt
		      + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
		ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],
				    jacInv_ext,
				    u_grad_trial_trace);
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],
			      u_ext);
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial_trace,
			       grad_u_ext);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		//
		//load the boundary values
		//
		double hK = elementDiameter[eN]/degree_polynomial;
		double epsHeaviside = epsFactHeaviside*hK;
		double dSu_ext = smoothedDerivativeSign(epsHeaviside,u_ext);
		//VRANS
		porosity_ext = 1.-ebqe_vos_ext[ebNE_kb];
		//
		//moving domain
		//
		double mesh_velocity[3];
		mesh_velocity[0] = xt_ext;
		mesh_velocity[1] = yt_ext;
		mesh_velocity[2] = zt_ext;
		//std::cout<<"ext J mesh_velocity"<<std::endl;
		for (int I=0;I<nSpace;I++)
		  df_ext[I] = porosity_ext*(ebqe_velocity_ext[ebNE_kb_nSpace+I]
					    - MOVING_DOMAIN*mesh_velocity[I]);
		//
		//calculate the numerical fluxes
		//
		exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
							 isFluxBoundary_u[ebNE_kb],
							 normal,
							 dSu_ext,
							 df_ext,//VRANS holds porosity
							 dflux_u_u_ext);
		//
		//calculate the flux jacobian
		//
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
		    register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		    fluxJacobian_u_u[j]=
		      ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,
								u_trial_trace_ref[ebN_local_kb_j]);
		  }//j
		//
		//update the global Jacobian from the flux Jacobian
		//
		for (int i=0;i<nDOF_test_element;i++)
		  {
		    register int eN_i = eN*nDOF_test_element+i;
		    //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
		    for (int j=0;j<nDOF_trial_element;j++)
		      {
			register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
			globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]]+=
			  fluxJacobian_u_u[j]*u_test_dS[i];
		      }//j
		  }//i
	      }//kb
	  }//ebNE
      }//computeJacobian for MCorr with CLSVOF

      void calculateMetricsAtEOS( //EOS=End Of Simulation
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
                                 int nElements_owned,
                                 int useMetrics,
                                 int* u_l2g,
                                 double* elementDiameter,
                                 double* nodeDiametersArray,
                                 double degree_polynomial,
                                 double epsFactHeaviside,
                                 double* u_dof,
                                 double* u0_dof,
                                 double* u_exact,
                                 int offset_u, int stride_u,
                                 double* global_I_err,
                                 double* global_sI_err,
                                 double* global_V,
                                 double* global_V0,
                                 double* global_sV,
                                 double* global_sV0,
                                 double* global_D_err,
				 double* global_L2_err,
				 double* global_L2Banded_err,
				 double* global_area_band,
				 double* global_sH_L2_err)
      {
        *global_I_err = 0.0;
        *global_sI_err = 0.0;
        *global_V = 0.0;
        *global_V0 = 0.0;
        *global_sV = 0.0;
        *global_sV0 = 0.0;
        *global_D_err = 0.0;
	*global_L2_err = 0.0;
	*global_L2Banded_err = 0.0;
	*global_area_band = 0.0;
	*global_sH_L2_err = 0.0;
        //////////////////////
        // ** LOOP IN CELLS //
        //////////////////////
        for(int eN=0;eN<nElements_global;eN++)
          {
            if (eN<nElements_owned) // just consider the locally owned cells
              {
                //declare local storage for local contributions and initialize
                double
		  cell_I_err = 0., cell_sI_err = 0.,
                  cell_V = 0., cell_V0 = 0., cell_sV = 0., cell_sV0 = 0.,
                  cell_D_err = 0.,
		  cell_L2_err = 0.,
		  cell_L2Banded_err = 0., cell_area_band = 0.,
		  cell_sH_L2_err = 0.;

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
                      dV,x,y,z,h_phi;
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
                    dV = fabs(jacDet)*dV_ref[k];
                    // get functions at quad points
                    ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],uh);
                    ck.valFromDOF(u0_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u0);
                    u = u_exact[eN_k];
                    // get gradients
                    ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                    ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_uh);

                    double epsHeaviside = epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial;
                    // compute (smoothed) heaviside functions //
                    double Hu0 = heaviside(u0);
                    double Hu = heaviside(u);
                    double Huh = heaviside(uh);
                    double sHu0 = smoothedHeaviside(epsHeaviside,u0);
                    double sHu = smoothedHeaviside(epsHeaviside,u);
                    double sHuh = smoothedHeaviside(epsHeaviside,uh);
		    //////////////////////////
                    // compute cell metrics //
		    //////////////////////////
		    // metrics on the interface
                    cell_I_err += fabs(Hu - Huh)*dV;
                    cell_sI_err += fabs(sHu - sHuh)*dV;
		    // L2 metrics on the level set
		    cell_L2_err += std::pow(u-uh,2)*dV;
		    if (fabs(uh) <= 2*epsHeaviside)
		      {
			cell_L2Banded_err += std::pow(u-uh,2)*dV;
			cell_area_band += dV;
		      }
		    // L2 metrics on the Heviside of the level set
		    cell_sH_L2_err += std::pow(sHu-sHuh,2)*dV;
		    // volume conservation
                    cell_V   += Huh*dV;
                    cell_V0  += Hu0*dV;
                    cell_sV  += sHuh*dV;
                    cell_sV0 += sHu0*dV;

                    double norm2_grad_uh = 0.;
                    for (int I=0; I<nSpace; I++)
                      norm2_grad_uh += grad_uh[I]*grad_uh[I];
                    cell_D_err += std::pow(std::sqrt(norm2_grad_uh) - 1, 2.)*dV;
                  }
                *global_V += cell_V;
                *global_V0 += cell_V0;
                *global_sV += cell_sV;
                *global_sV0 += cell_sV0;
                // metrics //
                *global_I_err    += cell_I_err;
                *global_sI_err += cell_sI_err;
                *global_D_err    += cell_D_err;
		*global_L2_err += cell_L2_err;
		*global_L2Banded_err += cell_L2Banded_err;
		*global_area_band += cell_area_band;
		*global_sH_L2_err += cell_sH_L2_err;
              }//elements
          }
        *global_D_err *= 0.5;
      }

      void calculateMetricsAtETS( // ETS=Every Time Step
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
                                 int nElements_owned,
                                 int useMetrics,
				 double* q_vos,
                                 int* u_l2g,
                                 double* elementDiameter,
                                 double* nodeDiametersArray,
                                 double degree_polynomial,
                                 double epsFactHeaviside,
                                 double* u_dof,
                                 double* u_dof_old,
                                 double* u0_dof,
                                 double* velocity,
                                 int offset_u, int stride_u,
                                 int numDOFs,
                                 double* R_vector,
                                 double* sR_vector,
                                 double* global_V,
                                 double* global_V0,
                                 double* global_sV,
                                 double* global_sV0,
                                 double* global_D_err)
      {
        *global_V = 0.0;
        *global_V0 = 0.0;
        *global_sV = 0.0;
        *global_sV0 = 0.0;
        *global_D_err = 0.0;
        //////////////////////////////////////////////
        // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
        //////////////////////////////////////////////
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double element_R[nDOF_test_element], element_sR[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                element_R[i] = 0.;
                element_sR[i] = 0.;
              }
            double
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
                  dV,x,y,z,h_phi,
		  porosity;
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
                dV = fabs(jacDet)*dV_ref[k];
                // get functions at quad points
                ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unp1);
                ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
                ck.valFromDOF(u0_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u0);
                // get gradients
                ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_unp1);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
                  }

		porosity = 1.0-q_vos[eN_k];

                double epsHeaviside = epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN])/degree_polynomial;
                // compute (smoothed) heaviside functions //
                double Hu0 = heaviside(u0);
                double Hunp1 = heaviside(unp1);
                double sHu0 = smoothedHeaviside(epsHeaviside,u0);
                double sHunp1 = smoothedHeaviside(epsHeaviside,unp1);

                // compute cell metrics //
                cell_V   += porosity*Hunp1*dV;
                cell_V0  += porosity*Hu0*dV;
                cell_sV  += porosity*sHunp1*dV;
                cell_sV0 += porosity*sHu0*dV;

                double norm2_grad_unp1 = 0.;
                for (int I=0; I<nSpace; I++)
                  norm2_grad_unp1 += grad_unp1[I]*grad_unp1[I];
                cell_D_err += std::pow(std::sqrt(norm2_grad_unp1) - 1, 2.)*dV;

                double Sunp1 = Sign(unp1);
                double Sun = Sign(un);
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
                    element_sR[i] += ((sSunp1-sSun)/dt*u_test_dV[i]
                                      + ck.Advection_weak(sFlux_np1,&u_grad_test_dV[i_nSpace]));
                  }
              }
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
                R_vector[gi] += element_R[i];
                sR_vector[gi] += element_sR[i];
              }
            if (eN<nElements_owned) // just consider the locally owned cells
              {
                *global_V += cell_V;
                *global_V0 += cell_V0;
                *global_sV += cell_sV;
                *global_sV0 += cell_sV0;
                *global_D_err    += cell_D_err;
              }
          }//elements
        *global_D_err *= 0.5;
      }

      void normalReconstruction(//element
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
                                double* u_dof,
                                int offset_u, int stride_u,
                                // PARAMETERS FOR EDGE VISCOSITY
                                int numDOFs,
				double* weighted_lumped_mass_matrix,
				// normal reconstruction via consistent mass matrix
				double* rhs_qx,
				double* rhs_qy,
				double* rhs_qz,
				int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				double* weighted_mass_matrix)
      {
        for (int i=0; i<numDOFs; i++)
          {
            weighted_lumped_mass_matrix[i]=0.;
	    rhs_qx[i]=0.;
	    rhs_qy[i]=0.;
	    rhs_qz[i]=0.;
          }
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double
              element_weighted_lumped_mass_matrix[nDOF_test_element],
              element_rhsx_normal_reconstruction[nDOF_test_element],
              element_rhsy_normal_reconstruction[nDOF_test_element],
              element_rhsz_normal_reconstruction[nDOF_test_element];
	    register double element_weighted_mass_matrix[nDOF_test_element][nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                element_weighted_lumped_mass_matrix[i]=0.0;
                element_rhsx_normal_reconstruction[i]=0.0;
                element_rhsy_normal_reconstruction[i]=0.0;
                element_rhsz_normal_reconstruction[i]=0.0;
		for (int j=0;j<nDOF_trial_element;j++)
		  element_weighted_mass_matrix[i][j]=0.0;
              }
            //loop over quadrature points and compute integrands
            for(int k=0;k<nQuadraturePoints_element;k++)
              {
                //compute indeces and declare local storage
                register int eN_k = eN*nQuadraturePoints_element+k,
                  eN_k_nSpace = eN_k*nSpace,
                  eN_nDOF_trial_element = eN*nDOF_trial_element;
                register double
                  //for mass matrix contributions
                  grad_u[nSpace],
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
                ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

                double rhsx = grad_u[0];
                double rhsy = grad_u[1];
                double rhsz = 0;
                if (nSpace==3)
                  rhsz = grad_u[2];

                double norm_grad_u = 0;
                for (int I=0;I<nSpace; I++)
                  norm_grad_u += grad_u[I]*grad_u[I];
                norm_grad_u = std::sqrt(norm_grad_u)+1E-10;

                for(int i=0;i<nDOF_test_element;i++)
                  {
                    element_weighted_lumped_mass_matrix[i] += norm_grad_u*u_test_dV[i];
                    element_rhsx_normal_reconstruction[i] += rhsx*u_test_dV[i];
                    element_rhsy_normal_reconstruction[i] += rhsy*u_test_dV[i];
                    element_rhsz_normal_reconstruction[i] += rhsz*u_test_dV[i];
		    for(int j=0;j<nDOF_trial_element;j++)
		      element_weighted_mass_matrix[i][j] +=
			norm_grad_u*u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i];
                  }
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index

                weighted_lumped_mass_matrix[gi] += element_weighted_lumped_mass_matrix[i];
		// rhs for reconstruction via consistent mass matrix
		rhs_qx[gi] += element_rhsx_normal_reconstruction[i];
		rhs_qy[gi] += element_rhsy_normal_reconstruction[i];
		rhs_qz[gi] += element_rhsz_normal_reconstruction[i];
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    weighted_mass_matrix[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]]
		      += element_weighted_mass_matrix[i][j];
		  }
              }//i
          }//elements
      }

      void calculateRhsL2Proj(//element
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
			      double he_for_disc_ICs,
			      double* u_dof,
			      int offset_u, int stride_u,
			      // PARAMETERS FOR EDGE VISCOSITY
			      int numDOFs,
			      double* rhs_l2_proj)
      {
        for (int i=0; i<numDOFs; i++)
	  rhs_l2_proj[i]=0.;
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double element_rhs_l2_proj[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
	      element_rhs_l2_proj[i]=0.0;

            //loop over quadrature points and compute integrands
            for(int k=0;k<nQuadraturePoints_element;k++)
              {
                //compute indeces and declare local storage
                register int eN_k = eN*nQuadraturePoints_element+k,
                  eN_k_nSpace = eN_k*nSpace,
                  eN_nDOF_trial_element = eN*nDOF_trial_element;
                register double
		  u,u_test_dV[nDOF_trial_element],
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
                ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      u);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

                for(int i=0;i<nDOF_test_element;i++)
		  element_rhs_l2_proj[i] += he_for_disc_ICs*u*u_test_dV[i];
		//element_rhs_l2_proj[i] += u*u_test_dV[i];
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		rhs_l2_proj[gi] += element_rhs_l2_proj[i];
              }//i
          }//elements
      }

      void calculateLumpedMassMatrix(//element
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
				     double* lumped_mass_matrix,
				     int offset_u, int stride_u)
      {
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double element_lumped_mass_matrix[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
	      element_lumped_mass_matrix[i]=0.0;
            //loop over quadrature points and compute integrands
            for(int k=0;k<nQuadraturePoints_element;k++)
              {
                //compute indeces and declare local storage
                register int eN_k = eN*nQuadraturePoints_element+k,
                  eN_k_nSpace = eN_k*nSpace,
                  eN_nDOF_trial_element = eN*nDOF_trial_element;
                register double
                  //for mass matrix contributions
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
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

                for(int i=0;i<nDOF_test_element;i++)
		  element_lumped_mass_matrix[i] += u_test_dV[i];
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
                lumped_mass_matrix[gi] += element_lumped_mass_matrix[i];
              }//i
          }//elements
      }

      void assembleSpinUpSystem(//element
				double* mesh_trial_ref,
				double* mesh_grad_trial_ref,
				double* mesh_dof,
				int* mesh_l2g,
				double* dV_ref,
				double* u_trial_ref,
				double* u_test_ref,
				//physics
				int nElements_global,
				int* u_l2g,
				double* uInitial,
				int offset_u, int stride_u,
				double* globalResidual,
				int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				double* globalMassMatrix)
      {
        for(int eN=0;eN<nElements_global;eN++)
          {
            register double
	      elementResidual_u[nDOF_test_element],
	      elementMassMatrix_u_u[nDOF_test_element][nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
	      {
		elementResidual_u[i]=0;
		for (int j=0;j<nDOF_trial_element;j++)
		  elementMassMatrix_u_u[i][j]=0.0;
	      }
            for  (int k=0;k<nQuadraturePoints_element;k++)
              {
                int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
                //declare local storage
                register double
                  jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                  u_test_dV[nDOF_test_element],
                  dV, x,y,z;
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
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref[k];
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

                //////////////////
                // LOOP ON DOFs //
                //////////////////
                for(int i=0;i<nDOF_test_element;i++)
                  {
		    elementResidual_u[i] += uInitial[eN_k]*u_test_dV[i];
                    for(int j=0;j<nDOF_trial_element;j++)
                      {
                        elementMassMatrix_u_u[i][j] +=
			  u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i];
                      }//j
                  }//i
              }//k
            //
            //load into element Jacobian into global Jacobian
            //
            for (int i=0;i<nDOF_test_element;i++)
              {
                int eN_i = eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index
		globalResidual[gi] += elementResidual_u[i];
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    int eN_i_j = eN_i*nDOF_trial_element+j;
                    globalMassMatrix[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
                      elementMassMatrix_u_u[i][j];
                  }//j
              }//i
          }//elements
      }

      void FCTStep(int NNZ, //number on non-zero entries on sparsity pattern
		   int numDOFs, //number of DOFs
		   double* lumped_mass_matrix, //lumped mass matrix (as vector)
		   double* soln,
		   double* solH, //DOFs of high order solution at tnp1
		   double* solL,
		   double* limited_solution,
		   int* csrRowIndeces_DofLoops, //csr row indeces
		   int* csrColumnOffsets_DofLoops, //csr column offsets
		   double* MassMatrix //mass matrix
		   )
      {
	Rpos.resize(numDOFs, 0.0);
	Rneg.resize(numDOFs, 0.0);
	FluxCorrectionMatrix.resize(NNZ, 0.0);
	//////////////////
	// LOOP in DOFs //
	//////////////////
	int ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    //read some vectors
	    double solHi = solH[i];
	    double solLi = solL[i];
	    double mi = lumped_mass_matrix[i];

	    double mini=-1.0, maxi=1.0; // global FCT
	    //double mini=1.0E10, maxi=-1.0E10;
	    double Pposi=0, Pnegi=0;
	    // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	    for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	      {
		int j = csrColumnOffsets_DofLoops[offset];
		// i-th row of flux correction matrix
		FluxCorrectionMatrix[ij] = ((i==j ? 1. : 0.)*mi - MassMatrix[ij]) * (solH[j]-solHi);

		//mini = fmin(mini,limited_solution[j]);
		//maxi = fmax(maxi,limited_solution[j]);

		///////////////////////
		// COMPUTE P VECTORS //
		///////////////////////
		Pposi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] > 0) ? 1. : 0.);
		Pnegi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] < 0) ? 1. : 0.);

		//update ij
		ij+=1;
	      }
	    ///////////////////////
	    // COMPUTE Q VECTORS //
	    ///////////////////////
	    double Qposi = mi*(maxi-solLi);
	    double Qnegi = mi*(mini-solLi);

	    ///////////////////////
	    // COMPUTE R VECTORS //
	    ///////////////////////
	    Rpos[i] = ((Pposi==0) ? 1. : std::min(1.0,Qposi/Pposi));
	    Rneg[i] = ((Pnegi==0) ? 1. : std::min(1.0,Qnegi/Pnegi));
	  } // i DOFs

	//////////////////////
	// COMPUTE LIMITERS //
	//////////////////////
	ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double ith_Limiter_times_FluxCorrectionMatrix = 0.;
	    double Rposi = Rpos[i], Rnegi = Rneg[i];
	    // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
	    for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	      {
		int j = csrColumnOffsets_DofLoops[offset];
		ith_Limiter_times_FluxCorrectionMatrix +=
		  ((FluxCorrectionMatrix[ij]>0) ? std::min(Rposi,Rneg[j]) : std::min(Rnegi,Rpos[j]))
		  * FluxCorrectionMatrix[ij];
		//ith_Limiter_times_FluxCorrectionMatrix += FluxCorrectionMatrix[ij];
		//update ij
		ij+=1;
	      }
	    limited_solution[i] = solL[i] + 1./lumped_mass_matrix[i]*ith_Limiter_times_FluxCorrectionMatrix;
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
