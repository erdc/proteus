#ifndef BlendedSpaces_H
#define BlendedSpaces_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"

#define ALPHA_VIA_LINEAR_SPACE 1
#define DO_CHECKS 0
#define uDOT_VIA_HIGH_ORDER_GALERKIN 1

namespace proteus
{
  inline double ENTROPY(const double& u)
  {
    return 1./2.*std::pow(u,2.);
  }
  inline double DENTROPY(const double& u)
  {
    return u;
  }
  class BlendedSpaces_base
  {
    //The base class defining the interface
  public:
    ///////////////////////
    // Auxiliary vectors //
    ///////////////////////
    std::valarray<double> lowOrderSolution;
    std::valarray<double> boundaryIntegralLowOrder;
    std::valarray<double> boundaryIntegral;
    std::valarray<double> umax;
    std::valarray<double> umin;
    // for entropy viscosity
    std::valarray<double> EntVisc;
    // for debugging
    std::valarray<double> highOrderAdvection;
    std::valarray<double> consistentTimeDerivative;
    std::valarray<double> fluxCorrection;
    std::valarray<double> wij, wji;
    virtual ~BlendedSpaces_base(){}
    virtual void calculateResidual(//element
                                   double dt,
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
                                   int* u_l2g,
                                   int* r_l2g,
                                   double* elementDiameter,
                                   double* u_dof,
                                   double* u_dof_old,
                                   double* velocity,
                                   double* q_u,
				   double* q_grad_u,
                                   int offset_u, int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isFluxBoundary_u,
                                   double* ebqe_bc_flux_u_ext,
                                   double* ebqe_phi,double epsFact,
                                   double* ebqe_u,
                                   double* ebqe_flux,
                                   // PARAMETERS FOR EDGE VISCOSITY
                                   int numDOFs,
                                   int NNZ,
                                   int* rowptr,
                                   int* colind,
                                   int* csrRowIndeces_CellLoops,
                                   int* csrColumnOffsets_CellLoops,
                                   int* csrColumnOffsets_eb_CellLoops,
				   // FOR BLENDING SPACES
				   double* force,
				   double* alpha_dof,
				   double* aux_test_ref,
				   double* aux_grad_test_ref,
				   double* aux_test_trace_ref,
				   double* dLow,
				   double epsilon,
				   // Type of problem to solve
				   int PROBLEM_TYPE,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   // FOR highOrderLim
				   int GET_POINT_VALUES,
				   double* flux_qij,
				   double* element_flux_qij,
				   double* element_MC,
				   double* vVector,
				   double* element_flux_i,
				   double* intBernMat,
				   double* edge_based_cfl,
				   // velocity at dofs
				   double* u_vel_dofs,
				   double* v_vel_dofs,
				   // lumped mass matrices
				   double* QH_ML,
				   // inverse of dissipative mass matrix
				   double* inv_element_ML_minus_MC,
				   // C-matrices
				   double* Cx,
				   double* Cy,
				   double* CTx,
				   double* CTy,
				   double* CTxElem,
				   double* CTyElem,
				   double* xGradRHS,
				   double* yGradRHS)=0;
    virtual void calculateResidualEntropyVisc(//element
                                   double dt,
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
                                   int* u_l2g,
                                   int* r_l2g,
                                   double* elementDiameter,
                                   double* u_dof,
                                   double* u_dof_old,
                                   double* velocity,
                                   double* q_u,
				   double* q_grad_u,
                                   int offset_u, int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isFluxBoundary_u,
                                   double* ebqe_bc_flux_u_ext,
                                   double* ebqe_phi,double epsFact,
                                   double* ebqe_u,
                                   double* ebqe_flux,
                                   // PARAMETERS FOR EDGE VISCOSITY
                                   int numDOFs,
                                   int NNZ,
                                   int* rowptr,
                                   int* colind,
                                   int* csrRowIndeces_CellLoops,
                                   int* csrColumnOffsets_CellLoops,
                                   int* csrColumnOffsets_eb_CellLoops,
				   // FOR BLENDING SPACES
				   double* force,
				   double* alpha_dof,
				   double* aux_test_ref,
				   double* aux_grad_test_ref,
				   double* aux_test_trace_ref,
				   double* dLow,
				   double epsilon,
				   // Type of problem to solve
				   int PROBLEM_TYPE,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   // FOR highOrderLim
				   int GET_POINT_VALUES,
				   double* flux_qij,
				   double* element_flux_qij,
				   double* element_MC,
				   double* vVector,
				   double* element_flux_i,
				   double* intBernMat,
				   double* edge_based_cfl,
				   // velocity at dofs
				   double* u_vel_dofs,
				   double* v_vel_dofs,
				   // lumped mass matrices
				   double* QH_ML,
				   // inverse of dissipative mass matrix
				   double* inv_element_ML_minus_MC,
				   // C-matrices
				   double* Cx,
				   double* Cy,
				   double* CTx,
				   double* CTy,
				   double* CTxElem,
				   double* CTyElem,
				   double* xGradRHS,
				   double* yGradRHS)=0;
    virtual void calculateRHSGradientReconstruction(//element
						    double dt,
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
						    int* u_l2g,
						    int* r_l2g,
						    double* elementDiameter,
						    double* u_dof,
						    double* u_dof_old,
						    double* velocity,
						    double* q_u,
						    double* q_grad_u,
						    int offset_u, int stride_u,
						    double* globalResidual,
						    int nExteriorElementBoundaries_global,
						    int* exteriorElementBoundariesArray,
						    int* elementBoundaryElementsArray,
						    int* elementBoundaryLocalElementBoundariesArray,
						    double* ebqe_velocity_ext,
						    int* isDOFBoundary_u,
						    double* ebqe_bc_u_ext,
						    int* isFluxBoundary_u,
						    double* ebqe_bc_flux_u_ext,
						    double* ebqe_phi,double epsFact,
						    double* ebqe_u,
						    double* ebqe_flux,
						    // PARAMETERS FOR EDGE VISCOSITY
						    int numDOFs,
						    int NNZ,
						    int* rowptr,
						    int* colind,
						    int* csrRowIndeces_CellLoops,
						    int* csrColumnOffsets_CellLoops,
						    int* csrColumnOffsets_eb_CellLoops,
						    // FOR BLENDING SPACES
						    double* force,
						    double* alpha_dof,
						    double* aux_test_ref,
						    double* aux_grad_test_ref,
						    double* aux_test_trace_ref,
						    double* dLow,
						    double epsilon,
						    // Type of problem to solve
						    int PROBLEM_TYPE,
						    // AUX QUANTITIES OF INTEREST
						    double* quantDOFs,
						    // FOR highOrderLim
						    int GET_POINT_VALUES,
						    double* flux_qij,
						    double* element_flux_qij,
						    double* element_MC,
						    double* vVector,
						    double* element_flux_i,
						    double* intBernMat,
						    double* edge_based_cfl,
						    // velocity at dofs
						    double* u_vel_dofs,
						    double* v_vel_dofs,
						    // lumped mass matrices
						    double* QH_ML,
						    // inverse of dissipative mass matrix
						    double* inv_element_ML_minus_MC,
						    // C-matrices
						    double* Cx,
						    double* Cy,
						    double* CTx,
						    double* CTy,
						    double* CTxElem,
						    double* CTyElem,
						    double* xGradRHS,
						    double* yGradRHS)=0;    
    virtual void calculateJacobian(//element
                                   double dt,
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
                                   int* u_l2g,
                                   int* r_l2g,
                                   double* elementDiameter,
                                   double* u_dof,
                                   double* velocity,
                                   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                                   double* globalJacobian,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isFluxBoundary_u,
                                   double* ebqe_bc_flux_u_ext,
                                   int* csrColumnOffsets_eb_u_u,
				   // FOR BLENDING SPACES
				   int numDOFs,
				   int* rowptr,
				   int* colind,
				   double* alpha_dof,
				   double* aux_test_ref,
				   double* aux_grad_test_ref,
				   double* dLow,
				   double epsilon,
				   // Type of problem to solve
				   int PROBLEM_TYPE)=0;
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
                                       int* u_l2g,
                                       double* elementDiameter,
				       double* meshSize,
                                       double* nodeDiametersArray,
                                       double epsFactHeaviside,
                                       double* q_uh,
				       double* q_grad_uh,
                                       double* u_exact,
				       double* gradx_u_exact,
				       double* gracy_u_exact,
                                       int offset_u, int stride_u,
                                       double* global_L2,
                                       double* global_H1,
                                       double* global_L2_Omega1,
                                       double* global_H1_Omega1,
                                       double* global_Omega1,
                                       double* global_L2_Omega2,
                                       double* global_H1_Omega2,
				       double* global_Omega2,
				       double* global_L2_sH,
				       double* global_L2_1msH)=0;
  };

  template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class BlendedSpaces : public BlendedSpaces_base
    {
    public:
      const int nDOF_test_X_trial_element;
      CompKernelType ck;
    BlendedSpaces():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
	ck()
	  {}

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

      inline double Heaviside(double u)
      {
	double H;
	if (u > 0.)
	  H=1.0;
	else if (u < 0.)
	  H=0.0;
	else
	  H=0.5;
	return H;
      }

      inline void Mult(const double mat[nSpace*nSpace],
		       const double vec[nSpace],
		       double *mat_times_vector)
      {
	for (int I=0; I<nSpace; I++)
	  {
	    mat_times_vector[I] = 0.;
	    for (int J=0; J<nSpace; J++)
	      mat_times_vector[I] += mat[I*nSpace+J] * vec[J];
	  }
      }

      inline void MULT(const double A[nSpace*nSpace],
		       const double B[nSpace*nSpace],
		       double *mat)
      {
	for (int I=0; I<nSpace; I++)
	  {
	    for (int J=0; J<nSpace; J++)
	      {
		mat[I*nSpace+J] = 0;
		for (int K=0; K<nSpace; K++)
		  {
		    mat[I*nSpace+J] += A[I*nSpace+K] * B[K*nSpace+J];
		  }
	      }
	  }
      }

      inline void calculateDTensor(double *DTensor)
      {
	double k1 = 1000.0;
	double k2 = 1.0;
	double theta = M_PI/6;
	double visc[nSpace*nSpace], RLeft[nSpace*nSpace], RRight[nSpace*nSpace], aux[nSpace*nSpace];

	visc[0] = k1;
	visc[1] = 0;
	visc[2] = 0;
	visc[3] = k2;

	RLeft[0] = cos(-theta);
	RLeft[1] = sin(-theta);
	RLeft[2] = -sin(-theta);
	RLeft[3] = cos(-theta);

	RRight[0] = cos(theta);
	RRight[1] = sin(theta);
	RRight[2] = -sin(theta);
	RRight[3] = cos(theta);

	MULT(visc,RRight,aux);
	MULT(RLeft,aux,DTensor);
      }

      void calculateResidual(//element
			     double dt,
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
			     int* u_l2g,
			     int* r_l2g,
			     double* elementDiameter,
			     double* u_dof,
			     double* u_dof_old,
			     double* velocity,
			     double* q_u,
			     double* q_grad_u,
			     int offset_u, int stride_u,
			     double* globalResidual,
			     int nExteriorElementBoundaries_global,
			     int* exteriorElementBoundariesArray,
			     int* elementBoundaryElementsArray,
			     int* elementBoundaryLocalElementBoundariesArray,
			     double* ebqe_velocity_ext,
			     int* isDOFBoundary_u,
			     double* ebqe_bc_u_ext,
			     int* isFluxBoundary_u,
			     double* ebqe_bc_flux_u_ext,
			     double* ebqe_phi,double epsFact,
			     double* ebqe_u,
			     double* ebqe_flux,
			     // PARAMETERS FOR EDGE VISCOSITY
			     int numDOFs,
			     int NNZ,
			     int* rowptr,
			     int* colind,
			     int* csrRowIndeces_CellLoops,
			     int* csrColumnOffsets_CellLoops,
			     int* csrColumnOffsets_eb_CellLoops,
			     // FOR BLENDING SPACES
			     double* force,
			     double* alpha_dof,
			     double* aux_test_ref,
			     double* aux_grad_test_ref,
			     double* aux_test_trace_ref,
			     double* dLow,
			     double epsilon,
			     // Type of problem to solve
			     int PROBLEM_TYPE,
			     // AUX QUANTITIES OF INTEREST			     
			     double* quantDOFs,
			     // For highOrderLim
			     int GET_POINT_VALUES,
			     double* flux_qij,
			     double* element_flux_qij,
			     double* element_MC,
			     double* vVector,
			     double* element_flux_i,
			     double* intBernMat,
			     double* edge_based_cfl,
			     // velocity at dofs
			     double* u_vel_dofs,
			     double* v_vel_dofs,
			     // lumped mass matrices
			     double* QH_ML,
			     // inverse of dissipative mass matrix
			     double* inv_element_ML_minus_MC,
			     // C-matrices
			     double* Cx,
			     double* Cy,
			     double* CTx,
			     double* CTy,
			     double* CTxElem,
			     double* CTyElem,
			     double* xGradRHS,
			     double* yGradRHS)
      {
	//
	//loop over elements to compute volume integrals and load them into element and global res.
	//
	//eN is the element index
	//eN_k is the quadrature point index for a scalar
	//eN_k_nSpace is the quadrature point index for a vector
	//eN_i is the element test function index
	//eN_j is the element trial function index
	//eN_k_j is the quadrature point index for a trial function
	//eN_k_i is the quadrature point index for a trial function

	////////////////////////////////
	// Zero out auxiliary vectors //
	////////////////////////////////
	lowOrderSolution.resize(numDOFs,0.0);
	boundaryIntegralLowOrder.resize(numDOFs,0.0);
	boundaryIntegral.resize(numDOFs,0.0);
	umax.resize(numDOFs,0.0);
	umin.resize(numDOFs,0.0);

	// for debugging
	highOrderAdvection.resize(numDOFs,0.0);
	fluxCorrection.resize(numDOFs,0.0);
	
	register double uDot[numDOFs];
	for (int i=0; i<numDOFs; i++)
	  uDot[i]=0.0;
	
	wij.resize(NNZ,0.0);
	wji.resize(NNZ,0.0);
	for (int i=0; i<NNZ; i++)
	  {
	    wij[i] = 0.0;
	    wji[i] = 0.0;
	  }

	///////////////////
	// BOUNDARY TERM //
	///////////////////
	// * Compute boundary integrals 
	for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
	  {
	    register int ebN = exteriorElementBoundariesArray[ebNE];
	    register int eN  = elementBoundaryElementsArray[ebN*2+0],
	      ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	      eN_nDOF_trial_element = eN*nDOF_trial_element;
	    register double elementBoundaryFluxLowOrder[nDOF_test_element],
	      elementBoundaryFluxHighOrder[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementBoundaryFluxLowOrder[i]=0.0;
		elementBoundaryFluxHighOrder[i]=0.0;
	      }
	    // loop on quad points
	    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
	      {
		register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		  ebNE_kb_nSpace = ebNE_kb*nSpace,
		  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
		register double
		  uExt=0.,
		  jac_ext[nSpace*nSpace],
		  jacDet_ext,
		  jacInv_ext[nSpace*nSpace],
		  boundaryJac[nSpace*(nSpace-1)],
		  metricTensor[(nSpace-1)*(nSpace-1)],
		  metricTensorDetSqrt,
		  dS,
		  u_test_dS[nDOF_test_element],
		  normal[nSpace],x_ext,y_ext,z_ext;
		// calculate mappings
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
		dS = ((1.0)*metricTensorDetSqrt )*dS_ref[kb];
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],
			      uExt);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		
		//calculate flow
		double flow = 0.;
		for (int I=0; I < nSpace; I++)
		  flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];

		for (int i=0;i<nDOF_test_element;i++)
		  {
		    int eN_i = eN*nDOF_test_element+i;
		    int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index
		    double uni = u_dof_old[gi];
		    double uiInlet = 0.0;
		    // low order part
		    double auxLowOrder = (flow >= 0 ? 0. : 1.) * (uni-uiInlet)*flow*u_test_dS[i];
		    elementBoundaryFluxLowOrder[i] += auxLowOrder;
		    // high order part
		    double auxHighOrder = (flow >= 0 ? 0. : 1.)*(uExt-uni)*flow*u_test_dS[i];
		    elementBoundaryFluxHighOrder[i] += auxHighOrder;
		  }
	      }//kb
	    // save to the corresponding vectors
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index
		boundaryIntegralLowOrder[gi] += elementBoundaryFluxLowOrder[i];		
		element_flux_i[eN_i] = elementBoundaryFluxHighOrder[i]; 
		fluxCorrection[gi] += elementBoundaryFluxHighOrder[i];
	      }
	  }//ebNE
	///////////////////////////
	// END OF BOUNDARY TERMS //
	///////////////////////////
	
	////////////////////////
	// FIRST LOOP IN DOFs //
	////////////////////////
	// * Compute the low order solution
	// * Compute (lagged) edge_based_cfl
	// * Compute the dissipative part of the global flux_qij
	// * Compute uDot (which we might use depending on the high-order stabilization)
	// * Compute umax and umin
	int ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	    double u_veli = u_vel_dofs[i];
	    double v_veli = v_vel_dofs[i];
	    double ith_dissipative_term = 0;
	    double dLowii = 0.;
	    //int ii=0;
	    double ith_flux_term = 0;

	    // for the computation of the local bounds
	    double umaxi = u_dof_old[i];
	    double umini = u_dof_old[i];

	    double fxi = u_veli*solni;
	    double fyi = v_veli*solni;

	    // loop over the sparsity pattern of the i-th DOF
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
		double u_velj = u_vel_dofs[j];
		double v_velj = v_vel_dofs[j];
		double dij = 0;

		double fxj = u_velj*solnj;
		double fyj = v_velj*solnj;
		
		ith_flux_term += Cx[ij]*(u_velj*solnj) + Cy[ij]*(v_velj*solnj);
		
		if (i != j)
		  {
		    dij = fmax(fmax(fabs(Cx[ij]*u_veli + Cy[ij]*v_veli),
				    fabs(Cx[ij]*u_velj + Cy[ij]*v_velj)),
			       fmax(fabs(CTx[ij]*u_veli + CTy[ij]*v_veli),
				    fabs(CTx[ij]*u_velj + CTy[ij]*v_velj)));
		    dLowii -= dij;
		    ith_dissipative_term += dij*(solnj-solni);
		    // compute anti-dissipative term of the flux_qij
		    flux_qij[ij] = -dij*(solnj-solni);

		    // computation of the local bounds
		    umaxi = fmax(solnj,umaxi);
		    umini = fmin(solnj,umini);

		    // save dij
		    dLow[ij] = dij;
		  }
		else
		  {
		    flux_qij[ij] = 0.0;
		    dLow[ij] = 0.0; // Not true but irrelevant due to (solnj-solni)
		  }
		// compute wij elements
		wij[ij] = (2*dij*(solni+solnj)/2.0 
			   -(Cx[ij]*(fxj-fxi) + Cy[ij]*(fyj-fyi)));
		// compute wji elements
		wji[ij] = (2*dij*(solni+solnj)/2.0 
			   -(CTx[ij]*(fxi-fxj) + CTy[ij]*(fyi-fyj)));

		// update index
		ij+=1;
	      }
	    // computation of the local bounds
	    umax[i] = umaxi;
	    umin[i] = umini;
	    
	    double QH_mi = QH_ML[i];
	    // compute edge based cfl
	    edge_based_cfl[i] = 2*fabs(dLowii)/QH_mi;
	    
	    // compute low order solution //
	    lowOrderSolution[i] = solni - dt/QH_mi * (ith_flux_term
						      - ith_dissipative_term
						      - boundaryIntegralLowOrder[i]);
	    // compute uDot //
	    // Note: we don't always use this. It depends on the form of the high-order stabilization
	    uDot[i] = -1.0/QH_mi * (ith_flux_term
				    - ith_dissipative_term
				    - boundaryIntegralLowOrder[i]);
	  }      
	///////////////////////////////
	// END OF FIRST LOOP IN DOFs //
	///////////////////////////////
 
	/////////////////////////////
	// FIRST LOOP ON INTEGRALS // Compute uDot via high order galerkin (with lumped m.mat)
	/////////////////////////////
	// * Compute the high order flux term.
	//    Used if the high-order stabilization is based on uDot=-1/mi*high_order_flux
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double elementResidual_u[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementResidual_u[i]=0.0;
	      }//i
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  un=0.0,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z;
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
		//get the trial function gradients based on the blended functions
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial); // high-order space
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un); // from high-order space
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		
		// coefficient of steady advection - diffusion
		double velocityBeta[2];
		velocityBeta[0] = velocity[eN_k_nSpace+0];
		velocityBeta[1] = velocity[eN_k_nSpace+1];
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    elementResidual_u[i] += un*ck.NumericalDiffusion(1.0,
								     velocityBeta,
								     &u_grad_test_dV[i_nSpace]);
		  }//i
	      }
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		highOrderAdvection[offset_u+stride_u*r_l2g[eN_i]] += elementResidual_u[i];
	      }//i
	  }//elements
	////////////////////////////////////
	// END OF FIRST LOOP IN INTEGRALS //
	////////////////////////////////////
	if (uDOT_VIA_HIGH_ORDER_GALERKIN==1)
	  {	
	    for (int i=0; i<numDOFs; i++)
	      uDot[i] = 1.0/QH_ML[i] * (highOrderAdvection[i]
					+ boundaryIntegralLowOrder[i]);
	  }
	
	//////////////////////////////
	// SECOND LOOP ON INTEGRALS //
	//////////////////////////////
	// * Compute element based part of element_flux_i
	// * Compute element inegrals of fluxCorrection as a global vector (only for debugging)
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double elementFlux[nDOF_test_element], elementMass[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementFlux[i]=0.0;
		elementMass[i]=0.0;
	      }//i

	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  un=0.0,uhDot=0.0, 
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_test_dV[nDOF_trial_element],
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z;
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
		//get the trial function gradients based on the blended functions
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial); // high-order space
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un); // from high-order space
		ck.valFromDOF(uDot,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      uhDot); // from high-order space
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }

		//save solution
		//q_u[eN_k] = u;
		// save grad
		//for(int I=0;I<nSpace;I++)
		//q_grad_u[eN_k_nSpace+I]=grad_u[I];

		// coefficient of steady advection - diffusion
		double velocityBeta[2];
		velocityBeta[0] = velocity[eN_k_nSpace+0];
		velocityBeta[1] = velocity[eN_k_nSpace+1];
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    elementFlux[i] += -uhDot*u_test_dV[i]
		      + un*ck.NumericalDiffusion(1.0,velocityBeta,&u_grad_test_dV[i_nSpace]);
		    elementMass[i] += u_test_dV[i];
		  }//i
	      } //quad points
	    //
	    //load elements into global vectors
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i=eN*nDOF_test_element+i;
	      int gi = offset_u+stride_u*r_l2g[eN_i];

	      // flux correction (only for debugging)
	      fluxCorrection[gi] += elementFlux[i];

	      double solni = u_dof_old[gi];
	      double u_veli = u_vel_dofs[gi];
	      double v_veli = v_vel_dofs[gi];
	      
	      double qi = elementFlux[i] + elementMass[i]*uDot[gi];
	      
	      // for flux of low-order method
	      for(int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  int eN_j = eN*nDOF_test_element+j;
		  int gj = offset_u+stride_u*r_l2g[eN_j];		  
		  double solnj = u_dof_old[gj];
		  double u_velj = u_vel_dofs[gj];
		  double v_velj = v_vel_dofs[gj];
		  // low-order flux part of qi
		  qi -= (CTxElem[eN_i_j]*(u_velj*solnj) + CTyElem[eN_i_j]*(v_velj*solnj));
		}
	      element_flux_i[eN_i] += qi;
	    }
	  }//elements
	/////////////////////////////////////
	// END OF SECOND LOOP IN INTEGRALS //
	/////////////////////////////////////
	
	//////////////////////////////
	// compute element vector v //
	//////////////////////////////
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for(int j=0;j<nDOF_test_element;j++)
		  {
		    int eN_j = eN*nDOF_test_element+j;
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    vVector[eN_i] += inv_element_ML_minus_MC[eN_i_j]*element_flux_i[eN_j];
		  }
	      }
	    // substract mean
	    double mean = 0;
	    for(int i=0;i<nDOF_test_element;i++)
	      {
	    	int eN_i = eN*nDOF_test_element+i;
	    	mean += vVector[eN_i]/nDOF_test_element;
	      }
	    for(int i=0;i<nDOF_test_element;i++)
	      {
	    	int eN_i = eN*nDOF_test_element+i;
	    	vVector[eN_i] -= mean;
	      }
	  }
	////////////////////////////////////////
	// end of computation of the vector v //
	////////////////////////////////////////

	//////////////////////////////////////////////////////////
	// compute the element_flux_qij and the global flux_qij //
	//////////////////////////////////////////////////////////
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for(int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_j = eN*nDOF_test_element+j;
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    
		    element_flux_qij[eN_i_j] = element_MC[eN_i_j]*(vVector[eN_i]-vVector[eN_j]);
		    
		    int jacIndex = csrRowIndeces_CellLoops[eN_i]+csrColumnOffsets_CellLoops[eN_i_j];
		    flux_qij[jacIndex] += element_flux_qij[eN_i_j];
		  }
	      }
	  }
	///////////////////////////////////////////////////////
	// End of computation of element and global flux_qij //
	///////////////////////////////////////////////////////

	if (DO_CHECKS==1)
	  {
	    /////////////////////////////////////////
	    // ********** for debugging ********** //
	    /////////////////////////////////////////
	    // * Check that element_flux_i is massless; i.e., that sum_i(element_flux_i)=0
	    // * Check that sum_j(element_flux_qij) = element_flux_i
	    // * Check that sum_j(flux_qij) = fluxCorrection[i]
	    // * Check that sum_i(sum_j(flux_qij)) = 0
	    for(int eN=0;eN<nElements_global;eN++) //loop in cells
	      {
		double sumi_element_flux_i = 0;
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    double sumj_element_flux_qij=0;
		    int eN_i = eN*nDOF_test_element+i;
		    
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int eN_j = eN*nDOF_test_element+j;
			int eN_i_j = eN_i*nDOF_trial_element+j;
			sumj_element_flux_qij += element_flux_qij[eN_i_j];
		      }
		    // Check that sum_j(element_flux_qij) = element_flux_i
		    if (fabs(element_flux_i[eN_i]-sumj_element_flux_qij) > 1E-13)
		      {
			std::cout << "|element_flux_i - sum_j(element_flux_qij)| > 1E-13 ... "
				  << "\t" 
				  << fabs(element_flux_i[eN_i]-sumj_element_flux_qij)
				  << std::endl;
			abort();
		      }
		    sumi_element_flux_i += sumj_element_flux_qij;
		  }
		// Check that element_flux_i is massless; i.e., that sum_i(element_flux_i)=0
		if (fabs(sumi_element_flux_i)>1E-13)
		  {
		    std::cout << "|sum_i(element_flux_i)| > 1E-13 ... "
			      << "\t"
			      << fabs(sumi_element_flux_i)
			      << std::endl;
		    abort();
		  }
	      }
	    double sumi_sumj_flux_qij = 0;
	    ij=0;
	    for (int i=0; i<numDOFs; i++)
	      {
		double solni = u_dof_old[i]; // solution at time tn for the ith DOF
		double u_veli = u_vel_dofs[i];
		double v_veli = v_vel_dofs[i];
		
		double sumj_flux_qij=0;
		double ith_advection_fluxCorrection = 0.;
		double ith_dissipative_fluxCorrection = 0.;
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  {
		    int j = colind[offset];
		    double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
		    double u_velj = u_vel_dofs[j];
		    double v_velj = v_vel_dofs[j];
		    double dij = 0;
		    
		    sumj_flux_qij += flux_qij[ij];
		    
		    // compute ith advective and dissipative terms (to add to fluxCorrection[i])
		    ith_advection_fluxCorrection += CTx[ij]*(u_velj*solnj) + CTy[ij]*(v_velj*solnj);
		    dij = fmax(fmax(fabs(Cx[ij]*u_veli + Cy[ij]*v_veli),
				    fabs(Cx[ij]*u_velj + Cy[ij]*v_velj)),
			       fmax(fabs(CTx[ij]*u_veli + CTy[ij]*v_veli),
				    fabs(CTx[ij]*u_velj + CTy[ij]*v_velj)));
		    ith_dissipative_fluxCorrection += dij*(solnj-solni);
		    
		    // update index
		    ij+=1;
		  }
		double mi = QH_ML[i];
		fluxCorrection[i] += (mi*uDot[i]
				      - ith_advection_fluxCorrection
				      - ith_dissipative_fluxCorrection);
		// * Check that sum_j(flux_qij) = fluxCorrection[i]
		if (fabs(sumj_flux_qij - fluxCorrection[i])>1E-10)
		  {
		    std::cout << "|sumj_flux_qij - fluxCorrection[i]| > 1E-10 ..."
			      << "\t" 
			      << fabs(sumj_flux_qij - fluxCorrection[i])
			      << std::endl;
		    abort();
		  }
		sumi_sumj_flux_qij += sumj_flux_qij;
	      }
	    // * Check that sum_i(sum_j(flux_qij)) = 0
	    if (fabs(sumi_sumj_flux_qij)>1E-13)
	      {
		std::cout << "|sum_i(sum_j(flux_qij))| > 1E-13 ... "
			  << "\t"
			  << fabs(sumi_sumj_flux_qij)
			  << std::endl;
		abort();
	      }
	    //////////////////////////////////////////
	    // *** END OF Section for debugging *** //
	    //////////////////////////////////////////
	  }
	
	///////////////////////
	// LAST LOOP IN DOFs //
	///////////////////////
	// * Compute the final solution and get the solution at DOFs (for visualization)
	ij = 0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double ith_galerkin_fluxCorrection = 0.;
	    double ith_limited_fluxCorrection = 0.;
	    double ith_flux_term = 0.;
	    double solni=u_dof_old[i];
	    
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		double fij = flux_qij[ij];
		double dij = dLow[ij];
		
		double fStarij = 0.0;
		if (i!=j)
		  {
		    if (fij > 0)
		      fStarij = fmin(fij,fmin(2*dij*umax[i] - wij[ij], wji[ij] - 2*dij*umin[j]));
		    else
		      fStarij = fmax(fij,fmax(2*dij*umin[i] - wij[ij], wji[ij] - 2*dij*umax[j]));
		  }

		// compute Galerkin and limited flux
		ith_galerkin_fluxCorrection += flux_qij[ij];
		ith_limited_fluxCorrection += fStarij;

		// compute high-order solution with entropy viscosity //
		double solnj = u_dof_old[j];
		double u_velj = u_vel_dofs[j];
		double v_velj = v_vel_dofs[j];
		ith_flux_term += Cx[ij]*(u_velj*solnj) + Cy[ij]*(v_velj*solnj);
		
		ij+=1;
	      }
	    double mi = QH_ML[i];

	    // COMPUTE SOLUTION //
	    globalResidual[i] = lowOrderSolution[i] + dt/mi*fluxCorrection[i]; // for debugging
	    //globalResidual[i] = lowOrderSolution[i];

	    
	    //globalResidual[i] = lowOrderSolution[i] + dt/mi*ith_galerkin_fluxCorrection;
	    //globalResidual[i] = lowOrderSolution[i] + dt/mi*ith_limited_fluxCorrection;
	  }
	//////////////////////////////
	// END OF LAST LOOP IN DOFs //
	//////////////////////////////

	///////////////////////////////////////////////
	// CONVERT BERNSTEIN DOFs TO SOLUTION VALUES //
	///////////////////////////////////////////////
	if (GET_POINT_VALUES==1)
	  {
	    ij=0;
	    for (int i=0; i<numDOFs; i++)
	      {
		quantDOFs[i] = 0;
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  {
		    int j = colind[offset];
		    quantDOFs[i] += intBernMat[ij]*globalResidual[j];
		    ij+=1;
		  }
		//std::cout << quantDOFs[i]
		//	  << "\t"
		//	  << globalResidual[i]
		//	  << std::endl;
	      }
	  }
      }

      void calculateResidualEntropyVisc(//element
					double dt,
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
					int* u_l2g,
					int* r_l2g,
					double* elementDiameter,
					double* u_dof,
					double* u_dof_old,
					double* velocity,
					double* q_u,
					double* q_grad_u,
					int offset_u, int stride_u,
					double* globalResidual,
					int nExteriorElementBoundaries_global,
					int* exteriorElementBoundariesArray,
					int* elementBoundaryElementsArray,
					int* elementBoundaryLocalElementBoundariesArray,
					double* ebqe_velocity_ext,
					int* isDOFBoundary_u,
					double* ebqe_bc_u_ext,
					int* isFluxBoundary_u,
					double* ebqe_bc_flux_u_ext,
					double* ebqe_phi,double epsFact,
					double* ebqe_u,
					double* ebqe_flux,
					// PARAMETERS FOR EDGE VISCOSITY
					int numDOFs,
					int NNZ,
					int* rowptr,
					int* colind,
					int* csrRowIndeces_CellLoops,
					int* csrColumnOffsets_CellLoops,
					int* csrColumnOffsets_eb_CellLoops,
					// FOR BLENDING SPACES
					double* force,
					double* alpha_dof,
					double* aux_test_ref,
					double* aux_grad_test_ref,
					double* aux_test_trace_ref,
					double* dLow,
					double epsilon,
					// Type of problem to solve
					int PROBLEM_TYPE,
					// AUX QUANTITIES OF INTEREST			     
					double* quantDOFs,
					// For highOrderLim
					int GET_POINT_VALUES,
					double* flux_qij,
					double* element_flux_qij,
					double* element_MC,
					double* vVector,
					double* element_flux_i,
					double* intBernMat,
					double* edge_based_cfl,
					// velocity at dofs
					double* u_vel_dofs,
					double* v_vel_dofs,
					// lumped mass matrices
					double* QH_ML,
					// inverse of dissipative mass matrix
					double* inv_element_ML_minus_MC,
					// C-matrices
					double* Cx,
					double* Cy,
					double* CTx,
					double* CTy,
					double* CTxElem,
					double* CTyElem,
					double* xGradRHS,
					double* yGradRHS)
      {
	////////////////////////////////
	// Zero out auxiliary vectors //
	////////////////////////////////
	boundaryIntegral.resize(numDOFs,0.0);

	// for entropy viscosity
	EntVisc.resize(numDOFs,0.0);
	highOrderAdvection.resize(numDOFs,0.0);
	consistentTimeDerivative.resize(numDOFs,0.0);
	    
	///////////////////
	// BOUNDARY TERM //
	///////////////////
	// * Compute boundary integrals 
	for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
	  {
	    register int ebN = exteriorElementBoundariesArray[ebNE];
	    register int eN  = elementBoundaryElementsArray[ebN*2+0],
	      ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	      eN_nDOF_trial_element = eN*nDOF_trial_element;
	    register double elementBoundaryFluxHighOrder[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementBoundaryFluxHighOrder[i]=0.0;
	    // loop on quad points
	    for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
	      {
		register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		  ebNE_kb_nSpace = ebNE_kb*nSpace,
		  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
		register double
		  uExt=0.,
		  jac_ext[nSpace*nSpace],
		  jacDet_ext,
		  jacInv_ext[nSpace*nSpace],
		  boundaryJac[nSpace*(nSpace-1)],
		  metricTensor[(nSpace-1)*(nSpace-1)],
		  metricTensorDetSqrt,
		  dS,
		  u_test_dS[nDOF_test_element],
		  normal[nSpace],x_ext,y_ext,z_ext;
		// calculate mappings
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
		dS = ((1.0)*metricTensorDetSqrt )*dS_ref[kb];
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element],
			      uExt);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		
		//calculate flow
		double flow = 0.;
		for (int I=0; I < nSpace; I++)
		  flow += normal[I]*ebqe_velocity_ext[ebNE_kb_nSpace+I];

		for (int i=0;i<nDOF_test_element;i++)
		  {
		    int eN_i = eN*nDOF_test_element+i;
		    // high order flux
		    double auxHighOrder = (flow >= 0 ? 1. : 0.)*uExt*flow*u_test_dS[i];
		    elementBoundaryFluxHighOrder[i] += auxHighOrder;
		  }
	      }//kb
	    // save to the corresponding vectors
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index
		boundaryIntegral[gi] += elementBoundaryFluxHighOrder[i];
	      }
	  }//ebNE
	///////////////////////////
	// END OF BOUNDARY TERMS //
	///////////////////////////
	
	////////////////////////
	// FIRST LOOP IN DOFs //
	////////////////////////
	// * Compute (lagged) edge_based_cfl
	// * Compute EntVisc
	// * Compute dLow[ij]
	int ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	    double u_veli = u_vel_dofs[i];
	    double v_veli = v_vel_dofs[i];
	    double dLowii = 0.;

	    // for entropy viscosity
	    double DEnti = DENTROPY(solni);
	    double DenEntViscPart1 = 0.;
	    double DenEntViscPart2 = 0.;
	    double ith_NumEntVisc = 0.;
	    double ith_DenEntVisc = 0.;
	    
	    // loop over the sparsity pattern of the i-th DOF
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
		double u_velj = u_vel_dofs[j];
		double v_velj = v_vel_dofs[j];
		double dij = 0;
		
		if (i != j)
		  {
		    dij = fmax(fmax(fabs(Cx[ij]*u_veli + Cy[ij]*v_veli),
				    fabs(Cx[ij]*u_velj + Cy[ij]*v_velj)),
			       fmax(fabs(CTx[ij]*u_veli + CTy[ij]*v_veli),
				    fabs(CTx[ij]*u_velj + CTy[ij]*v_velj)));
		    dLowii -= dij;
		    // save dij
		    dLow[ij] = dij;
		  }
		else
		  {
		    dLow[ij] = 0.0; // Not true but irrelevant due to use of ...(solnj-solni)
		  }		
		// For entropy viscosity production
		double fxj = u_velj*solnj;
		double fyj = v_velj*solnj;
		
		double x_EntFluxj = u_velj*ENTROPY(solnj) - u_veli*ENTROPY(solni);
		double y_EntFluxj = v_velj*ENTROPY(solnj) - v_veli*ENTROPY(solni);
		ith_NumEntVisc += Cx[ij]*(x_EntFluxj-DEnti*fxj) + Cy[ij]*(y_EntFluxj-DEnti*fyj);

		// aux parts to compute DenEntVisc
		DenEntViscPart1 += Cx[ij]*x_EntFluxj + Cy[ij]*y_EntFluxj;
		DenEntViscPart2 += Cx[ij]*fxj + Cy[ij]*fyj;

		//update index 
		ij+=1;
	      }
	    // compute DenEntVisc
	    ith_DenEntVisc = fabs(DenEntViscPart1) + fabs(DEnti)*fabs(DenEntViscPart2)+1E-15;
	    EntVisc[i] = fabs(ith_NumEntVisc)/ith_DenEntVisc;
	    
	    // compute edge based cfl
	    double QH_mi = QH_ML[i];
	    edge_based_cfl[i] = 2*fabs(dLowii)/QH_mi;
	  }
	///////////////////////////////
	// END OF FIRST LOOP IN DOFs //
	///////////////////////////////
 
	/////////////////////////////
	// FIRST LOOP ON INTEGRALS // Compute uDot via high order galerkin (with lumped m.mat)
	/////////////////////////////
	// * Compute the high order flux term.
	//    Used if the high-order stabilization is based on uDot=-1/mi*high_order_flux
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double elementTimeDerivative[nDOF_test_element],
	      elementFluxTerm[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementTimeDerivative[i]=0.0;
		elementFluxTerm[i]=0.0;
	      }
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  unp1=0.0,
		  un=0.0,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_test_dV[nDOF_test_element],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z;
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
		//get the trial function gradients based on the blended functions
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial); // high-order space
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un); // from high-order space
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      unp1); // from high-order space
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		
		// coefficient of steady advection - diffusion
		double velocityBeta[2];
		velocityBeta[0] = velocity[eN_k_nSpace+0];
		velocityBeta[1] = velocity[eN_k_nSpace+1];
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    elementTimeDerivative[i] += (unp1-un)*u_test_dV[i];
		    
		    register int i_nSpace=i*nSpace;
		    elementFluxTerm[i] += un*ck.NumericalDiffusion(1.0,
								    velocityBeta,
								    &u_grad_test_dV[i_nSpace]);
		  }//i
	      }
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		
		consistentTimeDerivative[gi] += elementTimeDerivative[i];
		highOrderAdvection[gi] += elementFluxTerm[i];
	      }//i
	  }//elements
	////////////////////////////////////
	// END OF FIRST LOOP IN INTEGRALS //
	////////////////////////////////////
 
	///////////////////////
	// LAST LOOP IN DOFs //
	///////////////////////
	// * Compute the final solution and get the solution at DOFs (for visualization)
	ij = 0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double ith_flux_term = 0.;
	    double ith_EntVisc_dissipative_term = 0.;
	    double solni=u_dof_old[i];
	    
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		double fij = flux_qij[ij];
		double dij = dLow[ij];
		
		// compute high-order solution with entropy viscosity //
		double solnj = u_dof_old[j];
		double u_velj = u_vel_dofs[j];
		double v_velj = v_vel_dofs[j];

		ith_flux_term += Cx[ij]*(u_velj*solnj) + Cy[ij]*(v_velj*solnj);
		//ith_EntVisc_dissipative_term += dLow[ij]*(solnj-solni);
		ith_EntVisc_dissipative_term += dLow[ij]*fmax(EntVisc[i],EntVisc[j])*(solnj-solni);
		
		ij+=1;
	      }
	    double mi = QH_ML[i];

	    //globalResidual[i] = solni - dt/mi*(ith_flux_term
	    //				       -ith_EntVisc_dissipative_term);
	    globalResidual[i] = solni - dt/mi*(-highOrderAdvection[i]
	    				       + boundaryIntegral[i]
	    				       - ith_EntVisc_dissipative_term);

	    //globalResidual[i] = consistentTimeDerivative[i] + dt*(-highOrderAdvection[i]
	    //							  + boundaryIntegral[i]
	    //							  - ith_EntVisc_dissipative_term);
	  }
	//////////////////////////////
	// END OF LAST LOOP IN DOFs //
	//////////////////////////////
      }
	    
      void calculateRHSGradientReconstruction(//element
					      double dt,
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
					      int* u_l2g,
					      int* r_l2g,
					      double* elementDiameter,
					      double* u_dof,
					      double* u_dof_old,
					      double* velocity,
					      double* q_u,
					      double* q_grad_u,
					      int offset_u, int stride_u,
					      double* globalResidual,
					      int nExteriorElementBoundaries_global,
					      int* exteriorElementBoundariesArray,
					      int* elementBoundaryElementsArray,
					      int* elementBoundaryLocalElementBoundariesArray,
					      double* ebqe_velocity_ext,
					      int* isDOFBoundary_u,
					      double* ebqe_bc_u_ext,
					      int* isFluxBoundary_u,
					      double* ebqe_bc_flux_u_ext,
					      double* ebqe_phi,double epsFact,
					      double* ebqe_u,
					      double* ebqe_flux,
					      // PARAMETERS FOR EDGE VISCOSITY
					      int numDOFs,
					      int NNZ,
					      int* rowptr,
					      int* colind,
					      int* csrRowIndeces_CellLoops,
					      int* csrColumnOffsets_CellLoops,
					      int* csrColumnOffsets_eb_CellLoops,
					      // FOR BLENDING SPACES
					      double* force,
					      double* alpha_dof,
					      double* aux_test_ref,
					      double* aux_grad_test_ref,
					      double* aux_test_trace_ref,
					      double* dLow,
					      double epsilon,
					      // Type of problem to solve
					      int PROBLEM_TYPE,
					      // AUX QUANTITIES OF INTEREST
					      double* quantDOFs,
					      // For highOrderLim
					      int GET_POINT_VALUES,
					      double* flux_qij,
					      double* element_flux_qij,
					      double* element_MC,
					      double* vVector,
					      double* element_flux_i,
					      double* intBernMat,
					      double* edge_based_cfl,
					      // velocity at dofs
					      double* u_vel_dofs,
					      double* v_vel_dofs,
					      // lumped mass matrices
					      double* QH_ML,
					      // inverse of dissipative mass matrix
					      double* inv_element_ML_minus_MC,
					      // C-matrices
					      double* Cx,
					      double* Cy,
					      double* CTx,
					      double* CTy,
					      double* CTxElem,
					      double* CTyElem,
					      // gradient reconstruction
					      double* xGradRHS,
					      double* yGradRHS)
      {
	///////////////////////
	// LOOP ON INTEGRALS //
	///////////////////////
	// * Compute grad(un)
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double
	      element_xGradRHS[nDOF_test_element],
	      element_yGradRHS[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		element_xGradRHS[i]=0.0;
		element_yGradRHS[i]=0.0;
	      }
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  grad_un[nSpace],
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_test_dV[nDOF_test_element],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  dV,x,y,z;
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
		//get the trial function gradients based on the blended functions
		ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial); // high-order space
		//get the solution based on the blended functions
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_un);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    element_xGradRHS[i] += grad_un[0]*u_test_dV[i];
		    element_yGradRHS[i] += grad_un[1]*u_test_dV[i];
		  }//i
	      }
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		xGradRHS[gi] += element_xGradRHS[i];
		yGradRHS[gi] += element_yGradRHS[i];
	      }//i
	  }//elements
	////////////////////////////////////
	// END OF FIRST LOOP IN INTEGRALS //
	////////////////////////////////////
      }
	    
      
      void calculateJacobian(//element
			     double dt,
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
			     int* u_l2g,
			     int* r_l2g,
			     double* elementDiameter,
			     double* u_dof,
			     double* velocity,
			     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			     double* globalJacobian,
			     int nExteriorElementBoundaries_global,
			     int* exteriorElementBoundariesArray,
			     int* elementBoundaryElementsArray,
			     int* elementBoundaryLocalElementBoundariesArray,
			     double* ebqe_velocity_ext,
			     int* isDOFBoundary_u,
			     double* ebqe_bc_u_ext,
			     int* isFluxBoundary_u,
			     double* ebqe_bc_flux_u_ext,
			     int* csrColumnOffsets_eb_u_u,
			     // FOR BLENDING SPACES
			     int numDOFs,
			     int* rowptr,
			     int* colind,
			     double* alpha_dof,
			     double* aux_test_ref,
			     double* aux_grad_test_ref,
			     double* dLow,
			     double epsilon,
			     // Type of problem to solve
			     int PROBLEM_TYPE)
      {
	//
	//loop over elements
	//
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  elementJacobian_u_u[i][j]=0.0;
		}
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

		//declare local storage
		register double
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  dV,
		  u_test_dV[nDOF_test_element],
		  x,y,z;
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

		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			elementJacobian_u_u[i][j] +=
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
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] +=
		      elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements

	// NO LOOP IN BOUNDARIES //
      }//computeJacobian

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
                                 int* u_l2g,
                                 double* elementDiameter,
				 double* meshSize,
                                 double* nodeDiametersArray,
                                 double epsFactHeaviside,
                                 double* q_uh,
				 double* q_grad_uh,
                                 double* u_exact,
				 double* gradx_u_exact,
				 double* grady_u_exact,
                                 int offset_u, int stride_u,
				 double* global_L2,
				 double* global_H1,
                                 double* global_L2_Omega1,
				 double* global_H1_Omega1,
                                 double* global_Omega1,
                                 double* global_L2_Omega2,
				 double* global_H1_Omega2,
                                 double* global_Omega2,
				 double* global_L2_sH,
				 double* global_L2_1msH)
      {
	*global_L2 = 0.0;
	*global_H1 = 0.0;
        *global_L2_Omega1 = 0.0;
	*global_H1_Omega1 = 0.0;
	*global_Omega1 = 0.0;
        *global_L2_Omega2 = 0.0;
	*global_H1_Omega2 = 0.0;
	*global_Omega2 = 0.0;
	*global_L2_sH = 0.0;
	*global_L2_1msH = 0.0;
        //////////////////////
        // ** LOOP IN CELLS //
        //////////////////////
        for(int eN=0;eN<nElements_global;eN++)
          {
            if (eN<nElements_owned) // just consider the locally owned cells
              {
                //declare local storage for local contributions and initialize
                double cell_L2 = 0., cell_H1 = 0.,
		  cell_L2_Omega1 = 0., cell_H1_Omega1 = 0., cell_Omega1 = 0.,
		  cell_L2_Omega2 = 0., cell_H1_Omega2 = 0., cell_Omega2 = 0.,
		  cell_L2_sH = 0., cell_L2_1msH = 0.;

                //loop over quadrature points and compute integrands
                for  (int k=0;k<nQuadraturePoints_element;k++)
                  {
                    //compute indeces and declare local storage
                    register int eN_k = eN*nQuadraturePoints_element+k,
                      eN_k_nSpace = eN_k*nSpace,
                      eN_nDOF_trial_element = eN*nDOF_trial_element;
                    register double
                      u, gradx_u, grady_u,
		      u_grad_trial[nDOF_trial_element*nSpace],
		      uh, grad_uh[nSpace],
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
		    // load exact solution and its gradient
		    u = u_exact[eN_k];
		    gradx_u = gradx_u_exact[eN_k];
		    grady_u = grady_u_exact[eN_k];
                    // get numerical solution
		    uh = q_uh[eN_k];
                    // get gradients of numerical solution
		    grad_uh[0] = q_grad_uh[eN_k_nSpace+0];
		    grad_uh[1] = q_grad_uh[eN_k_nSpace+1];
		    // Norms in whole domain //
		    cell_L2 += std::pow(u-uh,2)*dV;
		    cell_H1 += (std::pow(u-uh,2) +
		    		(std::pow(grad_uh[0] - gradx_u,2.0) +
		    		 std::pow(grad_uh[1] - grady_u,2.0)) )*dV;

		    // Norms in Omega 1 //
		    if (x < 0.5 - meshSize[0])
		      //if (x < 0.5 - 0.25)
		      {
			cell_L2_Omega1 += std::pow(u-uh,2)*dV;
			cell_H1_Omega1 += (std::pow(u-uh,2) +
					   (std::pow(grad_uh[0] - gradx_u,2) +
					    std::pow(grad_uh[1] - grady_u,2)) )*dV;
			cell_Omega1 += dV;
		      }
		    if (x > 0.5 + meshSize[0])
		      //if (x > 0.5 + 0.25)
		      {
			cell_L2_Omega2 += std::pow(u-uh,2)*dV;
			cell_H1_Omega2 += (std::pow(u-uh,2) +
					   (std::pow(grad_uh[0] - gradx_u,2) +
					    std::pow(grad_uh[1] - grady_u,2)) )*dV;
			cell_Omega2 += dV;
		      }
		    //double epsHeaviside = epsFactHeaviside*elementDiameter[eN]/2.0;
		    double epsHeaviside = 0.1;
		    double sH = smoothedHeaviside(epsHeaviside,x-0.75);
		    //double epsHeaviside = meshSize[0];
		    //double sH = smoothedHeaviside(epsHeaviside,x-0.5-4*epsHeaviside);
		    //double sH = Heaviside(x-0.75);
		    cell_L2_sH += std::pow((u-uh)*sH,2)*dV;
		    cell_L2_1msH += std::pow(sH,2)*dV;
		    //cell_L2_1msH += std::pow(u-uh,2)*(1.0-sH)*dV;
                  }
		*global_L2 += cell_L2;
		*global_H1 += cell_H1;
		*global_L2_Omega1 += cell_L2_Omega1;
		*global_H1_Omega1 += cell_H1_Omega1;
		*global_Omega1 += cell_Omega1;
		*global_L2_Omega2 += cell_L2_Omega2;
		*global_H1_Omega2 += cell_H1_Omega2;
		*global_Omega2 += cell_Omega2;
		*global_L2_sH += cell_L2_sH;
		*global_L2_1msH += cell_L2_1msH;
              }//elements
          }
      }
    };//BlendedSpaces

  inline BlendedSpaces_base* newBlendedSpaces(int nSpaceIn,
					      int nQuadraturePoints_elementIn,
					      int nDOF_mesh_trial_elementIn,
					      int nDOF_trial_elementIn,
					      int nDOF_test_elementIn,
					      int nQuadraturePoints_elementBoundaryIn,
					      int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<BlendedSpaces_base,BlendedSpaces,CompKernel>(nSpaceIn,
												     nQuadraturePoints_elementIn,
												     nDOF_mesh_trial_elementIn,
												     nDOF_trial_elementIn,
												     nDOF_test_elementIn,
												     nQuadraturePoints_elementBoundaryIn,
												     CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<BlendedSpaces_base,BlendedSpaces,CompKernel>(nSpaceIn,
												   nQuadraturePoints_elementIn,
												   nDOF_mesh_trial_elementIn,
												   nDOF_trial_elementIn,
												   nDOF_test_elementIn,
												   nQuadraturePoints_elementBoundaryIn,
												   CompKernelFlag);
  }
}//proteus
#endif
