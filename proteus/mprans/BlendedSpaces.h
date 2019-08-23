#ifndef BlendedSpaces_H
#define BlendedSpaces_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"

#define USE_Q1_STENCIL 1
#define USE_MACRO_CELL 1 // for computation of gamma indicator
#define ELEMENT_BASED_ENTROPY_RESIDUAL 0 // for entropy residual 
#define C_GAMMA 5.0
// METHOD 4
// 0: low-order
// 1: high-order non-limited
// 2: limiting without gamma indicator
// 3: limiting with gamma indicator based on DK and CL
// 4: limiting with gamma indicator based on MQL, DK and CK
namespace proteus
{
  // ENTROPY FUNCTION AND ITS DERIVATIVE //
  inline double ENTROPY(const double& soln)
  {
    return 1./2.*std::pow(soln,2.);
  }
  inline double DENTROPY(const double& soln)
  {
    return soln;
  }
  // ENTROPY FLUX //
  inline double xEntFlux(const double& u_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return u_vel*ENTROPY(soln);
    else if (PROBLEM_TYPE==1)
      return soln*soln*soln/3.0;
    else
      return soln*std::sin(soln)+std::cos(soln)-1;
  }
  inline double yEntFlux(const double& v_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return v_vel*ENTROPY(soln);
    else if (PROBLEM_TYPE==1)
      return soln*soln*soln/3.0;
    else
      return soln*std::cos(soln)-std::sin(soln);
  }
  //////////////////////
  // PROBLEM ORIENTED //
  //////////////////////
  // FLUX //
  inline double xFlux(const double& u_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return u_vel*soln;
    else if (PROBLEM_TYPE==1)
      return 0.5*soln*soln;
    else return std::sin(soln);
  }
  inline double yFlux(const double& v_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return v_vel*soln;
    else if (PROBLEM_TYPE==1)
      return 0.5*soln*soln;
    else
      return std::cos(soln);
  }
  inline double xFluxJacobian(const double& u_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return u_vel;
    else if (PROBLEM_TYPE==1)
      return soln;
    else
      return std::cos(soln);
  }
  inline double yFluxJacobian(const double& v_vel, const double& soln, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
      return v_vel;
    else if (PROBLEM_TYPE==1)
      return soln;
    else
      return -std::sin(soln);
  }

  inline double compute_dij(const double& solni, const double& solnj,
			    const double& u_veli, const double& v_veli,
			    const double& u_velj, const double& v_velj,
			    const double& Cx, const double& Cy,
			    const double& CTx, const double& CTy, const int& PROBLEM_TYPE)
  {
    if (PROBLEM_TYPE==0)
    return fmax(fmax(fabs(Cx*u_veli + Cy*v_veli),
    		     fabs(Cx*u_velj + Cy*v_velj)),
    		fmax(fabs(CTx*u_veli + CTy*v_veli),
    		     fabs(CTx*u_velj + CTy*v_velj)));
    else if (PROBLEM_TYPE==1)
      {
	double lambda = fmax(fabs(solni),fabs(solnj));
	return fmax(fabs(Cx+Cy), fabs(CTx+CTy))*lambda;
      }
    else
      {
	double lambda = 1;
	return fmax(fabs(Cx+Cy), fabs(CTx+CTy))*lambda;
      }
  }
			    
    
  class BlendedSpaces_base
  {
    //The base class defining the interface
  public:
    ///////////////////////
    // Auxiliary vectors //
    ///////////////////////
    // low order vectors //
    std::valarray<double> lowOrderSolution;
    std::valarray<double> lowOrderDissipativeTerm;
    std::valarray<double> lowOrderFluxTerm;
    std::valarray<double> lowOrderBoundaryIntegral;
    // high order vectors //
    std::valarray<double> highOrderBoundaryIntegral;
    // For limiting //
    std::valarray<double> umax;
    std::valarray<double> umin;
    std::valarray<double> wBarij, wBarji;
    std::valarray<double> fluxStar;
    // for smoothness indicator //
    std::valarray<double> element_He;
    std::valarray<double> global_hi;
    std::valarray<double> min_hiHe;
    // for Zalesak's FCT //
    std::valarray<double> Rpos, Rneg;

    std::valarray<double> weighted_lumped_mass_matrix;
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
				   double* u_hess_trial_ref,
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
				   // Type of problem to solve
				   int PROBLEM_TYPE,
				   int ONE_DIM_PROBLEM,
				   int METHOD,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   double* quantDOFs4,
				   double* quantDOFs5,
				   // FOR highOrderLim
				   double* q_uInitial,
				   double* ebqe_uInlet,
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
				   // high order stab
				   double umaxG,
				   double uminG,
				   double* uHDot,
				   double* EntVisc,
				   // for smoothness indicator
				   double* gamma_dof,
				   int* first_adjacent_dof_to_middle_dof,
				   int* second_adjacent_dof_to_middle_dof,
				   double* is_dof_external,
				   double* is_dof_internal,
				   double* den_hi,
				   // C-matrices
				   double* Cx,
				   double* Cy,
				   double* CTx,
				   double* CTy,
				   double* PrCx,
				   double* PrCy,
				   double* PrCTx,
				   double* PrCTy,
				   double* CxElem,
				   double* CyElem,
				   double* CTxElem,
				   double* CTyElem,
				   double* PrCxElem,
				   double* PrCyElem,
				   double* PrCTxElem,
				   double* PrCTyElem,
				   double* dLowElem,
				   double* Q1_sparsity,
				   double* qNorm,
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
				   double* u_hess_trial_ref,
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
				   // Type of problem to solve
				   int PROBLEM_TYPE,
				   int ONE_DIM_PROBLEM,
				   int METHOD,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   double* quantDOFs4,
				   double* quantDOFs5,
				   // FOR highOrderLim
				   double* q_uInitial,
				   double* ebqe_uInlet,
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
				   // high order stab
				   double umaxG,
				   double uminG,
				   double* uHDot,
				   double* EntVisc,
				   // for smoothness indicator
				   double* gamma_dof,
				   int* first_adjacent_dof_to_middle_dof,
				   int* second_adjacent_dof_to_middle_dof,
				   double* is_dof_external,
				   double* is_dof_internal,
				   double* den_hi,
				   // C-matrices
				   double* Cx,
				   double* Cy,
				   double* CTx,
				   double* CTy,
				   double* PrCx,
				   double* PrCy,
				   double* PrCTx,
				   double* PrCTy,
				   double* CxElem,
				   double* CyElem,				   
				   double* CTxElem,
				   double* CTyElem,
				   double* PrCxElem,
				   double* PrCyElem,
				   double* PrCTxElem,
				   double* PrCTyElem,
				   double* dLowElem,
				   double* Q1_sparsity,
				   double* qNorm,
				   double* xGradRHS,
				   double* yGradRHS)=0;
    virtual void calculateResidualProjection(//element
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
				   double* u_hess_trial_ref,
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
				   // Type of problem to solve
				   int PROBLEM_TYPE,
				   int ONE_DIM_PROBLEM,
				   int METHOD,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   double* quantDOFs4,
				   double* quantDOFs5,
				   // FOR highOrderLim
				   double* q_uInitial,
				   double* ebqe_uInlet,
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
				   // high order stab
				   double umaxG,
				   double uminG,
				   double* uHDot,
				   double* EntVisc,
				   // for smoothness indicator
				   double* gamma_dof,
				   int* first_adjacent_dof_to_middle_dof,
				   int* second_adjacent_dof_to_middle_dof,
				   double* is_dof_external,
				   double* is_dof_internal,
				   double* den_hi,
				   // C-matrices
				   double* Cx,
				   double* Cy,
				   double* CTx,
				   double* CTy,
				   double* PrCx,
				   double* PrCy,
				   double* PrCTx,
				   double* PrCTy,
				   double* CxElem,
				   double* CyElem,				   
				   double* CTxElem,
				   double* CTyElem,
				   double* PrCxElem,
				   double* PrCyElem,
				   double* PrCTxElem,
				   double* PrCTyElem,
				   double* dLowElem,
				   double* Q1_sparsity,
				   double* qNorm,
				   double* xGradRHS,
				   double* yGradRHS)=0;    
    virtual void calculateMassMatrix(//element
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
				   double* dLow)=0;
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
                                       double* u_dof,
                                       double* u_exact,
				       double* gradx_u_exact,
				       double* gracy_u_exact,
                                       int offset_u, int stride_u,
                                       double* global_L1,
				       double* global_L2,
				       double* global_LInf,
                                       double* global_H1)=0;
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
			     double* u_hess_trial_ref,
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
			     // Type of problem to solve
			     int PROBLEM_TYPE,
			     int ONE_DIM_PROBLEM,
			     int METHOD,
			     // AUX QUANTITIES OF INTEREST			     
			     double* quantDOFs,
			     double* quantDOFs4,
			     double* quantDOFs5,
			     // For highOrderLim
			     double* q_uInitial,
			     double* ebqe_uInlet,
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
			     // high order stab
			     double umaxG,
			     double uminG,
			     double* uHDot,
			     double* EntVisc,
			     // for smoothness indicator
			     double* gamma_dof,
			     int* first_adjacent_dof_to_middle_dof,
			     int* second_adjacent_dof_to_middle_dof,
			     double* is_dof_external,
			     double* is_dof_internal,
			     double* den_hi,
			     // C-matrices
			     double* Cx,
			     double* Cy,
			     double* CTx,
			     double* CTy,
			     double* PrCx,
			     double* PrCy,
			     double* PrCTx,
			     double* PrCTy,
			     double* CxElem,
			     double* CyElem,
			     double* CTxElem,
			     double* CTyElem,
			     double* PrCxElem,
			     double* PrCyElem,
			     double* PrCTxElem,
			     double* PrCTyElem,
			     double* dLowElem,
			     double* Q1_sparsity,
			     double* qNorm,
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

	//////////////////////////////////////////////////////
	// ********** Zero out auxiliary vectors ********** //
	//////////////////////////////////////////////////////
	// low order vectors //
	lowOrderDissipativeTerm.resize(numDOFs,0.0);
	lowOrderFluxTerm.resize(numDOFs,0.0);	
	lowOrderBoundaryIntegral.resize(numDOFs,0.0);
	lowOrderSolution.resize(numDOFs,0.0);
	
	// limited flux correction //
	umax.resize(numDOFs,0.0);
	umin.resize(numDOFs,0.0);
	
	fluxStar.resize(numDOFs,0.0);

	wBarij.resize(nElements_global*nDOF_test_element*nDOF_trial_element,0.0);
	wBarji.resize(nElements_global*nDOF_test_element*nDOF_trial_element,0.0);
	
	// for std Zalesak's FCT //
	Rpos.resize(numDOFs,0.0);
	Rneg.resize(numDOFs,0.0);
	
	int ij=0;	
	/////////////////////////////////////////
	// ********** BOUNDARY TERM ********** //
	/////////////////////////////////////////
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
		double fluxJacobian[2];
		double uInlet = ebqe_uInlet[ebNE_kb];
		fluxJacobian[0] = xFluxJacobian(ebqe_velocity_ext[ebNE_kb_nSpace+0],uExt,PROBLEM_TYPE);
		fluxJacobian[1] = yFluxJacobian(ebqe_velocity_ext[ebNE_kb_nSpace+1],uExt,PROBLEM_TYPE);
		
		for (int I=0; I < nSpace; I++)
		  flow += normal[I]*fluxJacobian[I];

		for (int i=0;i<nDOF_test_element;i++)
		  {
		    int eN_i = eN*nDOF_test_element+i;
		    int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index
		    double uni = u_dof_old[gi];
		    //double uInlet = ebqe_uInlet[gi];
		    
		    // low order part
		    double auxLowOrder = (flow >= 0 ? 0. : 1.)*(uni-uInlet)*flow*u_test_dS[i];
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
		lowOrderBoundaryIntegral[gi] += elementBoundaryFluxLowOrder[i];
		element_flux_i[eN_i] = elementBoundaryFluxHighOrder[i]; 
	      }
	  }//ebNE
	/////////////////////////////////////////////////
	// ********** END OF BOUNDARY TERMS ********** //
	/////////////////////////////////////////////////	

	///////////////////////////////////////////////////
	// ********** FIRST LOOP ON INTEGRALS ********** // 
	///////////////////////////////////////////////////
	// * Compute first part of (lagged) edge_based_cfl = sum_el(dLowElemii)
	// * Compute dLowElem matrix
	// * Compute and save wBarij and wBarji
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double
	      elementMass[nDOF_test_element],
	      elementFluxCorrection[nDOF_test_element];
	      
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementMass[i] = 0.0;
		elementFluxCorrection[i] = 0.0;
	      }
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  un=0.0, uhDot=0.0,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_test_dV[nDOF_trial_element],
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
		ck.valFromDOF(uHDot,
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
		
		double flux[2];
		flux[0] = xFlux(velocity[eN_k_nSpace+0],un,PROBLEM_TYPE);
		flux[1] = yFlux(velocity[eN_k_nSpace+1],un,PROBLEM_TYPE);
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    elementMass[i] += u_test_dV[i];
		    elementFluxCorrection[i] += -uhDot*u_test_dV[i]
		      +ck.NumericalDiffusion(1.0,flux,&u_grad_test_dV[i_nSpace]);
		  }//i
	      }
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		
		// solution, velocity and fluxes at node i (within the element eN)
		double solni = u_dof_old[gi];
		double u_veli = u_vel_dofs[gi];
		double v_veli = v_vel_dofs[gi];
		    
		double fxi = xFlux(u_veli,solni,PROBLEM_TYPE);
		double fyi = yFlux(v_veli,solni,PROBLEM_TYPE);

		// first part of vector qi
		double qi = elementFluxCorrection[i] + elementMass[i]*uHDot[gi];
	      
		// for edge_based_cfl
		double dLowElemii = 0.;
		
		// loop in j
		for(int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    int eN_j = eN*nDOF_test_element+j;
		    int gj = offset_u+stride_u*r_l2g[eN_j];
		    
		    // solution, velocity and fluxes at node j (within the element eN)
		    double solnj = u_dof_old[gj];
		    double u_velj = u_vel_dofs[gj];
		    double v_velj = v_vel_dofs[gj];
		    
		    double fxj = xFlux(u_velj,solnj,PROBLEM_TYPE);
		    double fyj = yFlux(v_velj,solnj,PROBLEM_TYPE);
		    
		    // low-order flux part of qi
		    qi -= (CTxElem[eN_i_j]*fxj + CTyElem[eN_i_j]*fyj);
		    qi -= ((CxElem[eN_i_j]-PrCxElem[eN_i_j])*fxj +
			   (CyElem[eN_i_j]-PrCyElem[eN_i_j])*fyj);
		    
		    // compute low order flux term
		    lowOrderFluxTerm[gi] += (PrCxElem[eN_i_j]*fxj + PrCyElem[eN_i_j]*fyj);

		    int jacIndex=csrRowIndeces_CellLoops[eN_i]+csrColumnOffsets_CellLoops[eN_i_j];
		    double dLowElemij=0;
		    if (i != j)
		      {
			dLowElemij = compute_dij(solni,solnj,
						 u_veli,v_veli,
						 u_velj,v_velj,
						 PrCxElem[eN_i_j], PrCyElem[eN_i_j],
						 PrCTxElem[eN_i_j], PrCTyElem[eN_i_j],
						 PROBLEM_TYPE);
			dLowElemii -= dLowElemij;

			// save anti-dissipative flux into element_flux_qij
			element_flux_qij[eN_i_j] =
			  (fmax(EntVisc[gi],EntVisc[gj])-1.0)*dLowElemij*(solnj-solni);
			
			// compute low order dissipative term
			lowOrderDissipativeTerm[gi] += dLowElemij*(solnj-solni);
			
			// compute wBar states (wBar = 2*dij*uBar)
			wBarij[eN_i_j]=
			  (2.0*dLowElemij*(solnj+solni)/2.0
			   -(PrCxElem[eN_i_j]*(fxj-fxi) + PrCyElem[eN_i_j]*(fyj-fyi)));
			wBarji[eN_i_j]=
			  (2.0*dLowElemij*(solnj+solni)/2.0
			   -(PrCTxElem[eN_i_j]*(fxi-fxj) + PrCTyElem[eN_i_j]*(fyi-fyj)));

			// save low order matrix
			dLowElem[eN_i_j] = dLowElemij;
		      }
		    else
		      {
			dLowElem[eN_i_j] = 0; // not true but irrelevant since we know that fii=0
		      }
		  }//j
		element_flux_i[eN_i] += qi;
		// No need to compute diagonal entry of wBarij and wBarji since fStarij=0 when i=j
		// compute dLowii and store it in edge_based_cfl
		edge_based_cfl[gi] += dLowElemii;
	      }//i
	  }//elements
	////////////////////////////////////
	// END OF FIRST LOOP IN INTEGRALS //
	////////////////////////////////////
 
	////////////////////////
	// FIRST LOOP IN DOFs //
	////////////////////////
	// * Compute edge_based_cfl = 2|dLowii|/mi
	// * Compute umax and umin
	// * Compute the low order solution
	ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double mi = QH_ML[i];
	    // compute edge_based_cfl
	    edge_based_cfl[i] = 2*fabs(edge_based_cfl[i])/mi;

	    // solution at node i
	    double solni = u_dof_old[i];
	    
	    // for the computation of the local bounds
	    double umaxi = solni;
	    double umini = solni;

	    // loop over the sparsity pattern of the i-th DOF
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		// solution at node j 
		double solnj = u_dof_old[j];

		if (USE_Q1_STENCIL==1)
		  {
		    if (Q1_sparsity[ij]==1)
		      {
			// computation of local bounds
			umaxi = fmax(solnj,umaxi);
			umini = fmin(solnj,umini);
		      }
		  }
		else
		  {
		    umaxi = fmax(solnj,umaxi);
		    umini = fmin(solnj,umini);
		  }
		ij+=1;
	      }
	    // Save local bounds 
	    umax[i] = umaxi;
	    umin[i] = umini;
	    
	    // compute and save low order solution 
	    lowOrderSolution[i] = solni - dt/mi * (lowOrderFluxTerm[i] 
						   - lowOrderDissipativeTerm[i]
						   - lowOrderBoundaryIntegral[i]);
	  }
	///////////////////////////////
	// END OF FIRST LOOP IN DOFs //
	///////////////////////////////
	
	/////////////////////////////
	// SECOND LOOP ON ELEMENTS //
	/////////////////////////////
	// * Compute vVector = inv(element_ML_minus_MC)*element_flux_i
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    // compute element vector v //
	    double mean = 0;
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for(int j=0;j<nDOF_test_element;j++)
		  {
		    int eN_j = eN*nDOF_test_element+j;
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    vVector[eN_i] += inv_element_ML_minus_MC[eN_i_j]*element_flux_i[eN_j];
		    mean += inv_element_ML_minus_MC[eN_i_j]*element_flux_i[eN_j];
		  }
	      }
	    mean /= nDOF_test_element;
	    // substract mean
	    for(int i=0;i<nDOF_test_element;i++)
	      {
	    	int eN_i = eN*nDOF_test_element+i;
	    	vVector[eN_i] -= mean;
	      }
	    // end of computation of the vector v //
	  }//elements
	/////////////////////////////////////
	// END OF SECOND LOOP IN INTEGRALS //
	/////////////////////////////////////

	//////////////////////////////
	// compute element_flux_qij //
	//////////////////////////////
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		for(int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_j = eN*nDOF_test_element+j;
		    int eN_i_j = eN_i*nDOF_trial_element+j;

		    element_flux_qij[eN_i_j] += element_MC[eN_i_j]*(vVector[eN_i]-vVector[eN_j]);
		  }
	      }
	  }
	
	/////////////////////////
	// * compute elementFluxStar_ij and fluxStar_i = sum_el(sum_j(elementFluxStar_ij))
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		int eN_i = eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		
		// get bounds at node i
		double umini = umin[gi];
		double umaxi = umax[gi];
		
		for(int j=0;j<nDOF_trial_element;j++)
		  {
		    int eN_j = eN*nDOF_test_element+j;
		    int eN_i_j = eN_i*nDOF_trial_element+j;
		    int gj = offset_u+stride_u*r_l2g[eN_j];
		    
		    // get bounds at node j
		    double uminj = umin[gj];
		    double umaxj = umax[gj];
		    
		    // compute element flux star
		    double fij = element_flux_qij[eN_i_j];
		    double dij = dLowElem[eN_i_j];
		    double wij = wBarij[eN_i_j];
		    double wji = wBarji[eN_i_j];
		    
		    if (i!=j)
		      {
			double gammaij = fmax(gamma_dof[i], gamma_dof[j]);
			if (METHOD==4)
			  gammaij = fmax(qNorm[i], qNorm[j]);
			if (fij > 0)
			  {
			    if (METHOD==0)
			      fluxStar[gi] += 0; // low-order method
			    else if (METHOD==1)
			      fluxStar[gi] += fij; // high-order entropy viscosity
			    else if (METHOD==2)
				fluxStar[gi] += fmin(fij,fmin(2*dij*umaxi-wij,wji-2*dij*uminj));
			    else // METHOD==3 or 4
			      {
				double fijStarLoc=fmin(fij,fmin(2*dij*umaxi-wij,wji-2*dij*uminj));
				double fijStarGlob=fmin(fij,fmin(2*dij*umaxG-wij,wji-2*dij*uminG));
				fluxStar[gi] += fmin(fijStarGlob, fmax(fijStarLoc,gammaij*fij));
				//fluxStar[gi] += fmax(fijStarLoc,gammaij*fij);
			      }
			  }
			else
			  {
			    if (METHOD==0)
			      fluxStar[gi] = 0;
			    else if (METHOD==1)
			      fluxStar[gi] += fij;
			    else if (METHOD==2)
			      fluxStar[gi] += fmax(fij,fmax(2*dij*umini-wij,wji-2*dij*umaxj));
			    else // METHOD==3 or 4
			      {
				double fijStarLoc =fmax(fij,fmax(2*dij*umini-wij,wji-2*dij*umaxj));
				double fijStarGlob=fmax(fij,fmax(2*dij*uminG-wij,wji-2*dij*umaxG));
				fluxStar[gi] += fmax(fijStarGlob, fmin(fijStarLoc,gammaij*fij));
				//fluxStar[gi] += fmin(fijStarLoc,gammaij*fij);
			      }
			  }
		      }
		  } //j
	      } //i
	  } //element
	
	for (int i=0; i<numDOFs; i++)
	  {
	    double mi = QH_ML[i];
	    // COMPUTE SOLUTION //
	    //globalResidual[i] = lowOrderSolution[i];
	    globalResidual[i] = lowOrderSolution[i] + dt/mi*fluxStar[i];
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
					double* u_hess_trial_ref,
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
					// Type of problem to solve
					int PROBLEM_TYPE,
					int ONE_DIM_PROBLEM,
					int METHOD,
					// AUX QUANTITIES OF INTEREST			     
					double* quantDOFs,
					double* quantDOFs4,
					double* quantDOFs5,
					// For highOrderLim
					double* q_uInitial,
					double* ebqe_uInlet,
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
					// high order stab
					double umaxG,
					double uminG,
					double* uHDot,
					double* EntVisc,
					// for smoothness indicator
					double* gamma_dof,
					int* first_adjacent_dof_to_middle_dof,
					int* second_adjacent_dof_to_middle_dof,
					double* is_dof_external,
					double* is_dof_internal,
					double* den_hi,
					// C-matrices
					double* Cx,
					double* Cy,
					double* CTx,
					double* CTy,
					double* PrCx,
					double* PrCy,
					double* PrCTx,
					double* PrCTy,
					double* CxElem,
					double* CyElem,
					double* CTxElem,
					double* CTyElem,
					double* PrCxElem,
					double* PrCyElem,
					double* PrCTxElem,
					double* PrCTyElem,
					double* dLowElem,
					double* Q1_sparsity,
					double* qNorm,
					double* xGradRHS,
					double* yGradRHS)
      {
	////////////////////////////////
	// Zero out auxiliary vectors //
	////////////////////////////////
	//highOrderBoundaryIntegral.resize(numDOFs,0.0);

	// for entropy viscosity
	element_He.resize(nElements_global,0.0);
	global_hi.resize(numDOFs,0.0);
	min_hiHe.resize(numDOFs,1E100);

	weighted_lumped_mass_matrix.resize(numDOFs,0.0);

	for (int i=0; i<numDOFs; i++)
	  {
	    xGradRHS[i]=0.0;
	    yGradRHS[i]=0.0;
	    qNorm[i] = 0.0;
	  }
	///////////////////
	// BOUNDARY TERM //
	///////////////////
	// * Compute boundary integrals 
	/* for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) */
	/*   { */
	/*     register int ebN = exteriorElementBoundariesArray[ebNE]; */
	/*     register int eN  = elementBoundaryElementsArray[ebN*2+0], */
	/*       ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0], */
	/*       eN_nDOF_trial_element = eN*nDOF_trial_element; */
	/*     register double elementBoundaryFluxHighOrder[nDOF_test_element]; */
	/*     for (int i=0;i<nDOF_test_element;i++) */
	/*       elementBoundaryFluxHighOrder[i]=0.0; */
	/*     // loop on quad points */
	/*     for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
	/*       { */
	/* 	register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb, */
	/* 	  ebNE_kb_nSpace = ebNE_kb*nSpace, */
	/* 	  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb, */
	/* 	  ebN_local_kb_nSpace = ebN_local_kb*nSpace; */
	/* 	register double */
	/* 	  uExt=0., uInlet=0., */
	/* 	  jac_ext[nSpace*nSpace], */
	/* 	  jacDet_ext, */
	/* 	  jacInv_ext[nSpace*nSpace], */
	/* 	  boundaryJac[nSpace*(nSpace-1)], */
	/* 	  metricTensor[(nSpace-1)*(nSpace-1)], */
	/* 	  metricTensorDetSqrt, */
	/* 	  dS, */
	/* 	  u_test_dS[nDOF_test_element], */
	/* 	  normal[nSpace],x_ext,y_ext,z_ext; */
	/* 	// calculate mappings */
	/* 	ck.calculateMapping_elementBoundary(eN, */
	/* 					    ebN_local, */
	/* 					    kb, */
	/* 					    ebN_local_kb, */
	/* 					    mesh_dof, */
	/* 					    mesh_l2g, */
	/* 					    mesh_trial_trace_ref, */
	/* 					    mesh_grad_trial_trace_ref, */
	/* 					    boundaryJac_ref, */
	/* 					    jac_ext, */
	/* 					    jacDet_ext, */
	/* 					    jacInv_ext, */
	/* 					    boundaryJac, */
	/* 					    metricTensor, */
	/* 					    metricTensorDetSqrt, */
	/* 					    normal_ref, */
	/* 					    normal, */
	/* 					    x_ext,y_ext,z_ext); */
	/* 	dS = ((1.0)*metricTensorDetSqrt )*dS_ref[kb]; */
	/* 	ck.valFromDOF(u_dof_old, */
	/* 		      &u_l2g[eN_nDOF_trial_element], */
	/* 		      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element], */
	/* 		      uExt); */
	/* 	ck.valFromDOF(uInlet_dofs, */
	/* 		      &u_l2g[eN_nDOF_trial_element], */
	/* 		      &u_trial_trace_ref[ebN_local_kb*nDOF_test_element], */
	/* 		      uInlet); */
	/* 	//precalculate test function products with integration weights */
	/* 	for (int j=0;j<nDOF_trial_element;j++) */
	/* 	  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS; */
		
	/* 	//calculate flow */
	/* 	double flow = 0.; */
	/* 	double fluxJacobian[2]; */
	/* 	fluxJacobian[0] = xFluxJacobian(ebqe_velocity_ext[ebNE_kb_nSpace+0],uExt,PROBLEM_TYPE); */
	/* 	fluxJacobian[1] = yFluxJacobian(ebqe_velocity_ext[ebNE_kb_nSpace+1],uExt,PROBLEM_TYPE); */
		
	/* 	for (int I=0; I < nSpace; I++) */
	/* 	  flow += normal[I]*fluxJacobian[I]; */
		
	/* 	for (int i=0;i<nDOF_test_element;i++) */
	/* 	  { */
	/* 	    int eN_i = eN*nDOF_test_element+i; */
	/* 	    // high order flux */
	/* 	    double auxHighOrder = (flow >= 0 ? uExt : uInlet)*flow*u_test_dS[i]; */
	/* 	    elementBoundaryFluxHighOrder[i] += auxHighOrder; */
	/* 	  } */
	/*       }//kb */
	/*     // save to the corresponding vectors */
	/*     for (int i=0;i<nDOF_test_element;i++) */
	/*       { */
	/* 	int eN_i = eN*nDOF_test_element+i; */
	/* 	int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index */
	/* 	highOrderBoundaryIntegral[gi] += elementBoundaryFluxHighOrder[i]; */
	/*       } */
	/*   }//ebNE */
	///////////////////////////
	// END OF BOUNDARY TERMS //
	///////////////////////////
	
	/////////////////////////////
	// FIRST LOOP ON INTEGRALS // Compute uDot via high order galerkin (with lumped m.mat)
	/////////////////////////////
	int nSpace2 = nSpace*nSpace;
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    double det_hess_Ke=0, area_Ke=0;
	    register double
	      element_weighted_lumped_mass_matrix[nDOF_test_element],
	      element_rhs_qx[nDOF_test_element],
	      element_rhs_qy[nDOF_test_element],
	      elementTimeDerivative[nDOF_test_element],
	      elementFluxTerm[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      {
		elementTimeDerivative[i]=0.0;
		elementFluxTerm[i]=0.0;
		element_rhs_qx[i]=0.0;
		element_rhs_qy[i]=0.0;
		element_weighted_lumped_mass_matrix[i]=0.0;
	      }
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  rhsx=0.0, rhsy=0.0, norm_grad_un=0.0,
		  un=0.0, grad_un[nSpace], hess_un[nSpace2], det_hess_un=0.,
		  uDot=0.0,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_test_dV[nDOF_test_element],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_grad_test_dV[nDOF_test_element*nSpace],
		  u_hess_trial[nDOF_trial_element*nSpace2],
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
		ck.hessTrialFromRef(&u_hess_trial_ref[k*nDOF_trial_element*nSpace2],
				    jacInv,
				    u_hess_trial);
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof_old,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      un); // from high-order space
		ck.gradFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_un);
		ck.hessFromDOF(u_dof_old,
			       &u_l2g[eN_nDOF_trial_element],
			       u_hess_trial,
			       hess_un);
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      uDot); // from high-order space
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }

		rhsx = grad_un[0];
                rhsy = grad_un[1];
		norm_grad_un=0.0;
                for (int I=0;I<nSpace; I++)
                  norm_grad_un += grad_un[I]*grad_un[I];
                norm_grad_un = std::sqrt(norm_grad_un+1E-5);
		
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    element_weighted_lumped_mass_matrix[i] += norm_grad_un*u_test_dV[i];
                    element_rhs_qx[i] += rhsx*u_test_dV[i];
                    element_rhs_qy[i] += rhsy*u_test_dV[i];
		    
		    elementTimeDerivative[i] += uDot*u_test_dV[i]; 
		    register int i_nSpace=i*nSpace;
		    if (PROBLEM_TYPE==0)
		      elementFluxTerm[i] += (velocity[eN_k_nSpace+0]*grad_un[0]
					     +velocity[eN_k_nSpace+1]*grad_un[1])*u_test_dV[i];
		    else if (PROBLEM_TYPE==1)
		      elementFluxTerm[i] += un*(grad_un[0]+grad_un[1])*u_test_dV[i];
		    else // KPP
		      elementFluxTerm[i] += (std::cos(un)-std::sin(un))*u_test_dV[i];
		  }//i
		// To compute smoothness indicator based on 2nd derivative //
		
		if (ONE_DIM_PROBLEM==1)
		  det_hess_un = hess_un[0];
		else
		  det_hess_un = hess_un[0]*hess_un[3] - hess_un[2]*hess_un[1];
		det_hess_Ke += det_hess_un*dV;
		area_Ke += dV;
	      }//kb
	    //
	    //load element into global residual and save element residual
	    //
	    element_He[eN] = det_hess_Ke/area_Ke;
	    for(int i=0;i<nDOF_test_element;i++)
	      {		
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];

		weighted_lumped_mass_matrix[gi] += element_weighted_lumped_mass_matrix[i];
		xGradRHS[gi] += element_rhs_qx[i];
		yGradRHS[gi] += element_rhs_qy[i];
		
		// residual of uDot
		globalResidual[gi] += elementTimeDerivative[i] + elementFluxTerm[i];
		// global hi (to compute smoothness indicator)
		global_hi[gi] += element_He[eN];
	      } //i
	  } //elements

	for (int i=0; i<numDOFs; i++)
	  {
	    global_hi[i] /= den_hi[i];
	    xGradRHS[i] /= weighted_lumped_mass_matrix[i];
	    yGradRHS[i] /= weighted_lumped_mass_matrix[i];	    
	    qNorm[i] = std::sqrt(xGradRHS[i]*xGradRHS[i] + yGradRHS[i]*yGradRHS[i]);
	    qNorm[i] = 1.0 - (qNorm[i] < 0.99 ? 0.0 : 1.0);

	    /*
	    double solni = u_dof_old[i];
	    double num_betai = 0;
	    double den_betai = 0;
	    
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];
		double solnj = u_dof_old[j];

		num_betai += solni-solnj;
		den_betai += fabs(solni-solnj);
	      }
	    qNorm[i] = qNorm[i]*std::pow(fabs(num_betai)/(den_betai+1E-10),2.0);
	    */
	  }
	
	///////////////////
	// LOOP IN CELLS // to compute min(hi*He)
	///////////////////
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    double He = element_He[eN];
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		register int gi = offset_u+stride_u*r_l2g[eN_i];
		double hi = global_hi[gi];
		min_hiHe[gi] = fmin(min_hiHe[gi], hi*He);
	      }
	  }

	for (int i=0; i<numDOFs; i++)
	  {
	    quantDOFs4[i] = global_hi[i]*global_hi[i];
	    quantDOFs5[i] = min_hiHe[i];
	  }
	  
	// ***************************** //
	// ***** Entropy viscosity ***** //
	// ***************************** //
	int ij=0;
	for (int i=0; i<numDOFs; i++)
	  {
	    double solni = u_dof_old[i];
	    double u_veli = u_vel_dofs[i];
	    double v_veli = v_vel_dofs[i];
	    double DEnti = DENTROPY(solni);
	    double DenEntViscPart1 = 0.;
	    double DenEntViscPart2 = 0.;
	    double ith_NumEntVisc = 0.;
	    double ith_DenEntVisc = 0.;	      
	    // loop on the support of ith shape function
	    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
	      {
		int j = colind[offset];			
		// solution, velocity and fluxes at node j
		double solnj = u_dof_old[j];
		double u_velj = u_vel_dofs[j];
		double v_velj = v_vel_dofs[j];
		
		double fxj = xFlux(u_velj,solnj,PROBLEM_TYPE);
		double fyj = yFlux(v_velj,solnj,PROBLEM_TYPE);
		
		// For entropy viscosity //
		double x_EntFluxj=xEntFlux(u_velj,solnj,PROBLEM_TYPE);// - xEntFlux(u_veli,solni,PROBLEM_TYPE);
		double y_EntFluxj=yEntFlux(v_velj,solnj,PROBLEM_TYPE);// - yEntFlux(v_veli,solni,PROBLEM_TYPE);
		ith_NumEntVisc += (Cx[ij]*(x_EntFluxj-DEnti*fxj) +
				   Cy[ij]*(y_EntFluxj-DEnti*fyj));
		// aux parts to compute DenEntVisc
		DenEntViscPart1 += Cx[ij]*x_EntFluxj + Cy[ij]*y_EntFluxj;
		DenEntViscPart2 += Cx[ij]*fxj + Cy[ij]*fyj;
		
		ij+=1;
	      } //j
	    ith_DenEntVisc = (fabs(DenEntViscPart1) + fabs(DEnti)*fabs(DenEntViscPart2)+1E-15);
	    EntVisc[i] = std::pow(fabs(ith_NumEntVisc)/ith_DenEntVisc,1.0);
	    // ***** End of entropy viscosity ***** //
	  }
	
	// *************************** //
	// ***** Gamma indicator ***** //
	// *************************** //
	if (USE_MACRO_CELL==1)
	  {
	    for (int i=0; i<numDOFs; i++)
	      {
		// compute gamma in DOFs owned by big Q1 mesh
		if (is_dof_external[i] == 1) // if dof is in corner of macro cell 
		  {
		    double eps = 1E-10;
		    double hi2 = global_hi[i] * global_hi[i];
		    //gamma_dof[i] = fmax(0, fmin(hi2,C_GAMMA*min_hiHe[i]))/(hi2+eps);
		    gamma_dof[i] = (fmax(0, fmin(hi2, C_GAMMA*min_hiHe[i])) + eps)/(hi2+eps);
		    //gamma_dof[i] = hi2 == 0 ? 1. : fmax(0, fmin(hi2, C_GAMMA*min_hiHe[i]))/hi2;
		  }// i
		else
		  {
		    gamma_dof[i] = 0;
		  }
	      } //i
	    // make average to "interpolate" gamma from big Q1 mesh to finer mesh
	    for (int i=0; i<numDOFs; i++)
	      {
		if (is_dof_internal[i] == 1) // dof is internal
		  {
		    double num_external_dofs = 0; // this should be 4
		    // loop on its support
		    for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		      {
			int j = colind[offset];
			if (is_dof_external[j] == 1) // external j-dof
			  {
			    num_external_dofs += 1;
			    gamma_dof[i] += gamma_dof[j];
			  }
		      }
		    // normalize
		    gamma_dof[i] /= num_external_dofs;
		  }
		else if (is_dof_external[i] == 0) // this is a middle dof
		  {
		    int gj1 = first_adjacent_dof_to_middle_dof[i];
		    int gj2 = second_adjacent_dof_to_middle_dof[i];
		    gamma_dof[i] = 0.5*(gamma_dof[gj1] + gamma_dof[gj2]);
		  }
	      }
	  }
	else
	  {
	    // compute gamma in DOFs owned by big Q1 mesh
	    for (int i=0; i<numDOFs; i++)
	      {
		double eps = 1E-10;
		double hi2 = global_hi[i] * global_hi[i];
		//gamma_dof[i] = fmax(0, fmin(hi2, C_GAMMA*min_hiHe[i]))/(hi2+eps);
		gamma_dof[i] = (fmax(0, fmin(hi2, C_GAMMA*min_hiHe[i])) + eps)/(hi2+eps);
		//gamma_dof[i] = hi2 == 0 ? 1. : fmax(0, fmin(hi2, C_GAMMA*min_hiHe[i]))/hi2;
	      }
	  }
	//////////
	
	//////////////////////
	// ENTROPY RESIDUAL // (based on elements)
	//////////////////////
	if (ELEMENT_BASED_ENTROPY_RESIDUAL==1)
	  {
	    for(int eN=0;eN<nElements_global;eN++)
	      {
		for(int i=0;i<nDOF_test_element;i++)
		  {
		    int eN_i = eN*nDOF_test_element+i;
		    int gi = offset_u+stride_u*r_l2g[eN_i]; //global i-th index
		    
		    // for entropy viscosity //
		    double solni = u_dof_old[gi];
		    double u_veli = u_vel_dofs[gi];
		    double v_veli = v_vel_dofs[gi];
		    double DEnti = DENTROPY(solni);
		    double DenEntViscPart1 = 0.;
		    double DenEntViscPart2 = 0.;
		    double ith_NumEntVisc = 0.;
		    double ith_DenEntVisc = 0.;
		    ///////////////////////////
		    
		    // loop in j
		    for(int j=0;j<nDOF_trial_element;j++)
		      {
			int eN_i_j = eN_i*nDOF_trial_element+j;
			int eN_j = eN*nDOF_test_element+j;
			int gj = offset_u+stride_u*r_l2g[eN_j];
			
			// solution, velocity and fluxes at node j (within the element eN)
			double solnj = u_dof_old[gj];
			double u_velj = u_vel_dofs[gj];
			double v_velj = v_vel_dofs[gj];
			
			double fxj = xFlux(u_velj,solnj,PROBLEM_TYPE);
			double fyj = yFlux(v_velj,solnj,PROBLEM_TYPE);
			
			// For entropy viscosity //
			double x_EntFluxj = xEntFlux(u_velj,solnj,PROBLEM_TYPE);// - xEntFlux(u_veli,solni,PROBLEM_TYPE);
			double y_EntFluxj = yEntFlux(v_velj,solnj,PROBLEM_TYPE);// - yEntFlux(v_veli,solni,PROBLEM_TYPE);
			ith_NumEntVisc += (PrCxElem[eN_i_j]*(x_EntFluxj-DEnti*fxj) +
					   PrCyElem[eN_i_j]*(y_EntFluxj-DEnti*fyj));
			// aux parts to compute DenEntVisc
			DenEntViscPart1 += PrCxElem[eN_i_j]*x_EntFluxj + PrCyElem[eN_i_j]*y_EntFluxj;
			DenEntViscPart2 += PrCxElem[eN_i_j]*fxj + PrCyElem[eN_i_j]*fyj;
		      }//j
		    // compute DenEntVisc //
		    ith_DenEntVisc=(fabs(DenEntViscPart1) + fabs(DEnti)*fabs(DenEntViscPart2)+1E-15);
		    EntVisc[gi] = std::pow(fabs(ith_NumEntVisc)/ith_DenEntVisc,1.0);
		  }//i
	      }//elements
	  }
      }

      void calculateResidualProjection(//element
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
					double* u_hess_trial_ref,
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
					// Type of problem to solve
					int PROBLEM_TYPE,
					int ONE_DIM_PROBLEM,
					int METHOD,
					// AUX QUANTITIES OF INTEREST			     
					double* quantDOFs,
					double* quantDOFs4,
					double* quantDOFs5,
					// For highOrderLim
					double* q_uInitial,
					double* ebqe_uInlet,
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
					// high order stab
					double umaxG,
					double uminG,
					double* uHDot,
					double* EntVisc,
					// for smoothness indicator
					double* gamma_dof,
					int* first_adjacent_dof_to_middle_dof,
					int* second_adjacent_dof_to_middle_dof,
					double* is_dof_external,
					double* is_dof_internal,
					double* den_hi,
					// C-matrices
					double* Cx,
					double* Cy,
					double* CTx,
					double* CTy,
					double* PrCx,
					double* PrCy,
					double* PrCTx,
					double* PrCTy,
					double* CxElem,
					double* CyElem,
					double* CTxElem,
					double* CTyElem,
					double* PrCxElem,
					double* PrCyElem,
					double* PrCTxElem,
					double* PrCTyElem,
					double* dLowElem,
					double* Q1_sparsity,
					double* qNorm,
					double* xGradRHS,
					double* yGradRHS)
      {
	///////////////////////
	// LOOP ON INTEGRALS //
	///////////////////////
	for(int eN=0;eN<nElements_global;eN++) //loop in cells
	  {
	    //declare local storage for element residual and initialize
	    register double elementResidual[nDOF_test_element];
	    for (int i=0;i<nDOF_test_element;i++)
	      elementResidual[i]=0.0;
	    
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  uh=0.0, uInitial=0.0,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_test_dV[nDOF_test_element],
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
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      uh); // from high-order space
		uInitial = q_uInitial[eN_k];
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

		for(int i=0;i<nDOF_test_element;i++)
		  {
		    elementResidual[i] += (uh-uInitial)*u_test_dV[i]; 
		  }//i
	      }//kb
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {		
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		globalResidual[gi] += elementResidual[i];
	      } //i
	  } //elements
      }
      
      void calculateMassMatrix(//element
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
			       double* dLow)
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
                                 double* u_dof,
				 double* u_exact,
				 double* gradx_u_exact,
				 double* grady_u_exact,
                                 int offset_u, int stride_u,
				 double* global_L1,
				 double* global_L2,
				 double* global_LInf,
				 double* global_H1)
      {
	*global_L1 = 0.0;
	*global_L2 = 0.0;
	*global_LInf = 0.0;
	*global_H1 = 0.0;
        //////////////////////
        // ** LOOP IN CELLS //
        //////////////////////
	double error_LInf = 0.0;
        for(int eN=0;eN<nElements_global;eN++)
          {
            if (eN<nElements_owned) // just consider the locally owned cells
              {
                //declare local storage for local contributions and initialize
                double cell_L1 = 0., cell_L2 = 0., cell_H1 = 0.;

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
					u_grad_trial); // high-order space
		    // load exact solution and its gradient
		    u = u_exact[eN_k];
		    gradx_u = gradx_u_exact[eN_k];
		    grady_u = grady_u_exact[eN_k];
                    // get numerical solution
		    ck.valFromDOF(u_dof,
				  &u_l2g[eN_nDOF_trial_element],
				  &u_trial_ref[k*nDOF_trial_element],
				  uh); // from high-order space
		    ck.gradFromDOF(u_dof,
				   &u_l2g[eN_nDOF_trial_element],
				   u_grad_trial,
				   grad_uh);
		    // Norms in whole domain //
		    cell_L2 += std::pow(u-uh,2)*dV;
		    cell_H1 += (std::pow(u-uh,2) +
		    		(std::pow(grad_uh[0] - gradx_u,2.0) +
		    		 std::pow(grad_uh[1] - grady_u,2.0)) )*dV;
		    cell_L1 += fabs(u-uh)*dV;
		    error_LInf = fmax(fabs(u-uh), error_LInf);
                  }
		*global_L1 += cell_L1;
		*global_L2 += cell_L2;
		*global_H1 += cell_H1;
              }//elements
          }
	*global_LInf = error_LInf;
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
