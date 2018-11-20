#ifndef VOF_H
#define VOF_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define POWER_SMOOTHNESS_INDICATOR 2
#define IS_BETAij_ONE 0
#define GLOBAL_FCT 0

/////////////////////
//ENTROPY FUNCTION //
/////////////////////
// Power entropy //
#define entropy_power 2. // phiL and phiR are dummy variables
#define ENTROPY(phi,phiL,phiR) 1./entropy_power*std::pow(fabs(phi),entropy_power)
#define DENTROPY(phi,phiL,phiR) std::pow(fabs(phi),entropy_power-1.)*(phi>=0 ? 1 : -1)
// Log entropy //
// LOG ENTROPY FOR LEVEL SET FROM 0 to 1
#define ENTROPY_LOG(phi,phiL,phiR) std::log(fabs((phi-phiL)*(phiR-phi))+1E-14)
#define DENTROPY_LOG(phi,phiL,phiR) (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(fabs((phi-phiL)*(phiR-phi))+1E-14) 

namespace proteus
{
  class VOF_base
  {
    //The base class defining the interface
  public:
    virtual ~VOF_base(){}
    virtual void calculateResidual_blending_spaces(//element
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
						     int lag_shockCapturing,
						     double shockCapturingDiffusion,
						     double sc_uref, 
						     double sc_alpha,
						     //VRANS
						     const double* q_porosity,
						     const double* porosity_dof,
						     //
						     int* u_l2g, 
						     double* elementDiameter,
						     double* meshSize,
						     int degree_polynomial,
						     double* u_dof,
						     double* u_dof_old,
						     double* velocity,
						     double* q_m,
						     double* q_u,
						     double* q_grad_u,
						     double* q_m_betaBDF,
						     double* q_dV,
						     double* q_dV_last,
						     double* cfl,
						     double* edge_based_cfl,
						     double* q_numDiff_u, 
						     double* q_numDiff_u_last, 
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
						     // PARAMETERS FOR EDGE BASED STABILIZATION
						     double cE,
						     double cK,
						     // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
						     double uL, 
						     double uR, 
						     // PARAMETERS FOR EDGE VISCOSITY 
						     int numDOFs,
						     int NNZ,
						     int* csrRowIndeces_DofLoops,
						     int* csrColumnOffsets_DofLoops,
						     int* csrRowIndeces_CellLoops,
						     int* csrColumnOffsets_CellLoops,
						     int* csrColumnOffsets_eb_CellLoops,
						     // C matrices
						     double* Cx, 
						     double* Cy,
						     double* Cz,
						     double* CTx,
						     double* CTy,
						     double* CTz,
						     double* ML,
						     // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
						     int LUMPED_MASS_MATRIX, 
						     int STABILIZATION_TYPE,
						     int ENTROPY_TYPE,
						     // FOR FCT
						     double* low_order_solution,
						     double* dt_times_dH_minus_dL,
						     double* min_u_bc,
						     double* max_u_bc,
						     // FOR BLENDING SPACES
						     double* force,
						     double* uexact,
						     double* gradx_uexact,
						     double* grady_uexact,
						     double* alpha_value,
						     double* alpha_dof,
						     double* aux_test_ref,
						     double* aux_grad_test_ref,
						     // AUX QUANTITIES OF INTEREST
						     double* dLow,
						     double* quantDOFs,
						     double beta,
						     double epsilon)=0;    
    virtual void calculateJacobian_blending_spaces(//element
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
				   int lag_shockCapturing,/*mwf not used yet*/
				   double shockCapturingDiffusion,
				   //VRANS
				   const double* q_porosity,
				   //
				   int* u_l2g,
				   double* elementDiameter,
				   int degree_polynomial,
				   double* u_dof, 
				   double* velocity,
				   double* q_m_betaBDF, 
				   double* cfl,
				   double* q_numDiff_u_last, 
				   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				   int numDOFs,
				   int* csrRowIndeces_DofLoops,
				   int* csrColumnOffsets_DofLoops,
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
				   int LUMPED_MASS_MATRIX,
				   // FOR BLENDING SPACES
				   double* alpha_value,
				   double* alpha_dof,
				   double* aux_test_ref,
				   double* aux_grad_test_ref,
				   double* dLow,
				   double beta,
				   double epsilon
						   )=0;    
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
    virtual void getLumpedL2Projection(double* mesh_trial_ref,
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
				       double* q_alpha,
				       int offset_u, int stride_u,
				       int numDOFs,
				       double* lumpedL2Projection)=0;
  };

  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class VOF : public VOF_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    VOF():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}
    inline
    void evaluateCoefficients(const double v[nSpace],
			      const double& u,
			      const double& porosity, //VRANS specific
			      double& m,
			      double& dm,
			      double f[nSpace],
			      double df[nSpace])
    {
    m = porosity*u;
    dm= porosity;
    for (int I=0; I < nSpace; I++)
      {
	f[I] = v[I]*porosity*u;
	df[I] = v[I]*porosity;
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
	  //std::cout << mat_times_vector[I] << std::endl;
	}
    }

    inline double Dot(const double vec1[nSpace],
		      const double vec2[nSpace])
    {
      double dot = 0;
      for (int I=0; I<nSpace; I++)
	dot += vec1[I]*vec2[I];
      return dot;
    }    

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
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double   G[nSpace*nSpace],
					const double&  A0,
					const double   Ai[nSpace],
					double& tau_v,
					double& cfl)	
    {
      double v_d_Gv=0.0; 
      for(int I=0;I<nSpace;I++) 
         for (int J=0;J<nSpace;J++) 
           v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];     
    
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + 1.0e-8);    
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
	      std::cout<<"warning: VOF open boundary with no external trace, setting to zero for inflow"<<std::endl;
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

    void calculateResidual_blending_spaces(//element
			   double dt,
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   double* mesh_velocity_dof,
			   double MOVING_DOMAIN,
			   int* mesh_l2g,
			   double* dV_ref,
			   double* u_trial_ref, ///
			   double* u_grad_trial_ref, ///
			   double* u_test_ref, ///
			   double* u_grad_test_ref, ///
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
			   int lag_shockCapturing, /*mwf not used yet*/
			   double shockCapturingDiffusion,
			   double sc_uref, double sc_alpha,
			   //VRANS
			   const double* q_porosity,
			   const double* porosity_dof,
			   //
			   int* u_l2g, 
			   double* elementDiameter,
			   double* meshSize,
			   int degree_polynomial,
			   double* u_dof,
			   double* u_dof_old,
			   double* velocity,
			   double* q_m,
			   double* q_u,
			   double* q_grad_u,
			   double* q_m_betaBDF,
                           double* q_dV,
                           double* q_dV_last,
			   double* cfl,
			   double* edge_based_cfl,
			   double* q_numDiff_u, 
			   double* q_numDiff_u_last, 
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
			   // PARAMETERS FOR EDGE BASED STABILIZATION
			   double cE,
			   double cK,
			   // PARAMETERS FOR LOG BASED ENTROPY FUNCTION 
			   double uL, 
			   double uR,
			   // PARAMETERS FOR EDGE VISCOSITY 
			   int numDOFs,
			   int NNZ,
			   int* csrRowIndeces_DofLoops,
			   int* csrColumnOffsets_DofLoops,
			   int* csrRowIndeces_CellLoops,
			   int* csrColumnOffsets_CellLoops,
			   int* csrColumnOffsets_eb_CellLoops,
			   // C matrices
			   double* Cx, 
			   double* Cy, 
			   double* Cz, 
			   double* CTx,
			   double* CTy, 
			   double* CTz, 
			   double* ML,
			   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
			   int LUMPED_MASS_MATRIX, 
			   int STABILIZATION_TYPE,
			   int ENTROPY_TYPE,
			   // FOR FCT
			   double* low_order_solution,
			   double* dt_times_dH_minus_dL,
			   double* min_u_bc,
			   double* max_u_bc,
			   // FOR BLENDING SPACES
			   double* force,
			   double* uexact,
			   double* gradx_uexact,
			   double* grady_uexact,
			   double* alpha_value,
			   double* alpha_dof,
			   double* aux_test_ref,
			   double* aux_grad_test_ref,
			   // AUX QUANTITIES OF INTEREST
			   double* dLow,
			   double* quantDOFs,
			   double beta,
			   double epsilon)
    {
      double Ct_sge = 4.0;	  
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      //eN is the element index
      //eN_k is the quadrature point index for a scalar
      //eN_k_nSpace is the quadrature point index for a vector
      //eN_i is the element test function index
      //eN_j is the element trial function index
      //eN_k_j is the quadrature point index for a trial function
      //eN_k_i is the quadrature point index for a trial function

      register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
      for (int i=0; i<NNZ; i++)
	{
	  TransportMatrix[i] = 0.;
	  TransposeTransportMatrix[i] = 0.;
	}
      
      register double lumped_mass_matrix[numDOFs];
      for (int i=0; i<numDOFs; i++)
	lumped_mass_matrix[i] = 0;
      
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double elementResidual_u[nDOF_test_element];
	  register double  elementTransport[nDOF_test_element][nDOF_trial_element];
	  register double  elementTransposeTransport[nDOF_test_element][nDOF_trial_element];	  
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  elementTransport[i][j]=0.0;
		  elementTransposeTransport[i][j]=0.0;
		}	      
	    }//i
	  
	  
	  ///////////////////////////////////////////////////////
	  // COMPUTE BLENDING SHAPE FUNCTIONS FOR CURRENT CELL //
	  ///////////////////////////////////////////////////////
	  register double //alpha[nQuadraturePoints_element],
	    blended_test_ref[nQuadraturePoints_element*nDOF_trial_element],
	    blended_grad_test_ref[nQuadraturePoints_element*nDOF_trial_element*nSpace];
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      register int eN_k = eN*nQuadraturePoints_element+k, eN_k_nSpace = eN_k*nSpace;
	      int counter = k*nDOF_trial_element*nSpace;
	      int eN_nDOF_trial_element = eN*nDOF_trial_element;
	      // COMPUTE ALPHA AND ITS GRADIENT AT CURRENT QUAD POINT // 
	      double
		alpha=0, grad_alpha[nSpace], grad_alpha_referenceElement[nSpace],
		aux_grad_phi[nDOF_trial_element*nSpace],
		jac[nSpace*nSpace],jacDet,jacInv[nSpace*nSpace],x,y,z;
	      // compute jacInv
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,jacDet,jacInv,
					  x,y,z);
	      // COMPUTE ALPHA VIA LINEAR SPACE //
	      ck.gradTrialFromRef(&aux_grad_test_ref[k*nDOF_trial_element*nSpace],
				  jacInv,aux_grad_phi); //assume trial=test space	      
	      ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
			    &aux_test_ref[k*nDOF_trial_element],alpha);

	      // COMPUTE ALPHA VIA HIGH ORDER SPACE //
	      //ck.gradTrialFromRef(&u_grad_test_ref[k*nDOF_trial_element*nSpace],
	      //		  jacInv,aux_grad_phi); //assume trial=test space	    
	      //ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
	      //	    &u_test_ref[k*nDOF_trial_element],alpha);

	      // compute grad_alpha from DOFs
	      ck.gradFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
			     aux_grad_phi,grad_alpha); //This grad_alpha is NOT at the ref element
	      Mult(jac,grad_alpha,grad_alpha_referenceElement);
	      
	      // Loop in local shape functions
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  double phiH_i = u_test_ref[k*nDOF_trial_element+i];
		  double phiL_i = aux_test_ref[k*nDOF_trial_element+i];
		  blended_test_ref[k*nDOF_trial_element+i] = alpha*phiH_i + (1.0-alpha)*phiL_i;

		  for (int I=0;I<nSpace;I++)
		    {
		      double grad_phiH_i = u_grad_test_ref[counter];
		      double grad_phiL_i = aux_grad_test_ref[counter];
		      
		      blended_grad_test_ref[counter] = 
			alpha*grad_phiH_i + (1-alpha)*grad_phiL_i
			+ grad_alpha_referenceElement[I]*(phiH_i-phiL_i);
		      // update counter
		      counter++;
		    }
		}
	    } //k

	  ////////////////////////////////////////////////////////////////
	  // END OF COMPUTING BLENDING SHAPE FUNCTIONS FOR CURRENT CELL //
	  ////////////////////////////////////////////////////////////////
	  
	  //loop over quadrature points and compute integrands
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indeces and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double
		forceViaInterpolation=0.0, grad_forceViaInterpolation[nSpace],
		u=0.0,un=0.0,
		grad_u[nSpace],grad_un[nSpace],
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		u_test_dV[nDOF_trial_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		// FOR BLENDED SPACES
		// aux shape functions (linear on quad mesh)
		//aux_grad_phi[nDOF_trial_element*nSpace],
		//aux_phi_dV[nDOF_trial_element],
		//aux_grad_phi_dV[nDOF_test_element*nSpace],
		// end of aux shape functions
		// blended FE Space
		//blended_grad_phi[nDOF_trial_element*nSpace],
		//blended_phi_dV[nDOF_trial_element],
		//blended_grad_phi_dV[nDOF_test_element*nSpace],
		// end of blended FE space
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
	      
	      //std::cout << "k=" << k << "\t"
	      //	<< x << "," << y << std::endl;
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      //get the trial(or test) function gradients
	      ck.gradTrialFromRef(&blended_grad_test_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      //get the solution
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],
			    &blended_test_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],
			    &blended_test_ref[k*nDOF_trial_element],un);	      

	      // NORMAL SPACES
	      //ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
	      //		  jacInv,u_grad_trial);
	      //ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],
	      //	    &u_trial_ref[k*nDOF_trial_element],u);

	      // compute gradient of solution 
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],
			     u_grad_trial,grad_u);
	      ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],
			     u_grad_trial,grad_un);	      
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  //u_test_dV[j] = aux_test_ref[k*nDOF_trial_element+j]*dV;
		  u_test_dV[j] = blended_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		    }
		}
	      //save solution
	      q_u[eN_k] = u;
	      q_m[eN_k] = u;
	      // save grad
	      for(int I=0;I<nSpace;I++)
		q_grad_u[eN_k_nSpace+I]=grad_u[I];
	      if (q_dV_last[eN_k] <= -100)
                q_dV_last[eN_k] = dV;
              q_dV[eN_k] = dV;
	      //
	      //update element residual 
	      //
	      double velocityBeta[2];
	      velocityBeta[0] = 1.0;
	      velocityBeta[1] = 3.0;

	      double grad_uexact[2];
	      grad_uexact[0] = gradx_uexact[eN_k];
	      grad_uexact[1] = grady_uexact[eN_k];
	      for(int i=0;i<nDOF_test_element;i++) 
		{
		  register int i_nSpace=i*nSpace;
		  // poisson equation //
		  //elementResidual_u[i] +=
		  ////  u*u_test_dV[i] // poisson like equation
		  //ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
		  //- force[eN_k]*u_test_dV[i];
		  
		  // Ritz projection //
		  //elementResidual_u[i] +=
		  //ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
		  //-ck.NumericalDiffusion(1.0,grad_uexact,&u_grad_test_dV[i_nSpace]);

		  ////////////////////////////////
		  // STEADY ADVECTION-DIFFUSION //
		  ////////////////////////////////
		  elementResidual_u[i] +=
		    // ADVECTION
		    ck.NumericalDiffusion(1.0, 
					  velocityBeta,
					  grad_u)*u_test_dV[i]
		    // DISSIPATION
		    +epsilon*ck.NumericalDiffusion(1.0, 
		  				   grad_u,
		  				   &u_grad_test_dV[i_nSpace]);

		  /////////////////////////////////
		  // EXPLICIT TRANSPORT EQUATION //
		  /////////////////////////////////
		  //elementResidual_u[i] +=
		  // TIME DERIVATIVE 
		  //(u-un)/dt*u_test_dV[i]
		  // ADVECTION
		  //  ck.NumericalDiffusion(1.0, 
		  //			  velocityBeta,
		  //			  grad_un)*u_test_dV[i];
		  // j-th LOOP // To construct transport matrices

		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      elementTransport[i][j] += // int[(vel.grad_wj)*wi*dx]
			ck.Advection_strong(velocityBeta,
					    &u_grad_test_dV[j_nSpace])
			*blended_test_ref[k*nDOF_trial_element+i];
		      
		      elementTransposeTransport[i][j] += // int[(vel.grad_wi)*wj*dx]
			ck.Advection_strong(velocityBeta,
					    &u_grad_test_dV[i_nSpace])
			*blended_test_ref[k*nDOF_trial_element+j];
		    }
		  /////////////////////////
		}//i
	    }
	  //abort();
	  //
	  //load element into global residual and save element residual
	  //
	  
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	      // distribute transport matrices
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  TransportMatrix[csrRowIndeces_CellLoops[eN_i] +
				  csrColumnOffsets_CellLoops[eN_i_j]]
		    += elementTransport[i][j];
		  TransposeTransportMatrix[csrRowIndeces_CellLoops[eN_i] +
					   csrColumnOffsets_CellLoops[eN_i_j]]
		    += elementTransposeTransport[i][j];
		}//j	      
	    }//i
	}//elements

      int ij=0;    
      for (int i=0; i<numDOFs; i++)
	{
	  double solni = u_dof[i]; // solution at time tn for the ith DOF
	  //double solni = u_dof_old[i]; // solution at time tn for the ith DOF
	  double ith_dissipative_term = 0;
	  double ith_flux_term = 0;
	  double dLii = 0.;
	  int ii=0;
	  double alphai = alpha_dof[i];
	  
	  // loop over the sparsity pattern of the i-th DOF
	  for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
	    {
	      int j = csrColumnOffsets_DofLoops[offset];
	      double solnj = u_dof[j]; // solution at time tn for the jth DOF
	      //double solnj = u_dof_old[j]; // solution at time tn for the jth DOF
	      double alphaj = alpha_dof[j];
	      double dij;
	      
	      ith_flux_term += TransportMatrix[ij]*solnj;
	      
	      if (i != j)
		{
		  dij = fmax(fabs((1-alphai)*TransportMatrix[ij]),
		  	     fabs((1-alphaj)*TransposeTransportMatrix[ij]));
		  //dij = fmax(fabs(TransportMatrix[ij]),
		  //	     fabs(TransposeTransportMatrix[ij]));
		  dLow[ij] = dij;
		  //dLij = fmax(0.,fmax(psi[i]*TransportMatrix[ij], // Approach by S. Badia
		  //		  psi[j]*TransposeTransportMatrix[ij]));
		  //dissipative terms
		  ith_dissipative_term += dij*(solnj-solni);
		  dLii -= dij;
		}
	      else
		{
		  ii=ij; // save index for diagonal entry to update matrix later
		}
	      //update ij
	      ij+=1;
	    }
	  dLow[ii]=dLii;
	  // compute edge_based_cfl
	  // update residual
	  //globalResidual[i] += dt*(ith_flux_term - ith_dissipative_term);
	  globalResidual[i] += -ith_dissipative_term;
	}
      // NO LOOP IN BOUNDARIES //
    }
        
    void calculateJacobian_blending_spaces(//element
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
			   int lag_shockCapturing,/*mwf not used yet*/
			   double shockCapturingDiffusion,
			   //VRANS
			   const double* q_porosity,
			   //
			   int* u_l2g,
			   double* elementDiameter,
			   int degree_polynomial,
			   double* u_dof, 
			   double* velocity,
			   double* q_m_betaBDF, 
			   double* cfl,
			   double* q_numDiff_u_last, 
			   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			   int numDOFs,
			   int* csrRowIndeces_DofLoops,
			   int* csrColumnOffsets_DofLoops,
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
			   int LUMPED_MASS_MATRIX,
			   // FOR BLENDING SPACES
			   double* alpha_value,
			   double* alpha_dof,
			   double* aux_test_ref,
			   double* aux_grad_test_ref,
			   double* dLow,
			   double beta,
			   double epsilon)
    {
      //std::cout<<"ndjaco  address "<<q_numDiff_u_last<<std::endl;
      double Ct_sge = 4.0;
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //

      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_u_u[i][j]=0.0;
	      }
	  ///////////////////////////////////////////////////////
	  // COMPUTE BLENDING SHAPE FUNCTIONS FOR CURRENT CELL //
	  ///////////////////////////////////////////////////////
	  register double //alpha[nQuadraturePoints_element],
	    blended_test_ref[nQuadraturePoints_element*nDOF_trial_element],
	    blended_grad_test_ref[nQuadraturePoints_element*nDOF_trial_element*nSpace];
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      register int eN_k = eN*nQuadraturePoints_element+k, eN_k_nSpace = eN_k*nSpace;
	      int counter = k*nDOF_trial_element*nSpace;
	      int eN_nDOF_trial_element = eN*nDOF_trial_element;
	      // COMPUTE ALPHA AND ITS GRADIENT AT CURRENT QUAD POINT // 
	      double
		alpha=0, grad_alpha[nSpace], grad_alpha_referenceElement[nSpace],
		aux_grad_phi[nDOF_trial_element*nSpace],
		jac[nSpace*nSpace],jacDet,jacInv[nSpace*nSpace],x,y,z;
	      // compute jacInv
	      ck.calculateMapping_element(eN,
					  k,
					  mesh_dof,
					  mesh_l2g,
					  mesh_trial_ref,
					  mesh_grad_trial_ref,
					  jac,jacDet,jacInv,
					  x,y,z);
	      // COMPUTE ALPHA VIA LINEAR SPACE //
	      ck.gradTrialFromRef(&aux_grad_test_ref[k*nDOF_trial_element*nSpace],
				  jacInv,aux_grad_phi); //assume trial=test space
	      ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
			    &aux_test_ref[k*nDOF_trial_element],alpha);

	      // COMPUTE ALPHA VIA HIGH ORDER SPACE //
	      //ck.gradTrialFromRef(&u_grad_test_ref[k*nDOF_trial_element*nSpace],
	      //		  jacInv,aux_grad_phi); //assume trial=test space	    
	      //ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
	      //	    &u_test_ref[k*nDOF_trial_element],alpha);

	      // compute grad_alpha from DOFs
	      ck.gradFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
			     aux_grad_phi,grad_alpha); //This grad_alpha is NOT at the ref element
	      Mult(jac,grad_alpha,grad_alpha_referenceElement); // grad_alpha at the ref element
	      
	      // Loop in local shape functions
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  double phiH_i = u_test_ref[k*nDOF_trial_element+i];
		  double phiL_i = aux_test_ref[k*nDOF_trial_element+i];
		  blended_test_ref[k*nDOF_trial_element+i] = alpha*phiH_i + (1.0-alpha)*phiL_i;
 
		  for (int I=0;I<nSpace;I++)
		    {
		      double grad_phiH_i = u_grad_test_ref[counter];
		      double grad_phiL_i = aux_grad_test_ref[counter];
		      
		      blended_grad_test_ref[counter] =
			alpha*grad_phiH_i + (1-alpha)*grad_phiL_i
			+ grad_alpha_referenceElement[I]*(phiH_i-phiL_i);
		      // update counter
		      counter++;
		    }
		}
	    }
	  ////////////////////////////////////////////////////////////////
	  // END OF COMPUTING BLENDING SHAPE FUNCTIONS FOR CURRENT CELL //
	  ////////////////////////////////////////////////////////////////
	  
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
		u_grad_trial[nDOF_trial_element*nSpace],
		dV,
		u_test_dV[nDOF_test_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
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
	      // use blended spaces
	      ck.gradTrialFromRef(&blended_grad_test_ref[k*nDOF_trial_element*nSpace],
				  jacInv,u_grad_trial);
	      //get the trial function gradients
	      //ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],
	      //		  jacInv,u_grad_trial);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  u_test_dV[j] = blended_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
		}
	      double velocityBeta[2];
	      velocityBeta[0] = 1.0;
	      velocityBeta[1] = 3.0;
	      
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  for(int j=0;j<nDOF_trial_element;j++)
		    {
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      ////////////////////////////////
		      // poisson like: u-Delta u = f//
		      ////////////////////////////////
		      //elementJacobian_u_u[i][j] +=
		      ////blended_test_ref[k*nDOF_trial_element+j]*u_test_dV[i]
		      //ck.NumericalDiffusion(1.0,
		      //		      &u_grad_trial[j_nSpace],
		      //		      &u_grad_test_dV[i_nSpace]);		    
		      
		      ////////////////////////////////
		      // STEADY ADVECTION-DIFFUSION //
		      ////////////////////////////////
		      elementJacobian_u_u[i][j] +=
			ck.NumericalDiffusion(1.0,
					      velocityBeta,
					      &u_grad_trial[j_nSpace])*u_test_dV[i]
			+epsilon*ck.NumericalDiffusion(1.0,
						       &u_grad_trial[j_nSpace],
						       &u_grad_test_dV[i_nSpace]);
		      
		      ////////////////////////
		      // EXPLICIT TRANSPORT //
		      ////////////////////////
		      //elementJacobian_u_u[i][j] +=
		      //1./dt*blended_test_ref[k*nDOF_trial_element+j]*u_test_dV[i];
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

      int ij=0;
      for (int i=0; i<numDOFs; i++)
      {
	for (int offset=csrRowIndeces_DofLoops[i]; offset<csrRowIndeces_DofLoops[i+1]; offset++)
          { // First loop in j (sparsity pattern)
	    int j = csrColumnOffsets_DofLoops[offset];
	    globalJacobian[ij] -= dLow[ij];
	    //update ij
            ij++;
          }
      }
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
                                 int useMetrics,
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

      void getLumpedL2Projection(//element
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
				 double* q_alpha,
				 int offset_u, int stride_u,
				 // PARAMETERS FOR EDGE VISCOSITY
				 int numDOFs,
				 double* lumpedL2Projection)
      {
        register double
          lumped_mass_matrix[numDOFs];
        for (int i=0; i<numDOFs; i++)
          {
            lumpedL2Projection[i]=0.;
            lumped_mass_matrix[i]=0.;
          }
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for local contributions and initialize
            register double
              element_lumped_mass_matrix[nDOF_test_element],	      
              element_lumpedL2Projection[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                element_lumped_mass_matrix[i]=0.0;
                element_lumpedL2Projection[i]=0.0;
              }
            //loop over quadrature points and compute integrands
            for(int k=0;k<nQuadraturePoints_element;k++)
              {
                //compute indeces and declare local storage
                register int eN_k = eN*nQuadraturePoints_element+k,
                  eN_k_nSpace = eN_k*nSpace,
                  eN_nDOF_trial_element = eN*nDOF_trial_element;
                register double
                  alpha,
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

		alpha=q_alpha[eN_k];
                for(int i=0;i<nDOF_test_element;i++)
                  {
                    element_lumped_mass_matrix[i] += u_test_dV[i];
                    element_lumpedL2Projection[i] += alpha*u_test_dV[i];
                  }
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g[eN_i]; //global i-th index

                lumped_mass_matrix[gi] += element_lumped_mass_matrix[i];
                lumpedL2Projection[gi] += element_lumpedL2Projection[i];
              }//i
          }//elements
        // COMPUTE LUMPED L2 PROJECTION
        for (int i=0; i<numDOFs; i++)
          {
            double mi = lumped_mass_matrix[i];
            lumpedL2Projection[i] /= mi;
          }
      }      
  };//VOF

  inline VOF_base* newVOF(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<VOF_base,VOF,CompKernel>(nSpaceIn,
										 nQuadraturePoints_elementIn,
										 nDOF_mesh_trial_elementIn,
										 nDOF_trial_elementIn,
										 nDOF_test_elementIn,
										 nQuadraturePoints_elementBoundaryIn,
										 CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<VOF_base,VOF,CompKernel>(nSpaceIn,
									       nQuadraturePoints_elementIn,
									       nDOF_mesh_trial_elementIn,
									       nDOF_trial_elementIn,
									       nDOF_test_elementIn,
									       nQuadraturePoints_elementBoundaryIn,
									       CompKernelFlag);
  }
}//proteus
#endif
