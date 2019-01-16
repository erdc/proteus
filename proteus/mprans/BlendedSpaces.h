#ifndef BlendedSpaces_H
#define BlendedSpaces_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

#define ALPHA_VIA_LINEAR_SPACE 1
#define SIGN(z) (z == 0 ? 0 : (z>0 ? 1 : -1) )
#define USE_MACRO_CELL 1

namespace proteus
{
  class BlendedSpaces_base
  {
    //The base class defining the interface
  public:
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
				   double* dLow,
				   double epsilon,
				   // Type of problem to solve
				   int PROBLEM_TYPE,
				   int GALERKIN_SOLUTION,
				   double* galerkin_solution,
				   double* gamma_dof,
                                   // AUX QUANTITIES OF INTEREST
                                   double* quantDOFs,
				   double* quantDOFs2
				   )=0;
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
				   int GALERKIN_SOLUTION,
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
    virtual void calculateSmoothnessIndicator(//element
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
					      //physics
					      int nElements_global,
					      int* u_l2g,
					      int* r_l2g,
					      double* u_dof,
					      int offset_u, int stride_u,
					      double* is_dof_external,
					      double* is_dof_internal,
					      double* num_hi,
					      double* den_hi,
					      double* global_hi,
					      double* element_He,
					      // for loops in DOFs
					      double* gamma_dof,
					      int numDOFs,
					      int NNZ,
					      int* rowptr,
					      int* colind
					      )=0;
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
			     double* dLow,
			     double epsilon,
			     // Type of problem to solve
			     int PROBLEM_TYPE,
			     int GALERKIN_SOLUTION,
			     double* galerkin_solution,
			     double* gamma_dof,
			     // AUX QUANTITIES OF INTEREST
			     double* quantDOFs,
			     double* quantDOFs2)
      {
	register double sigma[numDOFs];
	double maxSigma = 0;
	
	if (GALERKIN_SOLUTION==1)
	  {
	    for (int i=0; i<numDOFs; i++)	      
	      alpha_dof[i] = 1.0;
	  }
	else
	  {
	    register double D2uG[numDOFs];
	    // first loop on DOFs: compute estimate of second order derivative
	    for (int i=0; i<numDOFs; i++)
	      {
		double uGi = galerkin_solution[i];
		D2uG[i] = 0;		
		// loop over the sparsity pattern of the i-th DOF
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  {
		    int j = colind[offset];
		    double uGj = galerkin_solution[j];
		    D2uG[i] += uGi-uGj;		    
		  }		
	      }
	    double maxSigma=0;
	    // second loop
	    for (int i=0; i<numDOFs; i++)
	      {
		int argmin = i;
		double sign;
		int same_sign = 1;
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  {
		    int j = colind[offset];
		    argmin = fabs(D2uG[j]) < fabs(D2uG[argmin]) ? j : argmin;
		    if (offset==rowptr[i])
		      {
			sign = SIGN(D2uG[j]);
		      }
		    else
		      {
			same_sign = sign == SIGN(D2uG[j]) ? 1 : 0;
			sign = SIGN(D2uG[j]);
		      }
		  }
		if (same_sign == 0)
		  sigma[i] = 0;
		else
		  sigma[i]=fabs(D2uG[argmin]);
		
		maxSigma = fmax(sigma[i],maxSigma);	       	    
		//quantDOFs[i] = sigma[i];
		
	      }
	    
	    // third loop to compute alpha
	    for (int i=0; i<numDOFs; i++)
	      {
		double c1=0.2;
		double c2=0.3;

		//double m = -1./((c2-c1)*maxSigma);
		//double b = c2/(c2-c1);
		//if (sigma[i] <= c1*maxSigma)
		//alpha_dof[i] = 1.;
		//else if (sigma[i] >= c2*maxSigma)
		//alpha_dof[i] = 0.;
		//else
		//alpha_dof[i] = m*sigma[i]+b;

		
		//double m = 1./((c2-c1)*maxSigma);
		//double b = -c1/(c2-c1);
		//if (sigma[i] <= c1*maxSigma)
		//alpha_dof[i] = 0.;
		//else if (sigma[i] >= c2*maxSigma)
		//alpha_dof[i] = 1.;
		//else
		//alpha_dof[i] = m*sigma[i]+b;

		//quantDOFs2[i] = galerkin_solution[i] < -1 ? 1.0 : 0.0;
		alpha_dof[i] = galerkin_solution[i] < -1 ? 0.0 : 1.0;
		quantDOFs[i] = alpha_dof[i];
		
		alpha_dof[i] = fmax(alpha_dof[i],gamma_dof[i]);
		//alpha_dof[i] = gamma_dof[i];
		quantDOFs2[i] = alpha_dof[i];		
	      }	    
	  }
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

	register double TransportMatrix[NNZ], TransposeTransportMatrix[NNZ];
	for (int i=0; i<NNZ; i++)
	  {
	    TransportMatrix[i] = 0.;
	    TransposeTransportMatrix[i] = 0.;
	  }

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

	    ///////////////////////////////////////////////////////////////////
	    // ***** COMPUTE BLENDING SHAPE FUNCTIONS FOR CURRENT CELL ***** //
	    //////////////////////////////////////////////////////////////////
	    register double
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
		if (ALPHA_VIA_LINEAR_SPACE==1)
		  {
		    // COMPUTE ALPHA VIA LINEAR SPACE //
		    ck.gradTrialFromRef(&aux_grad_test_ref[k*nDOF_trial_element*nSpace],
					jacInv,aux_grad_phi); //assume trial=test space
		    ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
				  &aux_test_ref[k*nDOF_trial_element],alpha);
		  }
		else
		  {
		    // COMPUTE ALPHA VIA HIGH ORDER SPACE //
		    ck.gradTrialFromRef(&u_grad_test_ref[k*nDOF_trial_element*nSpace],
					jacInv,aux_grad_phi); //assume trial=test space
		    ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
				  &u_test_ref[k*nDOF_trial_element],alpha);
		  }
		// compute grad_alpha from DOFs
		ck.gradFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
			       aux_grad_phi,grad_alpha); //This grad_alpha is NOT at the ref element
		Mult(jac,grad_alpha,grad_alpha_referenceElement); // send grad alpha to ref element

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
	    ////////////////////////////////////////////////////////////////////////////
	    // ***** END OF COMPUTING BLENDING SHAPE FUNCTIONS FOR CURRENT CELL ***** //
	    ////////////////////////////////////////////////////////////////////////////

	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  u=0.0,grad_u[nSpace],D_times_grad_u[nSpace],DTensor[nSpace*nSpace],
		  D_times_u_grad_trial_j[nSpace], D_times_u_grad_trial_i[nSpace],
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
		ck.gradTrialFromRef(&blended_grad_test_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &blended_test_ref[k*nDOF_trial_element],
			      u);
		//get the solution gradients
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_u);
		//multiply D times grad_u
		if (PROBLEM_TYPE==2)
		  {
		    calculateDTensor(DTensor);
		    Mult(DTensor,grad_u,D_times_grad_u);
		  }
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = blended_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;

		  }

		//save solution
		q_u[eN_k] = u;
		// save grad
		for(int I=0;I<nSpace;I++)
		  q_grad_u[eN_k_nSpace+I]=grad_u[I];

		// coefficient of steady advection - diffusion
		double velocityBeta[2];
		velocityBeta[0] = velocity[eN_k_nSpace+0];
		velocityBeta[1] = velocity[eN_k_nSpace+1];

		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;

		    if (PROBLEM_TYPE==0)
		      {
			// poisson equation //
			elementResidual_u[i] +=
			  ck.NumericalDiffusion(1.0,grad_u,&u_grad_test_dV[i_nSpace])
			  - force[eN_k]*u_test_dV[i];
		      }
		    else if (PROBLEM_TYPE==1)
		      {
			// STEADY ADVECTION-DIFFUSION //
			elementResidual_u[i] +=
			  // ADVECTION
			  ck.NumericalDiffusion(1.0,velocityBeta,grad_u)*u_test_dV[i]
			  // DISSIPATION
			  +epsilon*ck.NumericalDiffusion(1.0,
							 grad_u,
							 &u_grad_test_dV[i_nSpace]);
			// j-th LOOP // To construct transport matrices
			for(int j=0;j<nDOF_trial_element;j++)
			  {
			    int j_nSpace = j*nSpace;
			    elementTransport[i][j] += // int[(vel.grad_wj)*wi*dx]
			      ck.Advection_strong(velocityBeta,
						  &u_grad_test_dV[j_nSpace])
			      *blended_test_ref[k*nDOF_trial_element+i];

			    elementTransposeTransport[i][j] += // int[(vel.grad_wi)*wj*dx]
			      ck.Advection_strong(velocityBeta,
						  &u_grad_test_dV[i_nSpace])
			      *blended_test_ref[k*nDOF_trial_element+j];
			  }
		      }
		    else // PROBLEM_TYPE==2
		      {
			// poisson equation //
			elementResidual_u[i] += (ck.NumericalDiffusion(1.0,
								       D_times_grad_u,
								       &u_grad_test_dV[i_nSpace])
						 - force[eN_k]*u_test_dV[i]);
			Mult(DTensor,&u_grad_trial[i_nSpace],D_times_u_grad_trial_i);
			// j-th LOOP // To construct transport matrices
			for(int j=0;j<nDOF_trial_element;j++)
			  {
			    int j_nSpace = j*nSpace;
			    Mult(DTensor,&u_grad_trial[j_nSpace],D_times_u_grad_trial_j);
			    elementTransport[i][j] += // int[(D*grad_wj).grad_wi*dx]
			      ck.NumericalDiffusion(1.0,
						    D_times_u_grad_trial_j,
						    &u_grad_test_dV[i_nSpace]);
			    elementTransposeTransport[i][j] += // int[(D*grad_wi).grad_wj*dx]
			      ck.NumericalDiffusion(1.0,
						    D_times_u_grad_trial_i,
						    &u_grad_test_dV[j_nSpace]);
			  }
		      }
		  }//i
	      }
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;
		globalResidual[offset_u+stride_u*r_l2g[eN_i]] += elementResidual_u[i];
		if (PROBLEM_TYPE==1 || PROBLEM_TYPE==2)
		  {
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
		      }
		  }//j
	      }//i
	  }//elements

	if (GALERKIN_SOLUTION==0)
	if (PROBLEM_TYPE==1 or PROBLEM_TYPE==2)
	  {
	    // LOOP in DOFs //
	    int ij=0;
	    for (int i=0; i<numDOFs; i++)
	      {
		double solni = u_dof[i]; // solution at time tn for the ith DOF
		double ith_dissipative_term = 0;
		double dLii = 0.;
		int ii=0;
		double alphai = alpha_dof[i];

		// loop over the sparsity pattern of the i-th DOF
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  {
		    int j = colind[offset];
		    double solnj = u_dof[j]; // solution at time tn for the jth DOF
		    double alphaj = alpha_dof[j];
		    double dij;

		    if (i != j)
		      {
			if (PROBLEM_TYPE==1)
			  dij = fmax(fabs((1-alphai)*TransportMatrix[ij]),
				     fabs((1-alphaj)*TransposeTransportMatrix[ij]));
			else
			  dij = fmax(0,fmax((1-alphai)*TransportMatrix[ij],
					    (1-alphaj)*TransposeTransportMatrix[ij]));
			dLow[ij] = dij;
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
		dLow[ii] = dLii;
		// update residual
		globalResidual[i] += -ith_dissipative_term;
	      }
	  }

	// NO LOOP IN BOUNDARIES //
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
			     int GALERKIN_SOLUTION,
			     // Type of problem to solve
			     int PROBLEM_TYPE)
      {
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
		if (ALPHA_VIA_LINEAR_SPACE==1)
		  {
		    // COMPUTE ALPHA VIA LINEAR SPACE //
		    ck.gradTrialFromRef(&aux_grad_test_ref[k*nDOF_trial_element*nSpace],
					jacInv,aux_grad_phi); //assume trial=test space
		    ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
				  &aux_test_ref[k*nDOF_trial_element],alpha);
		  }
		else
		  {
		    // COMPUTE ALPHA VIA HIGH ORDER SPACE //
		    ck.gradTrialFromRef(&u_grad_test_ref[k*nDOF_trial_element*nSpace],
					jacInv,aux_grad_phi); //assume trial=test space
		    ck.valFromDOF(alpha_dof,&u_l2g[eN_nDOF_trial_element],
				  &u_test_ref[k*nDOF_trial_element],alpha);
		  }
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
		  D_times_u_grad_trial_j[nSpace], DTensor[nSpace*nSpace],
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
		//get the trial function gradients using the blended space
		ck.gradTrialFromRef(&blended_grad_test_ref[k*nDOF_trial_element*nSpace],
				    jacInv,
				    u_grad_trial);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  {
		    u_test_dV[j] = blended_test_ref[k*nDOF_trial_element+j]*dV;
		    for (int I=0;I<nSpace;I++)
		      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
		  }
		// get DTensor
		if (PROBLEM_TYPE==2)
		  calculateDTensor(DTensor);

		// coefficient of steady advection - diffusion
		double velocityBeta[2];
		velocityBeta[0] = velocity[eN_k_nSpace+0];
		velocityBeta[1] = velocity[eN_k_nSpace+1];

		for(int i=0;i<nDOF_test_element;i++)
		  {
		    for(int j=0;j<nDOF_trial_element;j++)
		      {

			int j_nSpace = j*nSpace;
			int i_nSpace = i*nSpace;

			if (PROBLEM_TYPE==0)
			  {
			    // poisson equation //
			    elementJacobian_u_u[i][j] +=
			      ck.NumericalDiffusion(1.0,
						    &u_grad_trial[j_nSpace],
						    &u_grad_test_dV[i_nSpace]);
			  }
			else if (PROBLEM_TYPE==1)
			  {
			    // STEADY ADVECTION-DIFFUSION //
			    elementJacobian_u_u[i][j] +=
			      ck.NumericalDiffusion(1.0,
						    velocityBeta,
						    &u_grad_trial[j_nSpace])*u_test_dV[i]
			      +epsilon*ck.NumericalDiffusion(1.0,
							     &u_grad_trial[j_nSpace],
							     &u_grad_test_dV[i_nSpace]);
			  }
			else // PROBLEM_TYPE==2
			  {
			    Mult(DTensor,&u_grad_trial[j_nSpace],D_times_u_grad_trial_j);
			    // poisson equation //
			    elementJacobian_u_u[i][j] +=
			      ck.NumericalDiffusion(1.0,
						    D_times_u_grad_trial_j,
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
		    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
		  }//j
	      }//i
	  }//elements

	if (GALERKIN_SOLUTION==0)
	if (PROBLEM_TYPE == 1 or PROBLEM_TYPE==2)
	  {
	    int ij=0;
	    for (int i=0; i<numDOFs; i++)
	      {
		for (int offset=rowptr[i]; offset<rowptr[i+1]; offset++)
		  { // First loop in j (sparsity pattern)
		    int j = colind[offset];
		    globalJacobian[ij] -= dLow[ij];
		    //update ij
		    ij++;
		  }
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

      void calculateSmoothnessIndicator(//element
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
					//physics
					int nElements_global,
					int* u_l2g,
					int* r_l2g,
					double* u_dof,
					int offset_u, int stride_u,
					double* is_dof_external,
					double* is_dof_internal,
					double* num_hi,
					double* den_hi,
					double* global_hi,
					double* element_He,
					// for loops in DOFs
					double* gamma_dof,
					int numDOFs,
					int NNZ,
					int* rowptr,
					int* colind)
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

	register int first_adjacent_dof_to_middle_dof[numDOFs];
	register int second_adjacent_dof_to_middle_dof[numDOFs];

	register double lumped_mass_matrix[numDOFs];
	for (int i=0; i<numDOFs; i++)
	  lumped_mass_matrix[i]=0;
	
	int nSpace2 = nSpace*nSpace;
	for(int eN=0;eN<nElements_global;eN++)
	  {
	    //declare local storage for element residual and initialize
	    register double element_mass_matrix[nDOF_test_element], element_hi[nDOF_test_element];
	    double det_hess_Ke=0, area_Ke=0;
	    double hess_u0=0, hess_u1=0, hess_u2=0, hess_u3=0;
	    for (int i=0;i<nDOF_test_element;i++)
	      {
	  	element_mass_matrix[i]=0.0;
	   	element_hi[i]=0.0;
	      }
	    //loop over quadrature points and compute integrands
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace,
		  eN_nDOF_trial_element = eN*nDOF_trial_element;
		register double
		  u=0.0,grad_u[nSpace],hess_u[nSpace2],
		  det_hess_u=0.,
		  jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  u_grad_trial[nDOF_trial_element*nSpace],
		  u_hess_trial[nDOF_trial_element*nSpace2],
		  u_test_dV[nDOF_trial_element],
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
				    u_grad_trial);
		ck.hessTrialFromRef(&u_hess_trial_ref[k*nDOF_trial_element*nSpace2],
				    jacInv,
				    u_hess_trial);
		//get the solution based on the blended functions
		ck.valFromDOF(u_dof,
			      &u_l2g[eN_nDOF_trial_element],
			      &u_trial_ref[k*nDOF_trial_element],
			      u);
		//get the solution gradients
		ck.gradFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_grad_trial,
			       grad_u);
		ck.hessFromDOF(u_dof,
			       &u_l2g[eN_nDOF_trial_element],
			       u_hess_trial,
			       hess_u);
		//precalculate test function products with integration weights
		for (int j=0;j<nDOF_trial_element;j++)
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;

		// compute the determinan of the hessian
		det_hess_u = fabs(hess_u[0]*hess_u[3] - hess_u[2]*hess_u[1]);
		
		hess_u0 += hess_u[0]*dV;
		hess_u1 += hess_u[1]*dV;
		hess_u2 += hess_u[2]*dV;
		hess_u3 += hess_u[3]*dV;

		/*
		double type=3;
		hess_u0=0;
		hess_u1=0;
		hess_u2=0;
		hess_u3=0;
		if (type==1)
		  {
		    hess_u0 += 0*dV;
		    hess_u1 += 0*dV;
		    hess_u2 += 0*dV;
		    hess_u3 += 0*dV;
		  }
		else if (type==2)
		  {
		    hess_u0 += 2*dV;
		    hess_u1 += 0*dV;
		    hess_u2 += 0*dV;
		    hess_u3 += 2*dV;
		  }
		else
		  {
		    double pi = M_PI;
		    hess_u0 += -4*pi*pi*sin(2*pi*x)*sin(2*pi*y)*dV;
		    hess_u1 += 4*pi*pi*cos(2*pi*x)*cos(2*pi*y)*dV;     
		    hess_u2 += 4*pi*pi*cos(2*pi*x)*cos(2*pi*y)*dV;     
		    hess_u3 += -4*pi*pi*sin(2*pi*x)*sin(2*pi*y)*dV;
		  }
		*/
		// to compute average of det of hessian in cell Ke
		det_hess_Ke += det_hess_u*dV;
		area_Ke += dV;

		for(int i=0;i<nDOF_test_element;i++)
		  {
		    register int i_nSpace=i*nSpace;
		    element_mass_matrix[i] += u_test_dV[i];
		    element_hi[i] += det_hess_u*u_test_dV[i];
		  }//i
	      }
	    // compute average of det of hessian in cell Ke
	    element_He[eN] = det_hess_Ke/area_Ke;

	    hess_u0 /= area_Ke;
	    hess_u1 /= area_Ke;
	    hess_u2 /= area_Ke;
	    hess_u3 /= area_Ke;

	    //element_He[eN] = fabs(hess_u0*hess_u3 - hess_u2*hess_u1);
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {		
		register int eN_i=eN*nDOF_test_element+i;
		int gi = offset_u+stride_u*r_l2g[eN_i];
		lumped_mass_matrix[offset_u+stride_u*r_l2g[eN_i]] += element_mass_matrix[i];
		global_hi[offset_u+stride_u*r_l2g[eN_i]] += element_hi[i];
		
		num_hi[offset_u+stride_u*r_l2g[eN_i]] += element_He[eN];
		den_hi[offset_u+stride_u*r_l2g[eN_i]] += 1;

		
		if (i<4)
		  is_dof_external[gi] = 1;
		else if (i==8)
		  is_dof_internal[gi] = 1;
		else //if (i>=4 && i<8) // middle dof
		  {
		    int gj1 = 0, gj2 = 0;
		    if (i==4)
		      {
			gj1 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+0];
			gj2 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+1];
		      }
		    else if (i==5)
		      {
			gj1 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+1];
			gj2 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+3];
		      }
		    else if (i==6)
		      {
			gj1 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+3];
			gj2 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+2];
		      }
		    else 
		      {
			gj1 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+2];
			gj2 = offset_u+stride_u*r_l2g[eN*nDOF_test_element+0];
		      }		    

		    first_adjacent_dof_to_middle_dof[gi]  = gj1;
		    second_adjacent_dof_to_middle_dof[gi] = gj2;
		  }
	      }//i
	  }//elements
	// finish the computation of global_hi
	for (int i=0; i<numDOFs; i++)
	  global_hi[i] = num_hi[i]/den_hi[i];
	
	register double min_hiHe[numDOFs];
	for (int i=0; i<numDOFs; i++)
	  {
	    min_hiHe[i] = 1E100;
	  }
	// Loop in cells to compute min(hi*He)
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

	if (USE_MACRO_CELL==1)
	  {
	    // compute gamma in DOFs owned by big Q1 mesh
	    for (int i=0; i<numDOFs; i++)
	      {
		if (is_dof_external[i] == 1) // if dof is in big Q1 mesh
		  {
		    double eps = 1E-10;
		    double C = 3.0;
		    double hi2 = global_hi[i] * global_hi[i];
		    gamma_dof[i] = fmax(0, fmin(hi2, C*min_hiHe[i]))/(hi2+eps);
		    //gamma_dof[i] = (fmax(0, fmin(hi2, C*min_hiHe[i])) + eps)/(hi2+eps);
		    //gamma_dof[i] = hi2 == 0 ? 1. : fmax(0, fmin(hi2, C*min_hiHe[i]))/hi2;
		  }
		else
		  {
		    gamma_dof[i] = 0;
		  }
	      }
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
		double C = 3.0;
		double hi2 = global_hi[i] * global_hi[i];
		gamma_dof[i] = fmax(0, fmin(hi2, C*min_hiHe[i]))/(hi2+eps);
		//gamma_dof[i] = (fmax(0, fmin(hi2, C*min_hiHe[i])) + eps)/(hi2+eps);
		//gamma_dof[i] = hi2 == 0 ? 1. : fmax(0, fmin(hi2, C*min_hiHe[i]))/hi2;
	      }
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
