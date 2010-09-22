#ifndef RDLSQ_H
#define RDLSQ_H
#include <cmath>
#include <iostream>

#define RDLSQ_SPACE_DIM 3
#if (RDLSQ_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 15
 #define nDOF_trial_element 10
 #define nDOF_test_element 10
 #define nDOF_test_X_trial_element 100
 #define nQuadraturePoints_elementBoundary 7
 #define RDLSQ_NAME RDLSQ3D
 #define RDLSQ_RES calculateResidual_RDLSQ3D
 #define RDLSQ_JAC calculateJacobian_RDLSQ3D
#else
 #if (RDLSQ_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 7
  #define nDOF_trial_element 6
  #define nDOF_test_element 6
  #define nDOF_test_X_trial_element 36
  #define nQuadraturePoints_elementBoundary 5
  #define RDLSQ_NAME RDLSQ2D
  #define RDLSQ_RES calculateResidual_RDLSQ2D
  #define RDLSQ_JAC calculateJacobian_RDLSQ2D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 1
  #define RDLSQ_NAME RDLSQ1D
  #define RDLSQ_RES calculateResidual_RDLSQ1D
  #define RDLSQ_JAC calculateJacobian_RDLSQ1D
 #endif
#endif
namespace RDLSQ_NAME
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
double valFromDOF_c(const double& dof,const double& v)
{
  return dof*v;
}

inline
void valFromDOF(int eN, int k, int j, double* u,double* dof,int* l2g,double* v)
{
  u[eN*nQuadraturePoints_element+
    k]
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+
			    j]],
		    v[eN*nQuadraturePoints_element*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
void valFromDOF_ext(int ebNE, int eN, int k, int j, double* u, double* dof, int* l2g, double* v)
{
  u[ebNE*nQuadraturePoints_elementBoundary+
    k] 
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		    v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
double gradFromDOF_c(const double& dof,const double& grad_v_I)
{
  return dof*grad_v_I;
}


inline
void gradFromDOF(int eN, int k, int j, double* grad_u,double* dof,int* l2g,double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[eN*nQuadraturePoints_element*nSpace+
	   k*nSpace+
	   I] += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+
				       j]],
			       grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
				      k*nDOF_trial_element*nSpace+
				      j*nSpace+
				      I]);
}

inline
void gradFromDOF_ext(int ebNE, int eN, int k, int j, double* grad_u, double* dof, int* l2g, double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[ebNE*nQuadraturePoints_elementBoundary*nSpace+
	   k*nSpace+I]
      += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		       grad_v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element*nSpace+
			      k*nDOF_trial_element*nSpace+
			      j*nSpace+I]);
}

inline
void evaluateCoefficients_c(const double& eps,
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
void evaluateCoefficients(int K,
			  const double& eps,
			  const double *u_levelSet,
			  const double *u,
			  const double *grad_u,
			  double* m,
			  double* dm,
			  double* H,
			  double* dH,
			  double* r)
{
  evaluateCoefficients_c(eps,
			 u_levelSet[K],
                         u[K],
                         &grad_u[K*nSpace],
                         m[K],
                         dm[K],
                         H[K],
                         &dH[K*nSpace],
			 r[K]);
}
  
inline
void backwardEuler_c(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =(m-m_old)/dt;
  dmt = dm/dt;
}

inline
void backwardEuler(int eN, int k, double dt, double* m_old, double* m, double* dm, double* mt, double* dmt)
{  
  backwardEuler_c(dt,
		  m_old[eN*nQuadraturePoints_element+
			k],
		  m[eN*nQuadraturePoints_element+
		    k],
		  dm[eN*nQuadraturePoints_element+
		     k],
		  mt[eN*nQuadraturePoints_element+
		     k], 
		  dmt[eN*nQuadraturePoints_element+
		      k]);
}

inline
void bdf_c(const double& alpha, const double& beta, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =alpha*m + beta;
  dmt = alpha*dm;
}
inline
void bdf(int eN, int k, double alpha, double* beta, double* m, double* dm, double* mt, double* dmt)
{  
  bdf_c(alpha,
	beta[eN*nQuadraturePoints_element+
	     k],
	m[eN*nQuadraturePoints_element+
	  k],
	dm[eN*nQuadraturePoints_element+
	   k],
	mt[eN*nQuadraturePoints_element+
	   k], 
	dmt[eN*nQuadraturePoints_element+
	    k]);
}


inline
double Mass_weak_c(const double& mt,
		   const double& w_dV)
{
  return mt*w_dV;
}

inline
void updateMass_weak(int eN,
                     int k,
                     int i,
                     double* mt,
                     double* w_dV,
                     double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i]  += Mass_weak_c(mt[eN*nQuadraturePoints_element + 
				      k],
				   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					k*nDOF_test_element + 
					i]);
}

inline
double MassJacobian_weak_c(const double& dmt,
			   const double& v,
			   const double& w_dV)
{
  return dmt*v*w_dV;
}

inline
void updateMassJacobian_weak(int eN,
			     int k,
			     int j,
			     int i,
			     double* dmt,
			     double* v,
			     double* w_dV,
			     double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += MassJacobian_weak_c(dmt[eN*nQuadraturePoints_element + 
						       k],
						   v[eN*nQuadraturePoints_element*nDOF_trial_element + 
						     k*nDOF_trial_element + 
						     j],
						   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							k*nDOF_test_element + 
							i]);
}

inline
double Mass_strong_c(const double& mt)
{
  return mt; 
}

inline
void updateMass_strong(int eN,
		       int k,
		       double* mt,
		       double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Mass_strong_c(mt[eN*nQuadraturePoints_element+
					 k]);
}

inline
double MassJacobian_strong_c(const double& dmt,
			     const double& v)
{
  return dmt*v;
}

inline
void updateMassJacobian_strong(int eN,
                               int k,
                               int j,
                               double* dmt,
                               double* v,
                               double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += MassJacobian_strong_c(dmt[eN*nQuadraturePoints_element+
						   k],
					       v[eN*nQuadraturePoints_element*nDOF_trial_element+
						 k*nDOF_trial_element + 
						 j]);
}

inline
double Mass_adjoint_c(const double& dmt,
		      const double& w_dV)
{
  return dmt*w_dV;
}

inline
void updateMass_adjoint(int eN,
			int k,
			int i,
			double* dmt,
			double* w_dV,
			double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Mass_adjoint_c(dmt[eN*nQuadraturePoints_element + 
				      k],
				  w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
				       k*nDOF_test_element + 
				       i]);
}

inline
double Hamiltonian_weak_c(const double& H,
			  const double& w_dV)
{
  return H*w_dV;
}

inline
void updateHamiltonian_weak(int eN,
                            int k,
                            int i,
                            double* H,
                            double* w_dV,
                            double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += Hamiltonian_weak_c(H[eN*nQuadraturePoints_element + 
								  k],
								w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
								     k*nDOF_test_element + 
								     i]);
}

inline
double HamiltonianJacobian_weak_c(const double dH[nSpace],
				  const double grad_v[nSpace],
				  const double& w_dV)
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_v[I]*w_dV;
  return tmp;
}

inline
void updateHamiltonianJacobian_weak(int eN,
				    int k,
				    int j,
				    int i,
				    double* dH,
				    double* grad_v,
				    double* w_dV,
				    double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element +
			 j]+= HamiltonianJacobian_weak_c(&dH[eN*nQuadraturePoints_element*nSpace + 
							     k*nSpace],
							 &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace + 
								 k*nDOF_trial_element*nSpace +
								 j*nSpace],
							 w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							      k*nDOF_test_element +
							      i]);
}

inline
double Hamiltonian_strong_c(const double dH[nSpace],
			    const double grad_u[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_u[I];
  return tmp;
}

inline
void updateHamiltonian_strong(int eN,
			      int k,
			      double* dH,
			      double* grad_u,
			      double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+k] += Hamiltonian_strong_c(&dH[eN*nQuadraturePoints_element*nSpace+
									      k*nSpace],
									  &grad_u[eN*nQuadraturePoints_element*nSpace+
										  k*nSpace]);
}

inline
double HamiltonianJacobian_strong_c(const double dH[nSpace],
				    const double grad_v[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_v[I];
  return tmp;
}

inline
void updateHamiltonianJacobian_strong(int eN,
                                      int k,
                                      int j,
                                      double* dH,
                                      double* grad_v,
                                      double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += HamiltonianJacobian_strong_c(&dH[eN*nQuadraturePoints_element*nSpace+
							  k*nSpace],
						      &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
							      k*nDOF_trial_element*nSpace + 
							      j*nSpace]);
}

inline
double Hamiltonian_adjoint_c(const double dH[nSpace],
			     const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= dH[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateHamiltonian_adjoint(int eN,
			       int k,
			       int i,
			       double* dH,
			       double* grad_w_dV,
			       double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Hamiltonian_adjoint_c(&dH[eN*nQuadraturePoints_element*nSpace + 
					     k*nSpace],
					 &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						    k*nDOF_test_element*nSpace + 
						    i*nSpace]);
}
inline
double Reaction_weak_c(const double& r,
		       const double& w_dV)
{
  return r*w_dV;
}

inline
void updateReaction_weak(int eN,
			 int k,
			 int i,
			 double* r,
			 double* w_dV,
			 double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i]  += Reaction_weak_c(r[eN*nQuadraturePoints_element + 
					 k],
				       w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					    k*nDOF_test_element + 
					    i]);
}

inline
double ReactionJacobian_weak_c(const double& dr,
			       const double& v,
			       const double& w_dV)
{
  return dr*v*w_dV;
}

inline
void updateReactionJacobian_weak(int eN,
				 int k,
				 int j,
				 int i,
				 double* dr,
				 double* v,
				 double* w_dV,
				 double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += ReactionJacobian_weak_c(dr[eN*nQuadraturePoints_element + 
							  k],
						       v[eN*nQuadraturePoints_element*nDOF_trial_element + 
							 k*nDOF_trial_element + 
							 j],
						       w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							    k*nDOF_test_element + 
							    i]);
}

inline
double Reaction_strong_c(const double& r)
{
  return r; 
}

inline
void updateReaction_strong(int eN,
			   int k,
			   double* r,
			   double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Reaction_strong_c(r[eN*nQuadraturePoints_element+
					    k]);
}

inline
double ReactionJacobian_strong_c(const double& dr,
				     const double& v)
{
  return dr*v;
}

inline
void updateReactionJacobian_strong(int eN,
				   int k,
				   int j,
				   double* dr,
				   double* v,
				   double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += ReactionJacobian_strong_c(dr[eN*nQuadraturePoints_element+
						      k],
						   v[eN*nQuadraturePoints_element*nDOF_trial_element+
						     k*nDOF_trial_element + 
						     j]);
}

inline
double Reaction_adjoint_c(const double& dr,
			  const double& w_dV)
{
  return dr*w_dV;
}

inline
void updateReaction_adjoint(int eN,
			    int k,
			    int i,
			    double* dr,
			    double* w_dV,
			    double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Reaction_adjoint_c(dr[eN*nQuadraturePoints_element + 
					 k],
				      w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					   k*nDOF_test_element + 
					   i]);
}

inline
void calculateSubgridError_tau_c(const double& elementDiameter,
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

inline
void calculateSubgridError_tau(int eN,
                               int k,
                               double* elementDiameter,
                               double* dmt,
                               double* dH,
                               double* cfl,
                               double* tau)
{
  calculateSubgridError_tau_c(elementDiameter[eN],
			      dmt[eN*nQuadraturePoints_element+k],
			      &dH[eN*nQuadraturePoints_element*nSpace+
				k*nSpace],
			      cfl[eN*nQuadraturePoints_element+k],
			      tau[eN*nQuadraturePoints_element+k]);
}

inline
void calculateSubgridError_tauRes_c(const double& tau,
				    const double& pdeResidual,
				    double& subgridError)
{
  subgridError = -tau*pdeResidual;
}

inline
void calculateSubgridErrorDerivatives_tauRes_c(const double& tau,
					       const double dpdeResidual_du[nDOF_trial_element],
					       double dsubgridError_du[nDOF_trial_element])
{
  for (int j=0;j<nDOF_trial_element;j++)
    {
      dsubgridError_du[j] = -tau*dpdeResidual_du[j];
     }
}

inline
void calculateSubgridError_tauRes(int eN,
                                  int k,
                                  double* tau,
                                  double* pdeResidual,
                                  double* dpdeResidual_du,
				  double* subgridError,
				  double* dsubgridError_du)
{
  calculateSubgridError_tauRes_c(tau[eN*nQuadraturePoints_element+k],
				 pdeResidual[eN*nQuadraturePoints_element+k],
				 subgridError[eN*nQuadraturePoints_element+k]);
  calculateSubgridErrorDerivatives_tauRes_c(tau[eN*nQuadraturePoints_element+k],
					    &dpdeResidual_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
					    &dsubgridError_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element]);
}

inline
void calculateNumericalDiffusion_c(const double& shockCapturingDiffusion,
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
void calculateNumericalDiffusion(int eN,
				 int k,
				 double shockCapturingDiffusion,
				 double* elementDiameter,
				 double* strong_residual,
				 double* grad_u,
				 double* numDiff)
{
  calculateNumericalDiffusion_c(shockCapturingDiffusion,
				elementDiameter[eN],
				strong_residual[eN*nQuadraturePoints_element+k],
				&grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace],
				numDiff[eN*nQuadraturePoints_element+k]);
}

inline
double SubgridError_c(const double& error,
		      const double& Lstar_w_dV)
{
  return -error*Lstar_w_dV;
}

inline
void updateSubgridError(int eN,
			int k,
			int i,
			double* error,
			double* Lstar_w_dV,
			double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] +=  SubgridError_c(error[eN*nQuadraturePoints_element + 
								   k],
							     Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
									k*nDOF_test_element + 
									i]);
}

inline
double SubgridErrorJacobian_c(const double& derror,
				  const double& Lstar_w_dV)
{
  return -derror*Lstar_w_dV;
}

inline
void updateSubgridErrorJacobian(int eN,
				int k,
				int j,
				int i,
				double* derror,
				double* Lstar_w_dV,
				double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += SubgridErrorJacobian_c(derror[eN*nQuadraturePoints_element*nDOF_trial_element+
							     k*nDOF_trial_element+
							     j],
						      Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
								 k*nDOF_test_element + 
								 i]);
}

inline
double NumericalDiffusion_c(const double& numDiff,
			    const double grad_u[nSpace],
			    const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateNumericalDiffusion(int eN,
			      int k,
			      int i,
			      double* numDiff,
			      double* grad_u,
			      double* grad_w_dV,
			      double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += NumericalDiffusion_c(numDiff[eN*nQuadraturePoints_element + 
									  k],
								  &grad_u[eN*nQuadraturePoints_element*nSpace + 
									  k*nSpace],
								  &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
									     k*nDOF_test_element*nSpace + 
									     i*nSpace]);
}

inline
double NumericalDiffusionJacobian_c(const double& numDiff,
					const double grad_v[nSpace],
					const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp += numDiff*grad_v[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateNumericalDiffusionJacobian(int eN,
				      int k,
				      int j,
				      int i,
				      double* numDiff,
				      double* grad_v,
				      double* grad_w_dV,
				      double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += NumericalDiffusionJacobian_c(numDiff[eN*nQuadraturePoints_element + 
								    k],
							    &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace + 
								    k*nDOF_trial_element*nSpace + 
								    j*nSpace],
							    &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
								       k*nDOF_test_element*nSpace + 
								       i*nSpace]);
}


inline
double ExteriorElementBoundaryFlux_c(const double& flux,
				     const double& w_dS)
{
  return flux*w_dS;
}

inline
void updateExteriorElementBoundaryFlux(int ebNE,
				       int eN_global,
				       int k,
                                       int i,
                                       double* flux,
                                       double* w_dS,
                                       double* residual)
{
  residual[eN_global*nDOF_test_element+
	   i] +=  ExteriorElementBoundaryFlux_c(flux[ebNE*nQuadraturePoints_elementBoundary+
						     k],
						w_dS[ebNE*nQuadraturePoints_elementBoundary*nDOF_test_element+
						     k*nDOF_test_element+
						     i]);
}

inline
double ExteriorNumericalAdvectiveFluxJacobian_c(const double& dflux_left,
						const double& v)
{
  return dflux_left*v;
}

inline
void updateExteriorNumericalAdvectiveFluxJacobian(int ebNE,
						  int k,
						  int j,
						  double* dflux_left,
						  double* v,
						  double* fluxJacobian)
{
  fluxJacobian[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j] 
    += 
    ExteriorNumericalAdvectiveFluxJacobian_c(dflux_left[ebNE*nQuadraturePoints_elementBoundary+k],
					     v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j]);
}

}//RDLSQ namespace

extern "C"
{
  //void calculateResidual_RDLSQ(int nElements_global,
  void RDLSQ_RES(int nElements_global,
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
  
  //void calculateJacobian_RDLSQ(int nElements_global,
  void RDLSQ_JAC(int nElements_global,
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
