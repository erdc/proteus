#ifndef MCorr_H
#define MCorr_H
#include <cmath>
#include <iostream>

#define MCorr_SPACE_DIM 3
#if (MCorr_SPACE_DIM == 3)
 #define nSpace 3
 #define nQuadraturePoints_element 5
 #define nDOF_trial_element 4
 #define nDOF_test_element 4
 #define nDOF_test_X_trial_element 16
 #define nQuadraturePoints_elementBoundary 4
 #define MCORR_NAME MCorr3D
 #define MCORR_RES calculateResidual_MCorr3D
 #define MCORR_JAC calculateJacobian_MCorr3D
#else
 #if (MCorr_SPACE_DIM == 2)
  #define nSpace 2
  #define nQuadraturePoints_element 6
  #define nDOF_trial_element 3
  #define nDOF_test_element 3
  #define nDOF_test_X_trial_element 9
  #define nQuadraturePoints_elementBoundary 4
  #define MCORR_NAME MCorr2D
  #define MCORR_RES calculateResidual_MCorr2D
  #define MCORR_JAC calculateJacobian_MCorr2D
 #else 
  #define nSpace 1
  #define nQuadraturePoints_element 3
  #define nDOF_trial_element 2
  #define nDOF_test_element 2
  #define nDOF_test_X_trial_element 4
  #define nQuadraturePoints_elementBoundary 1
  #define MCORR_NAME MCorr1D
  #define MCORR_RES calculateResidual_MCorr1D
  #define MCORR_JAC calculateJacobian_MCorr1D
 #endif
#endif

namespace MCORR_NAME // MCorr
{
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

inline double smoothedHeaviside_integral(double eps, double phi)
{
  double HI;
  if (phi > eps)
    {
      HI= phi - eps + 	0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  else if (phi < -eps)
    {
      HI=0.0;
    }
  else
    {
      HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  return HI;
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
void evaluateCoefficients_c(const double& epsHeaviside,
			    const double& epsDirac,
			    const double& phi,
			    const double& H,
			    const double& u,
			    double& r,
			    double& dr)
{
  r = smoothedHeaviside(epsHeaviside,u + phi) - H;
  dr = smoothedDirac(epsDirac,u + phi);
}

  
inline
void evaluateCoefficients(int K,
			  const double& epsHeaviside,
			  const double& epsDirac,
			  const double *phi,
			  const double *H,
			  const double *u,
			  double* r,
			  double* dr)
{
  evaluateCoefficients_c(epsHeaviside,
			 epsDirac,
                         phi[K],
			 H[K],
			 u[K],
                         r[K],
                         dr[K]);
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

}//MCorr namespace
extern "C"
{
  //void calculateResidual_MCorr(int nElements_global,
  void MCORR_RES (int nElements_global,
			 double epsHeaviside,
			 double epsDirac,
			 double epsDiffusion,
			 int* u_l2g, 
			 double* elementDiameter,
			 double* u_dof,
			 double* u_trial, 
			 double* u_grad_trial, 
			 double* u_test_dV, 
			 double* u_grad_test_dV, 
			 double* q_phi,
			 double* q_H,
			 double* q_u,
			 double* q_r,
			 int offset_u, int stride_u, 
			 double* globalResidual);
  //void calculateJacobian_MCorr(int nElements_global,
  void MCORR_JAC(int nElements_global,
			 double epsHeaviside,
			 double epsDirac,
			 double epsDiffusion,
			 int* u_l2g,
			 double* elementDiameter,
			 double* u_dof, 
			 double* u_trial, 
			 double* u_grad_trial, 
			 double* u_test_dV, 
			 double* u_grad_test_dV, 
			 double* q_phi,
			 double* q_H,
			 int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			 double* globalJacobian);
}//extern "C"
#endif
