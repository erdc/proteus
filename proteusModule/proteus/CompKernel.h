#ifndef COMPKERNEL_H
#define COMPKERNEL_H
#include <cmath>

template<const int NSPACE>
class CompKernel
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
	HI= phi - eps							\
	  + 0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) \
	  - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
      }
    else if (phi < -eps)
      {
	HI=0.0;
      }
    else
      {
	HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) \
	  - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
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

  inline double valFromDOF(const double& dof,const double& v)
  {
    return dof*v;
  }

  inline double gradFromDOF(const double& dof,const double& grad_v_I)
  {
    return dof*grad_v_I;
  }
  
  inline void backwardEuler(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
  {  
    mt =(m-m_old)/dt;
    dmt = dm/dt;
  }

  inline void bdf(const double& alpha, const double& beta, const double& m, const double& dm, double& mt, double& dmt)
  {  
    mt =alpha*m + beta;
    dmt = alpha*dm;
  }

  inline double Mass_weak(const double& mt, const double& w_dV)
  {
    return mt*w_dV;
  }

  inline double MassJacobian_weak(const double& dmt,
				  const double& v,
				  const double& w_dV)
  {
    return dmt*v*w_dV;
  }

  inline double Mass_strong(const double& mt)
  {
    return mt; 
  }

  inline double MassJacobian_strong(const double& dmt,
				    const double& v)
  {
    return dmt*v;
  }

  inline double Mass_adjoint(const double& dmt,
			     const double& w_dV)
  {
    return dmt*w_dV;
  }

  inline double Advection_weak(const double  f[NSPACE],
			       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= f[I]*grad_w_dV[I];
    return tmp;
  }

  inline double AdvectionJacobian_weak(const double df[NSPACE],
				       const double& v,
				       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= df[I]*v*grad_w_dV[I];
    return tmp;
  }

  inline double Advection_strong(const double df[NSPACE],
				 const double grad_u[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += df[I]*grad_u[I];
    return tmp;
  }

  inline double AdvectionJacobian_strong(const double df[NSPACE],
					 const double grad_v[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += df[I]*grad_v[I];
    return tmp;
  }

  inline double Advection_adjoint(const double df[NSPACE],
				  const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= df[I]*grad_w_dV[I];
    return tmp;
  }

  inline double Hamiltonian_weak(const double& H,
				 const double& w_dV)
  {
    return H*w_dV;
  }

  inline double HamiltonianJacobian_weak(const double dH[NSPACE],
					 const double grad_v[NSPACE],
					 const double& w_dV)
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_v[I]*w_dV;
    return tmp;
  }

  inline double Hamiltonian_strong(const double dH[NSPACE],
				   const double grad_u[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_u[I];
    return tmp;
  }

  inline double HamiltonianJacobian_strong(const double dH[NSPACE],
					   const double grad_v[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp += dH[I]*grad_v[I];
    return tmp;
  }

  inline double Hamiltonian_adjoint(const double dH[NSPACE],
				    const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      tmp -= dH[I]*grad_w_dV[I];
    return tmp;
  }

  inline double Diffusion_weak(int* rowptr,
			       int* colind,
			       double* a,
			       const double grad_phi[NSPACE],
			       const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for(int I=0;I<NSPACE;I++)
      for (int m=rowptr[I];m<rowptr[I+1];m++)
	tmp += a[m]*grad_phi[colind[m]]*grad_w_dV[I];
    return tmp;
  }

  inline double DiffusionJacobian_weak(int* rowptr,
				       int* colind,
				       double* a,
				       double* da,
				       const double grad_phi[NSPACE],
				       const double grad_w_dV[NSPACE],
				       const double& dphi,
				       const double& v,
				       const double grad_v[NSPACE])
  {
    double daProduct=0.0,dphiProduct=0.0;
    for (int I=0;I<NSPACE;I++)
      for (int m=rowptr[I];m<rowptr[I+1];m++)
	{
	  daProduct += da[m]*grad_phi[colind[m]]*grad_w_dV[I];
	  dphiProduct += a[m]*grad_v[colind[m]]*grad_w_dV[I];
	}
    return daProduct*v+dphiProduct*dphi;
  }

  inline double Reaction_weak(const double& r,
			      const double& w_dV)
  {
    return r*w_dV;
  }
  
  inline double ReactionJacobian_weak(const double& dr,
				      const double& v,
				      const double& w_dV)
  {
    return dr*v*w_dV;
  }

  inline double Reaction_strong(const double& r)
  {
    return r; 
  }

  inline double ReactionJacobian_strong(const double& dr,
					const double& v)
  {
    return dr*v;
  }

  inline double Reaction_adjoint(const double& dr,
				 const double& w_dV)
  {
    return dr*w_dV;
  }

  inline void calculateNumericalDiffusion(const double& shockCapturingDiffusion,
					  const double& elementDiameter,
					  const double& strong_residual,
					  const double grad_u[NSPACE],
					  double& numDiff)
  {
    double h,
      num,
      den,
      n_grad_u;
    h = elementDiameter;
    n_grad_u = 0.0;
    for (int I=0;I<NSPACE;I++)
      n_grad_u += grad_u[I]*grad_u[I];
    num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
    den = sqrt(n_grad_u) + 1.0e-8;
    numDiff = num/den;
  }

  inline double SubgridError(const double& error,
			     const double& Lstar_w_dV)
  {
    return -error*Lstar_w_dV;
  }

  inline double SubgridErrorJacobian(const double& derror,
				     const double& Lstar_w_dV)
  {
    return -derror*Lstar_w_dV;
  }

  inline double NumericalDiffusion(const double& numDiff,
				   const double grad_u[NSPACE],
				   const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for (int I=0;I<NSPACE;I++)
      tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
    return tmp;
  }

  inline double NumericalDiffusionJacobian(const double& numDiff,
					   const double grad_v[NSPACE],
					   const double grad_w_dV[NSPACE])
  {
    double tmp=0.0;
    for (int I=0;I<NSPACE;I++)
      tmp += numDiff*grad_v[I]*grad_w_dV[I];
    return tmp;
  }

 
  inline void exteriorNumericalDiffusiveFlux(const double& eps,
					     const double& phi,
					     int* rowptr,
					     int* colind,
					     const int& isDOFBoundary,
					     const int& isFluxBoundary,
					     const double n[NSPACE],
					     double* bc_a,
					     const double& bc_u,
					     const double& bc_flux,
					     double* a,
					     const double grad_phi[NSPACE],
					     const double& u,
					     const double& penalty,
					     double& flux)
  {
    double diffusiveVelocityComponent_I,penaltyFlux,max_a;
    if(isDOFBoundary == 1)
      {
	flux = 0.0;
	max_a=0.0;
	for(int I=0;I<NSPACE;I++)
	  {
	    diffusiveVelocityComponent_I=0.0;
	    for(int m=rowptr[I];m<rowptr[I+1];m++)
	      {
		diffusiveVelocityComponent_I -= a[m]*grad_phi[colind[m]];
		max_a = fmax(max_a,a[m]);
	      }
	    flux+= diffusiveVelocityComponent_I*n[I];
	  }
	penaltyFlux = penalty*(u-bc_u);
	flux += penaltyFlux;
	//contact line slip
	flux*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
      }
    else if(isFluxBoundary == 1)
      {
	flux = bc_flux;
      }
    else
      {
	std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
	flux = 0.0;
      }
  }

  inline double ExteriorElementBoundaryFlux(const double& flux,
					    const double& w_dS)
  {
    return flux*w_dS;
  }

  inline double ExteriorNumericalAdvectiveFluxJacobian(const double& dflux_left,
						       const double& v)
  {
    return dflux_left*v;
  }

  inline double ExteriorNumericalDiffusiveFluxJacobian(const double& eps,
						       const double& phi,
						       int* rowptr,
						       int* colind,
						       const int& isDOFBoundary,
						       const double n[NSPACE],
						       double* a,
						       const double& v,
						       const double grad_v[NSPACE],
						       const double& penalty)
  {
    double dvel_I,tmp=0.0;
    if(isDOFBoundary >= 1)
      {
	for(int I=0;I<NSPACE;I++)
	  {
	    dvel_I=0.0;
	    for(int m=rowptr[I];m<rowptr[I+1];m++)
	      {
		dvel_I -= a[m]*grad_v[colind[m]];
	      }
	    tmp += dvel_I*n[I];
	  }
	tmp +=penalty*v;
	//contact line slip
	tmp*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
      }
    return tmp;
  }
};
#endif
