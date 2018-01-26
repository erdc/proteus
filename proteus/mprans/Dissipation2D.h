#ifndef Dissipation2D_H
#define Dissipation2D_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
  class Dissipation2D_base
  {
    //The base class defining the interface
  public:
    virtual ~Dissipation2D_base(){}
    virtual void calculateResidual(//element
                                   double* mesh_trial_ref,
                                   double* mesh_grad_trial_ref,
                                   double* mesh_dof,
                                   double* meshVelocity_dof,
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
                                   //diffusion terms
                                   double nu_0,
                                   double nu_1,
                                   double sigma_e,
                                   double c_mu,
                                   double c_1,
                                   double c_2,
                                   double c_e,
                                   double rho_0,
                                   double rho_1,
                                   int dissipation_model_flag,
                                   //end diffusion
                                   double useMetrics,
                                   double alphaBDF,
                                   int lag_shockCapturing,
                                   double shockCapturingDiffusion,
                                   double sc_uref, double sc_alpha,
                                   int* u_l2g,
                                   double* elementDiameter,
                                   double* u_dof,double* u_dof_old,
                                   double* velocity,
                                   double* phi_ls, //level set variable
                                   double* q_kappa, //kinetic energy variable
                                   double* q_grad_kappa,
                                   double* q_porosity, //VRANS
                                   //velocity dof
                                   double * velocity_dof_u,
                                   double * velocity_dof_v,
                                   double * velocity_dof_w,
                                   //end velocity dof
                                   double* q_m,
                                   double* q_u,
                                   double* q_grad_u,
                                   double* q_m_betaBDF,
                                   double* cfl,
                                   double* q_numDiff_u,
                                   double* q_numDiff_u_last,
                                   double* ebqe_penalty_ext, //penalty
                                   int offset_u, int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isAdvectiveFluxBoundary_u,
                                   double* ebqe_bc_advectiveFlux_u_ext,
                                   int* isDiffusiveFluxBoundary_u,
                                   double* ebqe_bc_diffusiveFlux_u_ext,
                                   double* ebqe_phi,double epsFact,
                                   double* ebqe_kappa, //kinetic energy variable on boundary
                                   double* ebqe_porosity, //VRANS
                                   double* ebqe_u,
                                   double* ebqe_flux)=0;
    virtual void calculateJacobian(//element
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
                                   //diffusion
                                   double nu_0,
                                   double nu_1,
                                   double sigma_e,
                                   double c_mu,
                                   double c_1,
                                   double c_2,
                                   double c_e,
                                   double rho_0,
                                   double rho_1,
                                   int dissipation_model_flag,
                                   //end diffusion
                                   double useMetrics,
                                   double alphaBDF,
                                   int lag_shockCapturing,/*mwf not used yet*/
                                   double shockCapturingDiffusion,
                                   int* u_l2g,
                                   double* elementDiameter,
                                   double* u_dof, double* u_dof_old,
                                   double* velocity,
                                   double* phi_ls, //level set variable
                                   double* q_kappa, //kinetic energy
                                   double* q_grad_kappa,
                                   double* q_porosity,//VRANS
                                   //velocity dof
                                   double * velocity_dof_u,
                                   double * velocity_dof_v,
                                   double * velocity_dof_w,
                                   //end velocity dof
                                   double* q_m_betaBDF,
                                   double* cfl,
                                   double* q_numDiff_u_last,
                                   double* ebqe_penalty_ext, //penalty
                                   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                                   double* globalJacobian,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   double* ebqe_velocity_ext,
                                   int* isDOFBoundary_u,
                                   double* ebqe_bc_u_ext,
                                   int* isAdvectiveFluxBoundary_u,
                                   double* ebqe_bc_advectiveFlux_u_ext,
                                   int* isDiffusiveFluxBoundary_u,
                                   double* ebqe_bc_diffusiveFlux_u_ext,
                                   int* csrColumnOffsets_eb_u_u,
                                   double* ebqe_phi,double epsFact,
                                   double* ebqe_kappa,//kinetic energy on boundary
                                   double* ebqe_porosity)=0; //VRANS
  };

  template<class CompKernelType,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nQuadraturePoints_elementBoundary>
  class Dissipation2D : public Dissipation2D_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    Dissipation2D():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}

    inline
    void computeK_OmegaCoefficients(const double& div_eps,
                                    const double& k,
                                    const double& omega,
                                    const double grad_k[nSpace],
                                    const double grad_omega[nSpace],
                                    const double grad_vx[nSpace], //gradient of x component of velocity
                                    const double grad_vy[nSpace], //gradient of x component of velocity
                                    //const double grad_vz[nSpace], //gradient of x component of velocity
                                    double& inverse_sigma_k,
                                    double& inverse_sigma_omega,
                                    double& beta_star,
                                    double& beta,
                                    double& gamma)
    {
      //take these from NASA Langley Turbulence Model page
      //brute force just to see if I can figure it out
      //use inverse of sigma_k to match standard k-epsilon form
      inverse_sigma_k = 2.0; inverse_sigma_omega=2.0; gamma = 13.0/25.0;
      const double beta0_star = 0.09; const double beta0 = 9.0/125.0;
      double Omega[nSpace][nSpace] = {{0.,0.},
                                      {0.,0.}};

      double S[nSpace][nSpace] = {{0.,0.},
                                  {0.,0.}};

      //Omega_ij = (\pd{u_i}{x_j} - \pd{u_j}{x_i})/2
      Omega[0][1] = 0.5*(grad_vx[1]-grad_vy[0]);
      Omega[1][0] =-Omega[0][1];

      //S_ij = (\pd{u_i}{x_j} + \pd{u_j}{x_i})/2
      S[0][0] = grad_vx[0]; S[0][1] = 0.5*(grad_vx[1]+grad_vy[0]);
      S[1][0] = S[0][1];    S[1][1] = grad_vy[1];

      double chi_omega = 0.0;
      for (int i=0; i < nSpace; i++)
        for (int k=0; k < nSpace; k++)
          for (int j=0; j < nSpace; j++)
            chi_omega += Omega[i][j]*Omega[j][k]*S[k][i];
      if (fabs(omega) > div_eps)
        {
          chi_omega = fabs(chi_omega/(beta0_star*omega*beta0_star*omega*beta0_star*omega));

          const double f_beta = (1.0+70.0*chi_omega)/(1.0 + 80.0*chi_omega);
          beta = beta0*f_beta;
        }
      else
        {
          beta = beta0;
        }
      double chi_k = grad_k[0]*grad_omega[0] + grad_k[1]*grad_omega[1];
      double f_beta_star = 1.0;

      const double omega3 = omega*omega*omega;
      if (fabs(omega3) > div_eps)
        {
          chi_k = chi_k/omega3;
          f_beta_star = (1.0 + 680.0*chi_k*chi_k)/(1.0 + 400.0*chi_k*chi_k);
        }
      else if (chi_k > 0.0)
        f_beta_star = 680.0/400.0;

      beta_star = beta0_star*f_beta_star;
      //if (beta < 0.875*beta0 || beta > beta0)
      //        {
      //  std::cout<<"Kappa K-Omega coef problem k= "<<k<<" omega= "<<omega<<" beta= "<<beta<<" beta0= "<<beta0 <<" chi_omega= "<<chi_omega<<std::endl;
      //        }
      beta = fmax(0.875*beta0,fmin(beta,beta0));
      //if (beta_star < beta0_star || beta_star > (680.0+1.0e-4)/400.0*beta0_star)
      //{
      //  std::cout<<"Kappa K-Omega coef problem k= "<<k<<" omega= "<<omega<<" beta_star= "<<beta_star<<" beta0_star= "<<beta0_star <<" chi_k= "<<chi_k<<std::endl;
      //}
      beta_star = fmax(beta0_star,fmin(beta_star,(680.0/400.0)*beta0_star));
      //mwf hack
      //beta = beta0; beta_star = beta0_star;
    }

    //Try Lew, Buscaglia approximation
    inline
    void evaluateCoefficients(const double v[nSpace],
                              const double eps_mu,
                              const double phi,
                              const double nu_0,
                              const double nu_1,
                              const double sigma_e,
                              const double c_mu,
                              const double c_1,
                              const double c_2,
                              const double c_e,
                              const double grad_vx[nSpace], //gradient of x component of velocity
                              const double grad_vy[nSpace], //gradient of x component of velocity
                              //const double grad_vz[nSpace], //gradient of x component of velocity
                              const double& dissipation,
                              const double& dissipation_old,
                              const double& k,
                              const double& porosity,
                              int dissipation_model_flag,
                              const double grad_k[nSpace],
                              const double grad_dissipation_old[nSpace],
                              double& m,
                              double& dm,
                              double f[nSpace],
                              double df[nSpace],
                              double& a,
                              double& da_de,
                              double& r,
                              double& dr_de)
    {
      double nu_t=0.0,dnu_t_de=0.0,PiD4=0.0,disp=0.0,ddisp_de=0.0;
      double gamma_e=0.0,F_e=0.0, gamma_production=0.0,sigma_a=sigma_e,
        dgamma_e_d_dissipation=0.0, dF_e_d_dissipation=0.0;
      //either K-Epsilon or K-Omega
      const double isKEpsilon = (dissipation_model_flag>=2) ? 0.0 : 1.0;
      m = dissipation*porosity;
      dm = porosity;

      for (int I=0; I < nSpace; I++)
        {
          f[I] = v[I]*porosity*dissipation;
          df[I] = v[I]*porosity;
        }
      const double H_mu = smoothedHeaviside(eps_mu,phi);
      const double nu = (1.0-H_mu)*nu_0 + H_mu*nu_1;
      const double div_eps = 1.0e-2*fmin(nu_0,nu_1);
      //eddy viscosity
      nu_t     = isKEpsilon*c_mu*k*k/(fabs(dissipation_old)+div_eps)
        + (1.0-isKEpsilon)*k/(fabs(dissipation_old)+div_eps);

      dnu_t_de = 0.0;
      //if (nu_t > 1.e6*nu)
      //{
      //  std::cout<<"Dissipation2D WARNING isKEpsilon = "<<isKEpsilon<<" nu_t = " <<nu_t<<" nu= "<<nu<<" k= "<<k<<" dissipation= "<<dissipation<<std::endl;
      //}
      nu_t = fmax(nu_t,1.e-4*nu);
      //mwf hack
      nu_t     = fmin(nu_t,1.0e6*nu);

      //Production term
      PiD4 = 2.0*(grad_vx[0]*grad_vx[0] +
                  grad_vy[1]*grad_vy[1])
        +
        (grad_vx[1] + grad_vy[0])*(grad_vx[1] + grad_vy[0]);

      //K-Omega, 1998
      if (dissipation_model_flag==2)
        {
          //temporaries
          double sigma_k=1.,beta_star=1.,beta=1.0;
          computeK_OmegaCoefficients(div_eps,
                                     k,
                                     dissipation_old,
                                     grad_k,
                                     grad_dissipation_old,
                                     grad_vx,
                                     grad_vy,
                                     //grad_vz,
                                     sigma_k,
                                     sigma_a,
                                     beta_star,
                                     beta,
                                     gamma_production);
          //--full lagging of Gamma_e
          //dgamma_e_d_dissipation=0.0;
          //gamma_e=fmax(beta*dissipation_old,0.0);
          //--try to couple to k
          dgamma_e_d_dissipation = 0.0;
          gamma_e = fmax(beta*k/nu_t,0.0);
          //--quadratic nonlinearity
          //dgamma_e_d_dissipation = fmax(beta,0.0);
          //gamma_e = dgamma_e_d_dissipation*dissipation;

          //-- full lagging of production
          dF_e_d_dissipation=0.0;
          F_e = fmax(PiD4*gamma_production,0.0);
          //dF_e_d_dissipation = fmax(gamma_production*nu_t/(k+div_eps)*PiD4,0.0);
          //F_e = dF_e_d_dissipation*dissipation;
        }
      else if (dissipation_model_flag==3) //K-Omega 1988
        {
          sigma_a=2.0; //1/sigma_omega,
          gamma_production = 5.0/9.0;
          double beta = 3.0/40.0;

          dgamma_e_d_dissipation = 0.0;
          gamma_e = fmax(beta*k/nu_t,0.0);

          //-- full lagging of production
          dF_e_d_dissipation=0.0;
          F_e = fmax(PiD4*gamma_production,0.0);

        }
      else
        {
          //K-Epsilon
          gamma_e = fmax(c_2*dissipation_old/(k+div_eps),0.0);
          dgamma_e_d_dissipation = 0.0;
          F_e = fmax(c_1*k*PiD4,0.0);
          dF_e_d_dissipation=0.0;
          sigma_a = sigma_e;
        }

      a = porosity*(nu_t/sigma_a + nu);
      da_de = porosity*dnu_t_de/sigma_a;

      r = -porosity*F_e + porosity*gamma_e*dissipation;
      dr_de = -porosity*dF_e_d_dissipation + porosity*gamma_e + porosity*dgamma_e_d_dissipation;
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

      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv);
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
                                        const int& isAdvectiveFluxBoundary_u,
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
      else if (isAdvectiveFluxBoundary_u == 1)
        {
          flux = bc_flux_u;
          //std::cout<<"Flux boundary flux and flow"<<flux<<'\t'<<flow<<std::endl;
        }
      else
        {
          if (flow >= 0.0)
            {
              flux = u*flow;
            }
          else
            {
              ///std::cout<<"warning: Dissipation open boundary with no external trace, setting to zero for inflow flow = "<<flow <<" n= ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]"<<std::endl;
              flux = 0.0;
            }
          //std::cout<<"No BC boundary flux and flow "<<flux<<'\t'<<flow<<std::endl;

        }
      //flux = flow;
      //std::cout<<"flux error "<<flux-flow<<std::endl;
      //std::cout<<"flux in computationa"<<flux<<std::endl;
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

    inline
    void exteriorNumericalAdvectiveFluxDerivative(const int& isDOFBoundary_u,
                                                  const int& isAdvectiveFluxBoundary,
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
      else if (isAdvectiveFluxBoundary == 1)
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
    inline
    void exteriorNumericalDiffusiveFlux(const double& bc_flux,
                                        const int& isDOFBoundary,
                                        const int& isDiffusiveFluxBoundary,
                                        double n[nSpace],
                                        double bc_u,
                                        double a,
                                        double grad_psi[nSpace],
                                        double u,
                                        double penalty,
                                        double& flux)
    {
      double v_I;
      flux = 0.0;
      if (isDiffusiveFluxBoundary)
        {
          flux = bc_flux;
        }
      else if (isDOFBoundary)
        {
          flux = 0.0;
          for(int I=0;I<nSpace;I++)
            {
              v_I = -a*grad_psi[I];
              flux += v_I*n[I];
            }
          flux += penalty*(u-bc_u);
        }
      else
        {
          //std::cerr<<"warning, Dissipation2D diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
          flux = 0.0;
        }
    }
    inline
    void exteriorNumericalDiffusiveFluxDerivative(const int& isDOFBoundary,
                                                  const int& isDiffusiveFluxBoundary,
                                                  double n[nSpace],
                                                  double a,
                                                  double da,
                                                  double grad_psi[nSpace],
                                                  const double grad_v[nSpace],
                                                  double v,
                                                  double penalty,
                                                  double& fluxJacobian)
    {
      if (isDiffusiveFluxBoundary == 0 && isDOFBoundary == 1)
        {
          fluxJacobian = 0.0;
          for(int I=0;I<nSpace;I++)
            {
              fluxJacobian -= (a*grad_v[I] + da*v*grad_psi[I])*n[I];
            }
          fluxJacobian += penalty*v;
        }
      else
        fluxJacobian = 0.0;
    }

    void calculateResidual(//element
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
                           //diffusion terms
                           double nu_0,
                           double nu_1,
                           double sigma_e,
                           double c_mu,
                           double c_1,
                           double c_2,
                           double c_e,
                           double rho_0,
                           double rho_1,
                           int dissipation_model_flag,
                           //end diffusion
                           double useMetrics,
                           double alphaBDF,
                           int lag_shockCapturing, /*mwf not used yet*/
                           double shockCapturingDiffusion,
                           double sc_uref, double sc_alpha,
                           int* u_l2g,
                           double* elementDiameter,
                           double* u_dof,double* u_dof_old,
                           double* velocity,
                           double* phi_ls, //level set variable
                           double* q_kappa, //kinetic energy
                           double* q_grad_kappa,
                           double* q_porosity, //VRANS
                           //velocity dof
                           double * velocity_dof_u,
                           double * velocity_dof_v,
                           double * velocity_dof_w,
                           //end velocity dof
                           double* q_m,
                           double* q_u,
                           double* q_grad_u,
                           double* q_m_betaBDF,
                           double* cfl,
                           double* q_numDiff_u,
                           double* q_numDiff_u_last,
                           double* ebqe_penalty_ext, //penalty
                           int offset_u, int stride_u,
                           double* globalResidual,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           double* ebqe_velocity_ext,
                           int* isDOFBoundary_u,
                           double* ebqe_bc_u_ext,
                           int* isAdvectiveFluxBoundary_u,
                           double* ebqe_bc_advectiveFlux_u_ext,
                           int* isDiffusiveFluxBoundary_u,
                           double* ebqe_bc_diffusiveFlux_u_ext,
                           double* ebqe_phi,double epsFact,
                           double* ebqe_kappa, //kinetic energy on boundary
                           double* ebqe_porosity, //VRANS
                           double* ebqe_u,
                           double* ebqe_flux)
    {
      //cek should this be read in?
      double Ct_sge = 4.0;

      //loop over elements to compute volume integrals and load them into element and global residual
      //
      //eN is the element index
      //eN_k is the quadrature point index for a scalar
      //eN_k_nSpace is the quadrature point index for a vector
      //eN_i is the element test function index
      //eN_j is the element trial function index
      //eN_k_j is the quadrature point index for a trial function
      //eN_k_i is the quadrature point index for a trial function
      for(int eN=0;eN<nElements_global;eN++)
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
              register double u=0.0,u_old=0.0,
                grad_u[nSpace],grad_u_old[nSpace],
                m=0.0,dm=0.0,
                f[nSpace],df[nSpace],df_minus_da_grad_u[nSpace],
                m_t=0.0,dm_t=0.0,
                a=0.0,da=0.0,
                r=0.0,dr=0.0,
                grad_vx[nSpace],grad_vy[nSpace],//grad_vz[nSpace],
                pdeResidual_u=0.0,
                Lstar_u[nDOF_test_element],
                subgridError_u=0.0,
                tau=0.0,tau0=0.0,tau1=0.0,
                numDiff0=0.0,numDiff1=0.0,
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                u_grad_trial[nDOF_trial_element*nSpace],
                u_test_dV[nDOF_trial_element],
                u_grad_test_dV[nDOF_test_element*nSpace],
                dV,x,y,z,xt,yt,zt,
                G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv;
              // //
              // //compute solution and gradients at quadrature points
              // //
              // u=0.0;
              // for (int I=0;I<nSpace;I++)
              //   {
              //     grad_u[I]=0.0;
              //   }
              // for (int j=0;j<nDOF_trial_element;j++)
              //   {
              //     int eN_j=eN*nDOF_trial_element+j;
              //     int eN_k_j=eN_k*nDOF_trial_element+j;
              //     int eN_k_j_nSpace = eN_k_j*nSpace;
              //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
              //     for (int I=0;I<nSpace;I++)
              //       {
              //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
              //       }
              //   }
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
              ck.calculateG(jacInv,G,G_dd_G,tr_G);
              //get the trial function gradients
              ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
              //get the solution
              ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
              ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u_old);
              //get the solution gradients
              ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
              ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_old);
              //
              //compute velocity production terms, ***assumes same spaces for velocity dofs and Dissipation2D!***
              ck.gradFromDOF(velocity_dof_u,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vx);
              ck.gradFromDOF(velocity_dof_v,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vy);
              //ck.gradFromDOF(velocity_dof_w,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vz);
              //

              //
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                  for (int I=0;I<nSpace;I++)
                    {
                      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                    }
                }
              //
              //calculate pde coefficients at quadrature points
              //
              evaluateCoefficients(&velocity[eN_k_nSpace],
                                   epsFact,
                                   phi_ls[eN_k],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx,
                                   grad_vy,
                                   //grad_vz,
                                   u,
                                   u_old,
                                   q_kappa[eN_k],
                                   q_porosity[eN_k],
                                   dissipation_model_flag,
                                   &q_grad_kappa[eN_k_nSpace],
                                   grad_u_old,
                                   m,
                                   dm,
                                   f,
                                   df,
                                   a,
                                   da,
                                   r,
                                   dr);
              //
              //moving mesh
              //
              f[0] -= MOVING_DOMAIN*m*xt;
              f[1] -= MOVING_DOMAIN*m*yt;
              //f[2] -= MOVING_DOMAIN*m*zt;
              df[0] -= MOVING_DOMAIN*dm*xt;
              df[1] -= MOVING_DOMAIN*dm*yt;
              //df[2] -= MOVING_DOMAIN*dm*zt;

              //combine df and da/du \grad u term for stabilization and jacobian calculations
              df_minus_da_grad_u[0] = df[0] - da*grad_u[0];
              df_minus_da_grad_u[1] = df[1] - da*grad_u[1];
              //df_minus_da_grad_u[2] = df[2] - da*grad_u[2];
              //
              //calculate time derivative at quadrature points
              //
              ck.bdf(alphaBDF,
                     q_m_betaBDF[eN_k],
                     m,
                     dm,
                     m_t,
                     dm_t);
              //
              //calculate subgrid error (strong residual and adjoint)
              //
              //calculate strong residual
              pdeResidual_u = ck.Mass_strong(m_t) + ck.Advection_strong(df_minus_da_grad_u,grad_u)
                + ck.Reaction_strong(r);
              //calculate adjoint
              for (int i=0;i<nDOF_test_element;i++)
                {
                  // register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
                  // Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
                  register int i_nSpace = i*nSpace;
                  Lstar_u[i]  = ck.Advection_adjoint(df_minus_da_grad_u,&u_grad_test_dV[i_nSpace]) +
                    ck.Reaction_adjoint(dr,u_test_dV[i]);
                }
              //calculate tau and tau*Res
              calculateSubgridError_tau(elementDiameter[eN],dm_t + dr,df_minus_da_grad_u,cfl[eN_k],tau0);
              calculateSubgridError_tau(Ct_sge,
                                        G,
                                        dm_t + dr,
                                        df_minus_da_grad_u,
                                        tau1,
                                        cfl[eN_k]);

              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

              subgridError_u = -tau*pdeResidual_u;
              //
              //calculate shock capturing diffusion
              //


              ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u,grad_u,numDiff0);
              ck.calculateNumericalDiffusion(shockCapturingDiffusion,sc_uref, sc_alpha,G,G_dd_G,pdeResidual_u,grad_u,numDiff1);
              q_numDiff_u[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0;
              //std::cout<<tau<<"   "<<q_numDiff_u[eN_k]<<std::endl;
              //
              //update element residual
              //
              for(int i=0;i<nDOF_test_element;i++)
                {
                  register int eN_k_i=eN_k*nDOF_test_element+i,
                    eN_k_i_nSpace = eN_k_i*nSpace,
                    i_nSpace=i*nSpace;

                  elementResidual_u[i] += ck.Mass_weak(m_t,u_test_dV[i]) +
                    ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) +
                    ck.SubgridError(subgridError_u,Lstar_u[i]) +
                    ck.NumericalDiffusion(a,grad_u,&u_grad_test_dV[i_nSpace]) + //scalar diffusion so steal numericalDiffusion approximation
                    ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[i_nSpace]) +
                    ck.Reaction_weak(r,u_test_dV[i]);

                }//i
              //
              //cek/ido todo, get rid of m, since u=m
              //save momentum for time history and velocity for subgrid error
              //save solution for other models
              //
              q_u[eN_k] = u;
              q_m[eN_k] = m;
              for (int I=0; I < nSpace; I++)
                q_grad_u[eN_k_nSpace+I] = grad_u[I];
            }
          //
          //load element into global residual and save element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;

              globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
            }//i
        }//elements
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
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
              register double u_ext=0.0,u_old_ext=0.0,
                grad_u_ext[nSpace],grad_vx_ext[nSpace],grad_vy_ext[nSpace],//grad_vz_ext[nSpace],
                grad_u_old_ext[nSpace],grad_kappa_ext_dummy[nSpace],
                m_ext=0.0,
                dm_ext=0.0,
                f_ext[nSpace],
                df_ext[nSpace],
                a_ext=0.0,da_ext=0.0,
                r_ext=0.0,dr_ext=0.0,
                flux_ext=0.0,
                diffusive_flux_ext=0.0,
                bc_u_ext=0.0,
                bc_grad_u_ext[nSpace],
                bc_m_ext=0.0,
                bc_dm_ext=0.0,
                bc_f_ext[nSpace],
                bc_df_ext[nSpace],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,
                u_test_dS[nDOF_test_element],
                u_grad_trial_trace[nDOF_trial_element*nSpace],
                normal[3],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                G[nSpace*nSpace],G_dd_G,tr_G;
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
              dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
              //get the metric tensor
              //cek todo use symmetry
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
              //solution and gradients
              ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
              ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_old_ext);
              ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
              ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_old_ext);

              //mwf hack, skip on boundary for now
              grad_kappa_ext_dummy[0] = 0.0; grad_kappa_ext_dummy[1] = 0.0; //grad_kappa_ext_dummy[2] = 0.0;

              //
              //compute velocity production terms, ***assumes same spaces for velocity dofs and Dissipation2D!***
              ck.gradFromDOF(velocity_dof_u,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vx_ext);
              ck.gradFromDOF(velocity_dof_v,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vy_ext);
              //ck.gradFromDOF(velocity_dof_w,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vz_ext);
              //

              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                }
              //
              //load the boundary values
              //
              bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
              //
              //calculate the pde coefficients using the solution and the boundary values for the solution
              //
              evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
                                   epsFact,
                                   ebqe_phi[ebNE_kb],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx_ext,
                                   grad_vy_ext,
                                   //grad_vz_ext,
                                   u_ext,
                                   u_old_ext,
                                   ebqe_kappa[ebNE_kb],
                                   ebqe_porosity[ebNE_kb],
                                   dissipation_model_flag,
                                   grad_kappa_ext_dummy,
                                   grad_u_old_ext,
                                   m_ext,
                                   dm_ext,
                                   f_ext,
                                   df_ext,
                                   a_ext,
                                   da_ext,
                                   r_ext,
                                   dr_ext);
              evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
                                   epsFact,
                                   ebqe_phi[ebNE_kb],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx_ext,
                                   grad_vy_ext,
                                   //grad_vz_ext,
                                   bc_u_ext,
                                   bc_u_ext,
                                   ebqe_kappa[ebNE_kb],
                                   ebqe_porosity[ebNE_kb],
                                   dissipation_model_flag,
                                   grad_kappa_ext_dummy,
                                   grad_u_old_ext,
                                   bc_m_ext,
                                   bc_dm_ext,
                                   bc_f_ext,
                                   bc_df_ext,
                                   a_ext,
                                   da_ext,
                                   r_ext,
                                   dr_ext);
              //
              //moving mesh
              //
              double velocity_ext[nSpace];
              velocity_ext[0] = ebqe_velocity_ext[ebNE_kb_nSpace+0] - MOVING_DOMAIN*xt_ext;
              velocity_ext[1] = ebqe_velocity_ext[ebNE_kb_nSpace+1] - MOVING_DOMAIN*yt_ext;
              //velocity_ext[2] = ebqe_velocity_ext[ebNE_kb_nSpace+2] - MOVING_DOMAIN*zt_ext;
              //
              //calculate the numerical fluxes
              //
              exteriorNumericalAdvectiveFlux(isDOFBoundary_u[ebNE_kb],
                                             isAdvectiveFluxBoundary_u[ebNE_kb],
                                             normal,
                                             bc_u_ext,
                                             ebqe_bc_advectiveFlux_u_ext[ebNE_kb],
                                             u_ext,//smoothedHeaviside(eps,ebqe_phi[ebNE_kb]),
                                             velocity_ext,
                                             flux_ext);
              //diffusive flux now as well
              //for now just apply flux boundary through advection term
              const double bc_diffusive_flux = ebqe_bc_diffusiveFlux_u_ext[ebNE_kb];
              exteriorNumericalDiffusiveFlux(bc_diffusive_flux,
                                             isDOFBoundary_u[ebNE_kb],
                                             isDiffusiveFluxBoundary_u[ebNE_kb],
                                             normal,
                                             bc_u_ext,
                                             a_ext,
                                             grad_u_ext,
                                             u_ext,
                                             ebqe_penalty_ext[ebNE_kb],//penalty,
                                             diffusive_flux_ext);
              //mwf debug
              //std::cout<<"Residual ebNE= "<<ebNE<<" kb= "<<kb <<" penalty= "<<ebqe_penalty_ext[ebNE_kb] <<std::endl;
              flux_ext += diffusive_flux_ext;
              ebqe_flux[ebNE_kb] = flux_ext;
              //save for other models? cek need to be consistent with numerical flux
              if(flux_ext >=0.0)
                ebqe_u[ebNE_kb] = u_ext;
              else
                ebqe_u[ebNE_kb] = bc_u_ext;
              //
              //update residuals
              //
              for (int i=0;i<nDOF_test_element;i++)
                {
                  int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

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
    }

    void calculateJacobian(//element
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
                           //diffusion terms
                           double nu_0,
                           double nu_1,
                           double sigma_e,
                           double c_mu,
                           double c_1,
                           double c_2,
                           double c_e,
                           double rho_0,
                           double rho_1,
                           int dissipation_model_flag,
                           //end diffusion
                           double useMetrics,
                           double alphaBDF,
                           int lag_shockCapturing,/*mwf not used yet*/
                           double shockCapturingDiffusion,
                           int* u_l2g,
                           double* elementDiameter,
                           double* u_dof, double* u_dof_old,
                           double* velocity,
                           double* phi_ls, //level set variable
                           double* q_kappa, //kinetic energy
                           double* q_grad_kappa,
                           double* q_porosity,//VRANS
                           //velocity dof
                           double * velocity_dof_u,
                           double * velocity_dof_v,
                           double * velocity_dof_w,
                           //end velocity dof
                           double* q_m_betaBDF,
                           double* cfl,
                           double* q_numDiff_u_last,
                           double* ebqe_penalty_ext, //penalty
                           int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                           double* globalJacobian,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           double* ebqe_velocity_ext,
                           int* isDOFBoundary_u,
                           double* ebqe_bc_u_ext,
                           int* isAdvectiveFluxBoundary_u,
                           double* ebqe_bc_advectiveFlux_u_ext,
                           int* isDiffusiveFluxBoundary_u,
                           double* ebqe_bc_diffusiveFlux_u_ext,
                           int* csrColumnOffsets_eb_u_u,
                           double* ebqe_phi,double epsFact,
                           double* ebqe_kappa, //kinetic energy on boundary
                           double* ebqe_porosity)//VRANS
    {
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
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

              //declare local storage
              register double u=0.0,u_old=0.0,
                grad_u[nSpace],grad_vx[nSpace],grad_vy[nSpace],//grad_vz[nSpace],
                grad_u_old[nSpace],
                m=0.0,dm=0.0,a=0.0,da=0.0,r=0.0,dr=0.0,
                f[nSpace],df[nSpace],df_minus_da_grad_u[nSpace],
                m_t=0.0,dm_t=0.0,
                dpdeResidual_u_u[nDOF_trial_element],
                Lstar_u[nDOF_test_element],
                dsubgridError_u_u[nDOF_trial_element],
                tau=0.0,tau0=0.0,tau1=0.0,
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                u_grad_trial[nDOF_trial_element*nSpace],
                dV,
                u_test_dV[nDOF_test_element],
                u_grad_test_dV[nDOF_test_element*nSpace],
                x,y,z,xt,yt,zt,
                G[nSpace*nSpace],G_dd_G,tr_G;
              //
              //calculate solution and gradients at quadrature points
              //
              // u=0.0;
              // for (int I=0;I<nSpace;I++)
              //   {
              //     grad_u[I]=0.0;
              //   }
              // for (int j=0;j<nDOF_trial_element;j++)
              //   {
              //     int eN_j=eN*nDOF_trial_element+j;
              //     int eN_k_j=eN_k*nDOF_trial_element+j;
              //     int eN_k_j_nSpace = eN_k_j*nSpace;

              //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
              //     for (int I=0;I<nSpace;I++)
              //       {
              //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
              //       }
              //   }
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
              ck.calculateG(jacInv,G,G_dd_G,tr_G);
              //get the trial function gradients
              ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
              //get the solution
              ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
              ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u_old);
              //get the solution gradients
              ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
              ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_old);
              //
              //compute velocity production terms, ***assumes same spaces for velocity dofs and Dissipation2D!***
              ck.gradFromDOF(velocity_dof_u,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vx);
              ck.gradFromDOF(velocity_dof_v,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vy);
              //ck.gradFromDOF(velocity_dof_w,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_vz);
              //

              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                  for (int I=0;I<nSpace;I++)
                    {
                      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                    }
                }
              //
              //calculate pde coefficients and derivatives at quadrature points
              //
              evaluateCoefficients(&velocity[eN_k_nSpace],
                                   epsFact,
                                   phi_ls[eN_k],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx,
                                   grad_vy,
                                   //grad_vz,
                                   u,
                                   u_old,
                                   q_kappa[eN_k],
                                   q_porosity[eN_k],
                                   dissipation_model_flag,
                                   &q_grad_kappa[eN_k_nSpace],
                                   grad_u_old,
                                   m,
                                   dm,
                                   f,
                                   df,
                                   a,
                                   da,
                                   r,
                                   dr);
              //
              //moving mesh
              //
              f[0] -= MOVING_DOMAIN*m*xt;
              f[1] -= MOVING_DOMAIN*m*yt;
              //f[2] -= MOVING_DOMAIN*m*zt;
              df[0] -= MOVING_DOMAIN*dm*xt;
              df[1] -= MOVING_DOMAIN*dm*yt;
              //df[2] -= MOVING_DOMAIN*dm*zt;

              //
              //account for nonlinearity in diffusion coefficient by piggy backing on advection routines
              //
              df_minus_da_grad_u[0] = df[0] - da*grad_u[0];
              df_minus_da_grad_u[1] = df[1] - da*grad_u[1];
              //df_minus_da_grad_u[2] = df[2] - da*grad_u[2];

              //
              //calculate time derivatives
              //
              ck.bdf(alphaBDF,
                     q_m_betaBDF[eN_k],
                     m,
                     dm,
                     m_t,
                     dm_t);
              //
              //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
              //
              //calculate the adjoint times the test functions
              for (int i=0;i<nDOF_test_element;i++)
                {
                  // int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
                  // Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
                  register int i_nSpace = i*nSpace;
                  Lstar_u[i]=ck.Advection_adjoint(df_minus_da_grad_u,&u_grad_test_dV[i_nSpace])
                    + ck.Reaction_adjoint(dr,u_test_dV[i]);
                }
              //calculate the Jacobian of strong residual
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  //int eN_k_j=eN_k*nDOF_trial_element+j;
                  //int eN_k_j_nSpace = eN_k_j*nSpace;
                  int j_nSpace = j*nSpace;
                  dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
                    ck.AdvectionJacobian_strong(df_minus_da_grad_u,&u_grad_trial[j_nSpace]) +
                    ck.ReactionJacobian_strong(dr,u_trial_ref[k*nDOF_trial_element+j]);
                }
              //tau and tau*Res
              calculateSubgridError_tau(elementDiameter[eN],
                                        dm_t + dr,
                                        df_minus_da_grad_u,
                                        cfl[eN_k],
                                        tau0);

              calculateSubgridError_tau(Ct_sge,
                                        G,
                                        dm_t + dr,
                                        df_minus_da_grad_u,
                                        tau1,
                                        cfl[eN_k]);
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

              for(int j=0;j<nDOF_trial_element;j++)
                dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
              for(int i=0;i<nDOF_test_element;i++)
                {
                  //int eN_k_i=eN_k*nDOF_test_element+i;
                  //int eN_k_i_nSpace=eN_k_i*nSpace;
                  for(int j=0;j<nDOF_trial_element;j++)
                    {
                      //int eN_k_j=eN_k*nDOF_trial_element+j;
                      //int eN_k_j_nSpace = eN_k_j*nSpace;
                      int j_nSpace = j*nSpace;
                      int i_nSpace = i*nSpace;
                      elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) +
                        ck.AdvectionJacobian_weak(df_minus_da_grad_u,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
                        ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) +
                        ck.NumericalDiffusionJacobian(a,&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]) + //steal numericalDiffusion for scalar term
                        ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]) +
                        ck.ReactionJacobian_weak(dr,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]);
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
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
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

              register double u_ext=0.0,u_old_ext=0.0,
                grad_u_ext[nSpace],
                grad_vx_ext[nSpace],grad_vy_ext[nSpace],//grad_vz_ext[nSpace],
                grad_u_old_ext[nSpace],grad_kappa_ext_dummy[nSpace],
                m_ext=0.0,
                dm_ext=0.0,
                f_ext[nSpace],
                df_ext[nSpace],
                dflux_u_u_ext=0.0,
                a_ext=0.0,da_ext=0.0,r_ext=0.0,dr_ext=0.0,
                bc_u_ext=0.0,
                //bc_grad_u_ext[nSpace],
                bc_m_ext=0.0,
                bc_dm_ext=0.0,
                bc_f_ext[nSpace],
                bc_df_ext[nSpace],
                fluxJacobian_u_u[nDOF_trial_element],
                diffusiveFluxJacobian_u_u[nDOF_trial_element],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,
                u_test_dS[nDOF_test_element],
                u_grad_trial_trace[nDOF_trial_element*nSpace],
                normal[3],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                G[nSpace*nSpace],G_dd_G,tr_G;
              //
              //calculate the solution and gradients at quadrature points
              //
              // u_ext=0.0;
              // for (int I=0;I<nSpace;I++)
              //   {
              //     grad_u_ext[I] = 0.0;
              //     bc_grad_u_ext[I] = 0.0;
              //   }
              // for (int j=0;j<nDOF_trial_element;j++)
              //   {
              //     register int eN_j = eN*nDOF_trial_element+j,
              //       ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
              //       ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
              //     u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]);

              //     for (int I=0;I<nSpace;I++)
              //       {
              //         grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]);
              //       }
              //   }
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
              dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
              //dS = metricTensorDetSqrt*dS_ref[kb];
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
              //solution and gradients
              ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
              ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_old_ext);
              ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
              ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_old_ext);

              //mwf hack, skip on boundary for now
              grad_kappa_ext_dummy[0] = 0.0; grad_kappa_ext_dummy[1] = 0.0; //grad_kappa_ext_dummy[2] = 0.0;

              //
              //compute velocity production terms, ***assumes same spaces for velocity dofs and Dissipation2D!***
              ck.gradFromDOF(velocity_dof_u,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vx_ext);
              ck.gradFromDOF(velocity_dof_v,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vy_ext);
              //ck.gradFromDOF(velocity_dof_w,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_vz_ext);
              //

              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                }
              //
              //load the boundary values
              //
              bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
              //
              //calculate the internal and external trace of the pde coefficients
              //
              evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
                                   epsFact,
                                   ebqe_phi[ebNE_kb],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx_ext,
                                   grad_vy_ext,
                                   //grad_vz_ext,
                                   u_ext,
                                   u_old_ext,
                                   ebqe_kappa[ebNE_kb],
                                   ebqe_porosity[ebNE_kb],
                                   dissipation_model_flag,
                                   grad_kappa_ext_dummy,
                                   grad_u_old_ext,
                                   m_ext,
                                   dm_ext,
                                   f_ext,
                                   df_ext,
                                   a_ext,
                                   da_ext,
                                   r_ext,
                                   dr_ext);
              evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
                                   epsFact,
                                   ebqe_phi[ebNE_kb],
                                   nu_0,
                                   nu_1,
                                   sigma_e,
                                   c_mu,
                                   c_1,
                                   c_2,
                                   c_e,
                                   grad_vx_ext,
                                   grad_vy_ext,
                                   //grad_vz_ext,
                                   bc_u_ext,
                                   bc_u_ext,
                                   ebqe_kappa[ebNE_kb],
                                   ebqe_porosity[ebNE_kb],
                                   dissipation_model_flag,
                                   grad_kappa_ext_dummy,
                                   grad_u_old_ext,
                                   bc_m_ext,
                                   bc_dm_ext,
                                   bc_f_ext,
                                   bc_df_ext,
                                   a_ext,
                                   da_ext,
                                   r_ext,
                                   dr_ext);
              //
              //moving domain
              //
              double velocity_ext[nSpace];
              velocity_ext[0] = ebqe_velocity_ext[ebNE_kb_nSpace+0] - MOVING_DOMAIN*xt_ext;
              velocity_ext[1] = ebqe_velocity_ext[ebNE_kb_nSpace+1] - MOVING_DOMAIN*yt_ext;
              //velocity_ext[2] = ebqe_velocity_ext[ebNE_kb_nSpace+2] - MOVING_DOMAIN*zt_ext;
              //
              //calculate the numerical fluxes
              //
              exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
                                                       isAdvectiveFluxBoundary_u[ebNE_kb],
                                                       normal,
                                                       velocity_ext,//ebqe_velocity_ext[ebNE_kb_nSpace],
                                                       dflux_u_u_ext);
              //
              //calculate the flux jacobian
              //
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
                  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
                  //diffusive flux
                  exteriorNumericalDiffusiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
                                                           isDiffusiveFluxBoundary_u[ebNE_kb],
                                                           normal,
                                                           a_ext,
                                                           da_ext,
                                                           grad_u_ext,
                                                           &u_grad_trial_trace[j*nSpace],
                                                           u_trial_trace_ref[ebN_local_kb_j],
                                                           ebqe_penalty_ext[ebNE_kb],//penalty,
                                                           diffusiveFluxJacobian_u_u[j]);
                  //mwf debug
                  //std::cout<<"Jacobian ebNE= "<<ebNE<<" kb= "<<kb <<" penalty= "<<ebqe_penalty_ext[ebNE_kb] <<std::endl;
                  fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
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

                      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS[i]
                        + diffusiveFluxJacobian_u_u[j]*u_test_dS[i];
                    }//j
                }//i
            }//kb
        }//ebNE
    }//computeJacobian
  };//Dissipation2D

  inline Dissipation2D_base* newDissipation2D(int nSpaceIn,
                                int nQuadraturePoints_elementIn,
                                int nDOF_mesh_trial_elementIn,
                                int nDOF_trial_elementIn,
                                int nDOF_test_elementIn,
                                int nQuadraturePoints_elementBoundaryIn,
                                int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization2D<Dissipation2D_base,Dissipation2D,CompKernel>(nSpaceIn,
                                                                                                   nQuadraturePoints_elementIn,
                                                                                                   nDOF_mesh_trial_elementIn,
                                                                                                   nDOF_trial_elementIn,
                                                                                                   nDOF_test_elementIn,
                                                                                                   nQuadraturePoints_elementBoundaryIn,
                                                                                                   CompKernelFlag);
  }
}//proteus
#endif
