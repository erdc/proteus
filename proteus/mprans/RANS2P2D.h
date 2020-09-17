#ifndef RANS2P2D_H
#define RANS2P2D_H
#include <valarray>
#include <cmath>
#include <iostream>
#include <set>
#include <map>
#include "CompKernel.h"
#include "MixedModelFactory.h"
#include "PyEmbeddedFunctions.h"
#include "equivalent_polynomials.h"
#include "ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"
#include "mpi.h"

namespace py = pybind11;

#define ZEROVEC {0.,0.}
const bool UPWIND_DIRICHLET=true;

const  double DM=0.0;//1-mesh conservation and divergence, 0 - weak div(v) only
const  double DM2=0.0;//1-point-wise mesh volume strong-residual, 0 - div(v) only
const  double DM3=1.0;//1-point-wise divergence, 0-point-wise rate of volume change
const double inertial_term=1.0;
namespace proteus
{
  template<int nSpace, int nP, int nQ, int nEBQ>
  using GeneralizedFunctions = equivalent_polynomials::GeneralizedFunctions_mix<nSpace, nP, nQ, nEBQ>;

  class RANS2P2D_base
  {
  public:
    virtual ~RANS2P2D_base(){}
    virtual void calculateResidual(arguments_dict& args) = 0;
    virtual void calculateJacobian(arguments_dict& args) = 0;
    virtual void calculateVelocityAverage(arguments_dict& args)=0;
    virtual void getTwoPhaseAdvectionOperator(arguments_dict& args) = 0;
    virtual void getTwoPhaseInvScaledLaplaceOperator(arguments_dict& args)=0;
    virtual void getTwoPhaseScaledMassOperator(arguments_dict& args)=0;
  };

  template<class CompKernelType,
           class CompKernelType_v,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nDOF_v_trial_element,
           int nDOF_v_test_element,
           int nQuadraturePoints_elementBoundary>
  class RANS2P2D : public RANS2P2D_base
  {
  public:
    std::set<int> ifem_boundaries, ifem_boundary_elements,
      cutfem_boundaries, cutfem_boundary_elements;
    std::valarray<bool> elementIsActive;

    const int nDOF_test_X_trial_element;
    const int nDOF_test_X_v_trial_element;
    const int nDOF_v_test_X_trial_element;
    const int nDOF_v_test_X_v_trial_element;
    CompKernelType ck;
    CompKernelType_v ck_v;
    GeneralizedFunctions<nSpace,3,nQuadraturePoints_element,nQuadraturePoints_elementBoundary> gf;
    GeneralizedFunctions<nSpace,3,nQuadraturePoints_element,nQuadraturePoints_elementBoundary> gf_p;
    GeneralizedFunctions<nSpace,3,nQuadraturePoints_element,nQuadraturePoints_elementBoundary> gf_s;
    RANS2P2D():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      nDOF_test_X_v_trial_element(nDOF_test_element*nDOF_v_trial_element),
      nDOF_v_test_X_trial_element(nDOF_v_test_element*nDOF_trial_element),
      nDOF_v_test_X_v_trial_element(nDOF_v_test_element*nDOF_v_trial_element),
      ck(),
      ck_v()
    {}

    inline
    void evaluateCoefficients(const double NONCONSERVATIVE_FORM,
                              const double sigma,
                              const double rho,
                              double nu,
                              const double h_e,
                              const double smagorinskyConstant,
                              const int turbulenceClosureModel,
                              const double g[nSpace],
                              const double useVF,
                              const double& vf,
                              const double& phi,
                              const double n[nSpace],
                              const double& kappa,
                              const double porosity,//VRANS specific
                              const double phi_solid,
                              const double p_old,
                              const double u_old,
                              const double v_old,
                              const double w_old,
                              const double grad_p_old[nSpace],
                              const double grad_u_old[nSpace],
                              const double grad_v_old[nSpace],
                              const double grad_w_old[nSpace],
                              const double& p,
                              const double grad_p[nSpace],
                              const double grad_u[nSpace],
                              const double grad_v[nSpace],
                              const double grad_w[nSpace],
                              const double& u,
                              const double& v,
                              const double& w,
                              const double LAG_LES,
                              double& eddy_viscosity,
                              double& eddy_viscosity_last,
                              double& mom_u_acc,
                              double& dmom_u_acc_u,
                              double& mom_v_acc,
                              double& dmom_v_acc_v,
                              double& mom_w_acc,
                              double& dmom_w_acc_w,
                              double mass_adv[nSpace],
                              double dmass_adv_u[nSpace],
                              double dmass_adv_v[nSpace],
                              double dmass_adv_w[nSpace],
                              double mom_u_adv[nSpace],
                              double dmom_u_adv_u[nSpace],
                              double dmom_u_adv_v[nSpace],
                              double dmom_u_adv_w[nSpace],
                              double mom_v_adv[nSpace],
                              double dmom_v_adv_u[nSpace],
                              double dmom_v_adv_v[nSpace],
                              double dmom_v_adv_w[nSpace],
                              double mom_w_adv[nSpace],
                              double dmom_w_adv_u[nSpace],
                              double dmom_w_adv_v[nSpace],
                              double dmom_w_adv_w[nSpace],
                              double mom_uu_diff_ten[nSpace],
                              double mom_vv_diff_ten[nSpace],
                              double mom_ww_diff_ten[nSpace],
                              double mom_uv_diff_ten[1],
                              double mom_uw_diff_ten[1],
                              double mom_vu_diff_ten[1],
                              double mom_vw_diff_ten[1],
                              double mom_wu_diff_ten[1],
                              double mom_wv_diff_ten[1],
                              double& mom_u_source,
                              double& mom_v_source,
                              double& mom_w_source,
                              double& mom_u_ham,
                              double dmom_u_ham_grad_p[nSpace],
                              double dmom_u_ham_grad_u[nSpace],
                              double& dmom_u_ham_u,
                              double& dmom_u_ham_v,
                              double& dmom_u_ham_w,
                              double& mom_v_ham,
                              double dmom_v_ham_grad_p[nSpace],
                              double dmom_v_ham_grad_v[nSpace],
                              double& dmom_v_ham_u,
                              double& dmom_v_ham_v,
                              double& dmom_v_ham_w,
                              double& mom_w_ham,
                              double dmom_w_ham_grad_p[nSpace],
                              double dmom_w_ham_grad_w[nSpace],
                              double& dmom_w_ham_u,
                              double& dmom_w_ham_v,
                              double& dmom_w_ham_w,
                              double forcex,
                              double forcey,
                              double forcez)
    {
      double mu,norm_n,nu_t;
      //calculate eddy viscosity
      switch (turbulenceClosureModel)
        {
          double norm_S;
        case 1:
          {
            norm_S = sqrt(2.0*(grad_u[0]*grad_u[0] + grad_v[1]*grad_v[1] +
                               0.5*(grad_u[1]+grad_v[0])*(grad_u[1]+grad_v[0])));
            nu_t = smagorinskyConstant*smagorinskyConstant*h_e*h_e*norm_S;
	    break;
          }
        case 2:
          {
            double re,cs=0.0;
            norm_S = sqrt(2.0*(grad_u[0]*grad_u[0] + grad_v[1]*grad_v[1] +
                               0.5*(grad_u[1]+grad_v[0])*(grad_u[1]+grad_v[0])));
            re = h_e*h_e*norm_S/nu;
            if (re > 1.0)
              cs=0.027*pow(10.0,-3.23*pow(re,-0.92));
            nu_t = cs*h_e*h_e*norm_S;
	    break;
          }
	default:
	  {
	    nu_t=0.0;
	  }
        }
      eddy_viscosity = nu_t;
      nu += (1.0-LAG_LES)*nu_t + LAG_LES*eddy_viscosity_last;
      mu  = rho*nu;
      if (NONCONSERVATIVE_FORM > 0.0)
        {
          //u momentum accumulation
          mom_u_acc=u;//trick for non-conservative form
          dmom_u_acc_u=rho*porosity;

          //v momentum accumulation
          mom_v_acc=v;
          dmom_v_acc_v=rho*porosity;

          //mass advective flux
          mass_adv[0]=porosity*u;
          mass_adv[1]=porosity*v;

          dmass_adv_u[0]=porosity;
          dmass_adv_u[1]=0.0;

          dmass_adv_v[0]=0.0;
          dmass_adv_v[1]=porosity;

          dmass_adv_w[0]=0.0;
          dmass_adv_w[1]=0.0;

          //u momentum advective flux
          mom_u_adv[0]=0.0;
          mom_u_adv[1]=0.0;

          dmom_u_adv_u[0]=0.0;
          dmom_u_adv_u[1]=0.0;

          dmom_u_adv_v[0]=0.0;
          dmom_u_adv_v[1]=0.0;

          //v momentum advective_flux
          mom_v_adv[0]=0.0;
          mom_v_adv[1]=0.0;

          dmom_v_adv_u[0]=0.0;
          dmom_v_adv_u[1]=0.0;

          dmom_v_adv_v[0]=0.0;
          dmom_v_adv_v[1]=0.0;

          //u momentum diffusion tensor
          mom_uu_diff_ten[0] = 2.0*porosity*mu;
          mom_uu_diff_ten[1] = porosity*mu;

          mom_uv_diff_ten[0]=porosity*mu;

          //v momentum diffusion tensor
          mom_vv_diff_ten[0] = porosity*mu;
          mom_vv_diff_ten[1] = 2.0*porosity*mu;

          mom_vu_diff_ten[0]=porosity*mu;

          //momentum sources
          norm_n = sqrt(n[0]*n[0]+n[1]*n[1]);
          mom_u_source = -porosity*rho*g[0];// - porosity*d_mu*sigma*kappa*n[0];
          mom_v_source = -porosity*rho*g[1];// - porosity*d_mu*sigma*kappa*n[1];

          //u momentum Hamiltonian (pressure)
          mom_u_ham = porosity*grad_p[0];
          dmom_u_ham_grad_p[0]=porosity;
          dmom_u_ham_grad_p[1]=0.0;

          //v momentum Hamiltonian (pressure)
          mom_v_ham = porosity*grad_p[1];
          dmom_v_ham_grad_p[0]=0.0;
          dmom_v_ham_grad_p[1] = porosity;
          
          //u momentum Hamiltonian (advection)
          mom_u_ham += inertial_term*rho * porosity * (u * grad_u[0] + v * grad_u[1]);
          dmom_u_ham_grad_u[0] = inertial_term*rho * porosity * u;
          dmom_u_ham_grad_u[1] = inertial_term*rho * porosity * v;
          dmom_u_ham_u = inertial_term*rho * porosity * grad_u[0];
          dmom_u_ham_v = inertial_term*rho * porosity * grad_u[1];
          
          //v momentum Hamiltonian (advection)
          mom_v_ham += inertial_term*rho * porosity * (u * grad_v[0] + v * grad_v[1]);
          dmom_v_ham_grad_v[0] = inertial_term*rho * porosity * u;
          dmom_v_ham_grad_v[1] = inertial_term*rho * porosity * v;
          dmom_v_ham_u = inertial_term*rho * porosity * grad_v[0];
          dmom_v_ham_v = inertial_term*rho * porosity * grad_v[1];
        }
      else
        {
          //u momentum accumulation
          mom_u_acc=porosity*u;
          dmom_u_acc_u=porosity;

          //v momentum accumulation
          mom_v_acc=porosity*v;
          dmom_v_acc_v=porosity;

          //mass advective flux
          mass_adv[0]=porosity*u;
          mass_adv[1]=porosity*v;

          dmass_adv_u[0]=porosity;
          dmass_adv_u[1]=0.0;

          dmass_adv_v[0]=0.0;
          dmass_adv_v[1]=porosity;

          //u momentum advective flux
          mom_u_adv[0]=inertial_term*porosity*u*u;
          mom_u_adv[1]=inertial_term*porosity*u*v;

          dmom_u_adv_u[0]=inertial_term*2.0*porosity*u;
          dmom_u_adv_u[1]=inertial_term*porosity*v;

          dmom_u_adv_v[0]=0.0;
          dmom_u_adv_v[1]=inertial_term*porosity*u;

          //v momentum advective_flux
          mom_v_adv[0]=inertial_term*porosity*v*u;
          mom_v_adv[1]=inertial_term*porosity*v*v;

          dmom_v_adv_u[0]=inertial_term*porosity*v;
          dmom_v_adv_u[1]=0.0;

          dmom_v_adv_v[0]=inertial_term*porosity*u;
          dmom_v_adv_v[1]=inertial_term*2.0*porosity*v;

          //u momentum diffusion tensor
          mom_uu_diff_ten[0] = 2.0*porosity*nu;
          mom_uu_diff_ten[1] = porosity*nu;

          mom_uv_diff_ten[0]=porosity*nu;

          //v momentum diffusion tensor
          mom_vv_diff_ten[0] = porosity*nu;
          mom_vv_diff_ten[1] = 2.0*porosity*nu;

          mom_vu_diff_ten[0]=porosity*nu;

          //momentum sources
          norm_n = sqrt(n[0]*n[0]+n[1]*n[1]);//+n[2]*n[2]);
          mom_u_source = -porosity*g[0];
          mom_v_source = -porosity*g[1];

          //u momentum Hamiltonian (pressure)
          mom_u_ham = porosity*grad_p[0]/rho;
          dmom_u_ham_grad_p[0]=porosity/rho;
          dmom_u_ham_grad_p[1]=0.0;

          //v momentum Hamiltonian (pressure)
          mom_v_ham = porosity*grad_p[1]/rho;
          dmom_v_ham_grad_p[0]=0.0;
          dmom_v_ham_grad_p[1]=porosity/rho;

          //u momentum Hamiltonian (advection)
          dmom_u_ham_grad_u[0]=0.0;
          dmom_u_ham_grad_u[1]=0.0;
          dmom_u_ham_u =0.0;
          dmom_u_ham_v =0.0;

          //v momentum Hamiltonian (advection)
          dmom_v_ham_grad_v[0]=0.0;
          dmom_v_ham_grad_v[1]=0.0;
          dmom_v_ham_u =0.0;
          dmom_v_ham_v =0.0;
        }
      mom_u_source -= forcex;
      mom_v_source -= forcey;
    }

    int get_distance_to_ball(int n_balls,const double* ball_center, const double* ball_radius, const double x, const double y, const double z, double& distance)
    {
      distance = 1e10;
      int index = -1;
      double d_ball_i;
      for (int i=0; i<n_balls; ++i)
        {
          d_ball_i = std::sqrt((ball_center[i*3+0]-x)*(ball_center[i*3+0]-x)
                               +(ball_center[i*3+1]-y)*(ball_center[i*3+1]-y)
                               ) - ball_radius[i];
          if(d_ball_i<distance)
            {
              distance = d_ball_i;
              index = i;
            }
        }
      return index;
    }

    void get_distance_to_ith_ball(int n_balls,const double* ball_center, const double* ball_radius,
                                  int I,
                                  const double x, const double y, const double z,
                                  double& distance)
    {
      distance = std::sqrt((ball_center[I*3+0]-x)*(ball_center[I*3+0]-x)
                           + (ball_center[I*3+1]-y)*(ball_center[I*3+1]-y)
                           ) - ball_radius[I];
    }
    void get_normal_to_ith_ball(int n_balls,const double* ball_center, const double* ball_radius,
                                int I,
                                const double x, const double y, const double z,
                                double& nx, double& ny)
    {
      double distance = std::sqrt((ball_center[I*3+0]-x)*(ball_center[I*3+0]-x)
                                  + (ball_center[I*3+1]-y)*(ball_center[I*3+1]-y)
                                  );
      if (distance > 1.0e-8)
        {
          nx = (x - ball_center[I*3+0])/distance;
          ny = (y - ball_center[I*3+1])/distance;
          assert(std::fabs(std::sqrt(nx*nx + ny*ny) - 1.0) < 1.0e-10);
        }
      else
        {
          nx = 1.0;
          ny = 0.0;
        }
    }
    void get_velocity_to_ith_ball(int n_balls,const double* ball_center, const double* ball_radius,
                                  const double* ball_velocity, const double* ball_angular_velocity,
                                  int I,
                                  const double x, const double y, const double z,
                                  double& vx, double& vy)
    {
      vx = ball_velocity[3*I + 0] - ball_angular_velocity[3*I + 2]*(y-ball_center[3*I + 1]);
      vy = ball_velocity[3*I + 1] + ball_angular_velocity[3*I + 2]*(x-ball_center[3*I + 0]);
    }
    inline void updateSolidParticleTerms(const double NONCONSERVATIVE_FORM,
                                         bool element_owned,
                                         const double particle_nitsche,
                                         const double dV,
                                         const int nParticles,
                                         const int sd_offset,
                                         double* particle_signed_distances,
                                         double* particle_signed_distance_normals,
                                         double* particle_velocities,
                                         double* particle_centroids,
                                         const int use_ball_as_particle,
                                         const double* ball_center,
                                         const double* ball_radius,
                                         const double* ball_velocity,
                                         const double* ball_angular_velocity,
                                         const double* ball_density,
                                         const double porosity, //VRANS specific
                                         const double penalty,
                                         const double alpha,
                                         const double beta,
                                         const double eps_rho,
                                         const double eps_mu,
                                         const double rho_0,
                                         const double nu_0,
                                         const double rho_1,
                                         const double nu_1,
                                         const double useVF,
                                         const double vf,
                                         const double phi,
                                         const double x,
                                         const double y,
                                         const double z,
                                         const double p,
                                         const double u,
                                         const double v,
                                         const double w,
                                         const double uStar,
                                         const double vStar,
                                         const double wStar,
                                         const double eps_s,
                                         const double grad_u[nSpace],
                                         const double grad_v[nSpace],
                                         const double grad_w[nSpace],
                                         double &mass_source,
                                         double &mom_u_source,
                                         double &mom_v_source,
                                         double &mom_w_source,
                                         double dmom_u_source[nSpace],
                                         double dmom_v_source[nSpace],
                                         double dmom_w_source[nSpace],
                                         double mom_u_adv[nSpace],
                                         double mom_v_adv[nSpace],
                                         double mom_w_adv[nSpace],
                                         double dmom_u_adv_u[nSpace],
                                         double dmom_v_adv_v[nSpace],
                                         double dmom_w_adv_w[nSpace],
                                         double &mom_u_ham,
                                         double dmom_u_ham_grad_u[nSpace],
                                         double dmom_u_ham_grad_v[nSpace],
                                         double &dmom_u_ham_u,
                                         double &dmom_u_ham_v,
                                         double &dmom_u_ham_w,
                                         double &mom_v_ham,
                                         double dmom_v_ham_grad_u[nSpace],
                                         double dmom_v_ham_grad_v[nSpace],
                                         double &dmom_v_ham_u,
                                         double &dmom_v_ham_v,
                                         double &dmom_v_ham_w,
                                         double &mom_w_ham,
                                         double dmom_w_ham_grad_w[nSpace],
                                         double &dmom_w_ham_u,
                                         double &dmom_w_ham_v,
                                         double &dmom_w_ham_w,
                                         double &mass_ham,
                                         double &dmass_ham_u,
                                         double &dmass_ham_v,
                                         double &dmass_ham_w,
                                         double *particle_netForces,
                                         double *particle_netMoments,
                                         double *particle_surfaceArea)
    {
      double C, rho, mu, nu, H_mu, ImH_mu, uc, duc_du, duc_dv, duc_dw, H_s, ImH_s, D_s, phi_s, u_s, v_s, w_s;
      double force_x, force_y, r_x, r_y, force_p_x, force_p_y, force_stress_x, force_stress_y;
      double phi_s_normal[nSpace]=ZEROVEC;
      double fluid_outward_normal[nSpace]=ZEROVEC;
      double vel[nSpace]=ZEROVEC;
      double center[nSpace]=ZEROVEC;
      H_mu = (1.0 - useVF) * gf.H(eps_mu, phi) + useVF * fmin(1.0, fmax(0.0, vf));
      ImH_mu = (1.0 - useVF) * gf.ImH(eps_mu, phi) + useVF * (1.0-fmin(1.0, fmax(0.0, vf)));
      nu = nu_0 * ImH_mu + nu_1 * H_mu;
      rho = rho_0 * ImH_mu + rho_1 * H_mu;
      mu = rho_0 * nu_0 * ImH_mu + rho_1 * nu_1 * H_mu;
      C = 0.0;
      for (int i = 0; i < nParticles; i++)
        {
          if(use_ball_as_particle==1)
            {
              get_distance_to_ith_ball(nParticles,ball_center,ball_radius,i,x,y,z,phi_s);
              get_velocity_to_ith_ball(nParticles,ball_center,ball_radius,
                                       ball_velocity,ball_angular_velocity,
                                       i,x,y,z,
                                       vel[0],vel[1]);
              center[0] = ball_center[3*i+0];
              center[1] = ball_center[3*i+1];
            }
          else
            {
              phi_s = particle_signed_distances[i * sd_offset];
              vel[0] = particle_velocities[i * sd_offset * 3 + 0];
              vel[1] = particle_velocities[i * sd_offset * 3 + 1];
              center[0] = particle_centroids[3*i+0];
              center[1] = particle_centroids[3*i+1];
            }
          for (int I=0;I<nSpace;I++)
            phi_s_normal[I] = particle_signed_distance_normals[I];
          assert(std::fabs(1.0-std::sqrt(phi_s_normal[0]*phi_s_normal[0] + phi_s_normal[1]*phi_s_normal[1])) < 1.0e-8);
          /* if (fabs(vel[0] - particle_velocities[i * sd_offset * 3 + 0])> 1.0e-12) */
          /*   std::cout<<"vel[0] "<<vel[0]<<'\t'<<particle_velocities[3*i+0]<<std::endl; */
          /* if(fabs(vel[1] - particle_velocities[i * sd_offset * 3 + 1])> 1.0e-12) */
          /*   std::cout<<"vel[1] "<<vel[1]<<'\t'<<particle_velocities[3*i+1]<<std::endl; */
          /* if(fabs(center[0] - particle_centroids[3*i+0])> 1.0e-12) */
          /*   std::cout<<"center[0] "<<center[0]<<'\t'<<particle_centroids[3*i+0]<<std::endl; */
          /* if(fabs(center[1] - particle_centroids[3*i+1])> 1.0e-12) */
          /*   std::cout<<"center[1] "<<center[1]<<'\t'<<particle_centroids[3*i+1]<<std::endl; */
          /* if(fabs(phi_s - particle_signed_distances[i * sd_offset]) > 1.0e-12) */
          /*   std::cout<<"phi_s "<<phi_s<<'\t'<<particle_signed_distances[i * sd_offset]<<std::endl; */
          /* if(fabs(phi_s_normal[0] - particle_signed_distance_normals[i * sd_offset * 3 + 0]) > 1.0e-12) */
          /*   std::cout<<"phi_s_normal[0] "<<phi_s_normal[0]<<'\t'<<particle_signed_distance_normals[i * sd_offset*3 + 0]<<std::endl; */
          /* if(fabs(phi_s_normal[1] - particle_signed_distance_normals[i * sd_offset * 3 + 1]) > 1.0e-12) */
          /*   std::cout<<"phi_s_normal[1] "<<phi_s_normal[1]<<'\t'<<particle_signed_distance_normals[i * sd_offset*3 + 1]<<std::endl; */
          
          fluid_outward_normal[0] = -phi_s_normal[0];
          fluid_outward_normal[1] = -phi_s_normal[1];
          assert(std::fabs(1.0-std::sqrt(fluid_outward_normal[0]*fluid_outward_normal[0] + fluid_outward_normal[1]*fluid_outward_normal[1])) < 1.0e-8);
          u_s = vel[0];
          v_s = vel[1];
          w_s = 0.;
          D_s = gf_s.D(eps_s, phi_s);

          double rel_vel_norm = sqrt((uStar - u_s) * (uStar - u_s) +
                                     (vStar - v_s) * (vStar - v_s) +
                                     (wStar - w_s) * (wStar - w_s));
          force_p_x = porosity * dV * D_s * p * fluid_outward_normal[0];
          force_p_y = porosity * dV * D_s * p * fluid_outward_normal[1];
          force_stress_x = porosity * dV * D_s * (-mu * (fluid_outward_normal[0] * 2* grad_u[0] + fluid_outward_normal[1] * (grad_u[1]+grad_v[0]))
                                                  +mu*penalty*(u-u_s));
          force_stress_y = porosity * dV * D_s * (-mu * (fluid_outward_normal[0] * (grad_u[1]+grad_v[0]) + fluid_outward_normal[1] * 2* grad_v[1])
                                                  +mu*penalty*(v-v_s));
          force_x = force_p_x + force_stress_x;
          force_y = force_p_y + force_stress_y;
          //always 3D for particle centroids
          r_x = x - center[0];
          r_y = y - center[1];

          if (element_owned)
            {
              particle_surfaceArea[i] += dV * D_s;
              particle_netForces[i * 3 + 0] += force_x;
              particle_netForces[i * 3 + 1] += force_y;
              particle_netForces[(i+  nParticles)*3+0]+= force_stress_x;
              particle_netForces[(i+2*nParticles)*3+0]+= force_p_x;
              particle_netForces[(i+  nParticles)*3+1]+= force_stress_y;
              particle_netForces[(i+2*nParticles)*3+1]+= force_p_y;
              particle_netMoments[i * 3 + 2] += (r_x * force_y - r_y * force_x);
            }
          
          
          mass_source += D_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
          
          if (NONCONSERVATIVE_FORM > 0.0)
            {
              //upwinded advective flux
              if (!UPWIND_DIRICHLET || (fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s) < 0.0)
                {
                  mom_u_source += rho*D_s*(u_s - u)*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  mom_v_source += rho*D_s*(v_s - v)*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  dmom_u_source[0] -= rho*D_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  dmom_v_source[1] -= rho*D_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                }

              //viscous flux
              mom_u_ham -= D_s * porosity * mu * (fluid_outward_normal[0] * 2* grad_u[0] + fluid_outward_normal[1] * (grad_u[1]+grad_v[0]));
              dmom_u_ham_grad_u[0] -= D_s * porosity * mu * 2 * fluid_outward_normal[0];
              dmom_u_ham_grad_u[1] -= D_s * porosity * mu * fluid_outward_normal[1];
              dmom_u_ham_grad_v[0] -= D_s * porosity * mu * fluid_outward_normal[1];
              
              mom_v_ham -= D_s * porosity * mu * (fluid_outward_normal[0] * (grad_u[1]+grad_v[0]) + fluid_outward_normal[1] * 2* grad_v[1]);
              dmom_v_ham_grad_u[1] -= D_s * porosity * mu * fluid_outward_normal[0];
              dmom_v_ham_grad_v[0] -= D_s * porosity * mu * fluid_outward_normal[0];
              dmom_v_ham_grad_v[1] -= D_s * porosity * mu * 2 * fluid_outward_normal[1];

              //Nitsche Dirichlet penalty
              mom_u_source += D_s*mu*penalty * (u - u_s);
              dmom_u_source[0] += D_s*mu*penalty;
              
              mom_v_source += D_s*mu*penalty * (v - v_s);
              dmom_v_source[1] += D_s*mu*penalty;

              //Nitsche adjoint consistency
              mom_u_adv[0] += D_s * porosity * mu * fluid_outward_normal[0] * (u - u_s);
              mom_u_adv[1] += D_s * porosity * mu * fluid_outward_normal[1] * (u - u_s);
              dmom_u_adv_u[0] += D_s * porosity * mu * fluid_outward_normal[0];
              dmom_u_adv_u[1] += D_s * porosity * mu * fluid_outward_normal[1];
              
              mom_v_adv[0] += D_s * porosity * mu * fluid_outward_normal[0] * (v - v_s);
              mom_v_adv[1] += D_s * porosity * mu * fluid_outward_normal[1] * (v - v_s);
              dmom_v_adv_v[0] += D_s * porosity * mu * fluid_outward_normal[0];
              dmom_v_adv_v[1] += D_s * porosity * mu * fluid_outward_normal[1];
            }
          else
            {
              //divided through by rho...
              //upwinded advective flux
              if (!UPWIND_DIRICHLET || (fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s) < 0.0)
                {
                  mom_u_source += D_s*u_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  mom_v_source += D_s*v_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                }
              else
                {
                  mom_u_source += D_s*u*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  dmom_u_source[0] += D_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  mom_v_source += D_s*v*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                  dmom_v_source[1] += D_s*(fluid_outward_normal[0]*u_s + fluid_outward_normal[1]*v_s);
                }

              //viscous flux
              mom_u_ham -= D_s * porosity * nu * (fluid_outward_normal[0] * 2* grad_u[0] + fluid_outward_normal[1] * (grad_u[1]+grad_v[0]));
              dmom_u_ham_grad_u[0] -= D_s * porosity * nu * 2 * fluid_outward_normal[0];
              dmom_u_ham_grad_u[1] -= D_s * porosity * nu * fluid_outward_normal[1];
              dmom_u_ham_grad_v[0] -= D_s * porosity * nu * fluid_outward_normal[1];
              
              mom_v_ham -= D_s * porosity * nu * (fluid_outward_normal[0] * (grad_u[1]+grad_v[0]) + fluid_outward_normal[1] * 2* grad_v[1]);
              dmom_v_ham_grad_u[1] -= D_s * porosity * nu * fluid_outward_normal[0];
              dmom_v_ham_grad_v[0] -= D_s * porosity * nu * fluid_outward_normal[0];
              dmom_v_ham_grad_v[1] -= D_s * porosity * nu * 2 * fluid_outward_normal[1];
              
              //Nitsche Dirichlet penalty
              mom_u_source += D_s*nu*penalty * (u - u_s);
              dmom_u_source[0] += D_s*nu*penalty;
              
              mom_v_source += D_s*nu*penalty * (v - v_s);
              dmom_v_source[1] += D_s*nu*penalty;

              //Nitsche adjoint consistency
              mom_u_adv[0] += D_s * porosity * nu * fluid_outward_normal[0] * (u - u_s);
              mom_u_adv[1] += D_s * porosity * nu * fluid_outward_normal[1] * (u - u_s);
              dmom_u_adv_u[0] += D_s * porosity * nu * fluid_outward_normal[0];
              dmom_u_adv_u[1] += D_s * porosity * nu * fluid_outward_normal[1];
              
              mom_v_adv[0] += D_s * porosity * nu * fluid_outward_normal[0] * (v - v_s);
              mom_v_adv[1] += D_s * porosity * nu * fluid_outward_normal[1] * (v - v_s);
              dmom_v_adv_v[0] += D_s * porosity * nu * fluid_outward_normal[0];
              dmom_v_adv_v[1] += D_s * porosity * nu * fluid_outward_normal[1];
            }
        }
    }
    //VRANS specific
    inline
    void updateDarcyForchheimerTerms_Ergun(const double NONCONSERVATIVE_FORM,
                                           /* const double linearDragFactor, */
                                           /* const double nonlinearDragFactor, */
                                           /* const double porosity, */
                                           /* const double meanGrainSize, */
                                           const double alpha,
                                           const double beta,
                                           const double eps_rho,
                                           const double eps_mu,
                                           const double rho_0,
                                           const double nu_0,
                                           const double rho_1,
                                           const double nu_1,
                                           const double useVF,
                                           const double vf,
                                           const double phi,
                                           const double u,
                                           const double v,
                                           const double w,
                                           const double uStar,
                                           const double vStar,
                                           const double wStar,
                                           const double eps_porous,
                                           const double phi_porous,
                                           const double u_porous,
                                           const double v_porous,
                                           const double w_porous,
                                           double& mom_u_source,
                                           double& mom_v_source,
                                           double& mom_w_source,
                                           double dmom_u_source[nSpace],
                                           double dmom_v_source[nSpace],
                                           double dmom_w_source[nSpace])
    {
      double rho,mu,nu,H_mu,ImH_mu, uc,duc_du,duc_dv,duc_dw,viscosity,H_porous;
      H_mu = (1.0-useVF)*gf.H(eps_mu,phi)+useVF*fmin(1.0,fmax(0.0,vf));
      ImH_mu = (1.0-useVF)*gf.ImH(eps_mu,phi)+useVF*(1.0-fmin(1.0,fmax(0.0,vf)));
      nu  = nu_0*ImH_mu+nu_1*H_mu;
      rho  = rho_0*ImH_mu+rho_1*H_mu;
      mu  = rho_0*nu_0*ImH_mu+rho_1*nu_1*H_mu;
      if (NONCONSERVATIVE_FORM > 0.0)
        {
          viscosity = mu;
        }
      else
        {
          viscosity = nu;
        }
      double x = fmax(0.0, fmin( 1.0, 0.5+phi_porous/(2.0*eps_porous)));//0 at phi_porous = -eps, 1 at phi_porous=eps

      // Relaxation function, Jacobsen et al. 2011, Mayer et al 1998
      H_porous = (exp(pow(x,3.5)) - 1.)/ (exp(1.) - 1.);

      //implicit
      /* uc = sqrt(u*u+v*v*+w*w);  */
      /* duc_du = u/(uc+1.0e-12); */
      /* duc_dv = v/(uc+1.0e-12); */
      /* duc_dw = w/(uc+1.0e-12); */
      //semi-implicit quadratic term
      uc = sqrt(uStar*uStar+vStar*vStar);
      duc_du = 0.0;
      duc_dv = 0.0;

      mom_u_source += H_porous*viscosity*(alpha + beta*uc)*(u-u_porous);
      mom_v_source += H_porous*viscosity*(alpha + beta*uc)*(v-v_porous);

      dmom_u_source[0] = H_porous*viscosity*(alpha + beta*uc + beta*duc_du*(u-u_porous));
      dmom_u_source[1] = H_porous*viscosity*beta*duc_dv*(u-u_porous);

      dmom_v_source[0] = H_porous*viscosity*beta*duc_du*(v-v_porous);
      dmom_v_source[1] = H_porous*viscosity*(alpha + beta*uc + beta*duc_dv*(v-v_porous));
    }

    inline
    void updateTurbulenceClosure(const double NONCONSERVATIVE_FORM,
                                 const int turbulenceClosureModel,
                                 const double eps_rho,
                                 const double eps_mu,
                                 const double rho_0,
                                 const double nu_0,
                                 const double rho_1,
                                 const double nu_1,
                                 const double useVF,
                                 const double vf,
                                 const double phi,
                                 const double porosity,
                                 const double eddy_visc_coef_0,
                                 const double turb_var_0, //k for k-eps or k-omega
                                 const double turb_var_1, //epsilon for k-epsilon, omega for k-omega
                                 const double turb_grad_0[nSpace],//grad k for k-eps,k-omega
                                 double& eddy_viscosity,
                                 double mom_uu_diff_ten[nSpace],
                                 double mom_vv_diff_ten[nSpace],
                                 double mom_ww_diff_ten[nSpace],
                                 double mom_uv_diff_ten[1],
                                 double mom_uw_diff_ten[1],
                                 double mom_vu_diff_ten[1],
                                 double mom_vw_diff_ten[1],
                                 double mom_wu_diff_ten[1],
                                 double mom_wv_diff_ten[1],
                                 double& mom_u_source,
                                 double& mom_v_source,
                                 double& mom_w_source)
    {
      /****
           eddy_visc_coef
           <= 2  LES (do nothing)
           == 3  k-epsilon

      */
      assert (turbulenceClosureModel >=3);
      double rho,nu,H_mu,ImH_mu, nu_t=0.0,nu_t_keps =0.0, nu_t_komega=0.0;
      double isKEpsilon = 1.0, dynamic_eddy_viscosity = 0.0;

      if (turbulenceClosureModel == 4)
        isKEpsilon = 0.0;
      H_mu = (1.0-useVF)*gf.H(eps_mu,phi)+useVF*fmin(1.0,fmax(0.0,vf));
      ImH_mu = (1.0-useVF)*gf.ImH(eps_mu,phi)+useVF*(1.0-fmin(1.0,fmax(0.0,vf)));
      nu  = nu_0*ImH_mu+nu_1*H_mu;
      rho  = rho_0*ImH_mu+rho_1*H_mu;

      const double twoThirds = 2.0/3.0; const double div_zero = 1.0e-2*fmin(nu_0,nu_1);
      mom_u_source += twoThirds*turb_grad_0[0];
      mom_v_source += twoThirds*turb_grad_0[1];

      //--- closure model specific ---
      //k-epsilon
      nu_t_keps = eddy_visc_coef_0*turb_var_0*turb_var_0/(fabs(turb_var_1) + div_zero);
      //k-omega
      nu_t_komega = turb_var_0/(fabs(turb_var_1) + div_zero);
      //
      nu_t = isKEpsilon*nu_t_keps + (1.0-isKEpsilon)*nu_t_komega;

      nu_t = fmax(nu_t,1.0e-4*nu); //limit according to Lew, Buscaglia etal 01
      //mwf hack
      nu_t     = fmin(nu_t,1.0e6*nu);
      eddy_viscosity = nu_t;
      if (NONCONSERVATIVE_FORM > 0.0)
        {
          dynamic_eddy_viscosity = nu_t*rho;
          //u momentum diffusion tensor
          mom_uu_diff_ten[0] += 2.0*porosity*dynamic_eddy_viscosity;
          mom_uu_diff_ten[1] += porosity*dynamic_eddy_viscosity;

          mom_uv_diff_ten[0] +=porosity*dynamic_eddy_viscosity;

          //v momentum diffusion tensor
          mom_vv_diff_ten[0] += porosity*dynamic_eddy_viscosity;
          mom_vv_diff_ten[1] += 2.0*porosity*dynamic_eddy_viscosity;

          mom_vu_diff_ten[0] += porosity*dynamic_eddy_viscosity;
        }
      else
        {
          //u momentum diffusion tensor
          mom_uu_diff_ten[0] += 2.0*porosity*eddy_viscosity;
          mom_uu_diff_ten[1] += porosity*eddy_viscosity;

          mom_uv_diff_ten[0]+=porosity*eddy_viscosity;

          //v momentum diffusion tensor
          mom_vv_diff_ten[0] += porosity*eddy_viscosity;
          mom_vv_diff_ten[1] += 2.0*porosity*eddy_viscosity;

          mom_vu_diff_ten[0]+=porosity*eddy_viscosity;
        }
    }

    inline
    void calculateSubgridError_tau(const double&  hFactor,
                                   const double& elementDiameter,
                                   const double& dmt,
                                   const double& dm,
                                   const double df[nSpace],
                                   const double& a,
                                   const double&  pfac,
                                   double& tau_v,
                                   double& tau_p,
                                   double& cfl)
    {
      double h,oneByAbsdt,density,viscosity,nrm_df;
      h = hFactor*elementDiameter;
      density = dm;
      viscosity =  a;
      nrm_df=0.0;
      for(int I=0;I<nSpace;I++)
        nrm_df+=df[I]*df[I];
      nrm_df = sqrt(nrm_df);
      cfl = nrm_df/(h*density);//this is really cfl/dt, but that's what we want to know, the step controller expect this
      oneByAbsdt =  fabs(dmt);
      tau_v = 1.0/(4.0*viscosity/(h*h) + inertial_term*(2.0*nrm_df/h + oneByAbsdt));
      tau_p = (4.0*viscosity + inertial_term*(2.0*nrm_df*h + oneByAbsdt*h*h))/pfac;
    }

    inline
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double&  Cd_sge,
                                        const double   G[nSpace*nSpace],
                                        const double&  G_dd_G,
                                        const double&  tr_G,
                                        const double&  A0,
                                        const double   Ai[nSpace],
                                        const double&  Kij,
                                        const double&  pfac,
                                        double& tau_v,
                                        double& tau_p,
                                        double& q_cfl)
    {
      double v_d_Gv=0.0;
      for(int I=0;I<nSpace;I++)
        for (int J=0;J<nSpace;J++)
          v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];
      tau_v = 1.0/sqrt(inertial_term*(Ct_sge*A0*A0 + v_d_Gv + 1.0e-12) + Cd_sge*Kij*Kij*G_dd_G);
      tau_p = 1.0/(pfac*tr_G*tau_v);
    }

    inline
    void calculateSubgridError_tauRes(const double& tau_p,
                                      const double& tau_v,
                                      const double& pdeResidualP,
                                      const double& pdeResidualU,
                                      const double& pdeResidualV,
                                      const double& pdeResidualW,
                                      double& subgridErrorP,
                                      double& subgridErrorU,
                                      double& subgridErrorV,
                                      double& subgridErrorW)
    {
      /* GLS pressure */
      subgridErrorP = -tau_p*pdeResidualP;
      /* GLS momentum */
      subgridErrorU = -tau_v*pdeResidualU;
      subgridErrorV = -tau_v*pdeResidualV;
    }

    inline
    void calculateSubgridErrorDerivatives_tauRes(const double& tau_p,
                                                 const double& tau_v,
                                                 const double dpdeResidualP_du[nDOF_v_trial_element],
                                                 const double dpdeResidualP_dv[nDOF_v_trial_element],
                                                 const double dpdeResidualP_dw[nDOF_v_trial_element],
                                                 const double dpdeResidualU_dp[nDOF_trial_element],
                                                 const double dpdeResidualU_du[nDOF_v_trial_element],
                                                 const double dpdeResidualV_dp[nDOF_trial_element],
                                                 const double dpdeResidualV_dv[nDOF_v_trial_element],
                                                 const double dpdeResidualW_dp[nDOF_trial_element],
                                                 const double dpdeResidualW_dw[nDOF_v_trial_element],
                                                 double dsubgridErrorP_du[nDOF_v_trial_element],
                                                 double dsubgridErrorP_dv[nDOF_v_trial_element],
                                                 double dsubgridErrorP_dw[nDOF_v_trial_element],
                                                 double dsubgridErrorU_dp[nDOF_trial_element],
                                                 double dsubgridErrorU_du[nDOF_v_trial_element],
                                                 double dsubgridErrorV_dp[nDOF_trial_element],
                                                 double dsubgridErrorV_dv[nDOF_v_trial_element],
                                                 double dsubgridErrorW_dp[nDOF_trial_element],
                                                 double dsubgridErrorW_dw[nDOF_v_trial_element])
    {
      for (int j=0;j<nDOF_v_trial_element;j++)
        {
          /* GLS pressure */
          dsubgridErrorP_du[j] = -tau_p*dpdeResidualP_du[j];
          dsubgridErrorP_dv[j] = -tau_p*dpdeResidualP_dv[j];
          /* GLS  momentum*/
          /* u */
          dsubgridErrorU_du[j] = -tau_v*dpdeResidualU_du[j];
          /* v */
          dsubgridErrorV_dv[j] = -tau_v*dpdeResidualV_dv[j];
        }
      for (int j=0;j<nDOF_trial_element;j++)
        {
          /* GLS  momentum*/
          /* u */
          dsubgridErrorU_dp[j] = -tau_v*dpdeResidualU_dp[j];
          /* v */
          dsubgridErrorV_dp[j] = -tau_v*dpdeResidualV_dp[j];
        }
    }

    inline
    void exteriorNumericalAdvectiveFlux(const double NONCONSERVATIVE_FORM,
                                        const int& isDOFBoundary_p,
                                        const int& isDOFBoundary_u,
                                        const int& isDOFBoundary_v,
                                        const int& isDOFBoundary_w,
                                        const int& isFluxBoundary_p,
                                        const int& isFluxBoundary_u,
                                        const int& isFluxBoundary_v,
                                        const int& isFluxBoundary_w,
                                        const double& oneByRho,
                                        const double& bc_oneByRho,
                                        const double n[nSpace],
                                        const double& bc_p,
                                        const double& bc_u,
                                        const double& bc_v,
                                        const double bc_f_mass[nSpace],
                                        const double bc_f_umom[nSpace],
                                        const double bc_f_vmom[nSpace],
                                        const double bc_f_wmom[nSpace],
                                        const double& bc_flux_mass,
                                        const double& bc_flux_umom,
                                        const double& bc_flux_vmom,
                                        const double& bc_flux_wmom,
                                        const double& p,
                                        const double& u,
                                        const double& v,
                                        const double f_mass[nSpace],
                                        const double f_umom[nSpace],
                                        const double f_vmom[nSpace],
                                        const double f_wmom[nSpace],
                                        const double df_mass_du[nSpace],
                                        const double df_mass_dv[nSpace],
                                        const double df_mass_dw[nSpace],
                                        const double df_umom_dp[nSpace],
                                        const double dham_grad[nSpace],
                                        const double df_umom_du[nSpace],
                                        const double df_umom_dv[nSpace],
                                        const double df_umom_dw[nSpace],
                                        const double df_vmom_dp[nSpace],
                                        const double df_vmom_du[nSpace],
                                        const double df_vmom_dv[nSpace],
                                        const double df_vmom_dw[nSpace],
                                        const double df_wmom_dp[nSpace],
                                        const double df_wmom_du[nSpace],
                                        const double df_wmom_dv[nSpace],
                                        const double df_wmom_dw[nSpace],
                                        double& flux_mass,
                                        double& flux_umom,
                                        double& flux_vmom,
                                        double& flux_wmom,
                                        double* velocity)
    {
      double flowSpeedNormal;
      flux_mass = 0.0;
      flux_umom = 0.0;
      flux_vmom = 0.0;
      if (NONCONSERVATIVE_FORM > 0.0)
        {
          flowSpeedNormal = n[0] * df_vmom_dv[0] + n[1] * df_umom_du[1]; //tricky, works for  moving and fixed domains
          flowSpeedNormal += n[0] * dham_grad[0] + n[1] * dham_grad[1];
        }
      else
        flowSpeedNormal = n[0] * df_vmom_dv[0] + n[1] * df_umom_du[1]; //tricky, works for  moving and fixed domains
      if (isDOFBoundary_u != 1)
        {
          flux_mass += n[0] * f_mass[0];
          velocity[0] = f_mass[0];
          if (flowSpeedNormal >= 0.0)
            {
              flux_umom += n[0] * f_umom[0];
              flux_vmom += n[0] * f_vmom[0];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  flux_umom += (0.0 - u) * flowSpeedNormal;
                }
            }
        }
      else
        {
          flux_mass += n[0] * f_mass[0];
          velocity[0] = f_mass[0];
          if (UPWIND_DIRICHLET && flowSpeedNormal >= 0.0)
            {
              flux_umom += n[0] * f_umom[0];
              flux_vmom += n[0] * f_vmom[0];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  flux_umom += (bc_u - u) * flowSpeedNormal;
                }
              else
                {
                  flux_umom += n[0] * bc_f_umom[0];
                  flux_vmom += n[0] * bc_f_vmom[0];
                }
            }
        }
      if (isDOFBoundary_v != 1)
        {
          flux_mass += n[1] * f_mass[1];
          velocity[1] = f_mass[1];
          if (flowSpeedNormal >= 0.0)
            {
              flux_umom += n[1] * f_umom[1];
              flux_vmom += n[1] * f_vmom[1];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  flux_vmom += (0.0 - v) * flowSpeedNormal;
                }
            }
        }
      else
        {
          flux_mass += n[1] * f_mass[1];
          velocity[1] = f_mass[1];
          if (UPWIND_DIRICHLET && flowSpeedNormal >= 0.0)
            {
              flux_umom += n[1] * f_umom[1];
              flux_vmom += n[1] * f_vmom[1];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  flux_vmom += (bc_v - v) * flowSpeedNormal;
                }
              else
                {
                  flux_umom += n[1] * bc_f_umom[1];
                  flux_vmom += n[1] * bc_f_vmom[1];
                }
            }
        }
      if (isDOFBoundary_p == 1)
        {
          if (NONCONSERVATIVE_FORM > 0.0)
            {
              flux_umom += n[0] * (bc_p - p);
              flux_vmom += n[1] * (bc_p - p);
            }
          else
            {
              flux_umom += n[0] * (bc_p * bc_oneByRho - p * oneByRho);
              flux_vmom += n[1] * (bc_p * bc_oneByRho - p * oneByRho);
            }
        }
      if (isFluxBoundary_p == 1)
        {
          velocity[0] += (bc_flux_mass - flux_mass) * n[0];
          velocity[1] += (bc_flux_mass - flux_mass) * n[1];
          flux_mass = bc_flux_mass;
        }
      if (isFluxBoundary_u == 1)
        {
          flux_umom = bc_flux_umom;
        }
      if (isFluxBoundary_v == 1)
        {
          flux_vmom = bc_flux_vmom;
        }
    }

    inline
    void exteriorNumericalAdvectiveFluxDerivatives(const double NONCONSERVATIVE_FORM,
                                                   const int& isDOFBoundary_p,
                                                   const int& isDOFBoundary_u,
                                                   const int& isDOFBoundary_v,
                                                   const int& isDOFBoundary_w,
                                                   const int& isFluxBoundary_p,
                                                   const int& isFluxBoundary_u,
                                                   const int& isFluxBoundary_v,
                                                   const int& isFluxBoundary_w,
                                                   const double& oneByRho,
                                                   const double n[nSpace],
                                                   const double& bc_p,
                                                   const double& bc_u,
                                                   const double& bc_v,
                                                   const double bc_f_mass[nSpace],
                                                   const double bc_f_umom[nSpace],
                                                   const double bc_f_vmom[nSpace],
                                                   const double bc_f_wmom[nSpace],
                                                   const double& bc_flux_mass,
                                                   const double& bc_flux_umom,
                                                   const double& bc_flux_vmom,
                                                   const double& bc_flux_wmom,
                                                   const double& p,
                                                   const double& u,
                                                   const double& v,
                                                   const double& dmom_u_acc_u,
                                                   const double f_mass[nSpace],
                                                   const double f_umom[nSpace],
                                                   const double f_vmom[nSpace],
                                                   const double f_wmom[nSpace],
                                                   const double df_mass_du[nSpace],
                                                   const double df_mass_dv[nSpace],
                                                   const double df_mass_dw[nSpace],
                                                   const double df_umom_dp[nSpace],
                                                   const double dham_grad[nSpace],
                                                   const double df_umom_du[nSpace],
                                                   const double df_umom_dv[nSpace],
                                                   const double df_umom_dw[nSpace],
                                                   const double df_vmom_dp[nSpace],
                                                   const double df_vmom_du[nSpace],
                                                   const double df_vmom_dv[nSpace],
                                                   const double df_vmom_dw[nSpace],
                                                   const double df_wmom_dp[nSpace],
                                                   const double df_wmom_du[nSpace],
                                                   const double df_wmom_dv[nSpace],
                                                   const double df_wmom_dw[nSpace],
                                                   double& dflux_mass_du,
                                                   double& dflux_mass_dv,
                                                   double& dflux_mass_dw,
                                                   double& dflux_umom_dp,
                                                   double& dflux_umom_du,
                                                   double& dflux_umom_dv,
                                                   double& dflux_umom_dw,
                                                   double& dflux_vmom_dp,
                                                   double& dflux_vmom_du,
                                                   double& dflux_vmom_dv,
                                                   double& dflux_vmom_dw,
                                                   double& dflux_wmom_dp,
                                                   double& dflux_wmom_du,
                                                   double& dflux_wmom_dv,
                                                   double& dflux_wmom_dw)
    {
      double flowSpeedNormal;
      dflux_mass_du = 0.0;
      dflux_mass_dv = 0.0;

      dflux_umom_dp = 0.0;
      dflux_umom_du = 0.0;
      dflux_umom_dv = 0.0;

      dflux_vmom_dp = 0.0;
      dflux_vmom_du = 0.0;
      dflux_vmom_dv = 0.0;

      dflux_wmom_dp = 0.0;
      dflux_wmom_du = 0.0;
      dflux_wmom_dv = 0.0;
      
      flowSpeedNormal=n[0]*df_vmom_dv[0]+n[1]*df_umom_du[1];//tricky, works for moving and fixed  domains
      flowSpeedNormal+=NONCONSERVATIVE_FORM*(n[0]*dham_grad[0]+n[1]*dham_grad[1]);
      if (isDOFBoundary_u != 1)
        {
          dflux_mass_du += n[0] * df_mass_du[0];
          if (flowSpeedNormal >= 0.0)
            {
              dflux_umom_du += n[0] * df_umom_du[0];
              dflux_umom_dv += n[0] * df_umom_dv[0];
              
              dflux_vmom_du += n[0] * df_vmom_du[0];
              dflux_vmom_dv += n[0] * df_vmom_dv[0];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  dflux_umom_du += dmom_u_acc_u * n[0] * (0.0 - u) - flowSpeedNormal;
                  dflux_umom_dv += dmom_u_acc_u * n[1] * (0.0 - u) ;
                }
            }
        }
      else
        {
          //cek still upwind the advection for Dirichlet?
          dflux_mass_du += n[0] * df_mass_du[0];
          if (UPWIND_DIRICHLET && flowSpeedNormal >= 0.0)
            {
              dflux_umom_du += n[0] * df_umom_du[0];
              dflux_umom_dv += n[0] * df_umom_dv[0];
              
              dflux_vmom_du += n[0] * df_vmom_du[0];
              dflux_vmom_dv += n[0] * df_vmom_dv[0];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  dflux_umom_du += dmom_u_acc_u * n[0] * (bc_u - u) - flowSpeedNormal;
                  dflux_umom_dv += dmom_u_acc_u * n[1] * (bc_u - u) ;
                }
              else
                {
                  if (isDOFBoundary_v != 1)
                    dflux_vmom_dv += n[0] * df_vmom_dv[0];
                }
            }
        }
      if (isDOFBoundary_v != 1)
        {
          dflux_mass_dv += n[1] * df_mass_dv[1];
          if (flowSpeedNormal >= 0.0)
            {
              dflux_umom_du += n[1] * df_umom_du[1];
              dflux_umom_dv += n[1] * df_umom_dv[1];
              
              dflux_vmom_du += n[1] * df_vmom_du[1];
              dflux_vmom_dv += n[1] * df_vmom_dv[1];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  dflux_vmom_du += dmom_u_acc_u * n[0] * (0.0 - v);
                  dflux_vmom_dv += dmom_u_acc_u * n[1] * (0.0 - v) - flowSpeedNormal;
                }
            }
        }
      else
        {
          //cek still upwind the advection for Dirichlet?
          dflux_mass_dv += n[1] * df_mass_dv[1];
          if (UPWIND_DIRICHLET && flowSpeedNormal >= 0.0)
            {
              dflux_umom_du += n[1] * df_umom_du[1];
              dflux_umom_dv += n[1] * df_umom_dv[1];
              
              dflux_vmom_du += n[1] * df_vmom_du[1];
              dflux_vmom_dv += n[1] * df_vmom_dv[1];
            }
          else
            {
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  dflux_vmom_du += dmom_u_acc_u * n[0] * (bc_v - v);
                  dflux_vmom_dv += dmom_u_acc_u * n[1] * (bc_v - v) - flowSpeedNormal;
                }
              else
                {
                  if (isDOFBoundary_u != 1)
                    dflux_umom_du += n[1] * df_umom_du[1];
                }
            }
        }
      if (isDOFBoundary_p == 1)
        {
          if (NONCONSERVATIVE_FORM > 0.0)
            {
              dflux_umom_dp = -n[0];
              dflux_vmom_dp = -n[1];
            }
          else
            {
              dflux_umom_dp = -n[0] * oneByRho;
              dflux_vmom_dp = -n[1] * oneByRho;
            }
        }
      if (isFluxBoundary_p == 1)
        {
          dflux_mass_du = 0.0;
          dflux_mass_dv = 0.0;
        }
      if (isFluxBoundary_u == 1)
        {
          dflux_umom_dp = 0.0;
          dflux_umom_du = 0.0;
          dflux_umom_dv = 0.0;
        }
      if (isFluxBoundary_v == 1)
        {
          dflux_vmom_dp = 0.0;
          dflux_vmom_du = 0.0;
          dflux_vmom_dv = 0.0;
        }
    }

    inline
    void exteriorNumericalDiffusiveFlux(const double& eps,
                                        const double& phi,
                                        int* rowptr,
                                        int* colind,
                                        const int& isDOFBoundary,
                                        const int& isFluxBoundary,
                                        const double n[nSpace],
                                        double* bc_a,
                                        const double& bc_u,
                                        const double& bc_flux,
                                        double* a,
                                        const double grad_potential[nSpace],
                                        const double& u,
                                        const double& penalty,
                                        double& flux)
    {
      double diffusiveVelocityComponent_I,penaltyFlux,max_a;
      if(isFluxBoundary == 1)
        {
          flux = bc_flux;
        }
      else if(isDOFBoundary == 1)
        {
          flux = 0.0;
          max_a=0.0;
          for(int I=0;I<nSpace;I++)
            {
              diffusiveVelocityComponent_I=0.0;
              for(int m=rowptr[I];m<rowptr[I+1];m++)
                {
                  diffusiveVelocityComponent_I -= a[m]*grad_potential[colind[m]];
                  max_a = fmax(max_a,a[m]);
                }
              flux+= diffusiveVelocityComponent_I*n[I];
            }
          penaltyFlux = max_a*penalty*(u-bc_u);
          flux += penaltyFlux;
          //cek: need to investigate this issue more
          //contact line slip
          //flux*=(gf.D(eps,0) - gf.D(eps,phi))/gf.D(eps,0);
        }
      else
        {
          std::cerr<<"RANS2P2D: warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
          flux = 0.0;
        }
    }


    inline
    double ExteriorNumericalDiffusiveFluxJacobian(const double& eps,
                                                  const double& phi,
                                                  int* rowptr,
                                                  int* colind,
                                                  const int& isDOFBoundary,
                                                  const int& isFluxBoundary,
                                                  const double n[nSpace],
                                                  double* a,
                                                  const double& v,
                                                  const double grad_v[nSpace],
                                                  const double& penalty)
    {
      double dvel_I,tmp=0.0,max_a=0.0;
      if(isFluxBoundary==0 && isDOFBoundary==1)
        {
          for(int I=0;I<nSpace;I++)
            {
              dvel_I=0.0;
              for(int m=rowptr[I];m<rowptr[I+1];m++)
                {
                  dvel_I -= a[m]*grad_v[colind[m]];
                  max_a = fmax(max_a,a[m]);
                }
              tmp += dvel_I*n[I];
            }
          tmp +=max_a*penalty*v;
          //cek: need to investigate this issue more
          //contact line slip
          //tmp*=(gf.D(eps,0) - gf.D(eps,phi))/gf.D(eps,0);
        }
      return tmp;
    }

    void calculateResidual(arguments_dict& args)
    {
      double NONCONSERVATIVE_FORM = args.m_dscalar["NONCONSERVATIVE_FORM"];
      double MOMENTUM_SGE = args.m_dscalar["MOMENTUM_SGE"];
      double PRESSURE_SGE = args.m_dscalar["PRESSURE_SGE"];
      double VELOCITY_SGE = args.m_dscalar["VELOCITY_SGE"];
      double PRESSURE_PROJECTION_STABILIZATION = args.m_dscalar["PRESSURE_PROJECTION_STABILIZATION"];
      xt::pyarray<double>& numerical_viscosity = args.m_darray["numerical_viscosity"];
      xt::pyarray<double>& mesh_trial_ref = args.m_darray["mesh_trial_ref"];
      xt::pyarray<double>& mesh_grad_trial_ref = args.m_darray["mesh_grad_trial_ref"];
      xt::pyarray<double>& mesh_dof = args.m_darray["mesh_dof"];
      xt::pyarray<double>& mesh_velocity_dof = args.m_darray["mesh_velocity_dof"];
      double MOVING_DOMAIN = args.m_dscalar["MOVING_DOMAIN"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& x_ref = args.m_darray["x_ref"];
      xt::pyarray<double>& dV_ref = args.m_darray["dV_ref"];
      xt::pyarray<double>& p_trial_ref = args.m_darray["p_trial_ref"];
      xt::pyarray<double>& p_grad_trial_ref = args.m_darray["p_grad_trial_ref"];
      xt::pyarray<double>& p_test_ref = args.m_darray["p_test_ref"];
      xt::pyarray<double>& p_grad_test_ref = args.m_darray["p_grad_test_ref"];
      xt::pyarray<double>& vel_trial_ref = args.m_darray["vel_trial_ref"];
      xt::pyarray<double>& vel_grad_trial_ref = args.m_darray["vel_grad_trial_ref"];
      xt::pyarray<double>& vel_test_ref = args.m_darray["vel_test_ref"];
      xt::pyarray<double>& vel_grad_test_ref = args.m_darray["vel_grad_test_ref"];
      xt::pyarray<double>& mesh_trial_trace_ref = args.m_darray["mesh_trial_trace_ref"];
      xt::pyarray<double>& mesh_grad_trial_trace_ref = args.m_darray["mesh_grad_trial_trace_ref"];
      xt::pyarray<double>& xb_ref = args.m_darray["xb_ref"];
      xt::pyarray<double>& dS_ref = args.m_darray["dS_ref"];
      xt::pyarray<double>& p_trial_trace_ref = args.m_darray["p_trial_trace_ref"];
      xt::pyarray<double>& p_grad_trial_trace_ref = args.m_darray["p_grad_trial_trace_ref"];
      xt::pyarray<double>& p_test_trace_ref = args.m_darray["p_test_trace_ref"];
      xt::pyarray<double>& p_grad_test_trace_ref = args.m_darray["p_grad_test_trace_ref"];
      xt::pyarray<double>& vel_trial_trace_ref = args.m_darray["vel_trial_trace_ref"];
      xt::pyarray<double>& vel_grad_trial_trace_ref = args.m_darray["vel_grad_trial_trace_ref"];
      xt::pyarray<double>& vel_test_trace_ref = args.m_darray["vel_test_trace_ref"];
      xt::pyarray<double>& vel_grad_test_trace_ref = args.m_darray["vel_grad_test_trace_ref"];
      xt::pyarray<double>& normal_ref = args.m_darray["normal_ref"];
      xt::pyarray<double>& boundaryJac_ref = args.m_darray["boundaryJac_ref"];
      double eb_adjoint_sigma = args.m_dscalar["eb_adjoint_sigma"];
      xt::pyarray<double>& elementDiameter = args.m_darray["elementDiameter"];
      xt::pyarray<double>& elementBoundaryDiameter = args.m_darray["elementBoundaryDiameter"];
      xt::pyarray<double>& nodeDiametersArray = args.m_darray["nodeDiametersArray"];
      double hFactor = args.m_dscalar["hFactor"];
      int nElements_global = args.m_iscalar["nElements_global"];
      int nElementBoundaries_owned = args.m_iscalar["nElementBoundaries_owned"];
      double useRBLES = args.m_dscalar["useRBLES"];
      double useMetrics = args.m_dscalar["useMetrics"];
      double alphaBDF = args.m_dscalar["alphaBDF"];
      double epsFact_rho = args.m_dscalar["epsFact_rho"];
      double epsFact_mu = args.m_dscalar["epsFact_mu"];
      double sigma = args.m_dscalar["sigma"];
      double rho_0 = args.m_dscalar["rho_0"];
      double nu_0 = args.m_dscalar["nu_0"];
      double rho_1 = args.m_dscalar["rho_1"];
      double nu_1 = args.m_dscalar["nu_1"];
      double smagorinskyConstant = args.m_dscalar["smagorinskyConstant"];
      int turbulenceClosureModel = args.m_iscalar["turbulenceClosureModel"];
      double Ct_sge = args.m_dscalar["Ct_sge"];
      double Cd_sge = args.m_dscalar["Cd_sge"];
      double C_dc = args.m_dscalar["C_dc"];
      double C_b = args.m_dscalar["C_b"];
      const xt::pyarray<double>& eps_solid = args.m_darray["eps_solid"];
      xt::pyarray<double>& phi_solid = args.m_darray["phi_solid"];
      const xt::pyarray<double>& eps_porous = args.m_darray["eps_porous"];
      xt::pyarray<double>& phi_porous = args.m_darray["phi_porous"];
      const xt::pyarray<double>& q_velocity_porous = args.m_darray["q_velocity_porous"];
      const xt::pyarray<double>& q_porosity = args.m_darray["q_porosity"];
      const xt::pyarray<double>& q_dragAlpha = args.m_darray["q_dragAlpha"];
      const xt::pyarray<double>& q_dragBeta = args.m_darray["q_dragBeta"];
      const xt::pyarray<double>& q_mass_source = args.m_darray["q_mass_source"];
      const xt::pyarray<double>& q_turb_var_0 = args.m_darray["q_turb_var_0"];
      const xt::pyarray<double>& q_turb_var_1 = args.m_darray["q_turb_var_1"];
      const xt::pyarray<double>& q_turb_var_grad_0 = args.m_darray["q_turb_var_grad_0"];
      const double LAG_LES = args.m_dscalar["LAG_LES"];
      xt::pyarray<double> & q_eddy_viscosity = args.m_darray["q_eddy_viscosity"];
      xt::pyarray<double> & q_eddy_viscosity_last = args.m_darray["q_eddy_viscosity_last"];
      xt::pyarray<double> & ebqe_eddy_viscosity = args.m_darray["ebqe_eddy_viscosity"];
      xt::pyarray<double> & ebqe_eddy_viscosity_last = args.m_darray["ebqe_eddy_viscosity_last"];
      xt::pyarray<int>& p_l2g = args.m_iarray["p_l2g"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<int>& rp_l2g = args.m_iarray["rp_l2g"];
      xt::pyarray<int>& rvel_l2g = args.m_iarray["rvel_l2g"];
      xt::pyarray<double>& p_dof = args.m_darray["p_dof"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      xt::pyarray<double>& p_old_dof = args.m_darray["p_old_dof"];
      xt::pyarray<double>& u_old_dof = args.m_darray["u_old_dof"];
      xt::pyarray<double>& v_old_dof = args.m_darray["v_old_dof"];
      xt::pyarray<double>& w_old_dof = args.m_darray["w_old_dof"];
      xt::pyarray<double>& g = args.m_darray["g"];
      const double useVF = args.m_dscalar["useVF"];
      xt::pyarray<double>& q_rho = args.m_darray["q_rho"];
      xt::pyarray<double>& vf = args.m_darray["vf"];
      xt::pyarray<double>& phi = args.m_darray["phi"];
      xt::pyarray<double>& phi_nodes = args.m_darray["phi_nodes"];
      xt::pyarray<double>& normal_phi = args.m_darray["normal_phi"];
      xt::pyarray<double>& kappa_phi = args.m_darray["kappa_phi"];
      xt::pyarray<double>& q_mom_u_acc = args.m_darray["q_mom_u_acc"];
      xt::pyarray<double>& q_mom_v_acc = args.m_darray["q_mom_v_acc"];
      xt::pyarray<double>& q_mom_w_acc = args.m_darray["q_mom_w_acc"];
      xt::pyarray<double>& q_mass_adv = args.m_darray["q_mass_adv"];
      xt::pyarray<double>& q_mom_u_acc_beta_bdf = args.m_darray["q_mom_u_acc_beta_bdf"];
      xt::pyarray<double>& q_mom_v_acc_beta_bdf = args.m_darray["q_mom_v_acc_beta_bdf"];
      xt::pyarray<double>& q_mom_w_acc_beta_bdf = args.m_darray["q_mom_w_acc_beta_bdf"];
      xt::pyarray<double>& q_dV = args.m_darray["q_dV"];
      xt::pyarray<double>& q_dV_last = args.m_darray["q_dV_last"];
      xt::pyarray<double>& q_velocity_sge = args.m_darray["q_velocity_sge"];
      xt::pyarray<double>& q_cfl = args.m_darray["q_cfl"];
      xt::pyarray<double>& q_numDiff_u = args.m_darray["q_numDiff_u"];
      xt::pyarray<double>& q_numDiff_v = args.m_darray["q_numDiff_v"];
      xt::pyarray<double>& q_numDiff_w = args.m_darray["q_numDiff_w"];
      xt::pyarray<double>& q_numDiff_u_last = args.m_darray["q_numDiff_u_last"];
      xt::pyarray<double>& q_numDiff_v_last = args.m_darray["q_numDiff_v_last"];
      xt::pyarray<double>& q_numDiff_w_last = args.m_darray["q_numDiff_w_last"];
      xt::pyarray<int>& sdInfo_u_u_rowptr = args.m_iarray["sdInfo_u_u_rowptr"];
      xt::pyarray<int>& sdInfo_u_u_colind = args.m_iarray["sdInfo_u_u_colind"];
      xt::pyarray<int>& sdInfo_u_v_rowptr = args.m_iarray["sdInfo_u_v_rowptr"];
      xt::pyarray<int>& sdInfo_u_v_colind = args.m_iarray["sdInfo_u_v_colind"];
      xt::pyarray<int>& sdInfo_u_w_rowptr = args.m_iarray["sdInfo_u_w_rowptr"];
      xt::pyarray<int>& sdInfo_u_w_colind = args.m_iarray["sdInfo_u_w_colind"];
      xt::pyarray<int>& sdInfo_v_v_rowptr = args.m_iarray["sdInfo_v_v_rowptr"];
      xt::pyarray<int>& sdInfo_v_v_colind = args.m_iarray["sdInfo_v_v_colind"];
      xt::pyarray<int>& sdInfo_v_u_rowptr = args.m_iarray["sdInfo_v_u_rowptr"];
      xt::pyarray<int>& sdInfo_v_u_colind = args.m_iarray["sdInfo_v_u_colind"];
      xt::pyarray<int>& sdInfo_v_w_rowptr = args.m_iarray["sdInfo_v_w_rowptr"];
      xt::pyarray<int>& sdInfo_v_w_colind = args.m_iarray["sdInfo_v_w_colind"];
      xt::pyarray<int>& sdInfo_w_w_rowptr = args.m_iarray["sdInfo_w_w_rowptr"];
      xt::pyarray<int>& sdInfo_w_w_colind = args.m_iarray["sdInfo_w_w_colind"];
      xt::pyarray<int>& sdInfo_w_u_rowptr = args.m_iarray["sdInfo_w_u_rowptr"];
      xt::pyarray<int>& sdInfo_w_u_colind = args.m_iarray["sdInfo_w_u_colind"];
      xt::pyarray<int>& sdInfo_w_v_rowptr = args.m_iarray["sdInfo_w_v_rowptr"];
      xt::pyarray<int>& sdInfo_w_v_colind = args.m_iarray["sdInfo_w_v_colind"];
      int offset_p = args.m_iscalar["offset_p"];
      int offset_u = args.m_iscalar["offset_u"];
      int offset_v = args.m_iscalar["offset_v"];
      int offset_w = args.m_iscalar["offset_w"];
      int stride_p = args.m_iscalar["stride_p"];
      int stride_u = args.m_iscalar["stride_u"];
      int stride_v = args.m_iscalar["stride_v"];
      int stride_w = args.m_iscalar["stride_w"];
      xt::pyarray<double>& globalResidual = args.m_darray["globalResidual"];
      int nExteriorElementBoundaries_global = args.m_iscalar["nExteriorElementBoundaries_global"];
      xt::pyarray<int>& exteriorElementBoundariesArray = args.m_iarray["exteriorElementBoundariesArray"];
      xt::pyarray<int>& elementBoundariesArray = args.m_iarray["elementBoundariesArray"];
      xt::pyarray<int>& elementBoundaryElementsArray = args.m_iarray["elementBoundaryElementsArray"];
      xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.m_iarray["elementBoundaryLocalElementBoundariesArray"];
      xt::pyarray<double>& ebqe_vf_ext = args.m_darray["ebqe_vf_ext"];
      xt::pyarray<double>& bc_ebqe_vf_ext = args.m_darray["bc_ebqe_vf_ext"];
      xt::pyarray<double>& ebqe_phi_ext = args.m_darray["ebqe_phi_ext"];
      xt::pyarray<double>& bc_ebqe_phi_ext = args.m_darray["bc_ebqe_phi_ext"];
      xt::pyarray<double>& ebqe_normal_phi_ext = args.m_darray["ebqe_normal_phi_ext"];
      xt::pyarray<double>& ebqe_kappa_phi_ext = args.m_darray["ebqe_kappa_phi_ext"];
      const xt::pyarray<double>& ebqe_porosity_ext = args.m_darray["ebqe_porosity_ext"];
      const xt::pyarray<double>& ebqe_turb_var_0 = args.m_darray["ebqe_turb_var_0"];
      const xt::pyarray<double>& ebqe_turb_var_1 = args.m_darray["ebqe_turb_var_1"];
      xt::pyarray<int>& isDOFBoundary_p = args.m_iarray["isDOFBoundary_p"];
      xt::pyarray<int>& isDOFBoundary_u = args.m_iarray["isDOFBoundary_u"];
      xt::pyarray<int>& isDOFBoundary_v = args.m_iarray["isDOFBoundary_v"];
      xt::pyarray<int>& isDOFBoundary_w = args.m_iarray["isDOFBoundary_w"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_p = args.m_iarray["isAdvectiveFluxBoundary_p"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_u = args.m_iarray["isAdvectiveFluxBoundary_u"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_v = args.m_iarray["isAdvectiveFluxBoundary_v"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_w = args.m_iarray["isAdvectiveFluxBoundary_w"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_u = args.m_iarray["isDiffusiveFluxBoundary_u"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_v = args.m_iarray["isDiffusiveFluxBoundary_v"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_w = args.m_iarray["isDiffusiveFluxBoundary_w"];
      xt::pyarray<double>& ebqe_bc_p_ext = args.m_darray["ebqe_bc_p_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mass_ext = args.m_darray["ebqe_bc_flux_mass_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_u_adv_ext = args.m_darray["ebqe_bc_flux_mom_u_adv_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_v_adv_ext = args.m_darray["ebqe_bc_flux_mom_v_adv_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_w_adv_ext = args.m_darray["ebqe_bc_flux_mom_w_adv_ext"];
      xt::pyarray<double>& ebqe_bc_u_ext = args.m_darray["ebqe_bc_u_ext"];
      xt::pyarray<double>& ebqe_bc_flux_u_diff_ext = args.m_darray["ebqe_bc_flux_u_diff_ext"];
      xt::pyarray<double>& ebqe_penalty_ext = args.m_darray["ebqe_penalty_ext"];
      xt::pyarray<double>& ebqe_bc_v_ext = args.m_darray["ebqe_bc_v_ext"];
      xt::pyarray<double>& ebqe_bc_flux_v_diff_ext = args.m_darray["ebqe_bc_flux_v_diff_ext"];
      xt::pyarray<double>& ebqe_bc_w_ext = args.m_darray["ebqe_bc_w_ext"];
      xt::pyarray<double>& ebqe_bc_flux_w_diff_ext = args.m_darray["ebqe_bc_flux_w_diff_ext"];
      xt::pyarray<double>& q_x = args.m_darray["q_x"];
      xt::pyarray<double>& q_u_0 = args.m_darray["q_u_0"];
      xt::pyarray<double>& q_u_1 = args.m_darray["q_u_1"];
      xt::pyarray<double>& q_u_2 = args.m_darray["q_u_2"];
      xt::pyarray<double>& q_u_3 = args.m_darray["q_u_3"];
      xt::pyarray<double>& q_velocity = args.m_darray["q_velocity"];
      xt::pyarray<double>& ebqe_velocity = args.m_darray["ebqe_velocity"];
      xt::pyarray<double>& flux = args.m_darray["flux"];
      xt::pyarray<double>& elementResidual_p_save = args.m_darray["elementResidual_p_save"];
      xt::pyarray<int>& elementFlags = args.m_iarray["elementFlags"];
      xt::pyarray<int>& boundaryFlags = args.m_iarray["boundaryFlags"];
      xt::pyarray<double>& barycenters = args.m_darray["barycenters"];
      xt::pyarray<double>& wettedAreas = args.m_darray["wettedAreas"];
      xt::pyarray<double>& netForces_p = args.m_darray["netForces_p"];
      xt::pyarray<double>& netForces_v = args.m_darray["netForces_v"];
      xt::pyarray<double>& netMoments = args.m_darray["netMoments"];
      xt::pyarray<double>& velocityError = args.m_darray["velocityError"];
      xt::pyarray<double>& velocityErrorNodal = args.m_darray["velocityErrorNodal"];
      xt::pyarray<double>& forcex = args.m_darray["forcex"];
      xt::pyarray<double>& forcey = args.m_darray["forcey"];
      xt::pyarray<double>& forcez = args.m_darray["forcez"];
      int     use_ball_as_particle = args.m_iscalar["use_ball_as_particle"];
      xt::pyarray<double>& ball_center = args.m_darray["ball_center"];
      xt::pyarray<double>& ball_radius = args.m_darray["ball_radius"];
      xt::pyarray<double>& ball_velocity = args.m_darray["ball_velocity"];
      xt::pyarray<double>& ball_angular_velocity = args.m_darray["ball_angular_velocity"];
      xt::pyarray<double>& ball_density = args.m_darray["ball_density"];
      xt::pyarray<double>& particle_signed_distances = args.m_darray["particle_signed_distances"];
      xt::pyarray<double>& particle_signed_distance_normals = args.m_darray["particle_signed_distance_normals"];
      xt::pyarray<double>& particle_velocities = args.m_darray["particle_velocities"];
      xt::pyarray<double>& particle_centroids = args.m_darray["particle_centroids"];
      xt::pyarray<double>& ebqe_phi_s = args.m_darray["ebqe_phi_s"];
      xt::pyarray<double>& ebq_global_grad_phi_s = args.m_darray["ebq_global_grad_phi_s"];
      xt::pyarray<double>& ebq_particle_velocity_s = args.m_darray["ebq_particle_velocity_s"];
      int nParticles = args.m_iscalar["nParticles"];
      xt::pyarray<double>& particle_netForces = args.m_darray["particle_netForces"];
      xt::pyarray<double>& particle_netMoments = args.m_darray["particle_netMoments"];
      xt::pyarray<double>& particle_surfaceArea = args.m_darray["particle_surfaceArea"];
      int nElements_owned = args.m_iscalar["nElements_owned"];
      double particle_nitsche = args.m_dscalar["particle_nitsche"];
      double particle_epsFact = args.m_dscalar["particle_epsFact"];
      double particle_alpha = args.m_dscalar["particle_alpha"];
      double particle_beta = args.m_dscalar["particle_beta"];
      double particle_penalty_constant = args.m_dscalar["particle_penalty_constant"];
      double ghost_penalty_constant = args.m_dscalar["ghost_penalty_constant"];
      xt::pyarray<double>& phi_solid_nodes = args.m_darray["phi_solid_nodes"];
      xt::pyarray<double>& distance_to_solids = args.m_darray["distance_to_solids"];
      bool useExact = args.m_iscalar["useExact"];
      xt::pyarray<double>& isActiveR = args.m_darray["isActiveR"];
      xt::pyarray<double>& isActiveDOF_p = args.m_darray["isActiveDOF_p"];
      xt::pyarray<double>& isActiveDOF_vel = args.m_darray["isActiveDOF_vel"];
      const bool normalize_pressure = args.m_iscalar["normalize_pressure"];
      xt::pyarray<double>& errors = args.m_darray["errors"];
      logEvent("Entered mprans calculateResidual",6);
      gf.useExact = false;//useExact;
      gf_p.useExact = false;//useExact;
      gf_s.useExact = useExact;
      ifem_boundaries.clear();
      ifem_boundary_elements.clear();
      cutfem_boundaries.clear();
      cutfem_boundary_elements.clear();
      elementIsActive.resize(nElements_global);
      const int nQuadraturePoints_global(nElements_global*nQuadraturePoints_element);
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      double p_dv=0.0,pa_dv=0.0,total_volume=0.0,total_surface_area=0.0,total_flux=0.0;
      double mesh_volume_conservation=0.0,
        mesh_volume_conservation_weak=0.0,
        mesh_volume_conservation_err_max=0.0,
        mesh_volume_conservation_err_max_weak=0.0,
        domain_volume=0.0,
        &p_L1=errors(0,0),&u_L1=errors(0,1),&v_L1=errors(0,2),&w_L1=errors(0,2),&velocity_L1=errors(0,4),
        &p_L2=errors(1,0),&u_L2=errors(1,1),&v_L2=errors(1,2),&w_L2=errors(1,2),&velocity_L2=errors(1,4),
        &p_LI=errors(2,0),&u_LI=errors(2,1),&v_LI=errors(2,2),&w_LI=errors(2,2),&velocity_LI=errors(2,4);
      p_L1=0.0; u_L1=0.0; v_L1=0.0; w_L1=0.0; velocity_L1=0.0;
      p_L2=0.0; u_L2=0.0; v_L2=0.0; w_L2=0.0; velocity_L2=0.0;
      p_LI=0.0; u_LI=0.0; v_LI=0.0; w_LI=0.0; velocity_LI=0.0;
      double globalConservationError=0.0;
      /* std::cout<<"Ball Info: center "<<ball_center[0]<<'\t'<<ball_center[1]<<std::endl */
      /*          <<"Ball Info: radius "<<ball_radius[0]<<std::endl */
      /*          <<"Ball Info: velocity "<<ball_velocity[0]<<'\t'<<ball_velocity[1]<<'\t'<<ball_velocity[2]<<std::endl */
      /*          <<"Ball Info: angular "<<ball_angular_velocity[0]<<ball_angular_velocity[1]<<ball_angular_velocity[2]<<std::endl; */
      for(int eN=0;eN<nElements_global;eN++)
        {
          //declare local storage for element residual and initialize
          register double elementResidual_p[nDOF_test_element],elementResidual_p_check[nDOF_test_element],elementResidual_mesh[nDOF_test_element],
            elementResidual_u[nDOF_v_test_element],
            elementResidual_v[nDOF_v_test_element],
            pelementResidual_u[nDOF_v_test_element],
            pelementResidual_v[nDOF_v_test_element],
            velocityErrorElement[nDOF_v_test_element],
            eps_rho,eps_mu;
          bool element_active=false;
          elementIsActive[eN]=false;
          const double* elementResidual_w(NULL);
          double mesh_volume_conservation_element=0.0,
            mesh_volume_conservation_element_weak=0.0;
          for (int i=0;i<nDOF_test_element;i++)
            {
              int eN_i = eN*nDOF_test_element+i;
              elementResidual_p_save.data()[eN_i]=0.0;
              elementResidual_mesh[i]=0.0;
              elementResidual_p[i]=0.0;
              elementResidual_p_check[i]=0.0;
            }
          for (int i=0;i<nDOF_v_test_element;i++)
            {
              elementResidual_u[i]=0.0;
              elementResidual_v[i]=0.0;
              pelementResidual_u[i]=0.0;
              pelementResidual_v[i]=0.0;
              velocityErrorElement[i]=0.0;
            }//i
          //Use for plotting result
          if(use_ball_as_particle==1 && nParticles > 0)
            {
              for (int I=0;I<nDOF_mesh_trial_element;I++)
                get_distance_to_ball(nParticles, ball_center.data(), ball_radius.data(),
                                     mesh_dof.data()[3*mesh_l2g.data()[eN*nDOF_mesh_trial_element+I]+0],
                                     mesh_dof.data()[3*mesh_l2g.data()[eN*nDOF_mesh_trial_element+I]+1],
                                     mesh_dof.data()[3*mesh_l2g.data()[eN*nDOF_mesh_trial_element+I]+2],
                                     phi_solid_nodes.data()[mesh_l2g.data()[eN*nDOF_mesh_trial_element+I]]);
            }
          else
            {
              //phi_solid_nodes is updated in PreStep
            }
          double element_phi[nDOF_mesh_trial_element], element_phi_s[nDOF_mesh_trial_element];
          for (int j=0;j<nDOF_mesh_trial_element;j++)
            {
              register int eN_j = eN*nDOF_mesh_trial_element+j;
              element_phi[j] = phi_nodes.data()[p_l2g.data()[eN_j]];
              element_phi_s[j] = phi_solid_nodes.data()[p_l2g.data()[eN_j]];
            }
          double element_nodes[nDOF_mesh_trial_element*3];
          for (int i=0;i<nDOF_mesh_trial_element;i++)
            {
              register int eN_i=eN*nDOF_mesh_trial_element+i;
              for(int I=0;I<3;I++)
                element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
            }//i
          int icase_s = gf_s.calculate(element_phi_s, element_nodes, x_ref.data(),false);
          if (icase_s == 0)
            {
              //only works for simplices
              for (int ebN_element=0;ebN_element < nDOF_mesh_trial_element; ebN_element++)
                {
                  const int ebN = elementBoundariesArray.data()[eN*nDOF_mesh_trial_element+ebN_element];
                  //internal and actually a cut edge
                  //if (elementBoundaryElementsArray[ebN*2+1] != -1 && element_phi_s[(ebN_element+1)%nDOF_mesh_trial_element]*element_phi_s[(ebN_element+2)%nDOF_mesh_trial_element] <= 0.0)
                  //  cutfem_boundaries.insert(ebN);
                  if (elementBoundaryElementsArray.data()[ebN*2+1] != -1 && (ebN < nElementBoundaries_owned))
                    cutfem_boundaries.insert(ebN);
                }
            }
#ifdef IFEM
          int icase_p = gf_p.calculate(element_phi, element_nodes, x_ref.data(), -rho_1*g.data()[1], -rho_0*g.data()[1],false,true);
          int icase = gf.calculate(element_phi, element_nodes, x_ref.data(), rho_1*nu_1, rho_0*nu_0,false,false);
#else
          int icase_p = gf_p.calculate(element_phi, element_nodes, x_ref.data(), 1.,1.,false,false);
          int icase = gf.calculate(element_phi, element_nodes, x_ref.data(), 1.,1.,false,false);
#endif
          if (icase == 0)
            {
              //only works for simplices
              for (int ebN_element=0;ebN_element < nDOF_mesh_trial_element; ebN_element++)
                {
                  const int ebN = elementBoundariesArray.data()[eN*nDOF_mesh_trial_element+ebN_element];
                  //if (elementBoundaryElementsArray.data()[ebN*2+1] != -1 && (ebN < nElementBoundaries_owned))
		  //  ifem_boundaries.insert(ebN);
                  ifem_boundaries.insert(ebN);
                }
            }
          //
          //loop over quadrature points and compute integrands
          //
          double numDiffMax=0.0;
          for(int fluid_phase=0;fluid_phase < 2 - abs(icase);fluid_phase++)
            {
              for(int k=0;k<nQuadraturePoints_element;k++)
                {
                  //compute indices and declare local storage
                  register int eN_k = eN*nQuadraturePoints_element+k,
                    eN_k_nSpace = eN_k*nSpace,
                    eN_k_3d = eN_k*3,
                    eN_nDOF_trial_element = eN*nDOF_trial_element,
                    eN_nDOF_v_trial_element = eN*nDOF_v_trial_element;
                  register double p=0.0,u=0.0,v=0.0,w=0.0,
                    grad_p[nSpace]=ZEROVEC,grad_u[nSpace]=ZEROVEC,grad_v[nSpace]=ZEROVEC,grad_w[nSpace]=ZEROVEC,
                    p_old=0.0,u_old=0.0,v_old=0.0,w_old=0.0,
                    grad_p_old[nSpace]=ZEROVEC,grad_u_old[nSpace]=ZEROVEC,grad_v_old[nSpace]=ZEROVEC,grad_w_old[nSpace]=ZEROVEC,
                    mom_u_acc=0.0,
                    dmom_u_acc_u=0.0,
                    mom_v_acc=0.0,
                    dmom_v_acc_v=0.0,
                    mom_w_acc=0.0,
                    dmom_w_acc_w=0.0,
                    mass_adv[nSpace]=ZEROVEC,
                    dmass_adv_u[nSpace]=ZEROVEC,
                    dmass_adv_v[nSpace]=ZEROVEC,
                    dmass_adv_w[nSpace]=ZEROVEC,
                    mass_ham=0.0,
                    dmass_ham_u=0.0,
                    dmass_ham_v=0.0,
                    dmass_ham_w=0.0,
                    mom_u_adv[nSpace]=ZEROVEC,
                    dmom_u_adv_u[nSpace]=ZEROVEC,
                    dmom_u_adv_v[nSpace]=ZEROVEC,
                    dmom_u_adv_w[nSpace]=ZEROVEC,
                    mom_v_adv[nSpace]=ZEROVEC,
                    dmom_v_adv_u[nSpace]=ZEROVEC,
                    dmom_v_adv_v[nSpace]=ZEROVEC,
                    dmom_v_adv_w[nSpace]=ZEROVEC,
                    mom_w_adv[nSpace]=ZEROVEC,
                    dmom_w_adv_u[nSpace]=ZEROVEC,
                    dmom_w_adv_v[nSpace]=ZEROVEC,
                    dmom_w_adv_w[nSpace]=ZEROVEC,
                    mom_uu_diff_ten[nSpace]=ZEROVEC,
                    mom_vv_diff_ten[nSpace]=ZEROVEC,
                    mom_ww_diff_ten[nSpace]=ZEROVEC,
                    mom_uv_diff_ten[1],
                    mom_uw_diff_ten[1],
                    mom_vu_diff_ten[1],
                    mom_vw_diff_ten[1],
                    mom_wu_diff_ten[1],
                    mom_wv_diff_ten[1],
                    mom_u_source=0.0,
                    mom_v_source=0.0,
                    mom_w_source=0.0,
                    mom_u_ham=0.0,
                    dmom_u_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_u[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_v[nSpace]=ZEROVEC,
                    dmom_u_ham_u=0.0,
                    dmom_u_ham_v=0.0,
                    dmom_u_ham_w=0.0,
                    mom_v_ham=0.0,
                    dmom_v_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_u[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_v[nSpace]=ZEROVEC,
                    dmom_v_ham_u=0.0,
                    dmom_v_ham_v=0.0,
                    dmom_v_ham_w=0.0,
                    mom_w_ham=0.0,
                    dmom_w_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_w_ham_grad_w[nSpace]=ZEROVEC,
                    dmom_w_ham_u=0.0,
                    dmom_w_ham_v=0.0,
                    dmom_w_ham_w=0.0,
                    mom_u_acc_t=0.0,
                    dmom_u_acc_u_t=0.0,
                    mom_v_acc_t=0.0,
                    dmom_v_acc_v_t=0.0,
                    mom_w_acc_t=0.0,
                    dmom_w_acc_w_t=0.0,
                    pdeResidual_p=0.0,
                    pdeResidual_u=0.0,
                    pdeResidual_v=0.0,
                    pdeResidual_w=0.0,
                    Lstar_u_p[nDOF_test_element],
                    Lstar_v_p[nDOF_test_element],
                    Lstar_w_p[nDOF_test_element],
                    Lstar_u_u[nDOF_v_test_element],
                    Lstar_v_v[nDOF_v_test_element],
                    Lstar_w_w[nDOF_v_test_element],
                    Lstar_p_u[nDOF_v_test_element],
                    Lstar_p_v[nDOF_v_test_element],
                    Lstar_p_w[nDOF_v_test_element],
                    subgridError_p=0.0,
                    subgridError_u=0.0,
                    subgridError_v=0.0,
                    subgridError_w=0.0,
                    tau_p=0.0,tau_p0=0.0,tau_p1=0.0,
                    tau_v=0.0,tau_v0=0.0,tau_v1=0.0,
                    jac[nSpace*nSpace],
                    jacDet,
                    jacInv[nSpace*nSpace],
                    p_trial[nDOF_trial_element], vel_trial[nDOF_v_trial_element],
                    p_grad_trial_ib[nDOF_trial_element*nSpace], vel_grad_trial_ib[nDOF_v_trial_element*nSpace],
                    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_v_trial_element*nSpace],
                    p_test_dV[nDOF_trial_element],vel_test_dV[nDOF_v_test_element],
                    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_v_test_element*nSpace],
                    dV,x,y,z,xt,yt,zt,
                    p_element_avg=0.0,
                    //
                    porosity,
                    //meanGrainSize,
                    mass_source,
                    dmom_u_source[nSpace]=ZEROVEC,
                    dmom_v_source[nSpace]=ZEROVEC,
                    dmom_w_source[nSpace]=ZEROVEC,
                    //
                    G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv,h_phi, dmom_adv_star[nSpace]=ZEROVEC,dmom_adv_sge[nSpace]=ZEROVEC,dmom_ham_grad_sge[nSpace]=ZEROVEC,
                    //embedded solid terms
                    mass_source_s=0.0,
                    mom_u_source_s=0.0,
                    mom_v_source_s=0.0,
                    mom_w_source_s=0.0,
                    dmom_u_source_s[nSpace]=ZEROVEC,
                    dmom_v_source_s[nSpace]=ZEROVEC,
                    dmom_w_source_s[nSpace]=ZEROVEC,
                    mom_u_adv_s[nSpace]=ZEROVEC,
                    mom_v_adv_s[nSpace]=ZEROVEC,
                    mom_w_adv_s[nSpace]=ZEROVEC,
                    dmom_u_adv_u_s[nSpace]=ZEROVEC,
                    dmom_v_adv_v_s[nSpace]=ZEROVEC,
                    dmom_w_adv_w_s[nSpace]=ZEROVEC,
                    mom_u_ham_s=0.0,
                    dmom_u_ham_grad_u_s[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_v_s[nSpace]=ZEROVEC,
                    dmom_u_ham_u_s=0.0,
                    dmom_u_ham_v_s=0.0,
                    dmom_u_ham_w_s=0.0,
                    mom_v_ham_s=0.0,
                    dmom_v_ham_grad_u_s[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_v_s[nSpace]=ZEROVEC,
                    dmom_v_ham_u_s=0.0,
                    dmom_v_ham_v_s=0.0,
                    dmom_v_ham_w_s=0.0,
                    mom_w_ham_s=0.0,
                    dmom_w_ham_grad_w_s[nSpace]=ZEROVEC,
                    dmom_w_ham_u_s=0.0,
                    dmom_w_ham_v_s=0.0,
                    dmom_w_ham_w_s=0.0,
                    mass_ham_s=0.0,
                    dmass_ham_u_s=0.0,
                    dmass_ham_v_s=0.0,
                    dmass_ham_w_s=0.0;
                  //get jacobian, etc for mapping reference element
                  gf_s.set_quad(k);
                  gf.set_quad(k);
                  gf_p.set_quad(k);
                  ck.calculateMapping_element(eN,
                                              k,
                                              mesh_dof.data(),
                                              mesh_l2g.data(),
                                              mesh_trial_ref.data(),
                                              mesh_grad_trial_ref.data(),
                                              jac,
                                              jacDet,
                                              jacInv,
                                              x,y,z);
                  ck.calculateH_element(eN,
                                        k,
                                        nodeDiametersArray.data(),
                                        mesh_l2g.data(),
                                        mesh_trial_ref.data(),
                                        h_phi);

                  ck.calculateMappingVelocity_element(eN,
                                                      k,
                                                      mesh_velocity_dof.data(),
                                                      mesh_l2g.data(),
                                                      mesh_trial_ref.data(),
                                                      xt,yt,zt);
                  //xt=0.0;yt=0.0;zt=0.0;
                  //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
                  //get the physical integration weight
                  dV = fabs(jacDet)*dV_ref.data()[k];
                  ck.calculateG(jacInv,G,G_dd_G,tr_G);
                  //ck.calculateGScale(G,&normal_phi.data()[eN_k_nSpace],h_phi);

                  eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                  eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);

                  //get the trial function gradients
                  ck.gradTrialFromRef(&p_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
                  ck_v.gradTrialFromRef(&vel_grad_trial_ref.data()[k*nDOF_v_trial_element*nSpace],jacInv,vel_grad_trial);
                  for (int i=0; i < nDOF_trial_element; i++)
                    {
                      p_trial[i] = p_trial_ref.data()[k*nDOF_trial_element + i];
                      p_grad_trial_ib[i*nSpace + 0] = p_grad_trial[i*nSpace+0];
                      p_grad_trial_ib[i*nSpace + 1] = p_grad_trial[i*nSpace+1];
                    }
                  for (int i=0; i < nDOF_v_trial_element; i++)
                    {
                      vel_trial[i] = vel_trial_ref.data()[k*nDOF_v_trial_element + i];
                      vel_grad_trial_ib[i*nSpace + 0] = vel_grad_trial[i*nSpace+0];
                      vel_grad_trial_ib[i*nSpace + 1] = vel_grad_trial[i*nSpace+1];
                    }
                  if (icase == 0)
                    {
#ifdef IFEMBASIS
                      for (int i=0; i < nDOF_trial_element; i++)
                        {
                          if (fluid_phase == 0)
                            {
                              if (not std::isnan(gf_p.VA(i)))
                                {
                                  p_trial[i] = gf_p.VA(i);
                                  p_grad_trial_ib[i*nSpace + 0] = gf_p.VA_x(i);
                                  p_grad_trial_ib[i*nSpace + 1] = gf_p.VA_y(i);
                                }
                            }
                          else
                            {
                              if (not std::isnan(gf_p.VB(i)))
                                {
                                  p_trial[i] = gf_p.VB(i);
                                  p_grad_trial_ib[i*nSpace + 0] = gf_p.VB_x(i);
                                  p_grad_trial_ib[i*nSpace + 1] = gf_p.VB_y(i);
                                }
                            }
                        }
                      if(nDOF_v_trial_element == nDOF_trial_element)
                        {
                          for (int vi=0; vi < nDOF_v_trial_element; vi++)
                            {
                              if (fluid_phase == 0)
                                {
                                  if (not std::isnan(gf.VA(vi)))
                                    {
                                      vel_trial[vi] = gf.VA(vi);
                                      vel_grad_trial_ib[vi*nSpace + 0] = gf.VA_x(vi);
                                      vel_grad_trial_ib[vi*nSpace + 1] = gf.VA_y(vi);
                                    }
                                }
                              else
                                {
                                  if (not std::isnan(gf.VB(vi)))
                                    {
                                      vel_trial[vi] = gf.VB(vi);
                                      vel_grad_trial_ib[vi*nSpace + 0] = gf.VB_x(vi);
                                      vel_grad_trial_ib[vi*nSpace + 1] = gf.VB_y(vi);
                                    }
                                }
                            }
                        }
#endif
#ifndef IFEM
                      bool prob=false;
                      for (int vi=0; vi < nDOF_v_trial_element; vi++)
                        {
                          //pressure
                          if (fabs(p_trial_ref.data()[k*nDOF_trial_element + vi] - p_trial[vi]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Trial "<<p_trial_ref.data()[k*nDOF_trial_element + vj]<<'\t'<<gf_p.VA(vj)<<'\t'<<gf_p.VB(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(p_grad_trial[vi*nSpace + 0] - p_grad_trial_ib[vi*nSpace+0]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Grad Trial x"<<p_grad_trial[vj*nSpace + 0]<<'\t'<<gf_p.VA_x(vj)<<'\t'<<gf_p.VB_x(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(p_grad_trial[vi*nSpace + 1] - p_grad_trial_ib[vi*nSpace+1]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Grad Trial y "<<p_grad_trial[vj*nSpace + 1]<<'\t'<<gf_p.VA_y(vj)<<'\t'<<gf_p.VB_y(vj)<<std::endl;
                              prob=true;
                            }
                          //velocity
                          if (fabs(vel_trial_ref.data()[k*nDOF_v_trial_element + vi] - vel_trial[vi]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Trial "<<vel_trial_ref.data()[k*nDOF_v_trial_element + vj]<<'\t'<<gf.VA(vj)<<'\t'<<gf.VB(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(vel_grad_trial[vi*nSpace + 0] - vel_grad_trial_ib[vi*nSpace+0]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Grad Trial x"<<vel_grad_trial[vj*nSpace + 0]<<'\t'<<gf.VA_x(vj)<<'\t'<<gf.VB_x(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(vel_grad_trial[vi*nSpace + 1] - vel_grad_trial_ib[vi*nSpace+1]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Grad Trial y "<<vel_grad_trial[vj*nSpace + 1]<<'\t'<<gf.VA_y(vj)<<'\t'<<gf.VB_y(vj)<<std::endl;
                              prob=true;
                            }
                          if (prob)
                            break;
                        }
                      assert(!prob);
#endif
                    }
                  //get the solution
                  ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_trial,p);
                  ck_v.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_trial,u);
                  ck_v.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_trial,v);
                  ck.valFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_trial,p_old);
                  ck_v.valFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_trial,u_old);
                  ck_v.valFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_trial,v_old);
                  //get the solution gradients
                  ck.gradFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_ib,grad_p);
                  ck_v.gradFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_u);
                  ck_v.gradFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_v);
                  ck.gradFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_ib,grad_p_old);
                  ck_v.gradFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_u_old);
                  ck_v.gradFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_v_old);
                  // calculate the average pressure value
                  if (PRESSURE_PROJECTION_STABILIZATION)
                    ck.DOFaverage(p_dof.data(), &p_l2g.data()[eN_nDOF_trial_element],p_element_avg);
                  //precalculate test function products with integration weights
#ifdef IFEMGALERKIN
                  for (int j=0;j<nDOF_test_element;j++)
                    {
                      p_test_dV[j] = p_trial[j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          p_grad_test_dV[j*nSpace+I]   = p_grad_trial_ib[j*nSpace+I]*dV;
                        }
                    }
                  //precalculate test function products with integration weights
                  for (int j=0;j<nDOF_v_test_element;j++)
                    {
                      vel_test_dV[j] = vel_trial[j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          vel_grad_test_dV[j*nSpace+I] = vel_grad_trial_ib[j*nSpace+I]*dV;
                        }
                    }
#else
                  for (int j=0;j<nDOF_test_element;j++)
                    {
                      p_test_dV[j] = p_test_ref.data()[k*nDOF_trial_element+j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//assume test_i = trial_i, not using ib basis here
                        }
                    }
                  //precalculate test function products with integration weights
                  for (int j=0;j<nDOF_v_test_element;j++)
                    {
                      vel_test_dV[j] = vel_test_ref.data()[k*nDOF_v_trial_element+j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//assume test_i = trial_i
                        }
                    }
#endif
                  //todo: extend this to higher-order meshes, for now assume mesh trial and p trial are same 
                  double div_mesh_velocity=0.0;
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int eN_j=eN*nDOF_trial_element+j;
                      div_mesh_velocity +=
                        mesh_velocity_dof.data()[mesh_l2g.data()[eN_j]*3+0]*p_grad_trial[j*nSpace+0] +
                        mesh_velocity_dof.data()[mesh_l2g.data()[eN_j]*3+1]*p_grad_trial[j*nSpace+1];
                    }
                  mesh_volume_conservation_element += (alphaBDF*(dV-q_dV_last.data()[eN_k])/dV - div_mesh_velocity)*dV;
                  div_mesh_velocity = DM3*div_mesh_velocity + (1.0-DM3)*alphaBDF*(dV-q_dV_last.data()[eN_k])/dV;
                  //VRANS
                  porosity      = q_porosity.data()[eN_k];
                  //
                  q_velocity.data()[eN_k_nSpace+0]=u;
                  q_velocity.data()[eN_k_nSpace+1]=v;
                  q_x.data()[eN_k_3d + 0] = x;
                  q_x.data()[eN_k_3d + 1] = y;
                  double ball_n[nSpace];
                  if (use_ball_as_particle == 1 && nParticles > 0)
                    {
                      int ball_index=get_distance_to_ball(nParticles, ball_center.data(), ball_radius.data(),x,y,z,distance_to_solids.data()[eN_k]);
                      get_normal_to_ith_ball(nParticles, ball_center.data(), ball_radius.data(),ball_index,x,y,z,ball_n[0],ball_n[1]);
                    }
                  else
                    {
                      //distance_to_solids is given in Prestep
                    }
                  if (nParticles > 0)
                    phi_solid.data()[eN_k] = distance_to_solids.data()[eN_k];
                  const double particle_eps  = particle_epsFact*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                  //
                  //calculate pde coefficients at quadrature points
                  //
                  const double H_s = gf_s.H(particle_eps,phi_solid.data()[eN_k]);
                  const double D_s = gf_s.D(particle_eps,phi_solid.data()[eN_k]);
                  if ( nParticles == 0 || H_s != 0.0 || D_s != 0.0)
                    {
                      element_active=true;
                      elementIsActive[eN]=true;
                    }
                  //save velocity at quadrature points for other models to use
                  double p_e = q_u_0.data()[eN_k] - p,
                    u_e = q_u_1.data()[eN_k] - u,
                    v_e = q_u_2.data()[eN_k] - v,
		    velocity_e=sqrt(u_e*u_e + v_e*v_e);
                  
                  
                  /* q_u_0.data()[eN_k] = p; */
                  /* q_u_1.data()[eN_k] = u; */
                  /* q_u_2.data()[eN_k] = v; */
                  /* q_u_3.data()[eN_k] = 0.0; */
                  
                  double rho,nu;
                  if (gf.useExact)
                    {
                      if (icase == 0)
                        {
                          if (fluid_phase == 0)
                            {
                              rho=rho_0;
                              nu=nu_0;
                            }
                          else
                            {
                              rho=rho_1;
                              nu=nu_1;
                            }
                        }
                      else if (icase == -1)
                        {
                          rho=rho_0;
                          nu=nu_0;
                        }
                      else if (icase == 1)
                        {
                          rho=rho_1;
                          nu=nu_1;
                        }
                      else
                        assert(false);
                    }
                  else
                    {
                      double H = (1.0-useVF)*gf.H(eps_rho,phi[eN_k]) + useVF*fmin(1.0,fmax(0.0,vf[eN_k]));
                      double ImH = (1.0-useVF)*gf.ImH(eps_rho,phi[eN_k]) + useVF*(1.0-fmin(1.0,fmax(0.0,vf[eN_k])));
                      
                      rho  = rho_0*ImH + rho_1*H;
                      nu  = nu_0*ImH + nu_1*H;
                    }
                  evaluateCoefficients(NONCONSERVATIVE_FORM,
                                       sigma,
                                       rho,
                                       nu,
                                       elementDiameter.data()[eN],
                                       smagorinskyConstant,
                                       turbulenceClosureModel,
                                       g.data(),
                                       useVF,
                                       vf.data()[eN_k],
                                       phi.data()[eN_k],
                                       &normal_phi.data()[eN_k_nSpace],
                                       kappa_phi.data()[eN_k],
                                       //VRANS
                                       porosity,
                                       phi_solid.data()[eN_k],//distance to solid
                                       p_old,
                                       u_old,
                                       v_old,
                                       w_old,
                                       grad_p_old,
                                       grad_u_old,
                                       grad_v_old,
                                       grad_w_old,
                                       //
                                       p,
                                       grad_p,
                                       grad_u,
                                       grad_v,
                                       grad_w,
                                       u,
                                       v,
                                       w,
                                       LAG_LES,
                                       q_eddy_viscosity.data()[eN_k],
                                       q_eddy_viscosity_last.data()[eN_k],
                                       mom_u_acc,
                                       dmom_u_acc_u,
                                       mom_v_acc,
                                       dmom_v_acc_v,
                                       mom_w_acc,
                                       dmom_w_acc_w,
                                       mass_adv,
                                       dmass_adv_u,
                                       dmass_adv_v,
                                       dmass_adv_w,
                                       mom_u_adv,
                                       dmom_u_adv_u,
                                       dmom_u_adv_v,
                                       dmom_u_adv_w,
                                       mom_v_adv,
                                       dmom_v_adv_u,
                                       dmom_v_adv_v,
                                       dmom_v_adv_w,
                                       mom_w_adv,
                                       dmom_w_adv_u,
                                       dmom_w_adv_v,
                                       dmom_w_adv_w,
                                       mom_uu_diff_ten,
                                       mom_vv_diff_ten,
                                       mom_ww_diff_ten,
                                       mom_uv_diff_ten,
                                       mom_uw_diff_ten,
                                       mom_vu_diff_ten,
                                       mom_vw_diff_ten,
                                       mom_wu_diff_ten,
                                       mom_wv_diff_ten,
                                       mom_u_source,
                                       mom_v_source,
                                       mom_w_source,
                                       mom_u_ham,
                                       dmom_u_ham_grad_p,
                                       dmom_u_ham_grad_u,
                                       dmom_u_ham_u,
                                       dmom_u_ham_v,
                                       dmom_u_ham_w,
                                       mom_v_ham,
                                       dmom_v_ham_grad_p,
                                       dmom_v_ham_grad_v,
                                       dmom_v_ham_u,
                                       dmom_v_ham_v,
                                       dmom_v_ham_w,
                                       mom_w_ham,
                                       dmom_w_ham_grad_p,
                                       dmom_w_ham_grad_w,
                                       dmom_w_ham_u,
                                       dmom_w_ham_v,
                                       dmom_w_ham_w,
                                       forcex.data()[eN_k],
                                       forcey.data()[eN_k],
                                       forcez.data()[eN_k]);
                  q_rho.data()[eN_k] = rho;
                  //VRANS
                  mass_source = q_mass_source.data()[eN_k];
                  //todo: decide if these should be lagged or not?
                  updateDarcyForchheimerTerms_Ergun(NONCONSERVATIVE_FORM,
                                                    /* linearDragFactor, */
                                                    /* nonlinearDragFactor, */
                                                    /* porosity, */
                                                    /* meanGrainSize, */
                                                    q_dragAlpha.data()[eN_k],
                                                    q_dragBeta.data()[eN_k],
                                                    eps_rho,
                                                    eps_mu,
                                                    rho_0,
                                                    nu_0,
                                                    rho_1,
                                                    nu_1,
                                                    useVF,
                                                    vf.data()[eN_k],
                                                    phi.data()[eN_k],
                                                    u,
                                                    v,
                                                    w,
                                                    q_velocity_sge.data()[eN_k_nSpace+0],
                                                    q_velocity_sge.data()[eN_k_nSpace+1],
                                                    q_velocity_sge.data()[eN_k_nSpace+1],//dummy entry for 2D
                                                    eps_porous.data()[elementFlags.data()[eN]],
                                                    phi_porous.data()[eN_k],
                                                    q_velocity_porous.data()[eN_k_nSpace+0],
                                                    q_velocity_porous.data()[eN_k_nSpace+1],
                                                    q_velocity_porous.data()[eN_k_nSpace+1],//dummy entry for 2D
                                                    mom_u_source,
                                                    mom_v_source,
                                                    mom_w_source,
                                                    dmom_u_source,
                                                    dmom_v_source,
                                                    dmom_w_source);
                  //Turbulence closure model
                  if (turbulenceClosureModel >= 3)
                    {
                      const double c_mu = 0.09;//mwf hack
                      updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                              turbulenceClosureModel,
                                              eps_rho,
                                              eps_mu,
                                              rho_0,
                                              nu_0,
                                              rho_1,
                                              nu_1,
                                              useVF,
                                              vf.data()[eN_k],
                                              phi.data()[eN_k],
                                              porosity,
                                              c_mu, //mwf hack
                                              q_turb_var_0.data()[eN_k],
                                              q_turb_var_1.data()[eN_k],
                                              &q_turb_var_grad_0.data()[eN_k_nSpace],
                                              q_eddy_viscosity.data()[eN_k],
                                              mom_uu_diff_ten,
                                              mom_vv_diff_ten,
                                              mom_ww_diff_ten,
                                              mom_uv_diff_ten,
                                              mom_uw_diff_ten,
                                              mom_vu_diff_ten,
                                              mom_vw_diff_ten,
                                              mom_wu_diff_ten,
                                              mom_wv_diff_ten,
                                              mom_u_source,
                                              mom_v_source,
                                              mom_w_source);
                    }
                  //
                  //moving mesh
                  //
                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_u_ham -= MOVING_DOMAIN*dmom_u_acc_u*(grad_u[0]*xt + grad_u[1]*yt);
                      dmom_u_ham_grad_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
                      dmom_u_ham_grad_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
                    }
                  else
                    {
                      mom_u_adv[0] -= MOVING_DOMAIN*mom_u_acc*xt;
                      mom_u_adv[1] -= MOVING_DOMAIN*mom_u_acc*yt;
                      dmom_u_adv_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
                      dmom_u_adv_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
                    }

                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_v_ham -= MOVING_DOMAIN*dmom_v_acc_v*(grad_v[0]*xt + grad_v[1]*yt);
                      dmom_v_ham_grad_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
                      dmom_v_ham_grad_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
                    }
                  else
                    {
                      mom_v_adv[0] -= MOVING_DOMAIN*mom_v_acc*xt;
                      mom_v_adv[1] -= MOVING_DOMAIN*mom_v_acc*yt;
                      dmom_v_adv_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
                      dmom_v_adv_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
                    }

                  //
                  //calculate time derivative at quadrature points
                  //
                  if (q_dV_last.data()[eN_k] <= -100)
                    q_dV_last.data()[eN_k] = dV;
                  q_dV.data()[eN_k] = dV;
                  ck.bdf(alphaBDF,
                         q_mom_u_acc_beta_bdf.data()[eN_k]*q_dV_last.data()[eN_k]/dV,
                         mom_u_acc,
                         dmom_u_acc_u,
                         mom_u_acc_t,
                         dmom_u_acc_u_t);
                  ck.bdf(alphaBDF,
                         q_mom_v_acc_beta_bdf.data()[eN_k]*q_dV_last.data()[eN_k]/dV,
                         mom_v_acc,
                         dmom_v_acc_v,
                         mom_v_acc_t,
                         dmom_v_acc_v_t);
                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_u_acc_t *= dmom_u_acc_u;
                      mom_v_acc_t *= dmom_v_acc_v;
                    }
                  //
                  //calculate subgrid error (strong residual and adjoint)
                  //
                  //calculate strong residual
                  pdeResidual_p = ck.Advection_strong(dmass_adv_u,grad_u) +
                    ck.Advection_strong(dmass_adv_v,grad_v) +
                    DM2*MOVING_DOMAIN*ck.Reaction_strong(alphaBDF*(dV-q_dV_last.data()[eN_k])/dV - div_mesh_velocity) +
                    ck.Reaction_strong(mass_source);

                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      dmom_adv_sge[0] = 0.0;
                      dmom_adv_sge[1] = 0.0;
                      dmom_ham_grad_sge[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
                      dmom_ham_grad_sge[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
                    }
                  else
                    {
                      dmom_adv_sge[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
                      dmom_adv_sge[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
                      dmom_ham_grad_sge[0] = 0.0;
                      dmom_ham_grad_sge[1] = 0.0;
                    }
                  double mv_tau[nSpace]=ZEROVEC;
                  mv_tau[0] = dmom_adv_sge[0] + dmom_ham_grad_sge[0];
                  mv_tau[1] = dmom_adv_sge[1] + dmom_ham_grad_sge[1];

                  pdeResidual_u = ck.Mass_strong(mom_u_acc_t) +
                    ck.Advection_strong(dmom_adv_sge,grad_u) +
                    ck.Hamiltonian_strong(dmom_ham_grad_sge,grad_u) +
                    ck.Hamiltonian_strong(dmom_u_ham_grad_p,grad_p) +
                    ck.Reaction_strong(mom_u_source) -
                    ck.Reaction_strong(dmom_u_acc_u*u*div_mesh_velocity);

                  pdeResidual_v = ck.Mass_strong(mom_v_acc_t) +
                    ck.Advection_strong(dmom_adv_sge,grad_v) +
                    ck.Hamiltonian_strong(dmom_ham_grad_sge,grad_v) +
                    ck.Hamiltonian_strong(dmom_v_ham_grad_p,grad_p) +
                    ck.Reaction_strong(mom_v_source) -
                    ck.Reaction_strong(dmom_v_acc_v*v*div_mesh_velocity);

                  //calculate tau and tau*Res
                  //add contributions from mass and source terms
                  double tmpR=dmom_u_acc_u_t + dmom_u_source[0];
                  calculateSubgridError_tau(hFactor,
                                            elementDiameter.data()[eN],
                                            tmpR,//dmom_u_acc_u_t,
                                            dmom_u_acc_u,
                                            mv_tau,//dmom_adv_sge,
                                            mom_uu_diff_ten[1],
                                            dmom_u_ham_grad_p[0],
                                            tau_v0,
                                            tau_p0,
                                            q_cfl.data()[eN_k]);

                  calculateSubgridError_tau(Ct_sge,Cd_sge,
                                            G,G_dd_G,tr_G,
                                            tmpR,//dmom_u_acc_u_t,
                                            mv_tau,//dmom_adv_sge,
                                            mom_uu_diff_ten[1],
                                            dmom_u_ham_grad_p[0],
                                            tau_v1,
                                            tau_p1,
                                            q_cfl.data()[eN_k]);

                  tau_v = useMetrics*tau_v1+(1.0-useMetrics)*tau_v0;
                  tau_p = useMetrics*tau_p1+(1.0-useMetrics)*tau_p0;

                  calculateSubgridError_tauRes(tau_p,
                                               tau_v,
                                               pdeResidual_p,
                                               pdeResidual_u,
                                               pdeResidual_v,
                                               pdeResidual_w,
                                               subgridError_p,
                                               subgridError_u,
                                               subgridError_v,
                                               subgridError_w);
                  // velocity used in adjoint (VMS or RBLES, with or without lagging the grid scale velocity)
                  dmom_adv_star[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt + useRBLES*subgridError_u);
                  dmom_adv_star[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt + useRBLES*subgridError_v);

                  mom_u_adv[0] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_u*q_velocity_sge.data()[eN_k_nSpace+0]);
                  mom_u_adv[1] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_v*q_velocity_sge.data()[eN_k_nSpace+0]);

                  mom_v_adv[0] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_u*q_velocity_sge.data()[eN_k_nSpace+1]);
                  mom_v_adv[1] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_v*q_velocity_sge.data()[eN_k_nSpace+1]);

                  // adjoint times the test functions
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
                      Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      //use the same advection adjoint for all three since we're approximating the linearized adjoint
                      Lstar_u_u[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);//cek COMP/INCOMP form have same adjoint
                      Lstar_v_v[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);//ditto
                      Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
                      Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]);

                      //VRANS account for drag terms, diagonal only here ... decide if need off diagonal terms too
                      Lstar_u_u[i]+=ck.Reaction_adjoint(dmom_u_source[0],vel_test_dV[i]);
                      Lstar_v_v[i]+=ck.Reaction_adjoint(dmom_v_source[1],vel_test_dV[i]);
                      //
                    }

                  norm_Rv = sqrt(pdeResidual_u*pdeResidual_u + pdeResidual_v*pdeResidual_v);
                  q_numDiff_u.data()[eN_k] = C_dc*norm_Rv*(useMetrics/sqrt(G_dd_G+1.0e-12)  +
                                                           (1.0-useMetrics)*hFactor*hFactor*elementDiameter.data()[eN]*elementDiameter.data()[eN]);
                  q_numDiff_v.data()[eN_k] = q_numDiff_u.data()[eN_k];
                  q_numDiff_w.data()[eN_k] = q_numDiff_u.data()[eN_k];
                  numDiffMax = std::fmax(q_numDiff_u.data()[eN_k], numDiffMax);
                  if(nParticles > 0)
                    {
                      //cek todo, this needs to be fixed for not exact
                      double level_set_normal[nSpace];
                      double sign=0.0;
                      if (gf_s.useExact)
                        {
                          double norm_exact=0.0,norm_cut=0.0;
                          if (use_ball_as_particle)
                            {
                              for (int I=0;I<nSpace;I++)
                                {
                                  sign += ball_n[I]*gf_s.get_normal()[I];
                                  level_set_normal[I] = gf_s.get_normal()[I];
                                  norm_cut += level_set_normal[I]*level_set_normal[I];
                                  norm_exact += ball_n[I]*ball_n[I];
                                }
                            }
                          else
                            {
                              for (int I=0;I<nSpace;I++)
                                {
                                  sign += particle_signed_distance_normals.data()[eN_k_3d+I]*gf_s.get_normal()[I];
                                  level_set_normal[I] = gf_s.get_normal()[I];
                                  norm_cut += level_set_normal[I]*level_set_normal[I];
                                  norm_exact += particle_signed_distance_normals.data()[eN_k_3d+I]*particle_signed_distance_normals.data()[eN_k_3d+I];
                                }
                            }
			  norm_cut = std::sqrt(norm_cut);
			  norm_exact = std::sqrt(norm_exact);
                          assert(std::fabs(1.0-norm_cut) < 1.0e-8);
                          assert(std::fabs(1.0-norm_exact) < 1.0e-8);
                          if (sign < 0.0)
                            for (int I=0;I<nSpace;I++)
                              level_set_normal[I]*=-1.0;
                          /* if(icase_s==0)// && (1.0-sign*sign) > 1.0e-3) */
                          /*   { */
                          /*     std::cout<<"phi normal and cut normal divergent "<<eN<<'\t'<<k<<std::endl; */
                          /*     for (int I=0;I<nSpace;I++) */
                          /*       std::cout<<level_set_normal[I]<<'\t'<<particle_signed_distance_normals[eN_k_3d+I]<<std::endl; */
                          /*   } */
                        }
		      else
			{
			  if (use_ball_as_particle)
			    for (int I=0;I<nSpace;I++)
			      level_set_normal[I] = ball_n[I];
			  else
			    for (int I=0;I<nSpace;I++)
			      level_set_normal[I] = particle_signed_distance_normals.data()[eN_k_3d+I];
			}
                      updateSolidParticleTerms(NONCONSERVATIVE_FORM,
                                               eN < nElements_owned,
                                               particle_nitsche,
                                               dV,
                                               nParticles,
                                               nQuadraturePoints_global,
                                               &particle_signed_distances.data()[eN_k],
                                               level_set_normal,
                                               &particle_velocities.data()[eN_k_3d],
                                               particle_centroids.data(),
                                               use_ball_as_particle,
                                               ball_center.data(),
                                               ball_radius.data(),
                                               ball_velocity.data(),
                                               ball_angular_velocity.data(),
                                               ball_density.data(),
                                               porosity,
                                               particle_penalty_constant/h_phi,//penalty,
                                               particle_alpha,
                                               particle_beta,
                                               eps_rho,
                                               eps_mu,
                                               rho_0,
                                               nu_0,
                                               rho_1,
                                               nu_1,
                                               useVF,
                                               vf.data()[eN_k],
                                               phi.data()[eN_k],
                                               x,
                                               y,
                                               z,
                                               p,
                                               u,
                                               v,
                                               w,
                                               q_velocity_sge.data()[eN_k_nSpace+0],
                                               q_velocity_sge.data()[eN_k_nSpace+1],
                                               q_velocity_sge.data()[eN_k_nSpace+1],//dummy entry for 2D
                                               particle_eps,
                                               grad_u,
                                               grad_v,
                                               grad_w,
                                               mass_source_s,
                                               mom_u_source_s,
                                               mom_v_source_s,
                                               mom_w_source_s,
                                               dmom_u_source_s,
                                               dmom_v_source_s,
                                               dmom_w_source_s,
                                               mom_u_adv_s,
                                               mom_v_adv_s,
                                               mom_w_adv_s,
                                               dmom_u_adv_u_s,
                                               dmom_v_adv_v_s,
                                               dmom_w_adv_w_s,
                                               mom_u_ham_s,
                                               dmom_u_ham_grad_u_s,
                                               dmom_u_ham_grad_v_s,
                                               dmom_u_ham_u_s,
                                               dmom_u_ham_v_s,
                                               dmom_u_ham_w_s,
                                               mom_v_ham_s,
                                               dmom_v_ham_grad_u_s,
                                               dmom_v_ham_grad_v_s,
                                               dmom_v_ham_u_s,
                                               dmom_v_ham_v_s,
                                               dmom_v_ham_w_s,
                                               mom_w_ham_s,
                                               dmom_w_ham_grad_w_s,
                                               dmom_w_ham_u_s,
                                               dmom_w_ham_v_s,
                                               dmom_w_ham_w_s,
                                               mass_ham_s,
                                               dmass_ham_u_s,
                                               dmass_ham_v_s,
                                               dmass_ham_w_s,
                                               particle_netForces.data(),
                                               particle_netMoments.data(),
                                               particle_surfaceArea.data());
                    }
                  //
                  //save momentum for time history and velocity for subgrid error
                  //
                  //cek this needs to go with the particle term updates if moved--or check particle_velocities[...] array
                  //cek on cut cells this is getting set twice. For now it's identical because of our formulations (neither includes density so it's either velocity or porosity*velocity--same for both phases
                  //cek but this won't be right when we use a modified basis because we'll need the whole phase velocity from its vasis. hmm. special backward euler class that tracks both?
                  //same situation with subgrid error velocity
                  if (element_active)
                    {
                      q_mom_u_acc.data()[eN_k] = mom_u_acc;
                      q_mom_v_acc.data()[eN_k] = mom_v_acc;
                      //subgrid error uses grid scale velocity
                      q_mass_adv.data()[eN_k_nSpace+0] = u;
                      q_mass_adv.data()[eN_k_nSpace+1] = v;
                    }
                  else//use the solid velocity
                    {
                      q_mom_u_acc.data()[eN_k] = particle_velocities.data()[eN_k_3d+0];
                      q_mom_v_acc.data()[eN_k] = particle_velocities.data()[eN_k_3d+1];
                      q_mass_adv.data()[eN_k_nSpace+0] = particle_velocities.data()[eN_k_3d+0];
                      q_mass_adv.data()[eN_k_nSpace+1] = particle_velocities.data()[eN_k_3d+1];
                    }
                  //
                  //update element residual
                  //
                  double mesh_vel[nSpace];
                  mesh_vel[0] = xt;
                  mesh_vel[1] = yt;
                  double H_f=1.0;
                  if (gf.useExact && icase == 0)
                    {
                      if (fluid_phase == 0)
                        H_f = gf.ImH(0.,0.);
                      else
                        H_f = gf.H(0.,0.);
                    }
                  else
                    H_f = 1.0;
                  if (icase == 0)
                    {
                      //std::cout<<"H_f "<<H_f<<" fluid_phase "<<fluid_phase<<" eN "<<eN<<std::endl;
                    }
                  else
                    {
                      assert(H_f == 1);
                    }

		  if ((eN < nElements_owned) && elementIsActive[eN])
		    {
		      domain_volume += H_s*dV*H_f;
		      p_L1 += fabs(p_e)*H_s*dV*H_f;
		      u_L1 += fabs(u_e)*H_s*dV*H_f;
		      v_L1 += fabs(v_e)*H_s*dV*H_f;
		      velocity_L1 += fabs(velocity_e)*H_s*dV*H_f;
		      
		      p_L2 += p_e*p_e*H_s*dV*H_f;
		      u_L2 += u_e*u_e*H_s*dV*H_f;
		      v_L2 += v_e*v_e*H_s*dV*H_f;
		      velocity_L2 += velocity_e*velocity_e*H_s*dV*H_f;
		      p_dv += p*H_s*H_f*dV;
		      pa_dv += q_u_0.data()[eN_k]*H_s*H_f*dV;
		      total_volume+=H_s*H_f*dV;
		      total_surface_area+=D_s*H_f*dV;
		      if (phi_solid.data()[eN_k] >= 0.0)
			{
			  p_LI = fmax(p_LI, fabs(p_e));
			  u_LI = fmax(u_LI, fabs(u_e));
			  v_LI = fmax(v_LI, fabs(v_e));
			  velocity_LI = fmax(velocity_LI, fabs(velocity_e));
			}
		    }
                  for(int i=0;i<nDOF_test_element;i++)
                    {
                      register int i_nSpace=i*nSpace;
                      elementResidual_mesh[i] += H_s*H_f*(ck.Reaction_weak(1.0,p_test_dV[i]) -
                                                          ck.Reaction_weak(1.0,p_test_dV[i]*q_dV_last.data()[eN_k]/dV) -
                                                          ck.Advection_weak(mesh_vel,&p_grad_test_dV[i_nSpace]));
                      elementResidual_p[i] += H_s*H_f*(ck.Advection_weak(mass_adv,&p_grad_test_dV[i_nSpace])
                                                       + ck.Hamiltonian_weak(mass_ham, p_test_dV[i])
                                                       + DM*MOVING_DOMAIN*(ck.Reaction_weak(alphaBDF*1.0,p_test_dV[i]) -
                                                                           ck.Reaction_weak(alphaBDF*1.0,p_test_dV[i]*q_dV_last.data()[eN_k]/dV) -
                                                                           ck.Advection_weak(mesh_vel,&p_grad_test_dV[i_nSpace])) +
                                                       ck.Reaction_weak(mass_source,p_test_dV[i]));
                      if (nDOF_test_element == nDOF_v_test_element)
                        {
                          elementResidual_p[i] +=
                            H_s*H_f*(PRESSURE_PROJECTION_STABILIZATION * ck.pressureProjection_weak(mom_uu_diff_ten[1], p, p_element_avg, p_test_ref.data()[k*nDOF_test_element+i], dV) +
                                     (1 - PRESSURE_PROJECTION_STABILIZATION) * ck.SubgridError(subgridError_u,Lstar_u_p[i]) +
                                     (1 - PRESSURE_PROJECTION_STABILIZATION) * ck.SubgridError(subgridError_v,Lstar_v_p[i]));
                        }
                      if (PRESSURE_PROJECTION_STABILIZATION==1. && mom_uu_diff_ten[1]==0.)
                        {
                          printf("Warning the Bochev-Dohrnmann-Gunzburger stabilization cannot be applied to inviscid fluids.");
                        }
                      if (nParticles > 0)//solid boundary terms
                        {
                          if (gf_s.D(0.,0.) == 0.0)
                            assert(mass_source_s == 0.0);
                          elementResidual_p[i] += H_f*(ck.Reaction_weak(mass_source_s,p_test_dV[i]));
                        }
                    }
                  for(int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int i_nSpace=i*nSpace;
                      elementResidual_u[i] += H_s*H_f*(ck.Mass_weak(mom_u_acc_t,vel_test_dV[i]) +
                                                       ck.Advection_weak(mom_u_adv,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Diffusion_weak(sdInfo_u_u_rowptr.data(),sdInfo_u_u_colind.data(),mom_uu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Diffusion_weak(sdInfo_u_v_rowptr.data(),sdInfo_u_v_colind.data(),mom_uv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Reaction_weak(mom_u_source+NONCONSERVATIVE_FORM*dmom_u_acc_u*u*div_mesh_velocity,vel_test_dV[i]) +
                                                       ck.Hamiltonian_weak(mom_u_ham,vel_test_dV[i]) +
                                                       MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridError(subgridError_u,Lstar_u_u[i]) +
                                                       ck.NumericalDiffusion(q_numDiff_u_last.data()[eN_k],grad_u,&vel_grad_test_dV[i_nSpace]));
                      elementResidual_v[i] += H_s*H_f*(ck.Mass_weak(mom_v_acc_t,vel_test_dV[i]) +
                                                       ck.Advection_weak(mom_v_adv,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Diffusion_weak(sdInfo_v_u_rowptr.data(),sdInfo_v_u_colind.data(),mom_vu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Diffusion_weak(sdInfo_v_v_rowptr.data(),sdInfo_v_v_colind.data(),mom_vv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Reaction_weak(mom_v_source+NONCONSERVATIVE_FORM*dmom_v_acc_v*v*div_mesh_velocity,vel_test_dV[i]) +
                                                       ck.Hamiltonian_weak(mom_v_ham,vel_test_dV[i]) +
                                                       MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridError(subgridError_v,Lstar_v_v[i]) +
                                                       ck.NumericalDiffusion(q_numDiff_v_last.data()[eN_k],grad_v,&vel_grad_test_dV[i_nSpace]));
                      elementResidual_u[i] +=  H_s*H_f*MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridError(subgridError_p,Lstar_p_u[i]);
                      elementResidual_v[i] +=  H_s*H_f*MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridError(subgridError_p,Lstar_p_v[i]);                      
                      if (nParticles > 0)//solid boundary terms
                        {
                          elementResidual_u[i] += H_f*(ck.Advection_weak(mom_u_adv_s,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Reaction_weak(mom_u_source_s,vel_test_dV[i]) +
                                                       ck.Hamiltonian_weak(mom_u_ham_s,vel_test_dV[i]));
                          elementResidual_v[i] += H_f*(ck.Advection_weak(mom_v_adv_s,&vel_grad_test_dV[i_nSpace]) +
                                                       ck.Reaction_weak(mom_v_source_s,vel_test_dV[i]) +
                                                       ck.Hamiltonian_weak(mom_v_ham_s,vel_test_dV[i]));
                        }
                    }//i
                  //estimate the numerical viscosity combining shock capturing and VMS/SUPG
                  numerical_viscosity.data()[eN_k] = q_numDiff_u_last.data()[eN_k] + MOMENTUM_SGE*VELOCITY_SGE*tau_v*(dmom_adv_star[0]*dmom_adv_star[0]+
                                                                                                                      dmom_adv_star[1]*dmom_adv_star[1]);
                  if (!elementIsActive[eN])
                    {
                      assert(std::fabs(gf_s.H(particle_eps,phi_solid.data()[eN_k])) == 0.0);
                      assert(std::fabs(gf_s.D(particle_eps,phi_solid.data()[eN_k])) == 0.0);
                    }
                }//k
            }//fluid_phase
#ifdef MAXNUMDIFF
          for(int k=0;k<nQuadraturePoints_element;k++)
            {
              //compute indices and declare local storage
              register int eN_k = eN*nQuadraturePoints_element+k;
              q_numDiff_u.data()[eN_k] = numDiffMax;
              q_numDiff_v.data()[eN_k] = numDiffMax;
              q_numDiff_w.data()[eN_k] = numDiffMax;
            }
#endif
          //
          //load element into global residual and save element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              elementResidual_p_save.data()[eN_i] +=  elementResidual_p[i];
              mesh_volume_conservation_element_weak += elementResidual_mesh[i];
              if (!elementIsActive[eN])
                {
                  assert(elementResidual_p[i]==0.0);
                }
              globalResidual.data()[offset_p+stride_p*rp_l2g.data()[eN_i]]+=elementResidual_p[i];
              if (element_active)
		{
		  isActiveR.data()[offset_p+stride_p*rp_l2g.data()[eN_i]] = 1.0;
		  isActiveDOF_p.data()[p_l2g.data()[eN_i]] = 1.0;
		}
	    }
          for(int i=0;i<nDOF_v_test_element;i++)
            {
              register int eN_i=eN*nDOF_v_test_element+i;
              if (!elementIsActive[eN])
                {
                  assert(elementResidual_u[i]==0.0);
                  assert(elementResidual_v[i]==0.0);
                }
              globalResidual.data()[offset_u+stride_u*rvel_l2g.data()[eN_i]]+=elementResidual_u[i];
              globalResidual.data()[offset_v+stride_v*rvel_l2g.data()[eN_i]]+=elementResidual_v[i];
              if (element_active)
                {
                  isActiveR.data()[offset_u+stride_u*rvel_l2g.data()[eN_i]] = 1.0;
                  isActiveR.data()[offset_v+stride_v*rvel_l2g.data()[eN_i]] = 1.0;
                  isActiveDOF_vel.data()[vel_l2g.data()[eN_i]] = 1.0;
                }
            }//i
          mesh_volume_conservation += mesh_volume_conservation_element;
          mesh_volume_conservation_weak += mesh_volume_conservation_element_weak;
          mesh_volume_conservation_err_max=fmax(mesh_volume_conservation_err_max,fabs(mesh_volume_conservation_element));
          mesh_volume_conservation_err_max_weak=fmax(mesh_volume_conservation_err_max_weak,fabs(mesh_volume_conservation_element_weak));
        }//elements
      //std::cout<<"p,u,v L2 error integrals (shoudl be non-negative) "<<p_L2<<'\t'<<u_L2<<'\t'<<v_L2<<'\t'<<"Flow Domain Volume = "<<domain_volume<<std::endl;
      for (std::set<int>::iterator it=cutfem_boundaries.begin(); it!=cutfem_boundaries.end(); )
        {
          if(elementIsActive[elementBoundaryElementsArray[(*it)*2+0]] && elementIsActive[elementBoundaryElementsArray[(*it)*2+1]])
            {
              std::map<int,double> DWp_Dn_jump, DW_Dn_jump;
              register double gamma_cutfem=ghost_penalty_constant,gamma_cutfem_p=ghost_penalty_constant,h_cutfem=elementBoundaryDiameter.data()[*it];
              int eN_nDOF_v_trial_element  = elementBoundaryElementsArray.data()[(*it)*2+0]*nDOF_v_trial_element;
              //See Massing Schott Wall 2018
              //cek todo modify for two-fluids: rho_0 != rho_1
              double norm_v=0.0;
              for (int i=0;i<nDOF_v_trial_element;i++)//MSW18 is just on face, but this is easier
                {
                  double u=u_old_dof.data()[vel_l2g.data()[eN_nDOF_v_trial_element+i]],
                    v=v_old_dof.data()[vel_l2g.data()[eN_nDOF_v_trial_element+i]];
                  norm_v=fmax(norm_v,sqrt(u*u+v*v));
                }
              double gamma_v_dim = rho_0*(nu_0 + norm_v*h_cutfem + alphaBDF*h_cutfem*h_cutfem);
              gamma_cutfem_p *= h_cutfem*h_cutfem/gamma_v_dim;
              if (NONCONSERVATIVE_FORM)
                gamma_cutfem*=gamma_v_dim;
              else
                gamma_cutfem*=(gamma_v_dim/rho_0);
              for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                {
                  register double Dp_Dn_jump=0.0, Du_Dn_jump=0.0, Dv_Dn_jump=0.0,dS;
                  for (int eN_side=0;eN_side < 2; eN_side++)
                    {
                      register int ebN = *it,
                        eN  = elementBoundaryElementsArray.data()[ebN*2+eN_side];
                      for (int i=0;i<nDOF_test_element;i++)
                        {
                          DWp_Dn_jump[rp_l2g.data()[eN*nDOF_test_element+i]] = 0.0;
                        }
                      for (int i=0;i<nDOF_v_test_element;i++)
                        {
                          DW_Dn_jump[rvel_l2g.data()[eN*nDOF_v_test_element+i]] = 0.0;
                        }
                    }
                  for (int eN_side=0;eN_side < 2; eN_side++)
                    {
                      register int ebN = *it,
                        eN  = elementBoundaryElementsArray[ebN*2+eN_side],
                        ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+eN_side],
                        eN_nDOF_trial_element = eN*nDOF_trial_element,
                        eN_nDOF_v_trial_element = eN*nDOF_v_trial_element,
                        ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                        ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                      register double p_int=0.0,
                        u_int=0.0,
                        v_int=0.0,
                        grad_p_int[nSpace]=ZEROVEC,
                        grad_u_int[nSpace]=ZEROVEC,
                        grad_v_int[nSpace]=ZEROVEC,
                        jac_int[nSpace*nSpace],
                        jacDet_int,
                        jacInv_int[nSpace*nSpace],
                        boundaryJac[nSpace*(nSpace-1)],
                        metricTensor[(nSpace-1)*(nSpace-1)],
                        metricTensorDetSqrt,
                        p_test_dS[nDOF_test_element],vel_test_dS[nDOF_v_test_element],
                        p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_v_trial_element*nSpace],
                        p_grad_test_dS[nDOF_trial_element*nSpace],vel_grad_test_dS[nDOF_v_trial_element*nSpace],
                        normal[nSpace],x_int,y_int,z_int,xt_int,yt_int,zt_int,integralScaling,
                        G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,
                        force_x,force_y,force_z,force_p_x,force_p_y,force_p_z,force_v_x,force_v_y,force_v_z,r_x,r_y,r_z;
                      //compute information about mapping from reference element to physical element
                      ck.calculateMapping_elementBoundary(eN,
                                                          ebN_local,
                                                          kb,
                                                          ebN_local_kb,
                                                          mesh_dof.data(),
                                                          mesh_l2g.data(),
                                                          mesh_trial_trace_ref.data(),
                                                          mesh_grad_trial_trace_ref.data(),
                                                          boundaryJac_ref.data(),
                                                          jac_int,
                                                          jacDet_int,
                                                          jacInv_int,
                                                          boundaryJac,
                                                          metricTensor,
                                                          metricTensorDetSqrt,
                                                          normal_ref.data(),
                                                          normal,
                                                          x_int,y_int,z_int);
                      //todo: check that physical coordinates match
                      ck.calculateMappingVelocity_elementBoundary(eN,
                                                                  ebN_local,
                                                                  kb,
                                                                  ebN_local_kb,
                                                                  mesh_velocity_dof.data(),
                                                                  mesh_l2g.data(),
                                                                  mesh_trial_trace_ref.data(),
                                                                  xt_int,yt_int,zt_int,
                                                                  normal,
                                                                  boundaryJac,
                                                                  metricTensor,
                                                                  integralScaling);
                      dS = metricTensorDetSqrt*dS_ref.data()[kb];
                      //compute shape and solution information
                      //shape
                      ck.gradTrialFromRef(&p_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_int,p_grad_trial_trace);
                      ck_v.gradTrialFromRef(&vel_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_v_trial_element],jacInv_int,vel_grad_trial_trace);
                      //solution and gradients
                      ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],p_int);
                      ck_v.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],u_int);
                      ck_v.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],v_int);
                      ck.gradFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_int);
                      ck_v.gradFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_u_int);
                      ck_v.gradFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_v_int);
                      for (int I=0;I<nSpace;I++)
                        {
                          Dp_Dn_jump += grad_p_int[I]*normal[I];
                          Du_Dn_jump += grad_u_int[I]*normal[I];
                          Dv_Dn_jump += grad_v_int[I]*normal[I];
                        }
                      for (int i=0;i<nDOF_test_element;i++)
                        {
                          for (int I=0;I<nSpace;I++)
                            DWp_Dn_jump[rp_l2g[eN_nDOF_trial_element+i]] += p_grad_trial_trace[i*nSpace+I]*normal[I];
                        }
                      for (int i=0;i<nDOF_v_test_element;i++)
                        {
                          for (int I=0;I<nSpace;I++)
                            DW_Dn_jump[rvel_l2g[eN_nDOF_v_trial_element+i]] += vel_grad_trial_trace[i*nSpace+I]*normal[I];
                        }
                    }//eN_side
                  for (std::map<int,double>::iterator W_it=DWp_Dn_jump.begin(); W_it!=DWp_Dn_jump.end(); ++W_it)
                    {
                      int i_global = W_it->first;
                      double DWp_Dn_jump_i = W_it->second;
                      globalResidual.data()[offset_p+stride_p*i_global]+=gamma_cutfem_p*h_cutfem*Dp_Dn_jump*DWp_Dn_jump_i*dS;
                    }
                  for (std::map<int,double>::iterator W_it=DW_Dn_jump.begin(); W_it!=DW_Dn_jump.end(); ++W_it)
                    {
                      int i_global = W_it->first;
                      double DW_Dn_jump_i = W_it->second;
                      globalResidual.data()[offset_u+stride_u*i_global]+=gamma_cutfem*h_cutfem*Du_Dn_jump*DW_Dn_jump_i*dS;
                      globalResidual.data()[offset_v+stride_v*i_global]+=gamma_cutfem*h_cutfem*Dv_Dn_jump*DW_Dn_jump_i*dS;
                    }//i
                }//kb
              ++it;
            }
          else
            {
              it = cutfem_boundaries.erase(it);
            }
        }//cutfem element boundaries
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
        {
          register int ebN = exteriorElementBoundariesArray.data()[ebNE],
            eN  = elementBoundaryElementsArray.data()[ebN*2+0],
            ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
            eN_nDOF_trial_element = eN*nDOF_trial_element,
            eN_nDOF_v_trial_element = eN*nDOF_v_trial_element;
	  if (boundaryFlags[ebN] < 1)
	    continue;
          register double elementResidual_mesh[nDOF_test_element],
            elementResidual_p[nDOF_test_element],
            elementResidual_u[nDOF_v_test_element],
            elementResidual_v[nDOF_v_test_element],
            eps_rho,eps_mu;
          const double* elementResidual_w(NULL);
          for (int i=0;i<nDOF_test_element;i++)
            {
              elementResidual_mesh[i]=0.0;
              elementResidual_p[i]=0.0;
            }
          for (int i=0;i<nDOF_v_test_element;i++)
            {
              elementResidual_u[i]=0.0;
              elementResidual_v[i]=0.0;
            }
          double element_phi[nDOF_mesh_trial_element], element_phi_s[nDOF_mesh_trial_element];
          for (int j=0;j<nDOF_mesh_trial_element;j++)
            {
              register int eN_j = eN*nDOF_mesh_trial_element+j;
              element_phi[j] = phi_nodes.data()[p_l2g.data()[eN_j]];
              element_phi_s[j] = phi_solid_nodes[p_l2g.data()[eN_j]];
            }
          double element_nodes[nDOF_mesh_trial_element*3];
          for (int i=0;i<nDOF_mesh_trial_element;i++)
            {
              register int eN_i=eN*nDOF_mesh_trial_element+i;
              for(int I=0;I<3;I++)
                element_nodes[i*3 + I] = mesh_dof[mesh_l2g.data()[eN_i]*3 + I];
            }//i
          double mesh_dof_ref[nDOF_mesh_trial_element*3]={0.,0.,0.,1.,0.,0.,0.,1.,0.};
          double xb_ref_calc[nQuadraturePoints_elementBoundary*3];
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              double x=0.0,y=0.0,z=0.0;
              for (int j=0;j<nDOF_mesh_trial_element;j++)
                {
                  int ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
                  int ebN_local_kb_j = ebN_local_kb*nDOF_mesh_trial_element+j;
                  x += mesh_dof_ref[j*3+0]*mesh_trial_trace_ref.data()[ebN_local_kb_j]; 
                  y += mesh_dof_ref[j*3+1]*mesh_trial_trace_ref.data()[ebN_local_kb_j]; 
                  z += mesh_dof_ref[j*3+2]*mesh_trial_trace_ref.data()[ebN_local_kb_j];
                }
              xb_ref_calc[3*kb+0] = x;
              xb_ref_calc[3*kb+1] = y;
              xb_ref_calc[3*kb+2] = z;
            }
          int icase_s = gf_s.calculate(element_phi_s, element_nodes, xb_ref_calc, true);
#ifdef IFEM
          int icase = gf.calculate(element_phi, element_nodes, xb_ref.data(), -rho_1*g.data()[1], -rho_0*g.data()[1],true,true);
#else
          int icase = gf.calculate(element_phi, element_nodes, xb_ref.data(), 1.0,1.0,true,false);
#endif
          //cek todo needs modification for twophase flow ibm
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                ebNE_kb_nSpace = ebNE_kb*nSpace,
                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                ebN_local_kb_nSpace = ebN_local_kb*nSpace;
              register double phi_s_ext=0.0,
                p_ext=0.0,
                u_ext=0.0,
                v_ext=0.0,
                w_ext=0.0,
                grad_p_ext[nSpace]=ZEROVEC,
                grad_u_ext[nSpace]=ZEROVEC,
                grad_v_ext[nSpace]=ZEROVEC,
                grad_w_ext[nSpace]=ZEROVEC,
                p_old=0.0,u_old=0.0,v_old=0.0,w_old=0.0,
                grad_p_old[nSpace]=ZEROVEC,grad_u_old[nSpace]=ZEROVEC,grad_v_old[nSpace]=ZEROVEC,grad_w_old[nSpace]=ZEROVEC,
                mom_u_acc_ext=0.0,
                dmom_u_acc_u_ext=0.0,
                mom_v_acc_ext=0.0,
                dmom_v_acc_v_ext=0.0,
                mom_w_acc_ext=0.0,
                dmom_w_acc_w_ext=0.0,
                mass_adv_ext[nSpace]=ZEROVEC,
                dmass_adv_u_ext[nSpace]=ZEROVEC,
                dmass_adv_v_ext[nSpace]=ZEROVEC,
                dmass_adv_w_ext[nSpace]=ZEROVEC,
                mom_u_adv_ext[nSpace]=ZEROVEC,
                dmom_u_adv_u_ext[nSpace]=ZEROVEC,
                dmom_u_adv_v_ext[nSpace]=ZEROVEC,
                dmom_u_adv_w_ext[nSpace]=ZEROVEC,
                mom_v_adv_ext[nSpace]=ZEROVEC,
                dmom_v_adv_u_ext[nSpace]=ZEROVEC,
                dmom_v_adv_v_ext[nSpace]=ZEROVEC,
                dmom_v_adv_w_ext[nSpace]=ZEROVEC,
                mom_w_adv_ext[nSpace]=ZEROVEC,
                dmom_w_adv_u_ext[nSpace]=ZEROVEC,
                dmom_w_adv_v_ext[nSpace]=ZEROVEC,
                dmom_w_adv_w_ext[nSpace]=ZEROVEC,
                mom_uu_diff_ten_ext[nSpace]=ZEROVEC,
                mom_vv_diff_ten_ext[nSpace]=ZEROVEC,
                mom_ww_diff_ten_ext[nSpace]=ZEROVEC,
                mom_uv_diff_ten_ext[1],
                mom_uw_diff_ten_ext[1],
                mom_vu_diff_ten_ext[1],
                mom_vw_diff_ten_ext[1],
                mom_wu_diff_ten_ext[1],
                mom_wv_diff_ten_ext[1],
                mom_u_source_ext=0.0,
                mom_v_source_ext=0.0,
                mom_w_source_ext=0.0,
                mom_u_ham_ext=0.0,
                dmom_u_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_u_ham_grad_u_ext[nSpace]=ZEROVEC,
                dmom_u_ham_u_ext=0.0,
                dmom_u_ham_v_ext=0.0,
                dmom_u_ham_w_ext=0.0,
                mom_v_ham_ext=0.0,
                dmom_v_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_v_ham_grad_v_ext[nSpace]=ZEROVEC,
                dmom_v_ham_u_ext=0.0,
                dmom_v_ham_v_ext=0.0,
                dmom_v_ham_w_ext=0.0,
                mom_w_ham_ext=0.0,
                dmom_w_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_w_ham_grad_w_ext[nSpace]=ZEROVEC,
                dmom_w_ham_u_ext=0.0,
                dmom_w_ham_v_ext=0.0,
                dmom_w_ham_w_ext=0.0,
                dmom_u_adv_p_ext[nSpace]=ZEROVEC,
                dmom_v_adv_p_ext[nSpace]=ZEROVEC,
                dmom_w_adv_p_ext[nSpace]=ZEROVEC,
                flux_mass_ext=0.0,
                flux_mom_u_adv_ext=0.0,
                flux_mom_v_adv_ext=0.0,
                flux_mom_w_adv_ext=0.0,
                flux_mom_uu_diff_ext=0.0,
                flux_mom_uv_diff_ext=0.0,
                flux_mom_uw_diff_ext=0.0,
                flux_mom_vu_diff_ext=0.0,
                flux_mom_vv_diff_ext=0.0,
                flux_mom_vw_diff_ext=0.0,
                flux_mom_wu_diff_ext=0.0,
                flux_mom_wv_diff_ext=0.0,
                flux_mom_ww_diff_ext=0.0,
                bc_p_ext=0.0,
                bc_u_ext=0.0,
                bc_v_ext=0.0,
                bc_w_ext=0.0,
                bc_mom_u_acc_ext=0.0,
                bc_dmom_u_acc_u_ext=0.0,
                bc_mom_v_acc_ext=0.0,
                bc_dmom_v_acc_v_ext=0.0,
                bc_mom_w_acc_ext=0.0,
                bc_dmom_w_acc_w_ext=0.0,
                bc_mass_adv_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_u_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_v_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_w_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_uu_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_vv_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_ww_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_uv_diff_ten_ext[1],
                bc_mom_uw_diff_ten_ext[1],
                bc_mom_vu_diff_ten_ext[1],
                bc_mom_vw_diff_ten_ext[1],
                bc_mom_wu_diff_ten_ext[1],
                bc_mom_wv_diff_ten_ext[1],
                bc_mom_u_source_ext=0.0,
                bc_mom_v_source_ext=0.0,
                bc_mom_w_source_ext=0.0,
                bc_mom_u_ham_ext=0.0,
                bc_dmom_u_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_u_ham_grad_u_ext[nSpace]=ZEROVEC,
                bc_dmom_u_ham_u_ext=0.0,
                bc_dmom_u_ham_v_ext=0.0,
                bc_dmom_u_ham_w_ext=0.0,
                bc_mom_v_ham_ext=0.0,
                bc_dmom_v_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_v_ham_grad_v_ext[nSpace]=ZEROVEC,
                bc_dmom_v_ham_u_ext=0.0,
                bc_dmom_v_ham_v_ext=0.0,
                bc_dmom_v_ham_w_ext=0.0,
                bc_mom_w_ham_ext=0.0,
                bc_dmom_w_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_w_ham_grad_w_ext[nSpace]=ZEROVEC,
                bc_dmom_w_ham_u_ext=0.0,
                bc_dmom_w_ham_v_ext=0.0,
                bc_dmom_w_ham_w_ext=0.0,
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,p_test_dS[nDOF_test_element],vel_test_dS[nDOF_v_test_element],
                p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_v_trial_element*nSpace],
                vel_grad_test_dS[nDOF_v_trial_element*nSpace],
                normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                //VRANS
                porosity_ext,
                //
                G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,
                force_x,force_y,force_z,force_p_x,force_p_y,force_p_z,force_v_x,force_v_y,force_v_z,r_x,r_y,r_z;
              //compute information about mapping from reference element to physical element
              gf_s.set_boundary_quad(kb);
              gf.set_boundary_quad(kb);
              ck.calculateMapping_elementBoundary(eN,
                                                  ebN_local,
                                                  kb,
                                                  ebN_local_kb,
                                                  mesh_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_trace_ref.data(),
                                                  mesh_grad_trial_trace_ref.data(),
                                                  boundaryJac_ref.data(),
                                                  jac_ext,
                                                  jacDet_ext,
                                                  jacInv_ext,
                                                  boundaryJac,
                                                  metricTensor,
                                                  metricTensorDetSqrt,
                                                  normal_ref.data(),
                                                  normal,
                                                  x_ext,y_ext,z_ext);
              ck.calculateMappingVelocity_elementBoundary(eN,
                                                          ebN_local,
                                                          kb,
                                                          ebN_local_kb,
                                                          mesh_velocity_dof.data(),
                                                          mesh_l2g.data(),
                                                          mesh_trial_trace_ref.data(),
                                                          xt_ext,yt_ext,zt_ext,
                                                          normal,
                                                          boundaryJac,
                                                          metricTensor,
                                                          integralScaling);
              //xt_ext=0.0;yt_ext=0.0;zt_ext=0.0;
              //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<zt_ext<<std::endl;
              //std::cout<<"x_ext "<<x_ext<<'\t'<<y_ext<<'\t'<<z_ext<<std::endl;
              //std::cout<<"integralScaling - metricTensorDetSrt ==============================="<<integralScaling-metricTensorDetSqrt<<std::endl;
              /* std::cout<<"metricTensorDetSqrt "<<metricTensorDetSqrt */
              /*             <<"dS_ref.data()[kb]"<<dS_ref.data()[kb]<<std::endl; */
              //dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];//cek need to test effect on accuracy
              dS = metricTensorDetSqrt*dS_ref.data()[kb];
              //get the metric tensor
              //cek todo use symmetry
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              ck.calculateGScale(G,&ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],h_phi);

              eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
              eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);

              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&p_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
              ck_v.gradTrialFromRef(&vel_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_v_trial_element],jacInv_ext,vel_grad_trial_trace);
              //solution and gradients
              ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],p_ext);
              ck_v.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],u_ext);
              ck_v.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],v_ext);
              ck.valFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],p_old);
              ck_v.valFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],u_old);
              ck_v.valFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],v_old);
              ck.gradFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
              ck_v.gradFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_u_ext);
              ck_v.gradFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_v_ext);
              ck.gradFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_old);
              ck_v.gradFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_u_old);
              ck_v.gradFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_v_old);
              ck.valFromDOF(phi_solid_nodes.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],phi_s_ext);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_test_element;j++)
                {
                  p_test_dS[j] = p_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                }
              for (int j=0;j<nDOF_v_test_element;j++)
                {
                  vel_test_dS[j] = vel_test_trace_ref.data()[ebN_local_kb*nDOF_v_test_element+j]*dS;
                  for (int I=0;I<nSpace;I++)
                    vel_grad_test_dS[j*nSpace+I] = vel_grad_trial_trace[j*nSpace+I]*dS;//assume test_j = trial_j
                }
              bc_p_ext = isDOFBoundary_p.data()[ebNE_kb]*ebqe_bc_p_ext.data()[ebNE_kb]+(1-isDOFBoundary_p.data()[ebNE_kb])*p_ext;
              //note, our convention is that bc values at moving boundaries are relative to boundary velocity so we add it here
              bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*(ebqe_bc_u_ext.data()[ebNE_kb] + MOVING_DOMAIN*xt_ext) + (1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;
              bc_v_ext = isDOFBoundary_v.data()[ebNE_kb]*(ebqe_bc_v_ext.data()[ebNE_kb] + MOVING_DOMAIN*yt_ext) + (1-isDOFBoundary_v.data()[ebNE_kb])*v_ext;
              //VRANS
              porosity_ext = ebqe_porosity_ext.data()[ebNE_kb];
              //
              //calculate the pde coefficients using the solution and the boundary values for the solution
              //
              double eddy_viscosity_ext(0.),bc_eddy_viscosity_ext(0.); //not interested in saving boundary eddy viscosity for now
              if (use_ball_as_particle == 1 && nParticles > 0)
                {
                  get_distance_to_ball(nParticles, ball_center.data(), ball_radius.data(),x_ext,y_ext,z_ext,ebqe_phi_s.data()[ebNE_kb]);
                }
              //else ebqe_phi_s.data()[ebNE_kb] is computed in Prestep
              const double particle_eps  = particle_epsFact*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);

              //cek needs to be fixed for two-phase ifem
              double H = (1.0-useVF)*gf.H(eps_rho,ebqe_phi_ext[ebNE_kb]) + useVF*fmin(1.0,fmax(0.0,ebqe_vf_ext[ebNE_kb]));
              double ImH = (1.0-useVF)*gf.ImH(eps_rho,ebqe_phi_ext[ebNE_kb]) + useVF*(1.0-fmin(1.0,fmax(0.0,ebqe_vf_ext[ebNE_kb])));
              double rho  = rho_0*ImH + rho_1*H;
              double nu  = nu_0*ImH + nu_1*H;
              //
              evaluateCoefficients(NONCONSERVATIVE_FORM,
                                   sigma,
                                   rho,
                                   nu,
                                   elementDiameter.data()[eN],
                                   smagorinskyConstant,
                                   turbulenceClosureModel,
                                   g.data(),
                                   useVF,
                                   ebqe_vf_ext.data()[ebNE_kb],
                                   ebqe_phi_ext.data()[ebNE_kb],
                                   &ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],
                                   ebqe_kappa_phi_ext.data()[ebNE_kb],
                                   //VRANS
                                   porosity_ext,
                                   //
                                   ebqe_phi_s.data()[ebNE_kb],
                                   p_old,
                                   u_old,
                                   v_old,
                                   w_old,
                                   grad_p_old,
                                   grad_u_old,
                                   grad_v_old,
                                   grad_w_old,
                                   p_ext,
                                   grad_p_ext,
                                   grad_u_ext,
                                   grad_v_ext,
                                   grad_w_ext,
                                   u_ext,
                                   v_ext,
                                   w_ext,
                                   LAG_LES,
                                   ebqe_eddy_viscosity.data()[ebNE_kb],
                                   ebqe_eddy_viscosity_last.data()[ebNE_kb],
                                   mom_u_acc_ext,
                                   dmom_u_acc_u_ext,
                                   mom_v_acc_ext,
                                   dmom_v_acc_v_ext,
                                   mom_w_acc_ext,
                                   dmom_w_acc_w_ext,
                                   mass_adv_ext,
                                   dmass_adv_u_ext,
                                   dmass_adv_v_ext,
                                   dmass_adv_w_ext,
                                   mom_u_adv_ext,
                                   dmom_u_adv_u_ext,
                                   dmom_u_adv_v_ext,
                                   dmom_u_adv_w_ext,
                                   mom_v_adv_ext,
                                   dmom_v_adv_u_ext,
                                   dmom_v_adv_v_ext,
                                   dmom_v_adv_w_ext,
                                   mom_w_adv_ext,
                                   dmom_w_adv_u_ext,
                                   dmom_w_adv_v_ext,
                                   dmom_w_adv_w_ext,
                                   mom_uu_diff_ten_ext,
                                   mom_vv_diff_ten_ext,
                                   mom_ww_diff_ten_ext,
                                   mom_uv_diff_ten_ext,
                                   mom_uw_diff_ten_ext,
                                   mom_vu_diff_ten_ext,
                                   mom_vw_diff_ten_ext,
                                   mom_wu_diff_ten_ext,
                                   mom_wv_diff_ten_ext,
                                   mom_u_source_ext,
                                   mom_v_source_ext,
                                   mom_w_source_ext,
                                   mom_u_ham_ext,
                                   dmom_u_ham_grad_p_ext,
                                   dmom_u_ham_grad_u_ext,
                                   dmom_u_ham_u_ext,
                                   dmom_u_ham_v_ext,
                                   dmom_u_ham_w_ext,
                                   mom_v_ham_ext,
                                   dmom_v_ham_grad_p_ext,
                                   dmom_v_ham_grad_v_ext,
                                   dmom_v_ham_u_ext,
                                   dmom_v_ham_v_ext,
                                   dmom_v_ham_w_ext,
                                   mom_w_ham_ext,
                                   dmom_w_ham_grad_p_ext,
                                   dmom_w_ham_grad_w_ext,
                                   dmom_w_ham_u_ext,
                                   dmom_w_ham_v_ext,
                                   dmom_w_ham_w_ext,
                                   0.0,
                                   0.0,
                                   0.0);
              //cek needs to be fixed for two-phase ifem
              H = (1.0-useVF)*gf.H(eps_rho,bc_ebqe_phi_ext[ebNE_kb]) + useVF*fmin(1.0,fmax(0.0,bc_ebqe_vf_ext[ebNE_kb]));
              ImH = (1.0-useVF)*gf.ImH(eps_rho,bc_ebqe_phi_ext[ebNE_kb]) + useVF*(1.0-fmin(1.0,fmax(0.0,bc_ebqe_vf_ext[ebNE_kb])));
              rho  = rho_0*ImH + rho_1*H;
              nu  = nu_0*ImH + nu_1*H;
              //
              evaluateCoefficients(NONCONSERVATIVE_FORM,
                                   sigma,
                                   rho,
                                   nu,
                                   elementDiameter.data()[eN],
                                   smagorinskyConstant,
                                   turbulenceClosureModel,
                                   g.data(),
                                   useVF,
                                   bc_ebqe_vf_ext.data()[ebNE_kb],
                                   bc_ebqe_phi_ext.data()[ebNE_kb],
                                   &ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],
                                   ebqe_kappa_phi_ext.data()[ebNE_kb],
                                   //VRANS
                                   porosity_ext,
                                   //
                                   ebqe_phi_s.data()[ebNE_kb],
                                   p_old,
                                   u_old,
                                   v_old,
                                   w_old,
                                   grad_p_old,
                                   grad_u_old,
                                   grad_v_old,
                                   grad_w_old,
                                   bc_p_ext,
                                   grad_p_ext,
                                   grad_u_ext,
                                   grad_v_ext,
                                   grad_w_ext,
                                   bc_u_ext,
                                   bc_v_ext,
                                   bc_w_ext,
                                   LAG_LES,
                                   bc_eddy_viscosity_ext,
                                   ebqe_eddy_viscosity_last.data()[ebNE_kb],
                                   bc_mom_u_acc_ext,
                                   bc_dmom_u_acc_u_ext,
                                   bc_mom_v_acc_ext,
                                   bc_dmom_v_acc_v_ext,
                                   bc_mom_w_acc_ext,
                                   bc_dmom_w_acc_w_ext,
                                   bc_mass_adv_ext,
                                   bc_dmass_adv_u_ext,
                                   bc_dmass_adv_v_ext,
                                   bc_dmass_adv_w_ext,
                                   bc_mom_u_adv_ext,
                                   bc_dmom_u_adv_u_ext,
                                   bc_dmom_u_adv_v_ext,
                                   bc_dmom_u_adv_w_ext,
                                   bc_mom_v_adv_ext,
                                   bc_dmom_v_adv_u_ext,
                                   bc_dmom_v_adv_v_ext,
                                   bc_dmom_v_adv_w_ext,
                                   bc_mom_w_adv_ext,
                                   bc_dmom_w_adv_u_ext,
                                   bc_dmom_w_adv_v_ext,
                                   bc_dmom_w_adv_w_ext,
                                   bc_mom_uu_diff_ten_ext,
                                   bc_mom_vv_diff_ten_ext,
                                   bc_mom_ww_diff_ten_ext,
                                   bc_mom_uv_diff_ten_ext,
                                   bc_mom_uw_diff_ten_ext,
                                   bc_mom_vu_diff_ten_ext,
                                   bc_mom_vw_diff_ten_ext,
                                   bc_mom_wu_diff_ten_ext,
                                   bc_mom_wv_diff_ten_ext,
                                   bc_mom_u_source_ext,
                                   bc_mom_v_source_ext,
                                   bc_mom_w_source_ext,
                                   bc_mom_u_ham_ext,
                                   bc_dmom_u_ham_grad_p_ext,
                                   bc_dmom_u_ham_grad_u_ext,
                                   bc_dmom_u_ham_u_ext,
                                   bc_dmom_u_ham_v_ext,
                                   bc_dmom_u_ham_w_ext,
                                   bc_mom_v_ham_ext,
                                   bc_dmom_v_ham_grad_p_ext,
                                   bc_dmom_v_ham_grad_v_ext,
                                   bc_dmom_v_ham_u_ext,
                                   bc_dmom_v_ham_v_ext,
                                   bc_dmom_v_ham_w_ext,
                                   bc_mom_w_ham_ext,
                                   bc_dmom_w_ham_grad_p_ext,
                                   bc_dmom_w_ham_grad_w_ext,
                                   bc_dmom_w_ham_u_ext,
                                   bc_dmom_w_ham_v_ext,
                                   bc_dmom_w_ham_w_ext,
                                   0.0,
                                   0.0,
                                   0.0);

              //Turbulence closure model
              if (turbulenceClosureModel >= 3)
                {
                  const double turb_var_grad_0_dummy[nSpace] = ZEROVEC;
                  const double c_mu = 0.09;//mwf hack
                  updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                          turbulenceClosureModel,
                                          eps_rho,
                                          eps_mu,
                                          rho_0,
                                          nu_0,
                                          rho_1,
                                          nu_1,
                                          useVF,
                                          ebqe_vf_ext.data()[ebNE_kb],
                                          ebqe_phi_ext.data()[ebNE_kb],
                                          porosity_ext,
                                          c_mu, //mwf hack
                                          ebqe_turb_var_0.data()[ebNE_kb],
                                          ebqe_turb_var_1.data()[ebNE_kb],
                                          turb_var_grad_0_dummy, //not needed
                                          ebqe_eddy_viscosity.data()[ebNE_kb],
                                          mom_uu_diff_ten_ext,
                                          mom_vv_diff_ten_ext,
                                          mom_ww_diff_ten_ext,
                                          mom_uv_diff_ten_ext,
                                          mom_uw_diff_ten_ext,
                                          mom_vu_diff_ten_ext,
                                          mom_vw_diff_ten_ext,
                                          mom_wu_diff_ten_ext,
                                          mom_wv_diff_ten_ext,
                                          mom_u_source_ext,
                                          mom_v_source_ext,
                                          mom_w_source_ext);

                  updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                          turbulenceClosureModel,
                                          eps_rho,
                                          eps_mu,
                                          rho_0,
                                          nu_0,
                                          rho_1,
                                          nu_1,
                                          useVF,
                                          bc_ebqe_vf_ext.data()[ebNE_kb],
                                          bc_ebqe_phi_ext.data()[ebNE_kb],
                                          porosity_ext,
                                          c_mu, //mwf hack
                                          ebqe_turb_var_0.data()[ebNE_kb],
                                          ebqe_turb_var_1.data()[ebNE_kb],
                                          turb_var_grad_0_dummy, //not needed
                                          bc_eddy_viscosity_ext,
                                          bc_mom_uu_diff_ten_ext,
                                          bc_mom_vv_diff_ten_ext,
                                          bc_mom_ww_diff_ten_ext,
                                          bc_mom_uv_diff_ten_ext,
                                          bc_mom_uw_diff_ten_ext,
                                          bc_mom_vu_diff_ten_ext,
                                          bc_mom_vw_diff_ten_ext,
                                          bc_mom_wu_diff_ten_ext,
                                          bc_mom_wv_diff_ten_ext,
                                          bc_mom_u_source_ext,
                                          bc_mom_v_source_ext,
                                          bc_mom_w_source_ext);
                }


              //
              //moving domain
              //
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  mom_u_ham_ext -= MOVING_DOMAIN*dmom_u_acc_u_ext*(grad_u_ext[0]*xt_ext + grad_u_ext[1]*yt_ext);
                  dmom_u_ham_grad_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
                  dmom_u_ham_grad_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
                }
              else
                {
                  mom_u_adv_ext[0] -= MOVING_DOMAIN*mom_u_acc_ext*xt_ext;
                  mom_u_adv_ext[1] -= MOVING_DOMAIN*mom_u_acc_ext*yt_ext;
                  dmom_u_adv_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
                  dmom_u_adv_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
                }


              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  mom_v_ham_ext -= MOVING_DOMAIN*dmom_v_acc_v_ext*(grad_v_ext[0]*xt_ext + grad_v_ext[1]*yt_ext);
                  dmom_v_ham_grad_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
                  dmom_v_ham_grad_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
                }
              else
                {
                  mom_v_adv_ext[0] -= MOVING_DOMAIN*mom_v_acc_ext*xt_ext;
                  mom_v_adv_ext[1] -= MOVING_DOMAIN*mom_v_acc_ext*yt_ext;
                  dmom_v_adv_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
                  dmom_v_adv_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
                }

              //bc's
              if (NONCONSERVATIVE_FORM < 1.0)
                {
                  bc_mom_u_adv_ext[0] -= MOVING_DOMAIN*bc_mom_u_acc_ext*xt_ext;
                  bc_mom_u_adv_ext[1] -= MOVING_DOMAIN*bc_mom_u_acc_ext*yt_ext;

                  bc_mom_v_adv_ext[0] -= MOVING_DOMAIN*bc_mom_v_acc_ext*xt_ext;
                  bc_mom_v_adv_ext[1] -= MOVING_DOMAIN*bc_mom_v_acc_ext*yt_ext;
                }
              //
              //calculate the numerical fluxes
              //
              ck.calculateGScale(G,normal,h_penalty);
              penalty = useMetrics*C_b/h_penalty + (1.0-useMetrics)*ebqe_penalty_ext.data()[ebNE_kb];
              exteriorNumericalAdvectiveFlux(NONCONSERVATIVE_FORM,
                                             isDOFBoundary_p.data()[ebNE_kb],
                                             isDOFBoundary_u.data()[ebNE_kb],
                                             isDOFBoundary_v.data()[ebNE_kb],
                                             isDOFBoundary_w.data()[ebNE_kb],
                                             isAdvectiveFluxBoundary_p.data()[ebNE_kb],
                                             isAdvectiveFluxBoundary_u.data()[ebNE_kb],
                                             isAdvectiveFluxBoundary_v.data()[ebNE_kb],
                                             isAdvectiveFluxBoundary_w.data()[ebNE_kb],
                                             dmom_u_ham_grad_p_ext[0],//=1/rho,
                                             bc_dmom_u_ham_grad_p_ext[0],//=1/bc_rho,
                                             normal,
                                             bc_p_ext,
                                             bc_u_ext,
                                             bc_v_ext,
                                             bc_mass_adv_ext,
                                             bc_mom_u_adv_ext,
                                             bc_mom_v_adv_ext,
                                             bc_mom_w_adv_ext,
                                             ebqe_bc_flux_mass_ext.data()[ebNE_kb]+MOVING_DOMAIN*(xt_ext*normal[0]+yt_ext*normal[1]),//BC is relative mass flux
                                             ebqe_bc_flux_mom_u_adv_ext.data()[ebNE_kb],
                                             ebqe_bc_flux_mom_v_adv_ext.data()[ebNE_kb],
                                             ebqe_bc_flux_mom_w_adv_ext.data()[ebNE_kb],
                                             p_ext,
                                             u_ext,
                                             v_ext,
                                             mass_adv_ext,
                                             mom_u_adv_ext,
                                             mom_v_adv_ext,
                                             mom_w_adv_ext,
                                             dmass_adv_u_ext,
                                             dmass_adv_v_ext,
                                             dmass_adv_w_ext,
                                             dmom_u_adv_p_ext,
                                             dmom_u_ham_grad_u_ext,
                                             dmom_u_adv_u_ext,
                                             dmom_u_adv_v_ext,
                                             dmom_u_adv_w_ext,
                                             dmom_v_adv_p_ext,
                                             dmom_v_adv_u_ext,
                                             dmom_v_adv_v_ext,
                                             dmom_v_adv_w_ext,
                                             dmom_w_adv_p_ext,
                                             dmom_w_adv_u_ext,
                                             dmom_w_adv_v_ext,
                                             dmom_w_adv_w_ext,
                                             flux_mass_ext,
                                             flux_mom_u_adv_ext,
                                             flux_mom_v_adv_ext,
                                             flux_mom_w_adv_ext,
                                             &ebqe_velocity.data()[ebNE_kb_nSpace]);
              for (int I=0;I<nSpace;I++)
                ebqe_velocity.data()[ebNE_kb_nSpace+I]/=porosity_ext;
              exteriorNumericalDiffusiveFlux(eps_rho,
                                             ebqe_phi_ext.data()[ebNE_kb],
                                             sdInfo_u_u_rowptr.data(),
                                             sdInfo_u_u_colind.data(),
                                             isDOFBoundary_u.data()[ebNE_kb],
                                             isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                             normal,
                                             bc_mom_uu_diff_ten_ext,
                                             bc_u_ext,
                                             ebqe_bc_flux_u_diff_ext.data()[ebNE_kb],
                                             mom_uu_diff_ten_ext,
                                             grad_u_ext,
                                             u_ext,
                                             penalty,//ebqe_penalty_ext.data()[ebNE_kb],
                                             flux_mom_uu_diff_ext);
              exteriorNumericalDiffusiveFlux(eps_rho,
                                             ebqe_phi_ext.data()[ebNE_kb],
                                             sdInfo_u_v_rowptr.data(),
                                             sdInfo_u_v_colind.data(),
                                             isDOFBoundary_v.data()[ebNE_kb],
                                             isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                             normal,
                                             bc_mom_uv_diff_ten_ext,
                                             bc_v_ext,
                                             0.0,//assume all of the flux gets applied in diagonal component
                                             mom_uv_diff_ten_ext,
                                             grad_v_ext,
                                             v_ext,
                                             penalty,//ebqe_penalty_ext.data()[ebNE_kb],
                                             flux_mom_uv_diff_ext);
              exteriorNumericalDiffusiveFlux(eps_rho,
                                             ebqe_phi_ext.data()[ebNE_kb],
                                             sdInfo_v_u_rowptr.data(),
                                             sdInfo_v_u_colind.data(),
                                             isDOFBoundary_u.data()[ebNE_kb],
                                             isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                             normal,
                                             bc_mom_vu_diff_ten_ext,
                                             bc_u_ext,
                                             0.0,//see above
                                             mom_vu_diff_ten_ext,
                                             grad_u_ext,
                                             u_ext,
                                             penalty,//ebqe_penalty_ext.data()[ebNE_kb],
                                             flux_mom_vu_diff_ext);
              exteriorNumericalDiffusiveFlux(eps_rho,
                                             ebqe_phi_ext.data()[ebNE_kb],
                                             sdInfo_v_v_rowptr.data(),
                                             sdInfo_v_v_colind.data(),
                                             isDOFBoundary_v.data()[ebNE_kb],
                                             isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                             normal,
                                             bc_mom_vv_diff_ten_ext,
                                             bc_v_ext,
                                             ebqe_bc_flux_v_diff_ext.data()[ebNE_kb],
                                             mom_vv_diff_ten_ext,
                                             grad_v_ext,
                                             v_ext,
                                             penalty,//ebqe_penalty_ext.data()[ebNE_kb],
                                             flux_mom_vv_diff_ext);
              flux.data()[ebN*nQuadraturePoints_elementBoundary+kb] = flux_mass_ext;
              /* std::cout<<"external u,v,u_n " */
              /*             <<ebqe_velocity.data()[ebNE_kb_nSpace+0]<<'\t' */
              /*             <<ebqe_velocity.data()[ebNE_kb_nSpace+1]<<'\t' */
              /*             <<flux.data()[ebN*nQuadraturePoints_elementBoundary+kb]<<std::endl; */
              //
              //integrate the net force and moment on flagged boundaries
              //
              if (ebN < nElementBoundaries_owned)
                {
                  force_v_x = (flux_mom_u_adv_ext + flux_mom_uu_diff_ext + flux_mom_uv_diff_ext + flux_mom_uw_diff_ext)/dmom_u_ham_grad_p_ext[0];//same as *rho
                  force_v_y = (flux_mom_v_adv_ext + flux_mom_vu_diff_ext + flux_mom_vv_diff_ext + flux_mom_vw_diff_ext)/dmom_u_ham_grad_p_ext[0];

                  force_p_x = p_ext*normal[0];
                  force_p_y = p_ext*normal[1];

                  force_x = force_p_x + force_v_x;
                  force_y = force_p_y + force_v_y;

                  r_x = x_ext - barycenters.data()[3*boundaryFlags.data()[ebN]+0];
                  r_y = y_ext - barycenters.data()[3*boundaryFlags.data()[ebN]+1];

                  wettedAreas.data()[boundaryFlags.data()[ebN]] += dS*(1.0-ebqe_vf_ext.data()[ebNE_kb]);

                  netForces_p.data()[3*boundaryFlags.data()[ebN]+0] += force_p_x*dS;
                  netForces_p.data()[3*boundaryFlags.data()[ebN]+1] += force_p_y*dS;

                  netForces_v.data()[3*boundaryFlags.data()[ebN]+0] += force_v_x*dS;
                  netForces_v.data()[3*boundaryFlags.data()[ebN]+1] += force_v_y*dS;

                  netMoments.data()[3*boundaryFlags.data()[ebN]+2] += (r_x*force_y - r_y*force_x)*dS;
                }
              //
              //update residuals
              //
              const double H_s = gf_s.H(particle_eps, ebqe_phi_s.data()[ebNE_kb]);
              if (elementIsActive[eN])
                { //if boundary flag positive, then include flux contributions on interpart boundaries
                  total_flux += flux_mass_ext*dS;
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      elementResidual_mesh[i] -= H_s*ck.ExteriorElementBoundaryFlux(MOVING_DOMAIN*(xt_ext*normal[0]+yt_ext*normal[1]),p_test_dS[i]);
                      elementResidual_p[i] += H_s*ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);
                      elementResidual_p[i] -= H_s*DM*ck.ExteriorElementBoundaryFlux(MOVING_DOMAIN*(xt_ext*normal[0]+yt_ext*normal[1]),p_test_dS[i]);
                      globalConservationError += H_s*ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      elementResidual_u[i] += H_s*(ck.ExteriorElementBoundaryFlux(flux_mom_u_adv_ext,vel_test_dS[i])+
                                                   ck.ExteriorElementBoundaryFlux(flux_mom_uu_diff_ext,vel_test_dS[i])+
                                                   ck.ExteriorElementBoundaryFlux(flux_mom_uv_diff_ext,vel_test_dS[i])+
                                                   ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_u.data()[ebNE_kb],
                                                                                              isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                                                              eb_adjoint_sigma,
                                                                                              u_ext,
                                                                                              bc_u_ext,
                                                                                              normal,
                                                                                              sdInfo_u_u_rowptr.data(),
                                                                                              sdInfo_u_u_colind.data(),
                                                                                              mom_uu_diff_ten_ext,
                                                                                              &vel_grad_test_dS[i*nSpace])+
                                                   ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_v.data()[ebNE_kb],
                                                                                              isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                                                              eb_adjoint_sigma,
                                                                                              v_ext,
                                                                                              bc_v_ext,
                                                                                              normal,
                                                                                              sdInfo_u_v_rowptr.data(),
                                                                                              sdInfo_u_v_colind.data(),
                                                                                              mom_uv_diff_ten_ext,
                                                                                              &vel_grad_test_dS[i*nSpace]));
                      elementResidual_v[i] += H_s*(ck.ExteriorElementBoundaryFlux(flux_mom_v_adv_ext,vel_test_dS[i]) +
                                                   ck.ExteriorElementBoundaryFlux(flux_mom_vu_diff_ext,vel_test_dS[i])+
                                                   ck.ExteriorElementBoundaryFlux(flux_mom_vv_diff_ext,vel_test_dS[i])+
                                                   ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_u.data()[ebNE_kb],
                                                                                              isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                                                                              eb_adjoint_sigma,
                                                                                              u_ext,
                                                                                              bc_u_ext,
                                                                                              normal,
                                                                                              sdInfo_v_u_rowptr.data(),
                                                                                              sdInfo_v_u_colind.data(),
                                                                                              mom_vu_diff_ten_ext,
                                                                                              &vel_grad_test_dS[i*nSpace])+
                                                   ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_v.data()[ebNE_kb],
                                                                                              isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                                                                              eb_adjoint_sigma,
                                                                                              v_ext,
                                                                                              bc_v_ext,
                                                                                              normal,
                                                                                              sdInfo_v_v_rowptr.data(),
                                                                                              sdInfo_v_v_colind.data(),
                                                                                              mom_vv_diff_ten_ext,
                                                                                              &vel_grad_test_dS[i*nSpace]));
                    }//i
                }//if boundary flag positive
            }//kb
          //
          //update the element and global residual storage
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              int eN_i = eN*nDOF_test_element+i;

              elementResidual_p_save.data()[eN_i] +=  elementResidual_p[i];
              mesh_volume_conservation_weak += elementResidual_mesh[i];
              globalResidual.data()[offset_p+stride_p*rp_l2g.data()[eN_i]]+=elementResidual_p[i];
            }
          for (int i=0;i<nDOF_v_test_element;i++)
            {
              int eN_i = eN*nDOF_v_test_element+i;
              globalResidual.data()[offset_u+stride_u*rvel_l2g.data()[eN_i]]+=elementResidual_u[i];
              globalResidual.data()[offset_v+stride_v*rvel_l2g.data()[eN_i]]+=elementResidual_v[i];
            }//i
        }//ebNE
      
      if (normalize_pressure)
        {
	  double send[4]={pa_dv,p_dv,total_volume, total_surface_area}, recv[4]={0.,0.,0.,0.};
	  MPI_Allreduce(send, recv,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  pa_dv = recv[0];
	  p_dv = recv[1];
	  total_volume = recv[2];
	  total_surface_area = recv[3];
	  //std::cout<<"Domain Volume: "<<total_volume<<std::endl;
	  //std::cout<<"Domain Surface Area: "<<total_surface_area<<std::endl;
	  //cek hack
	  // 1. This forces the pressure average to match the average of the analytical solution (or zero of no analytical solution is given)
	  // 2. I'm manually figuring out how many pressure dof there are
          /* std::cout<<"mesh volume conservation = "<<mesh_volume_conservation<<std::endl; */
          /* std::cout<<"mesh volume conservation weak = "<<mesh_volume_conservation_weak<<std::endl; */
          /* std::cout<<"mesh volume conservation err max= "<<mesh_volume_conservation_err_max<<std::endl; */
          /* std::cout<<"mesh volume conservation err max weak = "<<mesh_volume_conservation_err_max_weak<<std::endl; */
          /* std::cout<<"Pressure Integral "<<p_dv<<std::endl */
          /*          <<"Analytical Pressure Integral "<<pa_dv<<std::endl */
          /*          <<"Total Boundary Flux "<<total_flux<<std::endl; */
          int nDOF_pressure=0;
          for(int eN=0;eN<nElements_global;eN++)
            {
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  int eN_i = eN*nDOF_test_element+i;
		  if (p_l2g.data()[eN_i] > nDOF_pressure)
		    nDOF_pressure=p_l2g.data()[eN_i];
		}
	    }
	  nDOF_pressure +=1;
	  assert(p_dof.shape(0) == nDOF_pressure);
          //std::cout<<"nDOF_pressure "<<nDOF_pressure<<std::endl;
          for (int I=0;I<nDOF_pressure;I++)
            p_dof.data()[I] += (pa_dv - p_dv)/total_volume;
          double p_dv_new=0.0, pa_dv_new=0.0;
          p_L1=0.0;
          p_L2=0.0;
          p_LI=0.0;
          for (int eN=0 ; eN < nElements_owned ; ++eN)
            {
              double element_phi[nDOF_mesh_trial_element], element_phi_s[nDOF_mesh_trial_element];
              for (int j=0;j<nDOF_mesh_trial_element;j++)
                {
                  register int eN_j = eN*nDOF_mesh_trial_element+j;
                  element_phi[j] = phi_nodes.data()[p_l2g.data()[eN_j]];
                  element_phi_s[j] = phi_solid_nodes.data()[p_l2g.data()[eN_j]];
                }
              double element_nodes[nDOF_mesh_trial_element*3];
              for (int i=0;i<nDOF_mesh_trial_element;i++)
                {
                  register int eN_i=eN*nDOF_mesh_trial_element+i;
                  for(int I=0;I<3;I++)
                    element_nodes[i*3 + I] = mesh_dof[mesh_l2g[eN_i]*3 + I];
                }//i
              int icase_s = gf_s.calculate(element_phi_s, element_nodes, x_ref.data(), false);
              for (int k=0 ; k < nQuadraturePoints_element ; ++k)
                {
                  int eN_k = eN*nQuadraturePoints_element + k;
                  int eN_nDOF_trial_element = eN*nDOF_trial_element;
                  
                  double jac[nSpace*nSpace];
                  double jacInv[nSpace*nSpace];
                  double p=0.0,pe=0.0;
                  double jacDet, x, y, z, dV, h_phi;
                  gf_s.set_quad(k);
                  double H_s = gf_s.H(0.,0.);
                  ck.calculateMapping_element(eN,
                                              k,
                                              mesh_dof.data(),
                                              mesh_l2g.data(),
                                              mesh_trial_ref.data(),
                                              mesh_grad_trial_ref.data(),
                                              jac,
                                              jacDet,
                                              jacInv,
                                              x,y,z);
                  dV = fabs(jacDet)*dV_ref.data()[k];
                  ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_ref.data()[k*nDOF_trial_element],p);
		  if (elementIsActive[eN])
		    {
		      p_dv_new += p*H_s*dV;
		      pa_dv_new += q_u_0.data()[eN_k]*H_s*dV;
		      pe = p-q_u_0.data()[eN_k];
		      p_L1 += fabs(pe)*H_s*dV;
		      p_L2 += pe*pe*H_s*dV;
		      if (fabs(pe) > p_LI)
			p_LI = fabs(pe);
		    }
                }
            }
          p_L2 = sqrt(p_L2);
          //        std::cout<<"Pressure Integral Shifted"<<p_dv_new<<std::endl
          //         <<"Analytical Pressure Integral 2 "<<pa_dv_new<<std::endl
          //         <<"Errors "<<p_L1<<'\t'<<p_L2<<'\t'<<p_LI<<std::endl;
        }
      assert(errors.shape(0)*errors.shape(1) == 15);
      MPI_Allreduce(MPI_IN_PLACE, errors.data(),(errors.shape(0)-1)*errors.shape(1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, errors.data()+(errors.shape(0)-1)*errors.shape(1),1*errors.shape(1),MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      assert(p_L2 >= 0.0);
      assert(u_L2 >= 0.0);
      assert(v_L2 >= 0.0);
      assert(velocity_L2 >= 0.0);
      p_L2 = sqrt(p_L2);
      u_L2 = sqrt(u_L2);
      v_L2 = sqrt(v_L2);
      velocity_L2 = sqrt(velocity_L2);
    }

    void calculateJacobian(arguments_dict& args)
    {
      double NONCONSERVATIVE_FORM = args.m_dscalar["NONCONSERVATIVE_FORM"];
      double MOMENTUM_SGE = args.m_dscalar["MOMENTUM_SGE"];
      double PRESSURE_SGE = args.m_dscalar["PRESSURE_SGE"];
      double VELOCITY_SGE = args.m_dscalar["VELOCITY_SGE"];
      double PRESSURE_PROJECTION_STABILIZATION = args.m_dscalar["PRESSURE_PROJECTION_STABILIZATION"];
      xt::pyarray<double>& mesh_trial_ref = args.m_darray["mesh_trial_ref"];
      xt::pyarray<double>& mesh_grad_trial_ref = args.m_darray["mesh_grad_trial_ref"];
      xt::pyarray<double>& mesh_dof = args.m_darray["mesh_dof"];
      xt::pyarray<double>& mesh_velocity_dof = args.m_darray["mesh_velocity_dof"];
      double MOVING_DOMAIN = args.m_dscalar["MOVING_DOMAIN"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& x_ref = args.m_darray["x_ref"];
      xt::pyarray<double>& dV_ref = args.m_darray["dV_ref"];
      xt::pyarray<double>& p_trial_ref = args.m_darray["p_trial_ref"];
      xt::pyarray<double>& p_grad_trial_ref = args.m_darray["p_grad_trial_ref"];
      xt::pyarray<double>& p_test_ref = args.m_darray["p_test_ref"];
      xt::pyarray<double>& p_grad_test_ref = args.m_darray["p_grad_test_ref"];
      xt::pyarray<double>& vel_trial_ref = args.m_darray["vel_trial_ref"];
      xt::pyarray<double>& vel_grad_trial_ref = args.m_darray["vel_grad_trial_ref"];
      xt::pyarray<double>& vel_test_ref = args.m_darray["vel_test_ref"];
      xt::pyarray<double>& vel_grad_test_ref = args.m_darray["vel_grad_test_ref"];
      xt::pyarray<double>& mesh_trial_trace_ref = args.m_darray["mesh_trial_trace_ref"];
      xt::pyarray<double>& mesh_grad_trial_trace_ref = args.m_darray["mesh_grad_trial_trace_ref"];
      xt::pyarray<double>& xb_ref = args.m_darray["xb_ref"];
      xt::pyarray<double>& dS_ref = args.m_darray["dS_ref"];
      xt::pyarray<double>& p_trial_trace_ref = args.m_darray["p_trial_trace_ref"];
      xt::pyarray<double>& p_grad_trial_trace_ref = args.m_darray["p_grad_trial_trace_ref"];
      xt::pyarray<double>& p_test_trace_ref = args.m_darray["p_test_trace_ref"];
      xt::pyarray<double>& p_grad_test_trace_ref = args.m_darray["p_grad_test_trace_ref"];
      xt::pyarray<double>& vel_trial_trace_ref = args.m_darray["vel_trial_trace_ref"];
      xt::pyarray<double>& vel_grad_trial_trace_ref = args.m_darray["vel_grad_trial_trace_ref"];
      xt::pyarray<double>& vel_test_trace_ref = args.m_darray["vel_test_trace_ref"];
      xt::pyarray<double>& vel_grad_test_trace_ref = args.m_darray["vel_grad_test_trace_ref"];
      xt::pyarray<double>& normal_ref = args.m_darray["normal_ref"];
      xt::pyarray<double>& boundaryJac_ref = args.m_darray["boundaryJac_ref"];
      double eb_adjoint_sigma = args.m_dscalar["eb_adjoint_sigma"];
      xt::pyarray<double>& elementDiameter = args.m_darray["elementDiameter"];
      xt::pyarray<double>& elementBoundaryDiameter = args.m_darray["elementBoundaryDiameter"];
      xt::pyarray<double>& nodeDiametersArray = args.m_darray["nodeDiametersArray"];
      double hFactor = args.m_dscalar["hFactor"];
      int nElements_global = args.m_iscalar["nElements_global"];
      double useRBLES = args.m_dscalar["useRBLES"];
      double useMetrics = args.m_dscalar["useMetrics"];
      double alphaBDF = args.m_dscalar["alphaBDF"];
      double epsFact_rho = args.m_dscalar["epsFact_rho"];
      double epsFact_mu = args.m_dscalar["epsFact_mu"];
      double sigma = args.m_dscalar["sigma"];
      double rho_0 = args.m_dscalar["rho_0"];
      double nu_0 = args.m_dscalar["nu_0"];
      double rho_1 = args.m_dscalar["rho_1"];
      double nu_1 = args.m_dscalar["nu_1"];
      double smagorinskyConstant = args.m_dscalar["smagorinskyConstant"];
      int turbulenceClosureModel = args.m_iscalar["turbulenceClosureModel"];
      double Ct_sge = args.m_dscalar["Ct_sge"];
      double Cd_sge = args.m_dscalar["Cd_sge"];
      double C_dg = args.m_dscalar["C_dg"];
      double C_b = args.m_dscalar["C_b"];
      const xt::pyarray<double>& eps_solid = args.m_darray["eps_solid"];
      const xt::pyarray<double>& phi_solid = args.m_darray["phi_solid"];
      const xt::pyarray<double>& eps_porous = args.m_darray["eps_porous"];
      const xt::pyarray<double>& phi_porous = args.m_darray["phi_porous"];
      const xt::pyarray<double>& q_velocity_porous = args.m_darray["q_velocity_porous"];
      const xt::pyarray<double>& q_porosity = args.m_darray["q_porosity"];
      const xt::pyarray<double>& q_dragAlpha = args.m_darray["q_dragAlpha"];
      const xt::pyarray<double>& q_dragBeta = args.m_darray["q_dragBeta"];
      const xt::pyarray<double>& q_mass_source = args.m_darray["q_mass_source"];
      const xt::pyarray<double>& q_turb_var_0 = args.m_darray["q_turb_var_0"];
      const xt::pyarray<double>& q_turb_var_1 = args.m_darray["q_turb_var_1"];
      const xt::pyarray<double>& q_turb_var_grad_0 = args.m_darray["q_turb_var_grad_0"];
      const double LAG_LES = args.m_dscalar["LAG_LES"];
      xt::pyarray<double> & q_eddy_viscosity_last = args.m_darray["q_eddy_viscosity_last"];
      xt::pyarray<double> & ebqe_eddy_viscosity_last = args.m_darray["ebqe_eddy_viscosity_last"];
      xt::pyarray<int>& p_l2g = args.m_iarray["p_l2g"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<double>& p_dof = args.m_darray["p_dof"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      xt::pyarray<double>& p_old_dof = args.m_darray["p_old_dof"];
      xt::pyarray<double>& u_old_dof = args.m_darray["u_old_dof"];
      xt::pyarray<double>& v_old_dof = args.m_darray["v_old_dof"];
      xt::pyarray<double>& w_old_dof = args.m_darray["w_old_dof"];
      xt::pyarray<double>& g = args.m_darray["g"];
      const double useVF = args.m_dscalar["useVF"];
      xt::pyarray<double>& vf = args.m_darray["vf"];
      xt::pyarray<double>& phi = args.m_darray["phi"];
      xt::pyarray<double>& phi_nodes = args.m_darray["phi_nodes"];
      xt::pyarray<double>& normal_phi = args.m_darray["normal_phi"];
      xt::pyarray<double>& kappa_phi = args.m_darray["kappa_phi"];
      xt::pyarray<double>& q_mom_u_acc_beta_bdf = args.m_darray["q_mom_u_acc_beta_bdf"];
      xt::pyarray<double>& q_mom_v_acc_beta_bdf = args.m_darray["q_mom_v_acc_beta_bdf"];
      xt::pyarray<double>& q_mom_w_acc_beta_bdf = args.m_darray["q_mom_w_acc_beta_bdf"];
      xt::pyarray<double>& q_dV = args.m_darray["q_dV"];
      xt::pyarray<double>& q_dV_last = args.m_darray["q_dV_last"];
      xt::pyarray<double>& q_velocity_sge = args.m_darray["q_velocity_sge"];
      xt::pyarray<double>& q_cfl = args.m_darray["q_cfl"];
      xt::pyarray<double>& q_numDiff_u_last = args.m_darray["q_numDiff_u_last"];
      xt::pyarray<double>& q_numDiff_v_last = args.m_darray["q_numDiff_v_last"];
      xt::pyarray<double>& q_numDiff_w_last = args.m_darray["q_numDiff_w_last"];
      xt::pyarray<int>& sdInfo_u_u_rowptr = args.m_iarray["sdInfo_u_u_rowptr"];
      xt::pyarray<int>& sdInfo_u_u_colind = args.m_iarray["sdInfo_u_u_colind"];
      xt::pyarray<int>& sdInfo_u_v_rowptr = args.m_iarray["sdInfo_u_v_rowptr"];
      xt::pyarray<int>& sdInfo_u_v_colind = args.m_iarray["sdInfo_u_v_colind"];
      xt::pyarray<int>& sdInfo_u_w_rowptr = args.m_iarray["sdInfo_u_w_rowptr"];
      xt::pyarray<int>& sdInfo_u_w_colind = args.m_iarray["sdInfo_u_w_colind"];
      xt::pyarray<int>& sdInfo_v_v_rowptr = args.m_iarray["sdInfo_v_v_rowptr"];
      xt::pyarray<int>& sdInfo_v_v_colind = args.m_iarray["sdInfo_v_v_colind"];
      xt::pyarray<int>& sdInfo_v_u_rowptr = args.m_iarray["sdInfo_v_u_rowptr"];
      xt::pyarray<int>& sdInfo_v_u_colind = args.m_iarray["sdInfo_v_u_colind"];
      xt::pyarray<int>& sdInfo_v_w_rowptr = args.m_iarray["sdInfo_v_w_rowptr"];
      xt::pyarray<int>& sdInfo_v_w_colind = args.m_iarray["sdInfo_v_w_colind"];
      xt::pyarray<int>& sdInfo_w_w_rowptr = args.m_iarray["sdInfo_w_w_rowptr"];
      xt::pyarray<int>& sdInfo_w_w_colind = args.m_iarray["sdInfo_w_w_colind"];
      xt::pyarray<int>& sdInfo_w_u_rowptr = args.m_iarray["sdInfo_w_u_rowptr"];
      xt::pyarray<int>& sdInfo_w_u_colind = args.m_iarray["sdInfo_w_u_colind"];
      xt::pyarray<int>& sdInfo_w_v_rowptr = args.m_iarray["sdInfo_w_v_rowptr"];
      xt::pyarray<int>& sdInfo_w_v_colind = args.m_iarray["sdInfo_w_v_colind"];
      xt::pyarray<int>& csrRowIndeces_p_p = args.m_iarray["csrRowIndeces_p_p"];
      xt::pyarray<int>& csrColumnOffsets_p_p = args.m_iarray["csrColumnOffsets_p_p"];
      xt::pyarray<int>& csrRowIndeces_p_u = args.m_iarray["csrRowIndeces_p_u"];
      xt::pyarray<int>& csrColumnOffsets_p_u = args.m_iarray["csrColumnOffsets_p_u"];
      xt::pyarray<int>& csrRowIndeces_p_v = args.m_iarray["csrRowIndeces_p_v"];
      xt::pyarray<int>& csrColumnOffsets_p_v = args.m_iarray["csrColumnOffsets_p_v"];
      xt::pyarray<int>& csrRowIndeces_p_w = args.m_iarray["csrRowIndeces_p_w"];
      xt::pyarray<int>& csrColumnOffsets_p_w = args.m_iarray["csrColumnOffsets_p_w"];
      xt::pyarray<int>& csrRowIndeces_u_p = args.m_iarray["csrRowIndeces_u_p"];
      xt::pyarray<int>& csrColumnOffsets_u_p = args.m_iarray["csrColumnOffsets_u_p"];
      xt::pyarray<int>& csrRowIndeces_u_u = args.m_iarray["csrRowIndeces_u_u"];
      xt::pyarray<int>& csrColumnOffsets_u_u = args.m_iarray["csrColumnOffsets_u_u"];
      xt::pyarray<int>& csrRowIndeces_u_v = args.m_iarray["csrRowIndeces_u_v"];
      xt::pyarray<int>& csrColumnOffsets_u_v = args.m_iarray["csrColumnOffsets_u_v"];
      xt::pyarray<int>& csrRowIndeces_u_w = args.m_iarray["csrRowIndeces_u_w"];
      xt::pyarray<int>& csrColumnOffsets_u_w = args.m_iarray["csrColumnOffsets_u_w"];
      xt::pyarray<int>& csrRowIndeces_v_p = args.m_iarray["csrRowIndeces_v_p"];
      xt::pyarray<int>& csrColumnOffsets_v_p = args.m_iarray["csrColumnOffsets_v_p"];
      xt::pyarray<int>& csrRowIndeces_v_u = args.m_iarray["csrRowIndeces_v_u"];
      xt::pyarray<int>& csrColumnOffsets_v_u = args.m_iarray["csrColumnOffsets_v_u"];
      xt::pyarray<int>& csrRowIndeces_v_v = args.m_iarray["csrRowIndeces_v_v"];
      xt::pyarray<int>& csrColumnOffsets_v_v = args.m_iarray["csrColumnOffsets_v_v"];
      xt::pyarray<int>& csrRowIndeces_v_w = args.m_iarray["csrRowIndeces_v_w"];
      xt::pyarray<int>& csrColumnOffsets_v_w = args.m_iarray["csrColumnOffsets_v_w"];
      xt::pyarray<int>& csrRowIndeces_w_p = args.m_iarray["csrRowIndeces_w_p"];
      xt::pyarray<int>& csrColumnOffsets_w_p = args.m_iarray["csrColumnOffsets_w_p"];
      xt::pyarray<int>& csrRowIndeces_w_u = args.m_iarray["csrRowIndeces_w_u"];
      xt::pyarray<int>& csrColumnOffsets_w_u = args.m_iarray["csrColumnOffsets_w_u"];
      xt::pyarray<int>& csrRowIndeces_w_v = args.m_iarray["csrRowIndeces_w_v"];
      xt::pyarray<int>& csrColumnOffsets_w_v = args.m_iarray["csrColumnOffsets_w_v"];
      xt::pyarray<int>& csrRowIndeces_w_w = args.m_iarray["csrRowIndeces_w_w"];
      xt::pyarray<int>& csrColumnOffsets_w_w = args.m_iarray["csrColumnOffsets_w_w"];
      xt::pyarray<double>& globalJacobian = args.m_darray["globalJacobian"];
      int nExteriorElementBoundaries_global = args.m_iscalar["nExteriorElementBoundaries_global"];
      xt::pyarray<int>& exteriorElementBoundariesArray = args.m_iarray["exteriorElementBoundariesArray"];
      xt::pyarray<int>& elementBoundaryElementsArray = args.m_iarray["elementBoundaryElementsArray"];
      xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.m_iarray["elementBoundaryLocalElementBoundariesArray"];
      xt::pyarray<double>& ebqe_vf_ext = args.m_darray["ebqe_vf_ext"];
      xt::pyarray<double>& bc_ebqe_vf_ext = args.m_darray["bc_ebqe_vf_ext"];
      xt::pyarray<double>& ebqe_phi_ext = args.m_darray["ebqe_phi_ext"];
      xt::pyarray<double>& bc_ebqe_phi_ext = args.m_darray["bc_ebqe_phi_ext"];
      xt::pyarray<double>& ebqe_normal_phi_ext = args.m_darray["ebqe_normal_phi_ext"];
      xt::pyarray<double>& ebqe_kappa_phi_ext = args.m_darray["ebqe_kappa_phi_ext"];
      const xt::pyarray<double>& ebqe_porosity_ext = args.m_darray["ebqe_porosity_ext"];
      const xt::pyarray<double>& ebqe_turb_var_0 = args.m_darray["ebqe_turb_var_0"];
      const xt::pyarray<double>& ebqe_turb_var_1 = args.m_darray["ebqe_turb_var_1"];
      xt::pyarray<int>& isDOFBoundary_p = args.m_iarray["isDOFBoundary_p"];
      xt::pyarray<int>& isDOFBoundary_u = args.m_iarray["isDOFBoundary_u"];
      xt::pyarray<int>& isDOFBoundary_v = args.m_iarray["isDOFBoundary_v"];
      xt::pyarray<int>& isDOFBoundary_w = args.m_iarray["isDOFBoundary_w"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_p = args.m_iarray["isAdvectiveFluxBoundary_p"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_u = args.m_iarray["isAdvectiveFluxBoundary_u"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_v = args.m_iarray["isAdvectiveFluxBoundary_v"];
      xt::pyarray<int>& isAdvectiveFluxBoundary_w = args.m_iarray["isAdvectiveFluxBoundary_w"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_u = args.m_iarray["isDiffusiveFluxBoundary_u"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_v = args.m_iarray["isDiffusiveFluxBoundary_v"];
      xt::pyarray<int>& isDiffusiveFluxBoundary_w = args.m_iarray["isDiffusiveFluxBoundary_w"];
      xt::pyarray<double>& ebqe_bc_p_ext = args.m_darray["ebqe_bc_p_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mass_ext = args.m_darray["ebqe_bc_flux_mass_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_u_adv_ext = args.m_darray["ebqe_bc_flux_mom_u_adv_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_v_adv_ext = args.m_darray["ebqe_bc_flux_mom_v_adv_ext"];
      xt::pyarray<double>& ebqe_bc_flux_mom_w_adv_ext = args.m_darray["ebqe_bc_flux_mom_w_adv_ext"];
      xt::pyarray<double>& ebqe_bc_u_ext = args.m_darray["ebqe_bc_u_ext"];
      xt::pyarray<double>& ebqe_bc_flux_u_diff_ext = args.m_darray["ebqe_bc_flux_u_diff_ext"];
      xt::pyarray<double>& ebqe_penalty_ext = args.m_darray["ebqe_penalty_ext"];
      xt::pyarray<double>& ebqe_bc_v_ext = args.m_darray["ebqe_bc_v_ext"];
      xt::pyarray<double>& ebqe_bc_flux_v_diff_ext = args.m_darray["ebqe_bc_flux_v_diff_ext"];
      xt::pyarray<double>& ebqe_bc_w_ext = args.m_darray["ebqe_bc_w_ext"];
      xt::pyarray<double>& ebqe_bc_flux_w_diff_ext = args.m_darray["ebqe_bc_flux_w_diff_ext"];
      xt::pyarray<int>& csrColumnOffsets_eb_p_p = args.m_iarray["csrColumnOffsets_eb_p_p"];
      xt::pyarray<int>& csrColumnOffsets_eb_p_u = args.m_iarray["csrColumnOffsets_eb_p_u"];
      xt::pyarray<int>& csrColumnOffsets_eb_p_v = args.m_iarray["csrColumnOffsets_eb_p_v"];
      xt::pyarray<int>& csrColumnOffsets_eb_p_w = args.m_iarray["csrColumnOffsets_eb_p_w"];
      xt::pyarray<int>& csrColumnOffsets_eb_u_p = args.m_iarray["csrColumnOffsets_eb_u_p"];
      xt::pyarray<int>& csrColumnOffsets_eb_u_u = args.m_iarray["csrColumnOffsets_eb_u_u"];
      xt::pyarray<int>& csrColumnOffsets_eb_u_v = args.m_iarray["csrColumnOffsets_eb_u_v"];
      xt::pyarray<int>& csrColumnOffsets_eb_u_w = args.m_iarray["csrColumnOffsets_eb_u_w"];
      xt::pyarray<int>& csrColumnOffsets_eb_v_p = args.m_iarray["csrColumnOffsets_eb_v_p"];
      xt::pyarray<int>& csrColumnOffsets_eb_v_u = args.m_iarray["csrColumnOffsets_eb_v_u"];
      xt::pyarray<int>& csrColumnOffsets_eb_v_v = args.m_iarray["csrColumnOffsets_eb_v_v"];
      xt::pyarray<int>& csrColumnOffsets_eb_v_w = args.m_iarray["csrColumnOffsets_eb_v_w"];
      xt::pyarray<int>& csrColumnOffsets_eb_w_p = args.m_iarray["csrColumnOffsets_eb_w_p"];
      xt::pyarray<int>& csrColumnOffsets_eb_w_u = args.m_iarray["csrColumnOffsets_eb_w_u"];
      xt::pyarray<int>& csrColumnOffsets_eb_w_v = args.m_iarray["csrColumnOffsets_eb_w_v"];
      xt::pyarray<int>& csrColumnOffsets_eb_w_w = args.m_iarray["csrColumnOffsets_eb_w_w"];
      xt::pyarray<int>& elementFlags = args.m_iarray["elementFlags"];
      xt::pyarray<int>& boundaryFlags = args.m_iarray["boundaryFlags"];
      int use_ball_as_particle = args.m_iscalar["use_ball_as_particle"];
      xt::pyarray<double>& ball_center = args.m_darray["ball_center"];
      xt::pyarray<double>& ball_radius = args.m_darray["ball_radius"];
      xt::pyarray<double>& ball_velocity = args.m_darray["ball_velocity"];
      xt::pyarray<double>& ball_angular_velocity = args.m_darray["ball_angular_velocity"];
      xt::pyarray<double>& ball_density = args.m_darray["ball_density"];
      xt::pyarray<double>& particle_signed_distances = args.m_darray["particle_signed_distances"];
      xt::pyarray<double>& particle_signed_distance_normals = args.m_darray["particle_signed_distance_normals"];
      xt::pyarray<double>& particle_velocities = args.m_darray["particle_velocities"];
      xt::pyarray<double>& particle_centroids = args.m_darray["particle_centroids"];
      xt::pyarray<double>& ebqe_phi_s = args.m_darray["ebqe_phi_s"];
      xt::pyarray<double>& ebq_global_grad_phi_s = args.m_darray["ebq_global_grad_phi_s"];
      xt::pyarray<double>& ebq_particle_velocity_s = args.m_darray["ebq_particle_velocity_s"];
      xt::pyarray<double>& phi_solid_nodes = args.m_darray["phi_solid_nodes"];
      xt::pyarray<double>& distance_to_solids = args.m_darray["distance_to_solids"];
      int nParticles = args.m_iscalar["nParticles"];
      int nElements_owned = args.m_iscalar["nElements_owned"];
      double particle_nitsche = args.m_dscalar["particle_nitsche"];
      double particle_epsFact = args.m_dscalar["particle_epsFact"];
      double particle_alpha = args.m_dscalar["particle_alpha"];
      double particle_beta = args.m_dscalar["particle_beta"];
      double particle_penalty_constant = args.m_dscalar["particle_penalty_constant"];
      double ghost_penalty_constant = args.m_dscalar["ghost_penalty_constant"];
      const bool useExact = args.m_iscalar["useExact"];
      const int nQuadraturePoints_global(nElements_global*nQuadraturePoints_element);
      std::valarray<double> particle_surfaceArea_tmp(nParticles), particle_netForces_tmp(nParticles*3*3), particle_netMoments_tmp(nParticles*3);
      gf.useExact = false;//useExact;
      gf_p.useExact = false;//useExact;
      gf_s.useExact = useExact;
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
        {
          register double eps_rho,eps_mu;

          register double  elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
            elementJacobian_p_u[nDOF_test_element][nDOF_v_trial_element],
            elementJacobian_p_v[nDOF_test_element][nDOF_v_trial_element],
            elementJacobian_p_w[nDOF_test_element][nDOF_v_trial_element],
            elementJacobian_u_p[nDOF_v_test_element][nDOF_trial_element],
            elementJacobian_u_u[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_u_v[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_u_w[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_v_p[nDOF_v_test_element][nDOF_trial_element],
            elementJacobian_v_u[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_v_v[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_v_w[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_w_p[nDOF_v_test_element][nDOF_trial_element],
            elementJacobian_w_u[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_w_v[nDOF_v_test_element][nDOF_v_trial_element],
            elementJacobian_w_w[nDOF_v_test_element][nDOF_v_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_trial_element;j++)
              {
                elementJacobian_p_p[i][j]=0.0;
              }
          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_v_trial_element;j++)
              {
                elementJacobian_p_u[i][j]=0.0;
                elementJacobian_p_v[i][j]=0.0;
                elementJacobian_p_w[i][j]=0.0;
                elementJacobian_u_p[j][i]=0.0;
                elementJacobian_v_p[j][i]=0.0;
                elementJacobian_w_p[j][i]=0.0;
              }
          for (int i=0;i<nDOF_v_test_element;i++)
            for (int j=0;j<nDOF_v_trial_element;j++)
              {
                elementJacobian_u_u[i][j]=0.0;
                elementJacobian_u_v[i][j]=0.0;
                elementJacobian_u_w[i][j]=0.0;
                elementJacobian_v_u[i][j]=0.0;
                elementJacobian_v_v[i][j]=0.0;
                elementJacobian_v_w[i][j]=0.0;
                elementJacobian_w_u[i][j]=0.0;
                elementJacobian_w_v[i][j]=0.0;
                elementJacobian_w_w[i][j]=0.0;
              }
          double element_phi[nDOF_mesh_trial_element], element_phi_s[nDOF_mesh_trial_element];
          for (int j=0;j<nDOF_mesh_trial_element;j++)
            {
              register int eN_j = eN*nDOF_mesh_trial_element+j;
              element_phi[j] = phi_nodes.data()[p_l2g.data()[eN_j]];
              element_phi_s[j] = phi_solid_nodes.data()[p_l2g.data()[eN_j]];
            }
          double element_nodes[nDOF_mesh_trial_element*3];
          for (int i=0;i<nDOF_mesh_trial_element;i++)
            {
              register int eN_i=eN*nDOF_mesh_trial_element+i;
              for(int I=0;I<3;I++)
                element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
            }//i
          int icase_s = gf_s.calculate(element_phi_s, element_nodes, x_ref.data(), false);
#ifdef IFEM
          int icase_p = gf_p.calculate(element_phi, element_nodes, x_ref.data(), -rho_1*g.data()[1], -rho_0*g.data()[1],false,true);
          int icase = gf.calculate(element_phi, element_nodes, x_ref.data(), rho_1*nu_1, rho_0*nu_0,false,false);
#else
          int icase_p = gf_p.calculate(element_phi, element_nodes, x_ref.data(), 1.,1.,false,false);
          int icase = gf.calculate(element_phi, element_nodes, x_ref.data(), 1.,1.,false,false);
#endif
          for (int fluid_phase=0;fluid_phase < 2 - abs(icase); fluid_phase++)
            {
              for  (int k=0;k<nQuadraturePoints_element;k++)
                {
                  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                    eN_k_nSpace = eN_k*nSpace,
                    eN_k_3d = eN_k*3,
                    eN_nDOF_trial_element = eN*nDOF_trial_element, //index to a vector at a quadrature point
                    eN_nDOF_v_trial_element = eN*nDOF_v_trial_element; //index to a vector at a quadrature point

                  //declare local storage
                  register double p=0.0,u=0.0,v=0.0,w=0.0,
                    grad_p[nSpace]=ZEROVEC,grad_u[nSpace]=ZEROVEC,grad_v[nSpace]=ZEROVEC,grad_w[nSpace]=ZEROVEC,
                    p_old=0.0,u_old=0.0,v_old=0.0,w_old=0.0,
                    grad_p_old[nSpace]=ZEROVEC,grad_u_old[nSpace]=ZEROVEC,grad_v_old[nSpace]=ZEROVEC,grad_w_old[nSpace]=ZEROVEC,
                    mom_u_acc=0.0,
                    dmom_u_acc_u=0.0,
                    mom_v_acc=0.0,
                    dmom_v_acc_v=0.0,
                    mom_w_acc=0.0,
                    dmom_w_acc_w=0.0,
                    mass_adv[nSpace]=ZEROVEC,
                    dmass_adv_u[nSpace]=ZEROVEC,
                    dmass_adv_v[nSpace]=ZEROVEC,
                    dmass_adv_w[nSpace]=ZEROVEC,
                    mass_ham=0.0,
                    dmass_ham_u=0.0,
                    dmass_ham_v=0.0,
                    dmass_ham_w=0.0,
                    mom_u_adv[nSpace]=ZEROVEC,
                    dmom_u_adv_u[nSpace]=ZEROVEC,
                    dmom_u_adv_v[nSpace]=ZEROVEC,
                    dmom_u_adv_w[nSpace]=ZEROVEC,
                    mom_v_adv[nSpace]=ZEROVEC,
                    dmom_v_adv_u[nSpace]=ZEROVEC,
                    dmom_v_adv_v[nSpace]=ZEROVEC,
                    dmom_v_adv_w[nSpace]=ZEROVEC,
                    mom_w_adv[nSpace]=ZEROVEC,
                    dmom_w_adv_u[nSpace]=ZEROVEC,
                    dmom_w_adv_v[nSpace]=ZEROVEC,
                    dmom_w_adv_w[nSpace]=ZEROVEC,
                    mom_uu_diff_ten[nSpace]=ZEROVEC,
                    mom_vv_diff_ten[nSpace]=ZEROVEC,
                    mom_ww_diff_ten[nSpace]=ZEROVEC,
                    mom_uv_diff_ten[1],
                    mom_uw_diff_ten[1],
                    mom_vu_diff_ten[1],
                    mom_vw_diff_ten[1],
                    mom_wu_diff_ten[1],
                    mom_wv_diff_ten[1],
                    mom_u_source=0.0,
                    mom_v_source=0.0,
                    mom_w_source=0.0,
                    mom_u_ham=0.0,
                    dmom_u_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_u[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_v[nSpace]=ZEROVEC,
                    dmom_u_ham_u=0.0,
                    dmom_u_ham_v=0.0,
                    dmom_u_ham_w=0.0,
                    mom_v_ham=0.0,
                    dmom_v_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_u[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_v[nSpace]=ZEROVEC,
                    dmom_v_ham_u=0.0,
                    dmom_v_ham_v=0.0,
                    dmom_v_ham_w=0.0,
                    mom_w_ham=0.0,
                    dmom_w_ham_grad_p[nSpace]=ZEROVEC,
                    dmom_w_ham_grad_w[nSpace]=ZEROVEC,
                    dmom_w_ham_u=0.0,
                    dmom_w_ham_v=0.0,
                    dmom_w_ham_w=0.0,
                    mom_u_acc_t=0.0,
                    dmom_u_acc_u_t=0.0,
                    mom_v_acc_t=0.0,
                    dmom_v_acc_v_t=0.0,
                    mom_w_acc_t=0.0,
                    dmom_w_acc_w_t=0.0,
                    pdeResidual_p=0.0,
                    pdeResidual_u=0.0,
                    pdeResidual_v=0.0,
                    pdeResidual_w=0.0,
                    dpdeResidual_p_u[nDOF_v_trial_element],dpdeResidual_p_v[nDOF_v_trial_element],dpdeResidual_p_w[nDOF_v_trial_element],
                    dpdeResidual_u_p[nDOF_trial_element],dpdeResidual_u_u[nDOF_v_trial_element],
                    dpdeResidual_v_p[nDOF_trial_element],dpdeResidual_v_v[nDOF_v_trial_element],
                    dpdeResidual_w_p[nDOF_trial_element],dpdeResidual_w_w[nDOF_v_trial_element],
                    Lstar_u_p[nDOF_test_element],
                    Lstar_v_p[nDOF_test_element],
                    Lstar_w_p[nDOF_test_element],
                    Lstar_u_u[nDOF_v_test_element],
                    Lstar_v_v[nDOF_v_test_element],
                    Lstar_w_w[nDOF_v_test_element],
                    Lstar_p_u[nDOF_v_test_element],
                    Lstar_p_v[nDOF_v_test_element],
                    Lstar_p_w[nDOF_v_test_element],
                    subgridError_p=0.0,
                    subgridError_u=0.0,
                    subgridError_v=0.0,
                    subgridError_w=0.0,
                    dsubgridError_p_u[nDOF_v_trial_element],
                    dsubgridError_p_v[nDOF_v_trial_element],
                    dsubgridError_p_w[nDOF_v_trial_element],
                    dsubgridError_u_p[nDOF_trial_element],
                    dsubgridError_u_u[nDOF_v_trial_element],
                    dsubgridError_v_p[nDOF_trial_element],
                    dsubgridError_v_v[nDOF_v_trial_element],
                    dsubgridError_w_p[nDOF_trial_element],
                    dsubgridError_w_w[nDOF_v_trial_element],
                    tau_p=0.0,tau_p0=0.0,tau_p1=0.0,
                    tau_v=0.0,tau_v0=0.0,tau_v1=0.0,
                    jac[nSpace*nSpace],
                    jacDet,
                    jacInv[nSpace*nSpace],
                    p_trial[nDOF_trial_element], vel_trial[nDOF_v_trial_element],
                    p_grad_trial_ib[nDOF_trial_element*nSpace], vel_grad_trial_ib[nDOF_v_trial_element*nSpace],
                    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_v_trial_element*nSpace],
                    dV,
                    p_test_dV[nDOF_test_element],vel_test_dV[nDOF_v_test_element],
                    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_v_test_element*nSpace],
                    x,y,z,xt,yt,zt,
                    //VRANS
                    porosity,
                    //meanGrainSize,
                    dmom_u_source[nSpace]=ZEROVEC,
                    dmom_v_source[nSpace]=ZEROVEC,
                    dmom_w_source[nSpace]=ZEROVEC,
                    mass_source,
                    //
                    G[nSpace*nSpace],G_dd_G,tr_G,h_phi, dmom_adv_star[nSpace]=ZEROVEC, dmom_adv_sge[nSpace]=ZEROVEC, dmom_ham_grad_sge[nSpace]=ZEROVEC,
                    //embedded solid terms
                    mass_source_s=0.0,
                    mom_u_source_s=0.0,
                    mom_v_source_s=0.0,
                    mom_w_source_s=0.0,
                    dmom_u_source_s[nSpace]=ZEROVEC,
                    dmom_v_source_s[nSpace]=ZEROVEC,
                    dmom_w_source_s[nSpace]=ZEROVEC,
                    mom_u_adv_s[nSpace]=ZEROVEC,
                    mom_v_adv_s[nSpace]=ZEROVEC,
                    mom_w_adv_s[nSpace]=ZEROVEC,
                    dmom_u_adv_u_s[nSpace]=ZEROVEC,
                    dmom_v_adv_v_s[nSpace]=ZEROVEC,
                    dmom_w_adv_w_s[nSpace]=ZEROVEC,
                    mom_u_ham_s=0.0,
                    dmom_u_ham_grad_u_s[nSpace]=ZEROVEC,
                    dmom_u_ham_grad_v_s[nSpace]=ZEROVEC,
                    dmom_u_ham_u_s=0.0,
                    dmom_u_ham_v_s=0.0,
                    dmom_u_ham_w_s=0.0,
                    mom_v_ham_s=0.0,
                    dmom_v_ham_grad_u_s[nSpace]=ZEROVEC,
                    dmom_v_ham_grad_v_s[nSpace]=ZEROVEC,
                    dmom_v_ham_u_s=0.0,
                    dmom_v_ham_v_s=0.0,
                    dmom_v_ham_w_s=0.0,
                    mom_w_ham_s=0.0,
                    dmom_w_ham_grad_w_s[nSpace]=ZEROVEC,
                    dmom_w_ham_u_s=0.0,
                    dmom_w_ham_v_s=0.0,
                    dmom_w_ham_w_s=0.0,
                    mass_ham_s=0.0,
                    dmass_ham_u_s=0.0,
                    dmass_ham_v_s=0.0,
                    dmass_ham_w_s=0.0;
                  //get jacobian, etc for mapping reference element
                  gf_s.set_quad(k);
                  gf.set_quad(k);
                  gf_p.set_quad(k);
                  ck.calculateMapping_element(eN,
                                              k,
                                              mesh_dof.data(),
                                              mesh_l2g.data(),
                                              mesh_trial_ref.data(),
                                              mesh_grad_trial_ref.data(),
                                              jac,
                                              jacDet,
                                              jacInv,
                                              x,y,z);
                  ck.calculateH_element(eN,
                                        k,
                                        nodeDiametersArray.data(),
                                        mesh_l2g.data(),
                                        mesh_trial_ref.data(),
                                        h_phi);
                  ck.calculateMappingVelocity_element(eN,
                                                      k,
                                                      mesh_velocity_dof.data(),
                                                      mesh_l2g.data(),
                                                      mesh_trial_ref.data(),
                                                      xt,yt,zt);
                  //xt=0.0;yt=0.0;zt=0.0;
                  //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
                  //get the physical integration weight
                  dV = fabs(jacDet)*dV_ref.data()[k];
                  ck.calculateG(jacInv,G,G_dd_G,tr_G);
                  //ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);

                  eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                  eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                  //get the trial function gradients
                  ck.gradTrialFromRef(&p_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
                  ck_v.gradTrialFromRef(&vel_grad_trial_ref.data()[k*nDOF_v_trial_element*nSpace],jacInv,vel_grad_trial);
                  for (int i=0; i < nDOF_trial_element; i++)
                    {
                      p_trial[i] = p_trial_ref.data()[k*nDOF_trial_element + i];
                      p_grad_trial_ib[i*nSpace + 0] = p_grad_trial[i*nSpace+0];
                      p_grad_trial_ib[i*nSpace + 1] = p_grad_trial[i*nSpace+1];
                    }
                  for (int i=0; i < nDOF_v_trial_element; i++)
                    {
                      vel_trial[i] = vel_trial_ref.data()[k*nDOF_v_trial_element + i];
                      vel_grad_trial_ib[i*nSpace + 0] = vel_grad_trial[i*nSpace+0];
                      vel_grad_trial_ib[i*nSpace + 1] = vel_grad_trial[i*nSpace+1];
                    }
                  if (icase == 0)
                    {
#ifdef IFEMBASIS
                      for (int i=0; i < nDOF_trial_element; i++)
                        {
                          if (fluid_phase == 0)
                            {
                              if (not std::isnan(gf_p.VA(i)))
                                {
                                  p_trial[i] = gf_p.VA(i);
                                  p_grad_trial_ib[i*nSpace + 0] = gf_p.VA_x(i);
                                  p_grad_trial_ib[i*nSpace + 1] = gf_p.VA_y(i);
                                }
                            }
                          else
                            {
                              if (not std::isnan(gf_p.VB(i)))
                                {
                                  p_trial[i] = gf_p.VB(i);
                                  p_grad_trial_ib[i*nSpace + 0] = gf_p.VB_x(i);
                                  p_grad_trial_ib[i*nSpace + 1] = gf_p.VB_y(i);
                                }
                            }
                        }
                      if(nDOF_v_trial_element == nDOF_trial_element)
                        {
                          for (int vi=0; vi < nDOF_v_trial_element; vi++)
                            {
                              if (fluid_phase == 0)
                                {
                                  if (not std::isnan(gf.VA(vi)))
                                    {
                                      vel_trial[vi] = gf.VA(vi);
                                      vel_grad_trial_ib[vi*nSpace + 0] = gf.VA_x(vi);
                                      vel_grad_trial_ib[vi*nSpace + 1] = gf.VA_y(vi);
                                    }
                                }
                              else
                                {
                                  if (not std::isnan(gf.VB(vi)))
                                    {
                                      vel_trial[vi] = gf.VB(vi);
                                      vel_grad_trial_ib[vi*nSpace + 0] = gf.VB_x(vi);
                                      vel_grad_trial_ib[vi*nSpace + 1] = gf.VB_y(vi);
                                    }
                                }
                            }
                        }
#endif
#ifndef IFEM
                      bool prob=false;
                      for (int vi=0; vi < nDOF_v_trial_element; vi++)
                        {
                          //pressure
                          if (fabs(p_trial_ref.data()[k*nDOF_trial_element + vi] - p_trial[vi]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Trial "<<p_trial_ref.data()[k*nDOF_trial_element + vj]<<'\t'<<gf_p.VA(vj)<<'\t'<<gf_p.VB(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(p_grad_trial[vi*nSpace + 0] - p_grad_trial_ib[vi*nSpace+0]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Grad Trial x"<<p_grad_trial[vj*nSpace + 0]<<'\t'<<gf_p.VA_x(vj)<<'\t'<<gf_p.VB_x(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(p_grad_trial[vi*nSpace + 1] - p_grad_trial_ib[vi*nSpace+1]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_trial_element; vj++)
                                std::cout<<"Grad Trial y "<<p_grad_trial[vj*nSpace + 1]<<'\t'<<gf_p.VA_y(vj)<<'\t'<<gf_p.VB_y(vj)<<std::endl;
                              prob=true;
                            }
                          //velocity
                          if (fabs(vel_trial_ref.data()[k*nDOF_v_trial_element + vi] - vel_trial[vi]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Trial "<<vel_trial_ref.data()[k*nDOF_v_trial_element + vj]<<'\t'<<gf.VA(vj)<<'\t'<<gf.VB(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(vel_grad_trial[vi*nSpace + 0] - vel_grad_trial_ib[vi*nSpace+0]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Grad Trial x"<<vel_grad_trial[vj*nSpace + 0]<<'\t'<<gf.VA_x(vj)<<'\t'<<gf.VB_x(vj)<<std::endl;
                              prob=true;
                            }
                          if (fabs(vel_grad_trial[vi*nSpace + 1] - vel_grad_trial_ib[vi*nSpace+1]) > 1.0e-8)
                            {
                              for (int vj=0; vj < nDOF_v_trial_element; vj++)
                                std::cout<<"Grad Trial y "<<vel_grad_trial[vj*nSpace + 1]<<'\t'<<gf.VA_y(vj)<<'\t'<<gf.VB_y(vj)<<std::endl;
                              prob=true;
                            }
                          if (prob)
                            break;
                        }
                      assert(!prob);
#endif
                    }
                  //get the solution
                  ck.valFromDOF(p_dof.data(),&p_l2g[eN_nDOF_trial_element],p_trial,p);
                  ck_v.valFromDOF(u_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_trial,u);
                  ck_v.valFromDOF(v_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_trial,v);
                  ck.valFromDOF(p_old_dof.data(),&p_l2g[eN_nDOF_trial_element],p_trial,p_old);
                  ck_v.valFromDOF(u_old_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_trial,u_old);
                  ck_v.valFromDOF(v_old_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_trial,v_old);
                  //get the solution gradients
                  ck.gradFromDOF(p_dof.data(),&p_l2g[eN_nDOF_trial_element],p_grad_trial_ib,grad_p);
                  ck_v.gradFromDOF(u_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_u);
                  ck_v.gradFromDOF(v_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_v);
                  ck.gradFromDOF(p_dof.data(),&p_l2g[eN_nDOF_trial_element],p_grad_trial_ib,grad_p_old);
                  ck_v.gradFromDOF(u_old_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_u_old);
                  ck_v.gradFromDOF(v_old_dof.data(),&vel_l2g[eN_nDOF_v_trial_element],vel_grad_trial_ib,grad_v_old);
                  //precalculate test function products with integration weights
#ifdef IFEMGALERKIN
                  for (int j=0;j<nDOF_test_element;j++)
                    {
                      p_test_dV[j] = p_trial[j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          p_grad_test_dV[j*nSpace+I]   = p_grad_trial_ib[j*nSpace+I]*dV;
                        }
                    }
                  for (int j=0;j<nDOF_v_test_element;j++)
                    {
                      vel_test_dV[j] = vel_trial[j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          vel_grad_test_dV[j*nSpace+I] = vel_grad_trial_ib[j*nSpace+I]*dV;
                        }
                    }
#else
                  for (int j=0;j<nDOF_test_element;j++)
                    {
                      p_test_dV[j] = p_test_ref.data()[k*nDOF_trial_element+j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//assume test_j == trial_j, ok for ifem
                        }
                    }
                  for (int j=0;j<nDOF_v_test_element;j++)
                    {
                      vel_test_dV[j] = vel_test_ref.data()[k*nDOF_v_trial_element+j]*dV;
                      for (int I=0;I<nSpace;I++)
                        {
                          vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//assume test_j == trial_j, ok for ifem
                        }
                    }
#endif
                  //needs to be fixed for higher-order meshes, assuming mesh trial is same as p trial
                  double div_mesh_velocity=0.0;
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int eN_j=eN*nDOF_trial_element+j;
                      div_mesh_velocity +=
                        mesh_velocity_dof.data()[mesh_l2g.data()[eN_j]*3+0]*p_grad_trial[j*nSpace+0] +
                        mesh_velocity_dof.data()[mesh_l2g.data()[eN_j]*3+1]*p_grad_trial[j*nSpace+1];
                    }
                  div_mesh_velocity = DM3*div_mesh_velocity + (1.0-DM3)*alphaBDF*(dV-q_dV_last.data()[eN_k])/dV;
                  //
                  //VRANS
                  porosity = q_porosity.data()[eN_k];
                  //
                  double ball_n[nSpace];
                  if (use_ball_as_particle == 1 && nParticles > 0)
                    {
                      int ball_index=get_distance_to_ball(nParticles, ball_center.data(), ball_radius.data(),x,y,z,distance_to_solids.data()[eN_k]);
                      get_normal_to_ith_ball(nParticles, ball_center.data(), ball_radius.data(),ball_index,x,y,z,ball_n[0],ball_n[1]);
                    }
                  else
                    {
                      //distance_to_solids is given in Prestep
                    }
                  //
                  //calculate pde coefficients and derivatives at quadrature points
                  //
                  double eddy_viscosity(0.);//not really interested in saving eddy_viscosity in jacobian
                  const double particle_eps  = particle_epsFact*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                  const double H_s = gf_s.H(particle_eps, phi_solid.data()[eN_k]);
                  double rho,nu;
                  if (gf.useExact)
                    {
                      if (icase == 0)
                        {
                          if (fluid_phase == 0)
                            {
                              rho=rho_0;
                              nu=nu_0;
                            }
                          else
                            {
                              rho=rho_1;
                              nu=nu_1;
                            }
                        }
                      else if (icase == -1)
                        {
                          rho=rho_0;
                          nu=nu_0;
                        }
                      else if (icase == 1)
                        {
                          rho=rho_1;
                          nu=nu_1;
                        }
                      else
                        assert(false);
                    }
                  else
                    {
                      double H = (1.0-useVF)*gf.H(eps_rho,phi[eN_k]) + useVF*fmin(1.0,fmax(0.0,vf[eN_k]));
                      double ImH = (1.0-useVF)*gf.ImH(eps_rho,phi[eN_k]) + useVF*(1.0-fmin(1.0,fmax(0.0,vf[eN_k])));
                      
                      rho  = rho_0*ImH + rho_1*H;
                      nu  = nu_0*ImH + nu_1*H;
                    }
                  evaluateCoefficients(NONCONSERVATIVE_FORM,
                                       sigma,
                                       rho,
                                       nu,
                                       elementDiameter.data()[eN],
                                       smagorinskyConstant,
                                       turbulenceClosureModel,
                                       g.data(),
                                       useVF,
                                       vf.data()[eN_k],
                                       phi.data()[eN_k],
                                       &normal_phi.data()[eN_k_nSpace],
                                       kappa_phi.data()[eN_k],
                                       //VRANS
                                       porosity,
                                       //
                                       phi_solid.data()[eN_k],//updated in get residual
                                       p_old,
                                       u_old,
                                       v_old,
                                       w_old,
                                       grad_p_old,
                                       grad_u_old,
                                       grad_v_old,
                                       grad_w_old,
                                       p,
                                       grad_p,
                                       grad_u,
                                       grad_v,
                                       grad_w,
                                       u,
                                       v,
                                       w,
                                       LAG_LES,
                                       eddy_viscosity,
                                       q_eddy_viscosity_last.data()[eN_k],
                                       mom_u_acc,
                                       dmom_u_acc_u,
                                       mom_v_acc,
                                       dmom_v_acc_v,
                                       mom_w_acc,
                                       dmom_w_acc_w,
                                       mass_adv,
                                       dmass_adv_u,
                                       dmass_adv_v,
                                       dmass_adv_w,
                                       mom_u_adv,
                                       dmom_u_adv_u,
                                       dmom_u_adv_v,
                                       dmom_u_adv_w,
                                       mom_v_adv,
                                       dmom_v_adv_u,
                                       dmom_v_adv_v,
                                       dmom_v_adv_w,
                                       mom_w_adv,
                                       dmom_w_adv_u,
                                       dmom_w_adv_v,
                                       dmom_w_adv_w,
                                       mom_uu_diff_ten,
                                       mom_vv_diff_ten,
                                       mom_ww_diff_ten,
                                       mom_uv_diff_ten,
                                       mom_uw_diff_ten,
                                       mom_vu_diff_ten,
                                       mom_vw_diff_ten,
                                       mom_wu_diff_ten,
                                       mom_wv_diff_ten,
                                       mom_u_source,
                                       mom_v_source,
                                       mom_w_source,
                                       mom_u_ham,
                                       dmom_u_ham_grad_p,
                                       dmom_u_ham_grad_u,
                                       dmom_u_ham_u,
                                       dmom_u_ham_v,
                                       dmom_u_ham_w,
                                       mom_v_ham,
                                       dmom_v_ham_grad_p,
                                       dmom_v_ham_grad_v,
                                       dmom_v_ham_u,
                                       dmom_v_ham_v,
                                       dmom_v_ham_w,
                                       mom_w_ham,
                                       dmom_w_ham_grad_p,
                                       dmom_w_ham_grad_w,
                                       dmom_w_ham_u,
                                       dmom_w_ham_v,
                                       dmom_w_ham_w,
                                       0.0,
                                       0.0,
                                       0.0);
                  mass_source = q_mass_source.data()[eN_k];
                  updateDarcyForchheimerTerms_Ergun(NONCONSERVATIVE_FORM,
                                                    /* linearDragFactor, */
                                                    /* nonlinearDragFactor, */
                                                    /* porosity, */
                                                    /* meanGrainSize, */
                                                    q_dragAlpha.data()[eN_k],
                                                    q_dragBeta.data()[eN_k],
                                                    eps_rho,
                                                    eps_mu,
                                                    rho_0,
                                                    nu_0,
                                                    rho_1,
                                                    nu_1,
                                                    useVF,
                                                    vf.data()[eN_k],
                                                    phi.data()[eN_k],
                                                    u,
                                                    v,
                                                    w,
                                                    q_velocity_sge.data()[eN_k_nSpace+0],
                                                    q_velocity_sge.data()[eN_k_nSpace+1],
                                                    q_velocity_sge.data()[eN_k_nSpace+1],//dummy entry for 2D
                                                    eps_porous.data()[elementFlags.data()[eN]],
                                                    phi_porous.data()[eN_k],
                                                    q_velocity_porous.data()[eN_k_nSpace+0],
                                                    q_velocity_porous.data()[eN_k_nSpace+1],
                                                    q_velocity_porous.data()[eN_k_nSpace+1],//dummy entry for 2D
                                                    mom_u_source,
                                                    mom_v_source,
                                                    mom_w_source,
                                                    dmom_u_source,
                                                    dmom_v_source,
                                                    dmom_w_source);

                  //Turbulence closure model
                  if (turbulenceClosureModel >= 3)
                    {
                      const double c_mu = 0.09;//mwf hack
                      updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                              turbulenceClosureModel,
                                              eps_rho,
                                              eps_mu,
                                              rho_0,
                                              nu_0,
                                              rho_1,
                                              nu_1,
                                              useVF,
                                              vf.data()[eN_k],
                                              phi.data()[eN_k],
                                              porosity,
                                              c_mu, //mwf hack
                                              q_turb_var_0.data()[eN_k],
                                              q_turb_var_1.data()[eN_k],
                                              &q_turb_var_grad_0.data()[eN_k_nSpace],
                                              eddy_viscosity,
                                              mom_uu_diff_ten,
                                              mom_vv_diff_ten,
                                              mom_ww_diff_ten,
                                              mom_uv_diff_ten,
                                              mom_uw_diff_ten,
                                              mom_vu_diff_ten,
                                              mom_vw_diff_ten,
                                              mom_wu_diff_ten,
                                              mom_wv_diff_ten,
                                              mom_u_source,
                                              mom_v_source,
                                              mom_w_source);

                    }
                  //
                  //
                  //moving mesh
                  //
                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_u_ham -= MOVING_DOMAIN*dmom_u_acc_u*(grad_u[0]*xt + grad_u[1]*yt);
                      dmom_u_ham_grad_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
                      dmom_u_ham_grad_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
                    }
                  else
                    {
                      mom_u_adv[0] -= MOVING_DOMAIN*mom_u_acc*xt;
                      mom_u_adv[1] -= MOVING_DOMAIN*mom_u_acc*yt;
                      dmom_u_adv_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
                      dmom_u_adv_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
                    }

                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_v_ham -= MOVING_DOMAIN*dmom_v_acc_v*(grad_v[0]*xt + grad_v[1]*yt);
                      dmom_v_ham_grad_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
                      dmom_v_ham_grad_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
                    }
                  else
                    {
                      mom_v_adv[0] -= MOVING_DOMAIN*mom_v_acc*xt;
                      mom_v_adv[1] -= MOVING_DOMAIN*mom_v_acc*yt;
                      dmom_v_adv_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
                      dmom_v_adv_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
                    }
                  //
                  //calculate time derivatives
                  //
                  ck.bdf(alphaBDF,
                         q_mom_u_acc_beta_bdf.data()[eN_k]*q_dV_last.data()[eN_k]/dV,
                         mom_u_acc,
                         dmom_u_acc_u,
                         mom_u_acc_t,
                         dmom_u_acc_u_t);
                  ck.bdf(alphaBDF,
                         q_mom_v_acc_beta_bdf.data()[eN_k]*q_dV_last.data()[eN_k]/dV,
                         mom_v_acc,
                         dmom_v_acc_v,
                         mom_v_acc_t,
                         dmom_v_acc_v_t);
                  if (NONCONSERVATIVE_FORM > 0.0)
                    {
                      mom_u_acc_t *= dmom_u_acc_u;
                      mom_v_acc_t *= dmom_v_acc_v;
                    }
                  //
                  //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
                  //
                  if (NONCONSERVATIVE_FORM > 0.0)
		    {
                      dmom_adv_sge[0] = 0.0;
                      dmom_adv_sge[1] = 0.0;
                      dmom_ham_grad_sge[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
                      dmom_ham_grad_sge[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
                    }
                  else
                    {
                      dmom_adv_sge[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
                      dmom_adv_sge[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
                      dmom_ham_grad_sge[0] = 0.0;
                      dmom_ham_grad_sge[1] = 0.0;
                    }
                  double mv_tau[nSpace]=ZEROVEC;
                  mv_tau[0] = dmom_adv_sge[0] + dmom_ham_grad_sge[0];
                  mv_tau[1] = dmom_adv_sge[1] + dmom_ham_grad_sge[1];
                  //
                  //calculate strong residual
                  //
                  pdeResidual_p = ck.Advection_strong(dmass_adv_u,grad_u) +
                    ck.Advection_strong(dmass_adv_v,grad_v) +
                    DM2*MOVING_DOMAIN*ck.Reaction_strong(alphaBDF*(dV-q_dV_last.data()[eN_k])/dV - div_mesh_velocity) +
                    ck.Reaction_strong(mass_source);

                  pdeResidual_u = ck.Mass_strong(mom_u_acc_t) +
                    ck.Advection_strong(dmom_adv_sge,grad_u) +
                    ck.Hamiltonian_strong(dmom_ham_grad_sge,grad_u) +
                    ck.Hamiltonian_strong(dmom_u_ham_grad_p,grad_p) +
                    ck.Reaction_strong(mom_u_source) -
                    ck.Reaction_strong(dmom_u_acc_u*u*div_mesh_velocity);

                  pdeResidual_v = ck.Mass_strong(mom_v_acc_t) +
                    ck.Advection_strong(dmom_adv_sge,grad_v) +
                    ck.Hamiltonian_strong(dmom_ham_grad_sge,grad_v) +
                    ck.Hamiltonian_strong(dmom_v_ham_grad_p,grad_p) +
                    ck.Reaction_strong(mom_v_source)  -
                    ck.Reaction_strong(dmom_v_acc_v*v*div_mesh_velocity);

                  //calculate the Jacobian of strong residual
                  for (int j=0;j<nDOF_v_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace;
                      dpdeResidual_p_u[j]=ck.AdvectionJacobian_strong(dmass_adv_u,&vel_grad_trial_ib[j_nSpace]);
                      dpdeResidual_p_v[j]=ck.AdvectionJacobian_strong(dmass_adv_v,&vel_grad_trial_ib[j_nSpace]);
                      dpdeResidual_u_u[j]=ck.MassJacobian_strong(dmom_u_acc_u_t,vel_trial[j]) +
                        ck.HamiltonianJacobian_strong(dmom_ham_grad_sge,&vel_grad_trial_ib[j_nSpace]) +
                        ck.AdvectionJacobian_strong(dmom_adv_sge,&vel_grad_trial_ib[j_nSpace]) -
                        ck.ReactionJacobian_strong(dmom_u_acc_u*div_mesh_velocity,vel_trial[j]);
                      dpdeResidual_v_v[j]=ck.MassJacobian_strong(dmom_v_acc_v_t,vel_trial[j]) +
                        ck.HamiltonianJacobian_strong(dmom_ham_grad_sge,&vel_grad_trial_ib[j_nSpace]) +
                        ck.AdvectionJacobian_strong(dmom_adv_sge,&vel_grad_trial_ib[j_nSpace]) -
                        ck.ReactionJacobian_strong(dmom_v_acc_v*div_mesh_velocity,vel_trial[j]);
                      //VRANS account for drag terms, diagonal only here ... decide if need off diagonal terms too
                      dpdeResidual_u_u[j]+= ck.ReactionJacobian_strong(dmom_u_source[0],vel_trial[j]);
                      dpdeResidual_v_v[j]+= ck.ReactionJacobian_strong(dmom_v_source[1],vel_trial[j]);
                    }
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace;
                      dpdeResidual_u_p[j]=ck.HamiltonianJacobian_strong(dmom_u_ham_grad_p,&p_grad_trial_ib[j_nSpace]);
                      dpdeResidual_v_p[j]=ck.HamiltonianJacobian_strong(dmom_v_ham_grad_p,&p_grad_trial_ib[j_nSpace]);
                    }
                  //calculate tau and tau*Res
                  //add contributions from mass and sourced terms
                  double tmpR=dmom_u_acc_u_t + dmom_u_source[0];
                  calculateSubgridError_tau(hFactor,
                                            elementDiameter.data()[eN],
                                            tmpR,//dmom_u_acc_u_t,
                                            dmom_u_acc_u,
                                            mv_tau,//dmom_adv_sge,
                                            mom_uu_diff_ten[1],
                                            dmom_u_ham_grad_p[0],
                                            tau_v0,
                                            tau_p0,
                                            q_cfl.data()[eN_k]);

                  calculateSubgridError_tau(Ct_sge,Cd_sge,
                                            G,G_dd_G,tr_G,
                                            tmpR,//dmom_u_acc_u_t,
                                            mv_tau,//dmom_adv_sge,
                                            mom_uu_diff_ten[1],
                                            dmom_u_ham_grad_p[0],
                                            tau_v1,
                                            tau_p1,
                                            q_cfl.data()[eN_k]);

                  tau_v = useMetrics*tau_v1+(1.0-useMetrics)*tau_v0;
                  tau_p = useMetrics*tau_p1+(1.0-useMetrics)*tau_p0;
                  
                  calculateSubgridError_tauRes(tau_p,
                                               tau_v,
                                               pdeResidual_p,
                                               pdeResidual_u,
                                               pdeResidual_v,
                                               pdeResidual_w,
                                               subgridError_p,
                                               subgridError_u,
                                               subgridError_v,
                                               subgridError_w);

                  calculateSubgridErrorDerivatives_tauRes(tau_p,
                                                          tau_v,
                                                          dpdeResidual_p_u,
                                                          dpdeResidual_p_v,
                                                          dpdeResidual_p_w,
                                                          dpdeResidual_u_p,
                                                          dpdeResidual_u_u,
                                                          dpdeResidual_v_p,
                                                          dpdeResidual_v_v,
                                                          dpdeResidual_w_p,
                                                          dpdeResidual_w_w,
                                                          dsubgridError_p_u,
                                                          dsubgridError_p_v,
                                                          dsubgridError_p_w,
                                                          dsubgridError_u_p,
                                                          dsubgridError_u_u,
                                                          dsubgridError_v_p,
                                                          dsubgridError_v_v,
                                                          dsubgridError_w_p,
                                                          dsubgridError_w_w);
                  // velocity used in adjoint (VMS or RBLES, with or without lagging the grid scale velocity)
                  dmom_adv_star[0] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+0] - MOVING_DOMAIN*xt + useRBLES*subgridError_u);
                  dmom_adv_star[1] = inertial_term*dmom_u_acc_u*(q_velocity_sge.data()[eN_k_nSpace+1] - MOVING_DOMAIN*yt + useRBLES*subgridError_v);

                  //calculate the adjoint times the test functions
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
                      Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
                    }
                  //calculate the adjoint times the test functions
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      Lstar_u_u[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
                      Lstar_v_v[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
                      Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
                      Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
                      //VRANS account for drag terms, diagonal only here ... decide if need off diagonal terms too
                      Lstar_u_u[i]+=ck.Reaction_adjoint(dmom_u_source[0],vel_test_dV[i]);
                      Lstar_v_v[i]+=ck.Reaction_adjoint(dmom_v_source[1],vel_test_dV[i]);
                    }

                  // Assumes non-lagged subgrid velocity
                  dmom_u_adv_u[0] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_u);
                  dmom_u_adv_u[1] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_v);

                  dmom_v_adv_v[0] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_u);
                  dmom_v_adv_v[1] += inertial_term*dmom_u_acc_u*(useRBLES*subgridError_v);

                  if(nParticles > 0)
                    {
                      //cek todo, this needs to be fixed for not exact
                      double level_set_normal[nSpace];
                      double sign=0.0;
                      if (gf_s.useExact)
                        {
                          double norm_exact=0.0,norm_cut=0.0;
                          if (use_ball_as_particle)
                            {
                              for (int I=0;I<nSpace;I++)
                                {
                                  sign += ball_n[I]*gf_s.get_normal()[I];
                                  level_set_normal[I] = gf_s.get_normal()[I];
                                  norm_cut += level_set_normal[I]*level_set_normal[I];
                                  norm_exact += ball_n[I]*ball_n[I];
                                }
                            }
                          else
                            {
                              for (int I=0;I<nSpace;I++)
                                {
                                  sign += particle_signed_distance_normals.data()[eN_k_3d+I]*gf_s.get_normal()[I];
                                  level_set_normal[I] = gf_s.get_normal()[I];
                                  norm_cut += level_set_normal[I]*level_set_normal[I];
                                  norm_exact += particle_signed_distance_normals.data()[eN_k_3d+I]*particle_signed_distance_normals.data()[eN_k_3d+I];
                                }
                            }
                          assert(std::fabs(1.0-norm_cut) < 1.0e-8);
                          assert(std::fabs(1.0-norm_exact) < 1.0e-8);
                          if (sign < 0.0)
                            for (int I=0;I<nSpace;I++)
                              level_set_normal[I]*=-1.0;
                          /* if(icase_s==0)// && (1.0-sign*sign) > 1.0e-3) */
                          /*   { */
                          /*     std::cout<<"phi normal and cut normal divergent "<<eN<<'\t'<<k<<std::endl; */
                          /*     for (int I=0;I<nSpace;I++) */
                          /*       std::cout<<level_set_normal[I]<<'\t'<<particle_signed_distance_normals[eN_k_3d+I]<<std::endl; */
                          /*   } */
                        }
		      else
			{
			  if (use_ball_as_particle)
			    for (int I=0;I<nSpace;I++)
			      level_set_normal[I] = ball_n[I];
			  else
			    for (int I=0;I<nSpace;I++)
			      level_set_normal[I] = particle_signed_distance_normals.data()[eN_k_3d+I];
			}
                      updateSolidParticleTerms(NONCONSERVATIVE_FORM,
                                               eN < nElements_owned,
                                               particle_nitsche,
                                               dV,
                                               nParticles,
                                               nQuadraturePoints_global,
                                               &particle_signed_distances.data()[eN_k],
                                               level_set_normal,
                                               &particle_velocities.data()[eN_k_3d],
                                               particle_centroids.data(),
                                               use_ball_as_particle,
                                               ball_center.data(),
                                               ball_radius.data(),
                                               ball_velocity.data(),
                                               ball_angular_velocity.data(),
                                               ball_density.data(),
                                               porosity,
                                               particle_penalty_constant/h_phi,//penalty,
                                               particle_alpha,
                                               particle_beta,
                                               eps_rho,
                                               eps_mu,
                                               rho_0,
                                               nu_0,
                                               rho_1,
                                               nu_1,
                                               useVF,
                                               vf.data()[eN_k],
                                               phi.data()[eN_k],
                                               x,
                                               y,
                                               z,
                                               p,
                                               u,
                                               v,
                                               w,
                                               q_velocity_sge.data()[eN_k_nSpace+0],
                                               q_velocity_sge.data()[eN_k_nSpace+1],
                                               q_velocity_sge.data()[eN_k_nSpace+1],//dummy entry for 2D
                                               particle_eps,
                                               grad_u,
                                               grad_v,
                                               grad_w,
                                               mass_source_s,
                                               mom_u_source_s,
                                               mom_v_source_s,
                                               mom_w_source_s,
                                               dmom_u_source_s,
                                               dmom_v_source_s,
                                               dmom_w_source_s,
                                               mom_u_adv_s,
                                               mom_v_adv_s,
                                               mom_w_adv_s,
                                               dmom_u_adv_u_s,
                                               dmom_v_adv_v_s,
                                               dmom_w_adv_w_s,
                                               mom_u_ham_s,
                                               dmom_u_ham_grad_u_s,
                                               dmom_u_ham_grad_v_s,
                                               dmom_u_ham_u_s,
                                               dmom_u_ham_v_s,
                                               dmom_u_ham_w_s,
                                               mom_v_ham_s,
                                               dmom_v_ham_grad_u_s,
                                               dmom_v_ham_grad_v_s,
                                               dmom_v_ham_u_s,
                                               dmom_v_ham_v_s,
                                               dmom_v_ham_w_s,
                                               mom_w_ham_s,
                                               dmom_w_ham_grad_w_s,
                                               dmom_w_ham_u_s,
                                               dmom_w_ham_v_s,
                                               dmom_w_ham_w_s,
                                               mass_ham_s,
                                               dmass_ham_u_s,
                                               dmass_ham_v_s,
                                               dmass_ham_w_s,
                                               &particle_netForces_tmp[0],
                                               &particle_netMoments_tmp[0],
                                               &particle_surfaceArea_tmp[0]);
                    }
                  //cek todo add RBLES terms consistent to residual modifications or ignore the partials w.r.t the additional RBLES terms
                  double H_f=1.0;
                  if (gf.useExact && icase == 0)
                    {
                      if (fluid_phase == 0)
                        H_f = gf.ImH(0.,0.);
                      else
                        H_f = gf.H(0.,0.);
                    }
                  else
                    {
                      assert(fluid_phase == 0);
                      H_f = 1.0;
                    }
                  for(int i=0;i<nDOF_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      for(int j=0;j<nDOF_trial_element;j++)
                        {
                          register int j_nSpace = j*nSpace;
                          if (nDOF_test_element == nDOF_v_trial_element)
                            {
                              elementJacobian_p_p[i][j] += H_s*H_f*((1-PRESSURE_PROJECTION_STABILIZATION)*ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_p[i]) +
                                                                    (1-PRESSURE_PROJECTION_STABILIZATION)*ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_p[i]) +
                                                                    PRESSURE_PROJECTION_STABILIZATION*ck.pressureProjection_weak(mom_uu_diff_ten[1], p_trial[j], 1./3., p_test_ref.data()[k*nDOF_test_element +i],dV));
                            }
                        }
                    }
                  for(int i=0;i<nDOF_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      for(int j=0;j<nDOF_v_trial_element;j++)
                        {
                          register int j_nSpace = j*nSpace;
                          elementJacobian_p_u[i][j] += H_s*H_f*(ck.AdvectionJacobian_weak(dmass_adv_u,vel_trial[j],&p_grad_test_dV[i_nSpace]) +
                                                                ck.MassJacobian_weak(dmass_ham_u,vel_trial[j],p_test_dV[i]));
                          elementJacobian_p_v[i][j] += H_s*H_f*(ck.AdvectionJacobian_weak(dmass_adv_v,vel_trial[j],&p_grad_test_dV[i_nSpace]) +
                                                                ck.MassJacobian_weak(dmass_ham_v,vel_trial[j],p_test_dV[i]));
                          if (nDOF_test_element == nDOF_v_trial_element)
                            {
                              elementJacobian_p_u[i][j] += H_s*H_f*(1-PRESSURE_PROJECTION_STABILIZATION)*ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_p[i]);
                              elementJacobian_p_v[i][j] += H_s*H_f*(1-PRESSURE_PROJECTION_STABILIZATION)*ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_p[i]);
                            }
                        }
                    }
                  for(int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      for(int j=0;j<nDOF_trial_element;j++)
                        {
                          register int j_nSpace = j*nSpace;
                          elementJacobian_u_p[i][j] += H_s*H_f*(ck.HamiltonianJacobian_weak(dmom_u_ham_grad_p,&p_grad_trial_ib[j_nSpace],vel_test_dV[i])+
                                                                MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_u[i]));
                          elementJacobian_v_p[i][j] += H_s*H_f*(ck.HamiltonianJacobian_weak(dmom_v_ham_grad_p,&p_grad_trial_ib[j_nSpace],vel_test_dV[i])+
                                                                MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_v[i]));
                        }
                    }
                  for(int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int i_nSpace = i*nSpace;
                      for(int j=0;j<nDOF_v_trial_element;j++)
                        {
                          register int j_nSpace = j*nSpace;
                          elementJacobian_u_u[i][j] += H_s*H_f*(ck.MassJacobian_weak(dmom_u_acc_u_t,vel_trial[j],vel_test_dV[i]) +
                                                                ck.MassJacobian_weak(dmom_u_ham_u,vel_trial[j],vel_test_dV[i]) + //cek hack for nonlinear hamiltonian
                                                                ck.HamiltonianJacobian_weak(dmom_u_ham_grad_u,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_u_adv_u,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.SimpleDiffusionJacobian_weak(sdInfo_u_u_rowptr.data(),sdInfo_u_u_colind.data(),mom_uu_diff_ten,&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_u_source[0]+NONCONSERVATIVE_FORM*dmom_u_acc_u*div_mesh_velocity,vel_trial[j],vel_test_dV[i]) +
                                                                MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_u[i]) +
                                                                MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_u[i]) +
                                                                ck.NumericalDiffusionJacobian(q_numDiff_u_last.data()[eN_k],&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]));
                          elementJacobian_u_v[i][j] += H_s*H_f*(ck.HamiltonianJacobian_weak(dmom_u_ham_grad_v,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_u_adv_v,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.MassJacobian_weak(dmom_u_ham_v,vel_trial[j],vel_test_dV[i]) + //cek hack for nonlinear hamiltonian
                                                                ck.SimpleDiffusionJacobian_weak(sdInfo_u_v_rowptr.data(),sdInfo_u_v_colind.data(),mom_uv_diff_ten,&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_u_source[1],vel_trial[j],vel_test_dV[i]) +
                                                                MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_u[i]));                          
                          elementJacobian_v_u[i][j] += H_s*H_f*(ck.HamiltonianJacobian_weak(dmom_v_ham_grad_u,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_v_adv_u,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.MassJacobian_weak(dmom_v_ham_u,vel_trial[j],vel_test_dV[i]) + //cek hack for nonlinear hamiltonian
                                                                ck.SimpleDiffusionJacobian_weak(sdInfo_v_u_rowptr.data(),sdInfo_v_u_colind.data(),mom_vu_diff_ten,&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_v_source[0],vel_trial[j],vel_test_dV[i]) +
                                                                MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_v[i]));
                          elementJacobian_v_v[i][j] += H_s*H_f*(ck.MassJacobian_weak(dmom_v_acc_v_t,vel_trial[j],vel_test_dV[i]) +
                                                                ck.MassJacobian_weak(dmom_v_ham_v,vel_trial[j],vel_test_dV[i]) + //cek hack for nonlinear hamiltonian
                                                                ck.HamiltonianJacobian_weak(dmom_v_ham_grad_v,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_v_adv_v,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.SimpleDiffusionJacobian_weak(sdInfo_v_v_rowptr.data(),sdInfo_v_v_colind.data(),mom_vv_diff_ten,&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_v_source[1]+NONCONSERVATIVE_FORM*dmom_v_acc_v*div_mesh_velocity,vel_trial[j],vel_test_dV[i]) +
                                                                MOMENTUM_SGE*PRESSURE_SGE*ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_v[i]) +
                                                                MOMENTUM_SGE*VELOCITY_SGE*ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_v[i]) +
                                                                ck.NumericalDiffusionJacobian(q_numDiff_v_last.data()[eN_k],&vel_grad_trial_ib[j_nSpace],&vel_grad_test_dV[i_nSpace]));
                        }//j
                    }//i
                  if (nParticles > 0)
                    {
                      for(int i=0;i<nDOF_v_test_element;i++)
                        {
                          register int i_nSpace = i*nSpace;
                          for(int j=0;j<nDOF_v_trial_element;j++)
                            {
                              register int j_nSpace = j*nSpace;
                              elementJacobian_u_u[i][j] += H_f*(ck.MassJacobian_weak(dmom_u_ham_u_s,vel_trial[j],vel_test_dV[i]) +
                                                                ck.HamiltonianJacobian_weak(dmom_u_ham_grad_u_s,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_u_adv_u_s,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_u_source_s[0],vel_trial[j],vel_test_dV[i]));
                              
                              elementJacobian_u_v[i][j] += H_f*(ck.HamiltonianJacobian_weak(dmom_u_ham_grad_v_s,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.ReactionJacobian_weak(dmom_u_source_s[1],vel_trial[j],vel_test_dV[i]));
                              
                              elementJacobian_v_u[i][j] += H_f*(ck.HamiltonianJacobian_weak(dmom_v_ham_grad_u_s,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.ReactionJacobian_weak(dmom_v_source_s[0],vel_trial[j],vel_test_dV[i]));
                              
                              elementJacobian_v_v[i][j] += H_f*(ck.MassJacobian_weak(dmom_v_ham_v_s,vel_trial[j],vel_test_dV[i]) +
                                                                ck.HamiltonianJacobian_weak(dmom_v_ham_grad_v_s,&vel_grad_trial_ib[j_nSpace],vel_test_dV[i]) +
                                                                ck.AdvectionJacobian_weak(dmom_v_adv_v_s,vel_trial[j],&vel_grad_test_dV[i_nSpace]) +
                                                                ck.ReactionJacobian_weak(dmom_v_source_s[1],vel_trial[j],vel_test_dV[i]));
                            }//j
                        }//i
                    }
                }//k
            }//fluid_phase
          //
          //load into element Jacobian into global Jacobian
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian.data()[csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_p_p.data()[eN_i_j]] += elementJacobian_p_p[i][j];
                }
            }
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_v_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_v_trial_element+j;
                  globalJacobian.data()[csrRowIndeces_p_u.data()[eN_i] + csrColumnOffsets_p_u.data()[eN_i_j]] += elementJacobian_p_u[i][j];
                  globalJacobian.data()[csrRowIndeces_p_v.data()[eN_i] + csrColumnOffsets_p_v.data()[eN_i_j]] += elementJacobian_p_v[i][j];
                }
            }
          for (int i=0;i<nDOF_v_test_element;i++)
            {
              register int eN_i = eN*nDOF_v_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian.data()[csrRowIndeces_u_p.data()[eN_i] + csrColumnOffsets_u_p.data()[eN_i_j]] += elementJacobian_u_p[i][j];
                  globalJacobian.data()[csrRowIndeces_v_p.data()[eN_i] + csrColumnOffsets_v_p.data()[eN_i_j]] += elementJacobian_v_p[i][j];
                }
            }
          for (int i=0;i<nDOF_v_test_element;i++)
            {
              register int eN_i = eN*nDOF_v_test_element+i;
              for (int j=0;j<nDOF_v_trial_element;j++)
                {
                  register int eN_i_j = eN_i*nDOF_v_trial_element+j;
                  globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i][j];
                  globalJacobian.data()[csrRowIndeces_u_v.data()[eN_i] + csrColumnOffsets_u_v.data()[eN_i_j]] += elementJacobian_u_v[i][j];

                  globalJacobian.data()[csrRowIndeces_v_u.data()[eN_i] + csrColumnOffsets_v_u.data()[eN_i_j]] += elementJacobian_v_u[i][j];
                  globalJacobian.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_v_v.data()[eN_i_j]] += elementJacobian_v_v[i][j];
                }//j
            }//i
        }//elements
      for (std::set<int>::iterator it=cutfem_boundaries.begin(); it!=cutfem_boundaries.end(); ++it)
        {
          std::map<int,double> DWp_Dn_jump,DW_Dn_jump;
          std::map<std::pair<int, int>, int> p_p_nz, u_u_nz, v_v_nz;
          register double gamma_cutfem=ghost_penalty_constant,gamma_cutfem_p=ghost_penalty_constant,h_cutfem=elementBoundaryDiameter.data()[*it];
          int eN_nDOF_v_trial_element  = elementBoundaryElementsArray.data()[(*it)*2+0]*nDOF_v_trial_element;
          //See Massing Schott Wall 2018
          //cek todo modify for two-fluids: rho_0 != rho_1
          double norm_v=0.0;
          for (int i=0;i<nDOF_v_trial_element;i++)//MSW18 is just on face, but this is easier
            {
              double u=u_old_dof.data()[vel_l2g.data()[eN_nDOF_v_trial_element+i]],
                v=v_old_dof.data()[vel_l2g.data()[eN_nDOF_v_trial_element+i]];
              norm_v=fmax(norm_v,sqrt(u*u+v*v));
            }
          double gamma_v_dim = rho_0*(nu_0 + norm_v*h_cutfem + alphaBDF*h_cutfem*h_cutfem);
          gamma_cutfem_p *= h_cutfem*h_cutfem/gamma_v_dim;
          if (NONCONSERVATIVE_FORM)
            gamma_cutfem*=gamma_v_dim;
          else
            gamma_cutfem*=(gamma_v_dim/rho_0);
          for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register double Dp_Dn_jump=0.0, Du_Dn_jump=0.0, Dv_Dn_jump=0.0,dS;
              for (int eN_side=0;eN_side < 2; eN_side++)
                {
                  register int ebN = *it,
                    eN  = elementBoundaryElementsArray.data()[ebN*2+eN_side];
                  for (int i=0;i<nDOF_test_element;i++)
                    DWp_Dn_jump[p_l2g.data()[eN*nDOF_test_element+i]] = 0.0;
                  for (int i=0;i<nDOF_v_test_element;i++)
                    DW_Dn_jump[vel_l2g.data()[eN*nDOF_v_test_element+i]] = 0.0;
                }
              for (int eN_side=0;eN_side < 2; eN_side++)
                {
                  register int ebN = *it,
                    eN  = elementBoundaryElementsArray.data()[ebN*2+eN_side],
                    ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+eN_side],
                    eN_nDOF_trial_element = eN*nDOF_trial_element,
                    eN_nDOF_v_trial_element = eN*nDOF_v_trial_element,
                    ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                    ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                  register double p_int=0.0,
                    u_int=0.0,
                    v_int=0.0,
                    grad_p_int[nSpace]=ZEROVEC,
                    grad_u_int[nSpace]=ZEROVEC,
                    grad_v_int[nSpace]=ZEROVEC,
                    jac_int[nSpace*nSpace],
                    jacDet_int,
                    jacInv_int[nSpace*nSpace],
                    boundaryJac[nSpace*(nSpace-1)],
                    metricTensor[(nSpace-1)*(nSpace-1)],
                    metricTensorDetSqrt,
                    p_test_dS[nDOF_test_element],vel_test_dS[nDOF_v_test_element],
                    p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_v_trial_element*nSpace],
                    p_grad_test_dS[nDOF_trial_element*nSpace],vel_grad_test_dS[nDOF_v_trial_element*nSpace],
                    normal[nSpace],x_int,y_int,z_int,xt_int,yt_int,zt_int,integralScaling,
                    G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty,
                    force_x,force_y,force_z,force_p_x,force_p_y,force_p_z,force_v_x,force_v_y,force_v_z,r_x,r_y,r_z;
                  //compute information about mapping from reference element to physical element
                  ck.calculateMapping_elementBoundary(eN,
                                                      ebN_local,
                                                      kb,
                                                      ebN_local_kb,
                                                      mesh_dof.data(),
                                                      mesh_l2g.data(),
                                                      mesh_trial_trace_ref.data(),
                                                      mesh_grad_trial_trace_ref.data(),
                                                      boundaryJac_ref.data(),
                                                      jac_int,
                                                      jacDet_int,
                                                      jacInv_int,
                                                      boundaryJac,
                                                      metricTensor,
                                                      metricTensorDetSqrt,
                                                      normal_ref.data(),
                                                      normal,
                                                      x_int,y_int,z_int);
                  //todo: check that physical coordinates match
                  ck.calculateMappingVelocity_elementBoundary(eN,
                                                              ebN_local,
                                                              kb,
                                                              ebN_local_kb,
                                                              mesh_velocity_dof.data(),
                                                              mesh_l2g.data(),
                                                              mesh_trial_trace_ref.data(),
                                                              xt_int,yt_int,zt_int,
                                                              normal,
                                                              boundaryJac,
                                                              metricTensor,
                                                              integralScaling);
                  dS = metricTensorDetSqrt*dS_ref.data()[kb];
                  //compute shape and solution information
                  //shape
                  ck.gradTrialFromRef(&p_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_int,p_grad_trial_trace);
                  ck_v.gradTrialFromRef(&vel_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_v_trial_element],jacInv_int,vel_grad_trial_trace);
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      int eN_i = eN*nDOF_test_element + i;
                      for (int I=0;I<nSpace;I++)
                        DWp_Dn_jump[p_l2g.data()[eN_i]] += p_grad_trial_trace[i*nSpace+I]*normal[I];
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      int eN_i = eN*nDOF_v_test_element + i;
                      for (int I=0;I<nSpace;I++)
                        DW_Dn_jump[vel_l2g.data()[eN_i]] += vel_grad_trial_trace[i*nSpace+I]*normal[I];
                    }
                }//eN_side
              for (int eN_side=0;eN_side < 2; eN_side++)
                {
                  register int ebN = *it,
                    eN  = elementBoundaryElementsArray.data()[ebN*2+eN_side];
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_test_element+i;
                      for (int eN_side2=0;eN_side2 < 2; eN_side2++)
                        {
                          register int eN2  = elementBoundaryElementsArray.data()[ebN*2+eN_side2];
                          for (int j=0;j<nDOF_test_element;j++)
                            {
                              int eN_i_j = eN_i*nDOF_test_element + j;
                              int eN2_j = eN2*nDOF_test_element + j;
                              register int ebN_i_j = ebN*4*nDOF_test_X_trial_element +
                                eN_side*2*nDOF_test_X_trial_element +
                                eN_side2*nDOF_test_X_trial_element +
                                i*nDOF_trial_element +
                                j;
                              std::pair<int,int> ij = std::make_pair(p_l2g.data()[eN_i], p_l2g.data()[eN2_j]);
                              if (p_p_nz.count(ij))
                                {
                                  assert(p_p_nz[ij] == csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_eb_p_p.data()[ebN_i_j]);
                                }
                              else
                                p_p_nz[ij] =  csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_eb_p_p.data()[ebN_i_j];
                            }
                        }
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_v_test_element+i;
                      for (int eN_side2=0;eN_side2 < 2; eN_side2++)
                        {
                          register int eN2  = elementBoundaryElementsArray.data()[ebN*2+eN_side2];
                          for (int j=0;j<nDOF_v_test_element;j++)
                            {
                              int eN_i_j = eN_i*nDOF_v_test_element + j;
                              int eN2_j = eN2*nDOF_v_test_element + j;
                              register int ebN_i_j = ebN*4*nDOF_v_test_X_v_trial_element +
                                eN_side*2*nDOF_v_test_X_v_trial_element +
                                eN_side2*nDOF_v_test_X_v_trial_element +
                                i*nDOF_v_trial_element +
                                j;
                              std::pair<int,int> ij = std::make_pair(vel_l2g.data()[eN_i], vel_l2g.data()[eN2_j]);
                              if (u_u_nz.count(ij))
                                {
                                  assert(u_u_nz[ij] == csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_eb_u_u.data()[ebN_i_j]);
                                }
                              else
                                u_u_nz[ij] =  csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_eb_u_u.data()[ebN_i_j];
                              if (v_v_nz.count(ij))
                                {
                                  assert(v_v_nz[ij] == csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_eb_v_v.data()[ebN_i_j]);
                                }
                              else
                                v_v_nz[ij] =  csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_eb_v_v.data()[ebN_i_j];
                            }
                        }
                    }
                }
              for (std::map<int,double>::iterator Wi_it=DWp_Dn_jump.begin(); Wi_it!=DWp_Dn_jump.end(); ++Wi_it)
                for (std::map<int,double>::iterator Wj_it=DWp_Dn_jump.begin(); Wj_it!=DWp_Dn_jump.end(); ++Wj_it)
                  {
                    int i_global = Wi_it->first,
                      j_global = Wj_it->first;
                    double DWp_Dn_jump_i = Wi_it->second,
                      DWp_Dn_jump_j = Wj_it->second;
                    std::pair<int,int> ij = std::make_pair(i_global, j_global);
                    globalJacobian.data()[p_p_nz.at(ij)] += gamma_cutfem_p*h_cutfem*DWp_Dn_jump_j*DWp_Dn_jump_i*dS;
                  }//i,j
              for (std::map<int,double>::iterator Wi_it=DW_Dn_jump.begin(); Wi_it!=DW_Dn_jump.end(); ++Wi_it)
                for (std::map<int,double>::iterator Wj_it=DW_Dn_jump.begin(); Wj_it!=DW_Dn_jump.end(); ++Wj_it)
                  {
                    int i_global = Wi_it->first,
                      j_global = Wj_it->first;
                    double DW_Dn_jump_i = Wi_it->second,
                      DW_Dn_jump_j = Wj_it->second;
                    std::pair<int,int> ij = std::make_pair(i_global, j_global);
                    globalJacobian.data()[u_u_nz.at(ij)] += gamma_cutfem*h_cutfem*DW_Dn_jump_j*DW_Dn_jump_i*dS;
                    globalJacobian.data()[v_v_nz.at(ij)] += gamma_cutfem*h_cutfem*DW_Dn_jump_j*DW_Dn_jump_i*dS;
                  }//i,j
            }//kb
        }//cutfem element boundaries
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
        {
          register int ebN = exteriorElementBoundariesArray.data()[ebNE],
            eN  = elementBoundaryElementsArray.data()[ebN*2+0],
            eN_nDOF_trial_element = eN*nDOF_trial_element,
            eN_nDOF_v_trial_element = eN*nDOF_v_trial_element,
            ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0];
	  if (boundaryFlags[ebN] < 1)
	    continue;
          register double eps_rho,eps_mu;
          double element_phi[nDOF_mesh_trial_element], element_phi_s[nDOF_mesh_trial_element];
          for (int j=0;j<nDOF_mesh_trial_element;j++)
            {
              register int eN_j = eN*nDOF_mesh_trial_element+j;
              element_phi[j] = phi_nodes.data()[p_l2g.data()[eN_j]];
              element_phi_s[j] = phi_solid_nodes.data()[p_l2g.data()[eN_j]];
            }
          double element_nodes[nDOF_mesh_trial_element*3];
          for (int i=0;i<nDOF_mesh_trial_element;i++)
            {
              register int eN_i=eN*nDOF_mesh_trial_element+i;
              for(int I=0;I<3;I++)
                element_nodes[i*3 + I] = mesh_dof[mesh_l2g.data()[eN_i]*3 + I];
            }//i
          double mesh_dof_ref[nDOF_mesh_trial_element*3]={0.,0.,0.,1.,0.,0.,0.,1.,0.};
          double xb_ref_calc[nQuadraturePoints_elementBoundary*3];
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              double x=0.0,y=0.0,z=0.0;
              for (int j=0;j<nDOF_mesh_trial_element;j++)
                {
                  int ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb;
                  int ebN_local_kb_j = ebN_local_kb*nDOF_mesh_trial_element+j;
                  x += mesh_dof_ref[j*3+0]*mesh_trial_trace_ref.data()[ebN_local_kb_j]; 
                  y += mesh_dof_ref[j*3+1]*mesh_trial_trace_ref.data()[ebN_local_kb_j]; 
                  z += mesh_dof_ref[j*3+2]*mesh_trial_trace_ref.data()[ebN_local_kb_j];
                }
              xb_ref_calc[3*kb+0] = x;
              xb_ref_calc[3*kb+1] = y;
              xb_ref_calc[3*kb+2] = z;
            }
          int icase_s = gf_s.calculate(element_phi_s, element_nodes, xb_ref_calc,true);
#ifdef IFEM
          int icase = gf.calculate(element_phi, element_nodes, xb_ref.data(), rho_1*nu_1, rho_0*nu_0,true,false);
#else
          int icase = gf.calculate(element_phi, element_nodes, xb_ref.data(), 1.0,1.0,true, false);
#endif
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                ebNE_kb_nSpace = ebNE_kb*nSpace,
                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                ebN_local_kb_nSpace = ebN_local_kb*nSpace;

              register double phi_s_ext=0.0,
                p_ext=0.0,
                u_ext=0.0,
                v_ext=0.0,
                w_ext=0.0,
                grad_p_ext[nSpace]=ZEROVEC,
                grad_u_ext[nSpace]=ZEROVEC,
                grad_v_ext[nSpace]=ZEROVEC,
                grad_w_ext[nSpace]=ZEROVEC,
                p_old=0.0,u_old=0.0,v_old=0.0,w_old=0.0,
                grad_p_old[nSpace]=ZEROVEC,grad_u_old[nSpace]=ZEROVEC,grad_v_old[nSpace]=ZEROVEC,grad_w_old[nSpace]=ZEROVEC,
                mom_u_acc_ext=0.0,
                dmom_u_acc_u_ext=0.0,
                mom_v_acc_ext=0.0,
                dmom_v_acc_v_ext=0.0,
                mom_w_acc_ext=0.0,
                dmom_w_acc_w_ext=0.0,
                mass_adv_ext[nSpace]=ZEROVEC,
                dmass_adv_u_ext[nSpace]=ZEROVEC,
                dmass_adv_v_ext[nSpace]=ZEROVEC,
                dmass_adv_w_ext[nSpace]=ZEROVEC,
                mom_u_adv_ext[nSpace]=ZEROVEC,
                dmom_u_adv_u_ext[nSpace]=ZEROVEC,
                dmom_u_adv_v_ext[nSpace]=ZEROVEC,
                dmom_u_adv_w_ext[nSpace]=ZEROVEC,
                mom_v_adv_ext[nSpace]=ZEROVEC,
                dmom_v_adv_u_ext[nSpace]=ZEROVEC,
                dmom_v_adv_v_ext[nSpace]=ZEROVEC,
                dmom_v_adv_w_ext[nSpace]=ZEROVEC,
                mom_w_adv_ext[nSpace]=ZEROVEC,
                dmom_w_adv_u_ext[nSpace]=ZEROVEC,
                dmom_w_adv_v_ext[nSpace]=ZEROVEC,
                dmom_w_adv_w_ext[nSpace]=ZEROVEC,
                mom_uu_diff_ten_ext[nSpace]=ZEROVEC,
                mom_vv_diff_ten_ext[nSpace]=ZEROVEC,
                mom_ww_diff_ten_ext[nSpace]=ZEROVEC,
                mom_uv_diff_ten_ext[1],
                mom_uw_diff_ten_ext[1],
                mom_vu_diff_ten_ext[1],
                mom_vw_diff_ten_ext[1],
                mom_wu_diff_ten_ext[1],
                mom_wv_diff_ten_ext[1],
                mom_u_source_ext=0.0,
                mom_v_source_ext=0.0,
                mom_w_source_ext=0.0,
                mom_u_ham_ext=0.0,
                dmom_u_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_u_ham_grad_u_ext[nSpace]=ZEROVEC,
                dmom_u_ham_u_ext=0.0,
                dmom_u_ham_v_ext=0.0,
                dmom_u_ham_w_ext=0.0,
                mom_v_ham_ext=0.0,
                dmom_v_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_v_ham_grad_v_ext[nSpace]=ZEROVEC,
                dmom_v_ham_u_ext=0.0,
                dmom_v_ham_v_ext=0.0,
                dmom_v_ham_w_ext=0.0,
                mom_w_ham_ext=0.0,
                dmom_w_ham_grad_p_ext[nSpace]=ZEROVEC,
                dmom_w_ham_grad_w_ext[nSpace]=ZEROVEC,
                dmom_w_ham_u_ext=0.0,
                dmom_w_ham_v_ext=0.0,
                dmom_w_ham_w_ext=0.0,
                dmom_u_adv_p_ext[nSpace]=ZEROVEC,
                dmom_v_adv_p_ext[nSpace]=ZEROVEC,
                dmom_w_adv_p_ext[nSpace]=ZEROVEC,
                dflux_mass_u_ext=0.0,
                dflux_mass_v_ext=0.0,
                dflux_mass_w_ext=0.0,
                dflux_mom_u_adv_p_ext=0.0,
                dflux_mom_u_adv_u_ext=0.0,
                dflux_mom_u_adv_v_ext=0.0,
                dflux_mom_u_adv_w_ext=0.0,
                dflux_mom_v_adv_p_ext=0.0,
                dflux_mom_v_adv_u_ext=0.0,
                dflux_mom_v_adv_v_ext=0.0,
                dflux_mom_v_adv_w_ext=0.0,
                dflux_mom_w_adv_p_ext=0.0,
                dflux_mom_w_adv_u_ext=0.0,
                dflux_mom_w_adv_v_ext=0.0,
                dflux_mom_w_adv_w_ext=0.0,
                bc_p_ext=0.0,
                bc_u_ext=0.0,
                bc_v_ext=0.0,
                bc_w_ext=0.0,
                bc_mom_u_acc_ext=0.0,
                bc_dmom_u_acc_u_ext=0.0,
                bc_mom_v_acc_ext=0.0,
                bc_dmom_v_acc_v_ext=0.0,
                bc_mom_w_acc_ext=0.0,
                bc_dmom_w_acc_w_ext=0.0,
                bc_mass_adv_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmass_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_u_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_u_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_v_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_v_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_w_adv_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_u_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_v_ext[nSpace]=ZEROVEC,
                bc_dmom_w_adv_w_ext[nSpace]=ZEROVEC,
                bc_mom_uu_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_vv_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_ww_diff_ten_ext[nSpace]=ZEROVEC,
                bc_mom_uv_diff_ten_ext[1],
                bc_mom_uw_diff_ten_ext[1],
                bc_mom_vu_diff_ten_ext[1],
                bc_mom_vw_diff_ten_ext[1],
                bc_mom_wu_diff_ten_ext[1],
                bc_mom_wv_diff_ten_ext[1],
                bc_mom_u_source_ext=0.0,
                bc_mom_v_source_ext=0.0,
                bc_mom_w_source_ext=0.0,
                bc_mom_u_ham_ext=0.0,
                bc_dmom_u_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_u_ham_grad_u_ext[nSpace]=ZEROVEC,
                bc_dmom_u_ham_u_ext=0.0,
                bc_dmom_u_ham_v_ext=0.0,
                bc_dmom_u_ham_w_ext=0.0,
                bc_mom_v_ham_ext=0.0,
                bc_dmom_v_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_v_ham_grad_v_ext[nSpace]=ZEROVEC,
                bc_dmom_v_ham_u_ext=0.0,
                bc_dmom_v_ham_v_ext=0.0,
                bc_dmom_v_ham_w_ext=0.0,
                bc_mom_w_ham_ext=0.0,
                bc_dmom_w_ham_grad_p_ext[nSpace]=ZEROVEC,
                bc_dmom_w_ham_grad_w_ext[nSpace]=ZEROVEC,
                bc_dmom_w_ham_u_ext=0.0,
                bc_dmom_w_ham_v_ext=0.0,
                bc_dmom_w_ham_w_ext=0.0,
                fluxJacobian_p_p[nDOF_trial_element],
                fluxJacobian_p_u[nDOF_v_trial_element],
                fluxJacobian_p_v[nDOF_v_trial_element],
                fluxJacobian_p_w[nDOF_v_trial_element],
                fluxJacobian_u_p[nDOF_trial_element],
                fluxJacobian_u_u[nDOF_v_trial_element],
                fluxJacobian_u_v[nDOF_v_trial_element],
                fluxJacobian_u_w[nDOF_v_trial_element],
                fluxJacobian_v_p[nDOF_trial_element],
                fluxJacobian_v_u[nDOF_v_trial_element],
                fluxJacobian_v_v[nDOF_v_trial_element],
                fluxJacobian_v_w[nDOF_v_trial_element],
                fluxJacobian_w_p[nDOF_trial_element],
                fluxJacobian_w_u[nDOF_v_trial_element],
                fluxJacobian_w_v[nDOF_v_trial_element],
                fluxJacobian_w_w[nDOF_v_trial_element],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                p_grad_trial_trace[nDOF_trial_element*nSpace],
                vel_grad_trial_trace[nDOF_v_trial_element*nSpace],
                dS,
                p_test_dS[nDOF_test_element],
                vel_test_dS[nDOF_v_test_element],
                normal[nSpace],
                x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                vel_grad_test_dS[nDOF_v_trial_element*nSpace],
                //VRANS
                porosity_ext,
                //
                G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty;
              gf_s.set_boundary_quad(kb);
              gf.set_boundary_quad(kb);
              ck.calculateMapping_elementBoundary(eN,
                                                  ebN_local,
                                                  kb,
                                                  ebN_local_kb,
                                                  mesh_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_trace_ref.data(),
                                                  mesh_grad_trial_trace_ref.data(),
                                                  boundaryJac_ref.data(),
                                                  jac_ext,
                                                  jacDet_ext,
                                                  jacInv_ext,
                                                  boundaryJac,
                                                  metricTensor,
                                                  metricTensorDetSqrt,
                                                  normal_ref.data(),
                                                  normal,
                                                  x_ext,y_ext,z_ext);
              ck.calculateMappingVelocity_elementBoundary(eN,
                                                          ebN_local,
                                                          kb,
                                                          ebN_local_kb,
                                                          mesh_velocity_dof.data(),
                                                          mesh_l2g.data(),
                                                          mesh_trial_trace_ref.data(),
                                                          xt_ext,yt_ext,zt_ext,
                                                          normal,
                                                          boundaryJac,
                                                          metricTensor,
                                                          integralScaling);
              //dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
              dS = metricTensorDetSqrt*dS_ref.data()[kb];
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              ck.calculateGScale(G,&ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],h_phi);

              eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
              eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);

              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&p_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
              ck_v.gradTrialFromRef(&vel_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_v_trial_element],jacInv_ext,vel_grad_trial_trace);
              //solution and gradients
              ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],p_ext);
              ck_v.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],u_ext);
              ck_v.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],v_ext);
              ck.valFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],p_old);
              ck_v.valFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],u_old);
              ck_v.valFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],&vel_trial_trace_ref.data()[ebN_local_kb*nDOF_v_test_element],v_old);
              ck.gradFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
              ck_v.gradFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_u_ext);
              ck_v.gradFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_v_ext);
              ck.gradFromDOF(p_old_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_old);
              ck_v.gradFromDOF(u_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_u_old);
              ck_v.gradFromDOF(v_old_dof.data(),&vel_l2g.data()[eN_nDOF_v_trial_element],vel_grad_trial_trace,grad_v_old);
              ck.valFromDOF(phi_solid_nodes.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],phi_s_ext);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_test_element;j++)
                {
                  p_test_dS[j] = p_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                }
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_v_test_element;j++)
                {
                  vel_test_dS[j] = vel_test_trace_ref.data()[ebN_local_kb*nDOF_v_test_element+j]*dS;
                  for (int I=0;I<nSpace;I++)
                    vel_grad_test_dS[j*nSpace+I] = vel_grad_trial_trace[j*nSpace+I]*dS;//assume test_j == trial_j
                }
              //
              //load the boundary values
              //
              bc_p_ext = isDOFBoundary_p.data()[ebNE_kb]*ebqe_bc_p_ext.data()[ebNE_kb]+(1-isDOFBoundary_p.data()[ebNE_kb])*p_ext;
              //bc values at moving boundaries are specified relative to boundary motion so we need to add it here
              bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*(ebqe_bc_u_ext.data()[ebNE_kb] + MOVING_DOMAIN*xt_ext) + (1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;
              bc_v_ext = isDOFBoundary_v.data()[ebNE_kb]*(ebqe_bc_v_ext.data()[ebNE_kb] + MOVING_DOMAIN*yt_ext) + (1-isDOFBoundary_v.data()[ebNE_kb])*v_ext;
              //VRANS
              porosity_ext = ebqe_porosity_ext.data()[ebNE_kb];
              //
              //calculate the internal and external trace of the pde coefficients
              //
              double eddy_viscosity_ext(0.),bc_eddy_viscosity_ext(0.);//not interested in saving boundary eddy viscosity for now
              if (use_ball_as_particle == 1 && nParticles > 0)
                {
                  get_distance_to_ball(nParticles, ball_center.data(), ball_radius.data(),x_ext,y_ext,z_ext,ebqe_phi_s.data()[ebNE_kb]);
                }
              //else distance_to_solids is updated in PreStep
              const double particle_eps  = particle_epsFact*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
              //cek needs to be fixed for two-phase ifem
              double H = (1.0-useVF)*gf.H(eps_rho,ebqe_phi_ext.data()[ebNE_kb]) + useVF*fmin(1.0,fmax(0.0,ebqe_vf_ext.data()[ebNE_kb]));
              double ImH = (1.0-useVF)*gf.ImH(eps_rho,ebqe_phi_ext.data()[ebNE_kb]) + useVF*(1.0-fmin(1.0,fmax(0.0,ebqe_vf_ext.data()[ebNE_kb])));
              double rho  = rho_0*ImH + rho_1*H;
              double nu  = nu_0*ImH + nu_1*H;
              //
              evaluateCoefficients(NONCONSERVATIVE_FORM,
                                   sigma,
                                   rho,
                                   nu,
                                   elementDiameter.data()[eN],
                                   smagorinskyConstant,
                                   turbulenceClosureModel,
                                   g.data(),
                                   useVF,
                                   ebqe_vf_ext.data()[ebNE_kb],
                                   ebqe_phi_ext.data()[ebNE_kb],
                                   &ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],
                                   ebqe_kappa_phi_ext.data()[ebNE_kb],
                                   //VRANS
                                   porosity_ext,
                                   //
                                   ebqe_phi_s.data()[ebNE_kb],
                                   p_old,
                                   u_old,
                                   v_old,
                                   w_old,
                                   grad_p_old,
                                   grad_u_old,
                                   grad_v_old,
                                   grad_w_old,
                                   p_ext,
                                   grad_p_ext,
                                   grad_u_ext,
                                   grad_v_ext,
                                   grad_w_ext,
                                   u_ext,
                                   v_ext,
                                   w_ext,
                                   LAG_LES,
                                   eddy_viscosity_ext,
                                   ebqe_eddy_viscosity_last.data()[ebNE_kb],
                                   mom_u_acc_ext,
                                   dmom_u_acc_u_ext,
                                   mom_v_acc_ext,
                                   dmom_v_acc_v_ext,
                                   mom_w_acc_ext,
                                   dmom_w_acc_w_ext,
                                   mass_adv_ext,
                                   dmass_adv_u_ext,
                                   dmass_adv_v_ext,
                                   dmass_adv_w_ext,
                                   mom_u_adv_ext,
                                   dmom_u_adv_u_ext,
                                   dmom_u_adv_v_ext,
                                   dmom_u_adv_w_ext,
                                   mom_v_adv_ext,
                                   dmom_v_adv_u_ext,
                                   dmom_v_adv_v_ext,
                                   dmom_v_adv_w_ext,
                                   mom_w_adv_ext,
                                   dmom_w_adv_u_ext,
                                   dmom_w_adv_v_ext,
                                   dmom_w_adv_w_ext,
                                   mom_uu_diff_ten_ext,
                                   mom_vv_diff_ten_ext,
                                   mom_ww_diff_ten_ext,
                                   mom_uv_diff_ten_ext,
                                   mom_uw_diff_ten_ext,
                                   mom_vu_diff_ten_ext,
                                   mom_vw_diff_ten_ext,
                                   mom_wu_diff_ten_ext,
                                   mom_wv_diff_ten_ext,
                                   mom_u_source_ext,
                                   mom_v_source_ext,
                                   mom_w_source_ext,
                                   mom_u_ham_ext,
                                   dmom_u_ham_grad_p_ext,
                                   dmom_u_ham_grad_u_ext,
                                   dmom_u_ham_u_ext,
                                   dmom_u_ham_v_ext,
                                   dmom_u_ham_w_ext,
                                   mom_v_ham_ext,
                                   dmom_v_ham_grad_p_ext,
                                   dmom_v_ham_grad_v_ext,
                                   dmom_v_ham_u_ext,
                                   dmom_v_ham_v_ext,
                                   dmom_v_ham_w_ext,
                                   mom_w_ham_ext,
                                   dmom_w_ham_grad_p_ext,
                                   dmom_w_ham_grad_w_ext,
                                   dmom_w_ham_u_ext,
                                   dmom_w_ham_v_ext,
                                   dmom_w_ham_w_ext,
                                   0.0,
                                   0.0,
                                   0.0);
              //cek needs to be fixed for two-phase ifem
              H = (1.0-useVF)*gf.H(eps_rho,bc_ebqe_phi_ext.data()[ebNE_kb]) + useVF*fmin(1.0,fmax(0.0,bc_ebqe_vf_ext.data()[ebNE_kb]));
              ImH = (1.0-useVF)*gf.ImH(eps_rho,bc_ebqe_phi_ext.data()[ebNE_kb]) + useVF*(1.0-fmin(1.0,fmax(0.0,bc_ebqe_vf_ext.data()[ebNE_kb])));
              rho  = rho_0*ImH + rho_1*H;
              nu  = nu_0*ImH + nu_1*H;
              //
              evaluateCoefficients(NONCONSERVATIVE_FORM,
                                   sigma,
                                   rho,
                                   nu,
                                   elementDiameter.data()[eN],
                                   smagorinskyConstant,
                                   turbulenceClosureModel,
                                   g.data(),
                                   useVF,
                                   bc_ebqe_vf_ext.data()[ebNE_kb],
                                   bc_ebqe_phi_ext.data()[ebNE_kb],
                                   &ebqe_normal_phi_ext.data()[ebNE_kb_nSpace],
                                   ebqe_kappa_phi_ext.data()[ebNE_kb],
                                   //VRANS
                                   porosity_ext,
                                   //
                                   ebqe_phi_s.data()[ebNE_kb],
                                   p_old,
                                   u_old,
                                   v_old,
                                   w_old,
                                   grad_p_old,
                                   grad_u_old,
                                   grad_v_old,
                                   grad_w_old,
                                   bc_p_ext,
                                   grad_p_ext,
                                   grad_u_ext,
                                   grad_v_ext,
                                   grad_w_ext,
                                   bc_u_ext,
                                   bc_v_ext,
                                   bc_w_ext,
                                   LAG_LES,
                                   bc_eddy_viscosity_ext,
                                   ebqe_eddy_viscosity_last.data()[ebNE_kb],
                                   bc_mom_u_acc_ext,
                                   bc_dmom_u_acc_u_ext,
                                   bc_mom_v_acc_ext,
                                   bc_dmom_v_acc_v_ext,
                                   bc_mom_w_acc_ext,
                                   bc_dmom_w_acc_w_ext,
                                   bc_mass_adv_ext,
                                   bc_dmass_adv_u_ext,
                                   bc_dmass_adv_v_ext,
                                   bc_dmass_adv_w_ext,
                                   bc_mom_u_adv_ext,
                                   bc_dmom_u_adv_u_ext,
                                   bc_dmom_u_adv_v_ext,
                                   bc_dmom_u_adv_w_ext,
                                   bc_mom_v_adv_ext,
                                   bc_dmom_v_adv_u_ext,
                                   bc_dmom_v_adv_v_ext,
                                   bc_dmom_v_adv_w_ext,
                                   bc_mom_w_adv_ext,
                                   bc_dmom_w_adv_u_ext,
                                   bc_dmom_w_adv_v_ext,
                                   bc_dmom_w_adv_w_ext,
                                   bc_mom_uu_diff_ten_ext,
                                   bc_mom_vv_diff_ten_ext,
                                   bc_mom_ww_diff_ten_ext,
                                   bc_mom_uv_diff_ten_ext,
                                   bc_mom_uw_diff_ten_ext,
                                   bc_mom_vu_diff_ten_ext,
                                   bc_mom_vw_diff_ten_ext,
                                   bc_mom_wu_diff_ten_ext,
                                   bc_mom_wv_diff_ten_ext,
                                   bc_mom_u_source_ext,
                                   bc_mom_v_source_ext,
                                   bc_mom_w_source_ext,
                                   bc_mom_u_ham_ext,
                                   bc_dmom_u_ham_grad_p_ext,
                                   bc_dmom_u_ham_grad_u_ext,
                                   bc_dmom_u_ham_u_ext,
                                   bc_dmom_u_ham_v_ext,
                                   bc_dmom_u_ham_w_ext,
                                   bc_mom_v_ham_ext,
                                   bc_dmom_v_ham_grad_p_ext,
                                   bc_dmom_v_ham_grad_v_ext,
                                   bc_dmom_v_ham_u_ext,
                                   bc_dmom_v_ham_v_ext,
                                   bc_dmom_v_ham_w_ext,
                                   bc_mom_w_ham_ext,
                                   bc_dmom_w_ham_grad_p_ext,
                                   bc_dmom_w_ham_grad_w_ext,
                                   bc_dmom_w_ham_u_ext,
                                   bc_dmom_w_ham_v_ext,
                                   bc_dmom_w_ham_w_ext,
                                   0.0,
                                   0.0,
                                   0.0);
              //Turbulence closure model
              if (turbulenceClosureModel >= 3)
                {
                  const double turb_var_grad_0_dummy[nSpace] = ZEROVEC;
                  const double c_mu = 0.09;//mwf hack
                  updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                          turbulenceClosureModel,
                                          eps_rho,
                                          eps_mu,
                                          rho_0,
                                          nu_0,
                                          rho_1,
                                          nu_1,
                                          useVF,
                                          ebqe_vf_ext.data()[ebNE_kb],
                                          ebqe_phi_ext.data()[ebNE_kb],
                                          porosity_ext,
                                          c_mu, //mwf hack
                                          ebqe_turb_var_0.data()[ebNE_kb],
                                          ebqe_turb_var_1.data()[ebNE_kb],
                                          turb_var_grad_0_dummy, //not needed
                                          eddy_viscosity_ext,
                                          mom_uu_diff_ten_ext,
                                          mom_vv_diff_ten_ext,
                                          mom_ww_diff_ten_ext,
                                          mom_uv_diff_ten_ext,
                                          mom_uw_diff_ten_ext,
                                          mom_vu_diff_ten_ext,
                                          mom_vw_diff_ten_ext,
                                          mom_wu_diff_ten_ext,
                                          mom_wv_diff_ten_ext,
                                          mom_u_source_ext,
                                          mom_v_source_ext,
                                          mom_w_source_ext);

                  updateTurbulenceClosure(NONCONSERVATIVE_FORM,
                                          turbulenceClosureModel,
                                          eps_rho,
                                          eps_mu,
                                          rho_0,
                                          nu_0,
                                          rho_1,
                                          nu_1,
                                          useVF,
                                          ebqe_vf_ext.data()[ebNE_kb],
                                          ebqe_phi_ext.data()[ebNE_kb],
                                          porosity_ext,
                                          c_mu, //mwf hack
                                          ebqe_turb_var_0.data()[ebNE_kb],
                                          ebqe_turb_var_1.data()[ebNE_kb],
                                          turb_var_grad_0_dummy, //not needed
                                          bc_eddy_viscosity_ext,
                                          bc_mom_uu_diff_ten_ext,
                                          bc_mom_vv_diff_ten_ext,
                                          bc_mom_ww_diff_ten_ext,
                                          bc_mom_uv_diff_ten_ext,
                                          bc_mom_uw_diff_ten_ext,
                                          bc_mom_vu_diff_ten_ext,
                                          bc_mom_vw_diff_ten_ext,
                                          bc_mom_wu_diff_ten_ext,
                                          bc_mom_wv_diff_ten_ext,
                                          bc_mom_u_source_ext,
                                          bc_mom_v_source_ext,
                                          bc_mom_w_source_ext);
                }
              //
              //moving domain
              //
              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  mom_u_ham_ext -= MOVING_DOMAIN*dmom_u_acc_u_ext*(grad_u_ext[0]*xt_ext +
                                                                   grad_u_ext[1]*yt_ext);
                  dmom_u_ham_grad_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
                  dmom_u_ham_grad_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
                }
              else
                {
                  mom_u_adv_ext[0] -= MOVING_DOMAIN*mom_u_acc_ext*xt_ext;
                  mom_u_adv_ext[1] -= MOVING_DOMAIN*mom_u_acc_ext*yt_ext;
                  dmom_u_adv_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
                  dmom_u_adv_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
                }

              if (NONCONSERVATIVE_FORM > 0.0)
                {
                  mom_v_ham_ext -= MOVING_DOMAIN*dmom_v_acc_v_ext*(grad_v_ext[0]*xt_ext +
                                                                   grad_v_ext[1]*yt_ext);
                  dmom_v_ham_grad_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
                  dmom_v_ham_grad_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
                }
              else
                {
                  mom_v_adv_ext[0] -= MOVING_DOMAIN*mom_v_acc_ext*xt_ext;
                  mom_v_adv_ext[1] -= MOVING_DOMAIN*mom_v_acc_ext*yt_ext;
                  dmom_v_adv_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
                  dmom_v_adv_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
                }

              //moving domain bc's
              if (NONCONSERVATIVE_FORM < 1.0)
                {
                  bc_mom_u_adv_ext[0] -= MOVING_DOMAIN*bc_mom_u_acc_ext*xt_ext;
                  bc_mom_u_adv_ext[1] -= MOVING_DOMAIN*bc_mom_u_acc_ext*yt_ext;

                  bc_mom_v_adv_ext[0] -= MOVING_DOMAIN*bc_mom_v_acc_ext*xt_ext;
                  bc_mom_v_adv_ext[1] -= MOVING_DOMAIN*bc_mom_v_acc_ext*yt_ext;
                }
              //
              //calculate the numerical fluxes
              //
              exteriorNumericalAdvectiveFluxDerivatives(NONCONSERVATIVE_FORM,
                                                        isDOFBoundary_p.data()[ebNE_kb],
                                                        isDOFBoundary_u.data()[ebNE_kb],
                                                        isDOFBoundary_v.data()[ebNE_kb],
                                                        isDOFBoundary_w.data()[ebNE_kb],
                                                        isAdvectiveFluxBoundary_p.data()[ebNE_kb],
                                                        isAdvectiveFluxBoundary_u.data()[ebNE_kb],
                                                        isAdvectiveFluxBoundary_v.data()[ebNE_kb],
                                                        isAdvectiveFluxBoundary_w.data()[ebNE_kb],
                                                        dmom_u_ham_grad_p_ext[0],//=1/rho
                                                        normal,
                                                        bc_p_ext,
                                                        bc_u_ext,
                                                        bc_v_ext,
                                                        bc_mass_adv_ext,
                                                        bc_mom_u_adv_ext,
                                                        bc_mom_v_adv_ext,
                                                        bc_mom_w_adv_ext,
                                                        ebqe_bc_flux_mass_ext.data()[ebNE_kb]+MOVING_DOMAIN*(xt_ext*normal[0]+yt_ext*normal[1]),//bc is relative mass  flux
                                                        ebqe_bc_flux_mom_u_adv_ext.data()[ebNE_kb],
                                                        ebqe_bc_flux_mom_v_adv_ext.data()[ebNE_kb],
                                                        ebqe_bc_flux_mom_w_adv_ext.data()[ebNE_kb],
                                                        p_ext,
                                                        u_ext,
                                                        v_ext,
                                                        dmom_u_acc_u_ext,
                                                        mass_adv_ext,
                                                        mom_u_adv_ext,
                                                        mom_v_adv_ext,
                                                        mom_w_adv_ext,
                                                        dmass_adv_u_ext,
                                                        dmass_adv_v_ext,
                                                        dmass_adv_w_ext,
                                                        dmom_u_adv_p_ext,
                                                        dmom_u_ham_grad_u_ext,
                                                        dmom_u_adv_u_ext,
                                                        dmom_u_adv_v_ext,
                                                        dmom_u_adv_w_ext,
                                                        dmom_v_adv_p_ext,
                                                        dmom_v_adv_u_ext,
                                                        dmom_v_adv_v_ext,
                                                        dmom_v_adv_w_ext,
                                                        dmom_w_adv_p_ext,
                                                        dmom_w_adv_u_ext,
                                                        dmom_w_adv_v_ext,
                                                        dmom_w_adv_w_ext,
                                                        dflux_mass_u_ext,
                                                        dflux_mass_v_ext,
                                                        dflux_mass_w_ext,
                                                        dflux_mom_u_adv_p_ext,
                                                        dflux_mom_u_adv_u_ext,
                                                        dflux_mom_u_adv_v_ext,
                                                        dflux_mom_u_adv_w_ext,
                                                        dflux_mom_v_adv_p_ext,
                                                        dflux_mom_v_adv_u_ext,
                                                        dflux_mom_v_adv_v_ext,
                                                        dflux_mom_v_adv_w_ext,
                                                        dflux_mom_w_adv_p_ext,
                                                        dflux_mom_w_adv_u_ext,
                                                        dflux_mom_w_adv_v_ext,
                                                        dflux_mom_w_adv_w_ext);
              //
              //calculate the flux jacobian
              //
              ck.calculateGScale(G,normal,h_penalty);
              penalty = useMetrics*C_b/h_penalty + (1.0-useMetrics)*ebqe_penalty_ext.data()[ebNE_kb];
              if (elementIsActive[eN])
                //                if(true)//boundaryFlags[ebN] > 0)
                { //if boundary flag positive, then include flux contributions on interpart boundaries
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
                      fluxJacobian_p_p[j]=0.0;
                      fluxJacobian_u_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_p_ext,p_trial_trace_ref.data()[ebN_local_kb_j]);
                      fluxJacobian_v_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_p_ext,p_trial_trace_ref.data()[ebN_local_kb_j]);
                    }
                  for (int j=0;j<nDOF_v_trial_element;j++)
                    {
                      register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_v_trial_element+j;
                      fluxJacobian_p_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_u_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]);
                      fluxJacobian_p_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_v_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]);
                      fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_u_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]) +
                        ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
                                                               ebqe_phi_ext.data()[ebNE_kb],
                                                               sdInfo_u_u_rowptr.data(),
                                                               sdInfo_u_u_colind.data(),
                                                               isDOFBoundary_u.data()[ebNE_kb],
                                                               isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                               normal,
                                                               mom_uu_diff_ten_ext,
                                                               vel_trial_trace_ref.data()[ebN_local_kb_j],
                                                               &vel_grad_trial_trace[j_nSpace],
                                                               penalty);//ebqe_penalty_ext.data()[ebNE_kb]);
                      fluxJacobian_u_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_v_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]) +
                        ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
                                                               ebqe_phi_ext.data()[ebNE_kb],
                                                               sdInfo_u_v_rowptr.data(),
                                                               sdInfo_u_v_colind.data(),
                                                               isDOFBoundary_v.data()[ebNE_kb],
                                                               isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                                               normal,
                                                               mom_uv_diff_ten_ext,
                                                               vel_trial_trace_ref.data()[ebN_local_kb_j],
                                                               &vel_grad_trial_trace[j_nSpace],
                                                               penalty);//ebqe_penalty_ext.data()[ebNE_kb]);

                      fluxJacobian_v_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_u_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]) +
                        ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
                                                               ebqe_phi_ext.data()[ebNE_kb],
                                                               sdInfo_v_u_rowptr.data(),
                                                               sdInfo_v_u_colind.data(),
                                                               isDOFBoundary_u.data()[ebNE_kb],
                                                               isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                               normal,
                                                               mom_vu_diff_ten_ext,
                                                               vel_trial_trace_ref.data()[ebN_local_kb_j],
                                                               &vel_grad_trial_trace[j_nSpace],
                                                               penalty);//ebqe_penalty_ext.data()[ebNE_kb]);
                      fluxJacobian_v_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_v_ext,vel_trial_trace_ref.data()[ebN_local_kb_j]) +
                        ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
                                                               ebqe_phi_ext.data()[ebNE_kb],
                                                               sdInfo_v_v_rowptr.data(),
                                                               sdInfo_v_v_colind.data(),
                                                               isDOFBoundary_v.data()[ebNE_kb],
                                                               isDiffusiveFluxBoundary_v.data()[ebNE_kb],
                                                               normal,
                                                               mom_vv_diff_ten_ext,
                                                               vel_trial_trace_ref.data()[ebN_local_kb_j],
                                                               &vel_grad_trial_trace[j_nSpace],
                                                               penalty);//ebqe_penalty_ext.data()[ebNE_kb]);
                    }//j
                }//if boundaryFlags.data()[ebN] positive
              //
              //update the global Jacobian from the flux Jacobian
              //
              const double H_s = gf_s.H(particle_eps, ebqe_phi_s[ebNE_kb]);
              if (elementIsActive[eN])
                {
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_test_element+i;
                      for (int j=0;j<nDOF_trial_element;j++)
                        {
                          register int eN_j = eN*nDOF_trial_element+j;
                          register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;

                          globalJacobian.data()[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_eb_p_p.data()[ebN_i_j]] += H_s*fluxJacobian_p_p[j]*p_test_dS[i];
                        }
                    }
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_test_element+i;
                      for (int j=0;j<nDOF_v_trial_element;j++)
                        {
                          register int eN_j = eN*nDOF_v_trial_element+j;
                          register int ebN_i_j = ebN*4*nDOF_test_X_v_trial_element + i*nDOF_v_trial_element + j,ebN_local_kb_j=ebN_local_kb*nDOF_v_trial_element+j;
                          globalJacobian.data()[csrRowIndeces_p_u.data()[eN_i] + csrColumnOffsets_eb_p_u.data()[ebN_i_j]] += H_s*fluxJacobian_p_u[j]*p_test_dS[i];
                          globalJacobian.data()[csrRowIndeces_p_v.data()[eN_i] + csrColumnOffsets_eb_p_v.data()[ebN_i_j]] += H_s*fluxJacobian_p_v[j]*p_test_dS[i];
                        }
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_v_test_element+i;
                      for (int j=0;j<nDOF_trial_element;j++)
                        {
                          register int ebN_i_j = ebN*4*nDOF_v_test_X_trial_element + i*nDOF_trial_element + j,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
                          globalJacobian.data()[csrRowIndeces_u_p.data()[eN_i] + csrColumnOffsets_eb_u_p.data()[ebN_i_j]] += H_s*fluxJacobian_u_p[j]*vel_test_dS[i];
                          globalJacobian.data()[csrRowIndeces_v_p.data()[eN_i] + csrColumnOffsets_eb_v_p.data()[ebN_i_j]] += H_s*fluxJacobian_v_p[j]*vel_test_dS[i];
                        }
                    }
                  for (int i=0;i<nDOF_v_test_element;i++)
                    {
                      register int eN_i = eN*nDOF_v_test_element+i;
                      for (int j=0;j<nDOF_v_trial_element;j++)
                        {
                          register int eN_j = eN*nDOF_v_trial_element+j;
                          register int ebN_i_j = ebN*4*nDOF_v_test_X_v_trial_element + i*nDOF_v_trial_element + j,ebN_local_kb_j=ebN_local_kb*nDOF_v_trial_element+j;
                          globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_eb_u_u.data()[ebN_i_j]] +=
			    H_s*(fluxJacobian_u_u[j]*vel_test_dS[i]+
				 ck.ExteriorElementBoundaryDiffusionAdjointJacobian(isDOFBoundary_u.data()[ebNE_kb],
										    isDiffusiveFluxBoundary_u.data()[ebNE_kb],
										    eb_adjoint_sigma,
										    vel_trial_trace_ref.data()[ebN_local_kb_j],
										    normal,
										    sdInfo_u_u_rowptr.data(),
										    sdInfo_u_u_colind.data(),
										    mom_uu_diff_ten_ext,
										    &vel_grad_test_dS[i*nSpace]));
                          globalJacobian.data()[csrRowIndeces_u_v.data()[eN_i] + csrColumnOffsets_eb_u_v.data()[ebN_i_j]] +=
			    H_s*(fluxJacobian_u_v[j]*vel_test_dS[i]+
				 ck.ExteriorElementBoundaryDiffusionAdjointJacobian(isDOFBoundary_v.data()[ebNE_kb],
										    isDiffusiveFluxBoundary_u.data()[ebNE_kb],
										    eb_adjoint_sigma,
										    vel_trial_trace_ref.data()[ebN_local_kb_j],
										    normal,
										    sdInfo_u_v_rowptr.data(),
										    sdInfo_u_v_colind.data(),
										    mom_uv_diff_ten_ext,
										    &vel_grad_test_dS[i*nSpace]));
			  globalJacobian.data()[csrRowIndeces_v_u.data()[eN_i] + csrColumnOffsets_eb_v_u.data()[ebN_i_j]] +=
			    H_s*(fluxJacobian_v_u[j]*vel_test_dS[i]+
				 ck.ExteriorElementBoundaryDiffusionAdjointJacobian(isDOFBoundary_u.data()[ebNE_kb],
										    isDiffusiveFluxBoundary_v.data()[ebNE_kb],
										    eb_adjoint_sigma,
										    vel_trial_trace_ref.data()[ebN_local_kb_j],
										    normal,
										    sdInfo_v_u_rowptr.data(),
										    sdInfo_v_u_colind.data(),
										    mom_vu_diff_ten_ext,
										    &vel_grad_test_dS[i*nSpace]));
                          globalJacobian.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_eb_v_v.data()[ebN_i_j]] +=
			    H_s*(fluxJacobian_v_v[j]*vel_test_dS[i]+
				 ck.ExteriorElementBoundaryDiffusionAdjointJacobian(isDOFBoundary_v.data()[ebNE_kb],
										    isDiffusiveFluxBoundary_v.data()[ebNE_kb],
										    eb_adjoint_sigma,
										    vel_trial_trace_ref.data()[ebN_local_kb_j],
										    normal,
										    sdInfo_v_v_rowptr.data(),
										    sdInfo_v_v_colind.data(),
										    mom_vv_diff_ten_ext,
										    &vel_grad_test_dS[i*nSpace]));
                        }//j
                    }//i
                }
            }//kb
        }//ebNE
    }//computeJacobian
    
    void calculateVelocityAverage(arguments_dict& args)
    {
      int nExteriorElementBoundaries_global = args.m_iscalar["nExteriorElementBoundaries_global"];
      xt::pyarray<int>& exteriorElementBoundariesArray = args.m_iarray["exteriorElementBoundariesArray"];
      int nInteriorElementBoundaries_global = args.m_iscalar["nInteriorElementBoundaries_global"];
      xt::pyarray<int>& interiorElementBoundariesArray = args.m_iarray["interiorElementBoundariesArray"];
      xt::pyarray<int>& elementBoundaryElementsArray = args.m_iarray["elementBoundaryElementsArray"];
      xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.m_iarray["elementBoundaryLocalElementBoundariesArray"];
      xt::pyarray<double>& mesh_dof = args.m_darray["mesh_dof"];
      xt::pyarray<double>& mesh_velocity_dof = args.m_darray["mesh_velocity_dof"];
      double MOVING_DOMAIN = args.m_dscalar["MOVING_DOMAIN"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& mesh_trial_trace_ref = args.m_darray["mesh_trial_trace_ref"];
      xt::pyarray<double>& mesh_grad_trial_trace_ref = args.m_darray["mesh_grad_trial_trace_ref"];
      xt::pyarray<double>& normal_ref = args.m_darray["normal_ref"];
      xt::pyarray<double>& boundaryJac_ref = args.m_darray["boundaryJac_ref"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      xt::pyarray<double>& vel_trial_trace_ref = args.m_darray["vel_trial_trace_ref"];
      xt::pyarray<double>& ebqe_velocity = args.m_darray["ebqe_velocity"];
      xt::pyarray<double>& velocityAverage = args.m_darray["velocityAverage"];
      xt::pyarray<int>& elementMaterialTypes = args.m_iarray["elementMaterialTypes"];
      xt::pyarray<double>& porosityTypes = args.m_darray["porosityTypes"];
      int permutations[nQuadraturePoints_elementBoundary];
      double xArray_left[nQuadraturePoints_elementBoundary*nSpace],
        xArray_right[nQuadraturePoints_elementBoundary*nSpace];
      for (int i=0;i<nQuadraturePoints_elementBoundary;i++)
        permutations[i]=i;//just to initialize
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
        {
          register int ebN = exteriorElementBoundariesArray.data()[ebNE];
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace,
                ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
              velocityAverage.data()[ebN_kb_nSpace+0]=ebqe_velocity.data()[ebNE_kb_nSpace+0];
              velocityAverage.data()[ebN_kb_nSpace+1]=ebqe_velocity.data()[ebNE_kb_nSpace+1];
            }//ebNE
        }
      for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++)
        {
          register int ebN = interiorElementBoundariesArray.data()[ebNI],
            left_eN_global   = elementBoundaryElementsArray.data()[ebN*2+0],
            left_ebN_element  = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
            right_eN_global  = elementBoundaryElementsArray.data()[ebN*2+1],
            right_ebN_element = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+1],
            left_eN_nDOF_trial_element = left_eN_global*nDOF_trial_element,
            right_eN_nDOF_trial_element = right_eN_global*nDOF_trial_element;
          double jac[nSpace*nSpace],
            jacDet,
            jacInv[nSpace*nSpace],
            boundaryJac[nSpace*(nSpace-1)],
            metricTensor[(nSpace-1)*(nSpace-1)],
            metricTensorDetSqrt,
            normal[nSpace],
            x,y,z,
            xt,yt,zt,integralScaling,
            left_porosity  = porosityTypes[elementMaterialTypes[left_eN_global]],
            right_porosity = porosityTypes[elementMaterialTypes[right_eN_global]];

          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              ck.calculateMapping_elementBoundary(left_eN_global,
                                                  left_ebN_element,
                                                  kb,
                                                  left_ebN_element*nQuadraturePoints_elementBoundary+kb,
                                                  mesh_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_trace_ref.data(),
                                                  mesh_grad_trial_trace_ref.data(),
                                                  boundaryJac_ref.data(),
                                                  jac,
                                                  jacDet,
                                                  jacInv,
                                                  boundaryJac,
                                                  metricTensor,
                                                  metricTensorDetSqrt,
                                                  normal_ref.data(),
                                                  normal,
                                                  x,y,z);
              xArray_left[kb*nSpace+0] = x;
              xArray_left[kb*nSpace+1] = y;
              ck.calculateMapping_elementBoundary(right_eN_global,
                                                  right_ebN_element,
                                                  kb,
                                                  right_ebN_element*nQuadraturePoints_elementBoundary+kb,
                                                  mesh_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_trace_ref.data(),
                                                  mesh_grad_trial_trace_ref.data(),
                                                  boundaryJac_ref.data(),
                                                  jac,
                                                  jacDet,
                                                  jacInv,
                                                  boundaryJac,
                                                  metricTensor,
                                                  metricTensorDetSqrt,
                                                  normal_ref.data(),
                                                  normal,
                                                  x,y,z);
              ck.calculateMappingVelocity_elementBoundary(left_eN_global,
                                                          left_ebN_element,
                                                          kb,
                                                          left_ebN_element*nQuadraturePoints_elementBoundary+kb,
                                                          mesh_velocity_dof.data(),
                                                          mesh_l2g.data(),
                                                          mesh_trial_trace_ref.data(),
                                                          xt,yt,zt,
                                                          normal,
                                                          boundaryJac,
                                                          metricTensor,
                                                          integralScaling);
              xArray_right[kb*nSpace+0] = x;
              xArray_right[kb*nSpace+1] = y;
            }
          for  (int kb_left=0;kb_left<nQuadraturePoints_elementBoundary;kb_left++)
            {
              double errorNormMin = 1.0;
              for  (int kb_right=0;kb_right<nQuadraturePoints_elementBoundary;kb_right++)
                {
                  double errorNorm=0.0;
                  for (int I=0;I<nSpace;I++)
                    {
                      errorNorm += fabs(xArray_left[kb_left*nSpace+I]
                                        -
                                        xArray_right[kb_right*nSpace+I]);
                    }
                  if (errorNorm < errorNormMin)
                    {
                      permutations[kb_right] = kb_left;
                      errorNormMin = errorNorm;
                    }
                }
            }
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
              register double u_left=0.0,
                v_left=0.0,
                w_left=0.0,
                u_right=0.0,
                v_right=0.0,
                w_right=0.0;
              register int left_kb = kb,
                right_kb = permutations[kb],
                left_ebN_element_kb_nDOF_test_element=(left_ebN_element*nQuadraturePoints_elementBoundary+left_kb)*nDOF_test_element,
                right_ebN_element_kb_nDOF_test_element=(right_ebN_element*nQuadraturePoints_elementBoundary+right_kb)*nDOF_test_element;
              //
              //calculate the velocity solution at quadrature points on left and right
              //
              ck.valFromDOF(u_dof.data(),&vel_l2g.data()[left_eN_nDOF_trial_element],&vel_trial_trace_ref.data()[left_ebN_element_kb_nDOF_test_element],u_left);
              ck.valFromDOF(v_dof.data(),&vel_l2g.data()[left_eN_nDOF_trial_element],&vel_trial_trace_ref.data()[left_ebN_element_kb_nDOF_test_element],v_left);
              //
              ck.valFromDOF(u_dof.data(),&vel_l2g.data()[right_eN_nDOF_trial_element],&vel_trial_trace_ref.data()[right_ebN_element_kb_nDOF_test_element],u_right);
              ck.valFromDOF(v_dof.data(),&vel_l2g.data()[right_eN_nDOF_trial_element],&vel_trial_trace_ref.data()[right_ebN_element_kb_nDOF_test_element],v_right);
              //
              velocityAverage.data()[ebN_kb_nSpace+0]=0.5*(left_porosity*u_left + right_porosity*u_right);
              velocityAverage.data()[ebN_kb_nSpace+1]=0.5*(left_porosity*v_left + right_porosity*v_right);
            }//ebNI
        }
    }

    inline
    void evaluateTPAdvectionCoefficients(const double eps_rho,
                                         const double rho_0,
                                         const double rho_1,
                                         const double useVF,
                                         const double& vf,
                                         const double& phi,
                                         const double& u,
                                         const double& v,
                                         double dmass_adv_p[nSpace],
                                         double dmom_u_adv_u[nSpace],
                                         double dmom_v_adv_v[nSpace])
    {
      double H_rho, ImH_rho, rho;

      H_rho = (1.0-useVF)*gf.H(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf));
      ImH_rho = (1.0-useVF)*gf.ImH(eps_rho,phi) + useVF*(1.0-fmin(1.0,fmax(0.0,vf)));

      rho = rho_0*ImH_rho + rho_1*H_rho;

      dmass_adv_p[0] = rho*u;
      dmass_adv_p[1] = rho*v;

      dmom_u_adv_u[0] = rho*u;
      dmom_u_adv_u[1] = rho*v;

      dmom_v_adv_v[0] = rho*u;
      dmom_v_adv_v[1] = rho*v;
    }
    inline
    void evaluateTPInvViscosityMassCoefficients(const int use_numerical_viscosity,
                                                const double numerical_viscosity,
                                                const double eps_rho,
                                                const double eps_mu,
                                                const double rho_0,
                                                double nu_0,
                                                const double rho_1,
                                                double nu_1,
                                                const double useVF,
                                                const double& vf,
                                                const double& phi,
                                                const double& p,
                                                const double& u,
                                                const double& v,
                                                double& mom_p_acc,
                                                double& dmom_p_acc_p,
                                                double& mom_u_acc,
                                                double& dmom_u_acc_u,
                                                double& mom_v_acc,
                                                double& dmom_v_acc_v)
    {
      // This should be split off into a seperate function
      double H_rho, ImH_rho, H_mu, ImH_mu, rho, nu, mu;

      H_rho = (1.0-useVF)*gf.H(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf));
      ImH_rho = (1.0-useVF)*gf.ImH(eps_rho,phi) + useVF*(1.0-fmin(1.0,fmax(0.0,vf)));
      H_mu = (1.0-useVF)*gf.H(eps_mu,phi) + useVF*fmin(1.0,fmax(0.0,vf));
      ImH_mu = (1.0-useVF)*gf.ImH(eps_mu,phi) + useVF*(1.0-fmin(1.0,fmax(0.0,vf)));

      rho = rho_0*ImH_rho + rho_1*H_rho;
      nu = nu_0*ImH_mu + nu_1*H_mu;

      mu = rho_0*nu_0*ImH_mu + rho_1*nu_1*H_mu + use_numerical_viscosity*numerical_viscosity;
      //mu = rho*nu;

      mom_p_acc = p / mu;
      dmom_p_acc_p = 1. / mu;

      mom_u_acc = u / mu;
      dmom_u_acc_u = 1. / mu;

      mom_v_acc = v / mu;
      dmom_v_acc_v = 1. / mu;
    }
    inline
    void evaluateTPDensityMassCoefficients(const double eps_rho,
                                           const double rho_0,
                                           const double rho_1,
                                           const double useVF,
                                           const double& vf,
                                           const double& phi,
                                           const double& p,
                                           const double& u,
                                           const double& v,
                                           double& mom_p_acc,
                                           double& dmom_p_acc_p,
                                           double& mom_u_acc,
                                           double& dmom_u_acc_u,
                                           double& mom_v_acc,
                                           double& dmom_v_acc_v)
    {
      double H_rho, ImH_rho, rho;

      H_rho = (1.0-useVF)*gf.H(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf));
      ImH_rho = (1.0-useVF)*gf.ImH(eps_rho,phi) + useVF*(1.0-fmin(1.0,fmax(0.0,vf)));

      rho = rho_0*ImH_rho + rho_1*H_rho;

      mom_p_acc = p * rho;
      dmom_p_acc_p = rho;

      mom_u_acc = u * rho;
      dmom_u_acc_u = rho;

      mom_v_acc = v * rho;
      dmom_v_acc_v = rho;
    }
    inline
    void evaluateTPInvDensityLaplaceCoefficients(const double eps_rho,
                                                 const double rho_0,
                                                 const double rho_1,
                                                 const double useVF,
                                                 const double& vf,
                                                 const double& phi,
                                                 double mom_p_diff_ten[nSpace],
                                                 double mom_u_diff_ten[nSpace],
                                                 double mom_v_diff_ten[nSpace])
    {
      double H_rho, ImH_rho, rho;

      H_rho = (1.0-useVF)*gf.H(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf));
      ImH_rho = (1.0-useVF)*gf.ImH(eps_rho,phi) + useVF*(1.0-fmin(1.0,fmax(0.0,vf)));

      rho = rho_0*ImH_rho + rho_1*H_rho;

      mom_p_diff_ten[0] = 1.0 / rho ;
      mom_p_diff_ten[1] = 1.0 / rho ;

      mom_u_diff_ten[0] = 1.0 / rho ;
      mom_u_diff_ten[1] = 1.0 / rho ;

      mom_v_diff_ten[0] = 1.0 / rho ;
      mom_v_diff_ten[1] = 1.0 / rho ;

    }

    void getTwoPhaseAdvectionOperator(arguments_dict& args)
    {
      xt::pyarray<double>& mesh_trial_ref = args.m_darray["mesh_trial_ref"];
      xt::pyarray<double>& mesh_grad_trial_ref = args.m_darray["mesh_grad_trial_ref"];
      xt::pyarray<double>& mesh_dof = args.m_darray["mesh_dof"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& dV_ref = args.m_darray["dV_ref"];
      xt::pyarray<double>& p_trial_ref = args.m_darray["p_trial_ref"];
      xt::pyarray<double>& p_grad_trial_ref = args.m_darray["p_grad_trial_ref"];
      xt::pyarray<double>& vel_trial_ref = args.m_darray["vel_trial_ref"];
      xt::pyarray<double>& vel_grad_trial_ref = args.m_darray["vel_grad_trial_ref"];
      xt::pyarray<double>& elementDiameter = args.m_darray["elementDiameter"];
      xt::pyarray<double>& nodeDiametersArray = args.m_darray["nodeDiametersArray"];
      int nElements_global = args.m_iscalar["nElements_global"];
      double useMetrics = args.m_dscalar["useMetrics"];
      double epsFact_rho = args.m_dscalar["epsFact_rho"];
      double epsFact_mu = args.m_dscalar["epsFact_mu"];
      double rho_0 = args.m_dscalar["rho_0"];
      double nu_0 = args.m_dscalar["nu_0"];
      double rho_1 = args.m_dscalar["rho_1"];
      double nu_1 = args.m_dscalar["nu_1"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      const double useVF = args.m_dscalar["useVF"];
      xt::pyarray<double> &vf = args.m_darray["&vf"];
      xt::pyarray<double> &phi = args.m_darray["&phi"];
      xt::pyarray<int>& csrRowIndeces_p_p = args.m_iarray["csrRowIndeces_p_p"];
      xt::pyarray<int>& csrColumnOffsets_p_p = args.m_iarray["csrColumnOffsets_p_p"];
      xt::pyarray<int>& csrRowIndeces_u_u = args.m_iarray["csrRowIndeces_u_u"];
      xt::pyarray<int>& csrColumnOffsets_u_u = args.m_iarray["csrColumnOffsets_u_u"];
      xt::pyarray<int>& csrRowIndeces_v_v = args.m_iarray["csrRowIndeces_v_v"];
      xt::pyarray<int>& csrColumnOffsets_v_v = args.m_iarray["csrColumnOffsets_v_v"];
      xt::pyarray<int>& csrRowIndeces_w_w = args.m_iarray["csrRowIndeces_w_w"];
      xt::pyarray<int>& csrColumnOffsets_w_w = args.m_iarray["csrColumnOffsets_w_w"];
      xt::pyarray<double>& advection_matrix = args.m_darray["advection_matrix"];
      gf.useExact = false;
      for (int eN=0 ; eN < nElements_global ; ++eN)
        {
          // local matrix allocations
          double eps_rho;

          double local_matrix_p_p[nDOF_test_element][nDOF_trial_element];
          double local_matrix_u_u[nDOF_test_element][nDOF_trial_element];
          double local_matrix_v_v[nDOF_test_element][nDOF_trial_element];

          // clear local matrix entries
          for (int i=0 ; i < nDOF_test_element ; ++i)
            for (int j=0 ; j < nDOF_trial_element ; ++j){
              local_matrix_p_p[i][j] = 0. ;
              local_matrix_u_u[i][j] = 0. ;
              local_matrix_v_v[i][j] = 0. ;
            }

          for (int k=0 ; k < nQuadraturePoints_element ; ++k){

            int eN_k = eN*nQuadraturePoints_element + k;
            int eN_nDOF_trial_element = eN*nDOF_trial_element;

            double jac[nSpace*nSpace];
            double jacInv[nSpace*nSpace];
            double u=0.0, v=0.0;
            double dmass_adv_p[nSpace], dmom_u_adv_u[nSpace], dmom_v_adv_v[nSpace];
            double p_grad_trial[nDOF_trial_element*nSpace],
              vel_grad_trial[nDOF_trial_element*nSpace];
            double p_test_dV[nDOF_test_element], vel_test_dV[nDOF_test_element];
            double p_grad_test_dV[nDOF_test_element*nSpace],
              vel_grad_test_dV[nDOF_test_element*nSpace];
            double jacDet, x, y, z, dV, h_phi;

            ck.calculateMapping_element(eN,
                                        k,
                                        mesh_dof.data(),
                                        mesh_l2g.data(),
                                        mesh_trial_ref.data(),
                                        mesh_grad_trial_ref.data(),
                                        jac,
                                        jacDet,
                                        jacInv,
                                        x,y,z);

            ck.calculateH_element(eN,
                                  k,
                                  nodeDiametersArray.data(),
                                  mesh_l2g.data(),
                                  mesh_trial_ref.data(),
                                  h_phi);

            dV = fabs(jacDet)*dV_ref.data()[k];

            eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);

            ck.gradTrialFromRef(&p_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
            ck.gradTrialFromRef(&vel_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);

            ck.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],&vel_trial_ref.data()[k*nDOF_trial_element],u);
            ck.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],&vel_trial_ref.data()[k*nDOF_trial_element],v);

            for (int j=0; j<nDOF_trial_element;++j)
              {
                p_test_dV[j] = p_trial_ref.data()[k*nDOF_trial_element+j]*dV;
                vel_test_dV[j] = vel_trial_ref.data()[k*nDOF_trial_element+j]*dV;
                for (int i=0; i<nSpace; ++i)
                  {
                    p_grad_test_dV[j*nSpace+i] = p_grad_trial[j*nSpace+i]*dV;
                    vel_grad_test_dV[j*nSpace+i] = vel_grad_trial[j*nSpace+i]*dV;
                  }
              }


            evaluateTPAdvectionCoefficients(eps_rho,
                                            rho_0,
                                            rho_1,
                                            useVF,
                                            vf.data()[eN_k],
                                            phi.data()[eN_k],
                                            u,
                                            v,
                                            dmass_adv_p,
                                            dmom_u_adv_u,
                                            dmom_v_adv_v);


            for(int i=0; i<nDOF_test_element;++i){
              int i_nSpace = i*nSpace;

              for(int j=0; j<nDOF_trial_element;++j){

                int j_nSpace = j*nSpace;

                local_matrix_p_p[i][j] -= ck.HamiltonianJacobian_weak(dmass_adv_p,&p_grad_test_dV[i_nSpace],p_trial_ref.data()[j]);
                //local_matrix_p_p[i][j] += ck.HamiltonianJacobian_weak(dmass_adv_p ,&p_grad_trial[j_nSpace]  ,p_test_dV[i]);
                local_matrix_u_u[i][j] += ck.HamiltonianJacobian_weak(dmom_u_adv_u,&vel_grad_trial[j_nSpace],vel_test_dV[i]);
                local_matrix_v_v[i][j] += ck.HamiltonianJacobian_weak(dmom_v_adv_v,&vel_grad_trial[j_nSpace],vel_test_dV[i]);
              }
            }


          }//k

          // Write local matrix information into global system
          for (int i=0 ; i < nDOF_test_element ; ++i)
            {
              int eN_i = eN*nDOF_test_element + i;
              for (int j=0 ; j < nDOF_trial_element ; ++j)
                {
                  int eN_i_j = eN_i*nDOF_trial_element + j;
                  advection_matrix.data()[csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_p_p.data()[eN_i_j]] += local_matrix_p_p[i][j] ;
                  advection_matrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += local_matrix_u_u[i][j] ;
                  advection_matrix.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_v_v.data()[eN_i_j]] += local_matrix_v_v[i][j] ;
                }
            }

        }//eN
    } // getTwoPhaseAdvectionOperator

    void getTwoPhaseInvScaledLaplaceOperator(arguments_dict& args)
    {
      xt::pyarray<double>& mesh_trial_ref = args.m_darray["mesh_trial_ref"];
      xt::pyarray<double>& mesh_grad_trial_ref = args.m_darray["mesh_grad_trial_ref"];
      xt::pyarray<double>& mesh_dof = args.m_darray["mesh_dof"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& dV_ref = args.m_darray["dV_ref"];
      xt::pyarray<double>& p_grad_trial_ref = args.m_darray["p_grad_trial_ref"];
      xt::pyarray<double>& vel_grad_trial_ref = args.m_darray["vel_grad_trial_ref"];
      xt::pyarray<double>& elementDiameter = args.m_darray["elementDiameter"];
      xt::pyarray<double>& nodeDiametersArray = args.m_darray["nodeDiametersArray"];
      int nElements_global = args.m_iscalar["nElements_global"];
      double useMetrics = args.m_dscalar["useMetrics"];
      double epsFact_rho = args.m_dscalar["epsFact_rho"];
      double epsFact_mu = args.m_dscalar["epsFact_mu"];
      double rho_0 = args.m_dscalar["rho_0"];
      double nu_0 = args.m_dscalar["nu_0"];
      double rho_1 = args.m_dscalar["rho_1"];
      double nu_1 = args.m_dscalar["nu_1"];
      xt::pyarray<int>& p_l2g = args.m_iarray["p_l2g"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<double>& p_dof = args.m_darray["p_dof"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      const double useVF = args.m_dscalar["useVF"];
      xt::pyarray<double>& vf = args.m_darray["vf"];
      xt::pyarray<double>& phi = args.m_darray["phi"];
      xt::pyarray<int>& sdInfo_p_p_rowptr = args.m_iarray["sdInfo_p_p_rowptr"];
      xt::pyarray<int>& sdInfo_p_p_colind = args.m_iarray["sdInfo_p_p_colind"];
      xt::pyarray<int>& sdInfo_u_u_rowptr = args.m_iarray["sdInfo_u_u_rowptr"];
      xt::pyarray<int>& sdInfo_u_u_colind = args.m_iarray["sdInfo_u_u_colind"];
      xt::pyarray<int>& sdInfo_v_v_rowptr = args.m_iarray["sdInfo_v_v_rowptr"];
      xt::pyarray<int>& sdInfo_v_v_colind = args.m_iarray["sdInfo_v_v_colind"];
      xt::pyarray<int>& sdInfo_w_w_rowptr = args.m_iarray["sdInfo_w_w_rowptr"];
      xt::pyarray<int>& sdInfo_w_w_colind = args.m_iarray["sdInfo_w_w_colind"];
      xt::pyarray<int>& csrRowIndeces_p_p = args.m_iarray["csrRowIndeces_p_p"];
      xt::pyarray<int>& csrColumnOffsets_p_p = args.m_iarray["csrColumnOffsets_p_p"];
      xt::pyarray<int>& csrRowIndeces_u_u = args.m_iarray["csrRowIndeces_u_u"];
      xt::pyarray<int>& csrColumnOffsets_u_u = args.m_iarray["csrColumnOffsets_u_u"];
      xt::pyarray<int>& csrRowIndeces_v_v = args.m_iarray["csrRowIndeces_v_v"];
      xt::pyarray<int>& csrColumnOffsets_v_v = args.m_iarray["csrColumnOffsets_v_v"];
      xt::pyarray<int>& csrRowIndeces_w_w = args.m_iarray["csrRowIndeces_w_w"];
      xt::pyarray<int>& csrColumnOffsets_w_w = args.m_iarray["csrColumnOffsets_w_w"];
      xt::pyarray<double>& laplace_matrix = args.m_darray["laplace_matrix"];
      gf.useExact = false;
      for (int eN=0 ; eN < nElements_global ; ++eN)
        {
          // local matrix allocations
          double eps_rho, eps_mu;

          double local_matrix_p_p[nDOF_test_element][nDOF_trial_element];
          double local_matrix_u_u[nDOF_test_element][nDOF_trial_element];
          double local_matrix_v_v[nDOF_test_element][nDOF_trial_element];

          // reset local matrix entries
          for (int i=0 ; i < nDOF_test_element ; ++i)
            for (int j=0 ; j < nDOF_trial_element ; ++j){
              // set local matrices to 0
              local_matrix_p_p[i][j] = 0.;
              local_matrix_u_u[i][j] = 0.;
              local_matrix_v_v[i][j] = 0.;
            }

          // Loop over quadrature points on element
          for (int k=0 ; k < nQuadraturePoints_element; ++k){

            int eN_k = eN*nQuadraturePoints_element + k;
            int eN_nDOF_trial_element = eN*nDOF_trial_element;

            double grad_p[nSpace], grad_u[nSpace], grad_v[nSpace];
            double jac[nSpace*nSpace];
            double jacInv[nSpace*nSpace];
            double mom_pp_diff_ten[nSpace];
            double mom_uu_diff_ten[nSpace];
            double mom_vv_diff_ten[nSpace];
            double p_grad_trial[nDOF_trial_element*nSpace],
              vel_grad_trial[nDOF_trial_element*nSpace];
            double p_grad_test_dV[nDOF_test_element*nSpace],
              vel_grad_test_dV[nDOF_test_element*nSpace];
            double jacDet, x, y, z, dV, h_phi;

            ck.calculateMapping_element(eN,
                                        k,
                                        mesh_dof.data(),
                                        mesh_l2g.data(),
                                        mesh_trial_ref.data(),
                                        mesh_grad_trial_ref.data(),
                                        jac,
                                        jacDet,
                                        jacInv,
                                        x,y,z);

            ck.calculateH_element(eN,
                                  k,
                                  nodeDiametersArray.data(),
                                  mesh_l2g.data(),
                                  mesh_trial_ref.data(),
                                  h_phi);

            dV = fabs(jacDet)*dV_ref.data()[k];

            eps_mu = epsFact_mu * (useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
            eps_rho = epsFact_rho * (useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);

            ck.gradTrialFromRef(&p_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
            ck.gradTrialFromRef(&vel_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);

            ck.gradFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],p_grad_trial,grad_p);
            ck.gradFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],vel_grad_trial,grad_u);
            ck.gradFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],vel_grad_trial,grad_v);

            for (int j=0; j<nDOF_trial_element;++j)
              for (int i=0; i<nSpace; ++i)
                {
                  p_grad_test_dV[j*nSpace+i] = p_grad_trial[j*nSpace+i]*dV;
                  vel_grad_test_dV[j*nSpace+i] = vel_grad_trial[j*nSpace+i]*dV;
                }

            evaluateTPInvDensityLaplaceCoefficients(eps_rho,
                                                    rho_0,
                                                    rho_1,
                                                    useVF,
                                                    vf.data()[eN_k],
                                                    phi.data()[eN_k],
                                                    mom_pp_diff_ten,
                                                    mom_uu_diff_ten,
                                                    mom_vv_diff_ten);

            // loop over test and weighted trial functions to evaluate local inner products
            for (int i=0 ; i < nDOF_test_element ; ++i)
              {
                int i_nSpace = i*nSpace ;
                for (int j=0; j < nDOF_trial_element ; ++j){
                  int j_nSpace = j*nSpace ;
                  /* local_matrix_p_p[i][j] += ck.SimpleDiffusionJacobian_weak(sdInfo_p_p_rowptr.data(), */
                  /*                                                        sdInfo_p_p_colind.data(), */
                  /*                                                        mom_pp_diff_ten, */
                  /*                                                        &p_grad_trial[j_nSpace], */
                  /*                                                        &p_grad_test_dV[i_nSpace]); */

                  /* local_matrix_u_u[i][j] += ck.SimpleDiffusionJacobian_weak(sdInfo_u_u_rowptr.data(), */
                  /*                                                        sdInfo_u_u_colind.data(), */
                  /*                                                        mom_uu_diff_ten, */
                  /*                                                        &vel_grad_trial[j_nSpace], */
                  /*                                                        &vel_grad_test_dV[i_nSpace]); */

                  /* local_matrix_v_v[i][j] += ck.SimpleDiffusionJacobian_weak(sdInfo_v_v_rowptr.data(), */
                  /*                                                        sdInfo_v_v_colind.data(), */
                  /*                                                        mom_vv_diff_ten, */
                  /*                                                        &vel_grad_trial[j_nSpace], */
                  /*                                                        &vel_grad_test_dV[i_nSpace]); */
                  local_matrix_p_p[i][j] += ck.NumericalDiffusionJacobian(mom_pp_diff_ten[0],
                                                                          &p_grad_trial[j_nSpace],
                                                                          &p_grad_test_dV[i_nSpace]);

                  local_matrix_u_u[i][j] += ck.NumericalDiffusionJacobian(mom_uu_diff_ten[0],
                                                                          &vel_grad_trial[j_nSpace],
                                                                          &vel_grad_test_dV[i_nSpace]);

                  local_matrix_v_v[i][j] += ck.NumericalDiffusionJacobian(mom_vv_diff_ten[0],
                                                                          &vel_grad_trial[j_nSpace],
                                                                          &vel_grad_test_dV[i_nSpace]);

                } // j
              } // i

          } // k

            // Write local matrix information into global system
          for (int i=0 ; i < nDOF_test_element ; ++i)
            {
              int eN_i = eN*nDOF_test_element + i;
              for (int j=0 ; j < nDOF_trial_element ; ++j)
                {
                  int eN_i_j = eN_i*nDOF_trial_element + j;
                  laplace_matrix.data()[csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_p_p.data()[eN_i_j]] += local_matrix_p_p[i][j] ;
                  laplace_matrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += local_matrix_u_u[i][j] ;
                  laplace_matrix.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_v_v.data()[eN_i_j]] += local_matrix_v_v[i][j] ;
                }
            }

        } // eN
    }

    void getTwoPhaseScaledMassOperator(arguments_dict& args)
    {
      int scale_type = args.m_iscalar["scale_type"];
      int use_numerical_viscosity = args.m_iscalar["use_numerical_viscosity"];
      int lumped = args.m_iscalar["lumped"];
      xt::pyarray<double> &mesh_trial_ref = args.m_darray["&mesh_trial_ref"];
      xt::pyarray<double> &mesh_grad_trial_ref = args.m_darray["&mesh_grad_trial_ref"];
      xt::pyarray<double> &mesh_dof = args.m_darray["&mesh_dof"];
      xt::pyarray<int>& mesh_l2g = args.m_iarray["mesh_l2g"];
      xt::pyarray<double>& dV_ref = args.m_darray["dV_ref"];
      xt::pyarray<double>& p_trial_ref = args.m_darray["p_trial_ref"];
      xt::pyarray<double>& p_test_ref = args.m_darray["p_test_ref"];
      xt::pyarray<double>& vel_trial_ref = args.m_darray["vel_trial_ref"];
      xt::pyarray<double>& vel_test_ref = args.m_darray["vel_test_ref"];
      xt::pyarray<double>& elementDiameter = args.m_darray["elementDiameter"];
      xt::pyarray<double>& nodeDiametersArray = args.m_darray["nodeDiametersArray"];
      xt::pyarray<double>& numerical_viscosity = args.m_darray["numerical_viscosity"];
      int nElements_global = args.m_iscalar["nElements_global"];
      double useMetrics = args.m_dscalar["useMetrics"];
      double epsFact_rho = args.m_dscalar["epsFact_rho"];
      double epsFact_mu = args.m_dscalar["epsFact_mu"];
      double rho_0 = args.m_dscalar["rho_0"];
      double nu_0 = args.m_dscalar["nu_0"];
      double rho_1 = args.m_dscalar["rho_1"];
      double nu_1 = args.m_dscalar["nu_1"];
      xt::pyarray<int>& p_l2g = args.m_iarray["p_l2g"];
      xt::pyarray<int>& vel_l2g = args.m_iarray["vel_l2g"];
      xt::pyarray<double>& p_dof = args.m_darray["p_dof"];
      xt::pyarray<double>& u_dof = args.m_darray["u_dof"];
      xt::pyarray<double>& v_dof = args.m_darray["v_dof"];
      xt::pyarray<double>& w_dof = args.m_darray["w_dof"];
      const double useVF = args.m_dscalar["useVF"];
      xt::pyarray<double>& vf = args.m_darray["vf"];
      xt::pyarray<double>& phi = args.m_darray["phi"];
      xt::pyarray<int>& csrRowIndeces_p_p = args.m_iarray["csrRowIndeces_p_p"];
      xt::pyarray<int>& csrColumnOffsets_p_p = args.m_iarray["csrColumnOffsets_p_p"];
      xt::pyarray<int>& csrRowIndeces_u_u = args.m_iarray["csrRowIndeces_u_u"];
      xt::pyarray<int>& csrColumnOffsets_u_u = args.m_iarray["csrColumnOffsets_u_u"];
      xt::pyarray<int>& csrRowIndeces_v_v = args.m_iarray["csrRowIndeces_v_v"];
      xt::pyarray<int>& csrColumnOffsets_v_v = args.m_iarray["csrColumnOffsets_v_v"];
      xt::pyarray<int>& csrRowIndeces_w_w = args.m_iarray["csrRowIndeces_w_w"];
      xt::pyarray<int>& csrColumnOffsets_w_w = args.m_iarray["csrColumnOffsets_w_w"];
      xt::pyarray<double>& mass_matrix = args.m_darray["mass_matrix"];
      // Step 1.1 - Initialize local matrix

      for (int eN=0 ; eN < nElements_global; ++eN){

        double local_matrix_p_p[nDOF_test_element][nDOF_trial_element];
        double local_matrix_u_u[nDOF_test_element][nDOF_trial_element];
        double local_matrix_v_v[nDOF_test_element][nDOF_trial_element];
        double eps_rho, eps_mu;

        // reset local matrix entries
        for (int i=0; i<nDOF_test_element; ++i)
          for (int j=0; j<nDOF_trial_element; ++j){
            local_matrix_p_p[i][j] = 0.0 ;
            local_matrix_u_u[i][j] = 0.0 ;
            local_matrix_v_v[i][j] = 0.0 ;
          }
        // Step 1.2 - Loop over quadrature points on element
        for (int k=0 ; k < nQuadraturePoints_element; ++k){

          int eN_k = eN*nQuadraturePoints_element+k;
          int eN_nDOF_trial_element = eN*nDOF_trial_element;
          // *** Local storage arrays ***
          double p = 0.0, u = 0.0, v= 0.0 ;
          double dV;
          double mom_p_acc = 0.0, dmom_p_acc_p = 0.0;
          double mom_u_acc = 0.0, dmom_u_acc_u = 0.0;
          double mom_v_acc = 0.0, dmom_v_acc_v = 0.0;
          double jac[nSpace*nSpace] ;
          double jacInv[nSpace*nSpace] ;
          double jacDet,x,y,z ;
          double p_test_dV[nDOF_test_element], vel_test_dV[nDOF_test_element];
          double h_phi;

          // Step 1.2.1 Calculate integration weights

          ck.calculateMapping_element(eN,
                                      k,
                                      mesh_dof.data(),
                                      mesh_l2g.data(),
                                      mesh_trial_ref.data(),
                                      mesh_grad_trial_ref.data(),
                                      jac,
                                      jacDet,
                                      jacInv,
                                      x,y,z);

          ck.calculateH_element(eN,
                                k,
                                nodeDiametersArray.data(),
                                mesh_l2g.data(),
                                mesh_trial_ref.data(),
                                h_phi);

          dV = fabs(jacDet)*dV_ref.data()[k];

          ck.valFromDOF(p_dof.data(),&p_l2g.data()[eN_nDOF_trial_element],&p_trial_ref.data()[k*nDOF_trial_element],p);
          ck.valFromDOF(u_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],&vel_trial_ref.data()[k*nDOF_trial_element],u);
          ck.valFromDOF(v_dof.data(),&vel_l2g.data()[eN_nDOF_trial_element],&vel_trial_ref.data()[k*nDOF_trial_element],v);

          eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
          eps_mu = epsFact_mu * (useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
          // Step 1.2.2 Scale test functions with integration weights.
          for (int j=0 ; j<nDOF_trial_element ; ++j){
            p_test_dV[j] = p_test_ref.data()[k*nDOF_trial_element + j]*dV;
            vel_test_dV[j] = vel_test_ref.data()[k*nDOF_trial_element + j] * dV;
          }

          // Step 1.2.2 Evaluate coefficients
          if (scale_type==0){
            evaluateTPInvViscosityMassCoefficients(use_numerical_viscosity,
                                                   numerical_viscosity.data()[eN_k],
                                                   eps_rho,
                                                   eps_mu,
                                                   rho_0,
                                                   nu_0,
                                                   rho_1,
                                                   nu_1,
                                                   useVF,
                                                   vf.data()[eN_k],
                                                   phi.data()[eN_k],
                                                   p,
                                                   u,
                                                   v,
                                                   mom_p_acc,
                                                   dmom_p_acc_p,
                                                   mom_u_acc,
                                                   dmom_u_acc_u,
                                                   mom_v_acc,
                                                   dmom_v_acc_v) ; }
          else if(scale_type==1){
            evaluateTPDensityMassCoefficients(eps_rho,
                                              rho_0,
                                              rho_1,
                                              useVF,
                                              vf.data()[eN_k],
                                              phi.data()[eN_k],
                                              p,
                                              u,
                                              v,
                                              mom_p_acc,
                                              dmom_p_acc_p,
                                              mom_u_acc,
                                              dmom_u_acc_u,
                                              mom_v_acc,
                                              dmom_v_acc_v) ;
          }

          // Step 1.2.3 Loop over test and weighted trial functions
          // to evaluate local inner product contrubtions
          for (int i=0 ; i < nDOF_test_element; ++i)
            {
              int i_nSpace = i*nSpace;
              for (int j=0 ; j < nDOF_trial_element; ++j)
                {
                  int j_nSpace = j*nSpace;
                  local_matrix_p_p[i][j] += ck.MassJacobian_weak(dmom_p_acc_p,
                                                                 p_trial_ref.data()[k*nDOF_trial_element+j],
                                                                 p_test_dV[i]) ;
                  local_matrix_u_u[i][j] += ck.MassJacobian_weak(dmom_u_acc_u,
                                                                 vel_trial_ref.data()[k*nDOF_trial_element+j],
                                                                 vel_test_dV[i]) ;
                  local_matrix_v_v[i][j] += ck.MassJacobian_weak(dmom_v_acc_v,
                                                                 vel_trial_ref.data()[k*nDOF_trial_element+j],
                                                                 vel_test_dV[i]) ;
                }//j
            }//i


        } // k

          // Step 1.3 - Write local matrix information into global system
        for (int i=0 ; i<nDOF_test_element; ++i)
          {
            int eN_i = eN*nDOF_test_element+i;
            int eN_i_i = eN_i*nDOF_trial_element + i;
            for (int j=0 ; j < nDOF_trial_element; ++j)
              {
                int eN_i_j = eN_i*nDOF_trial_element + j;
                if (lumped)
                  {
                    mass_matrix.data()[csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_p_p.data()[eN_i_i]] += local_matrix_p_p[i][j] ;
                    mass_matrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_i]] += local_matrix_u_u[i][j] ;
                    mass_matrix.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_v_v.data()[eN_i_i]] += local_matrix_v_v[i][j] ;
                  }
                else
                  {
                    mass_matrix.data()[csrRowIndeces_p_p.data()[eN_i] + csrColumnOffsets_p_p.data()[eN_i_j]] += local_matrix_p_p[i][j] ;
                    mass_matrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += local_matrix_u_u[i][j] ;
                    mass_matrix.data()[csrRowIndeces_v_v.data()[eN_i] + csrColumnOffsets_v_v.data()[eN_i_j]] += local_matrix_v_v[i][j] ;
                  }
              }
          }
      } // eN
    }
  };//RANS2P2D

  inline RANS2P2D_base* newRANS2P2D(int nSpaceIn,
                                    int nQuadraturePoints_elementIn,
                                    int nDOF_mesh_trial_elementIn,
                                    int nDOF_trial_elementIn,
                                    int nDOF_test_elementIn,
                                    int nDOF_v_trial_elementIn,
                                    int nDOF_v_test_elementIn,
                                    int nQuadraturePoints_elementBoundaryIn,
                                    int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization2D<RANS2P2D_base,RANS2P2D,CompKernel,CompKernel>(nSpaceIn,
                                                                                                    nQuadraturePoints_elementIn,
                                                                                                    nDOF_mesh_trial_elementIn,
                                                                                                    nDOF_trial_elementIn,
                                                                                                    nDOF_test_elementIn,
                                                                                                    nDOF_v_trial_elementIn,
                                                                                                    nDOF_v_test_elementIn,
                                                                                                    nQuadraturePoints_elementBoundaryIn,
                                                                                                    CompKernelFlag);
  }
}//proteus

#endif
