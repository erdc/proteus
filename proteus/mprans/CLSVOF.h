#ifndef CLSVOF_H
#define CLSVOF_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"

namespace py = pybind11;

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
    virtual void calculateResidual(arguments_dict& args)=0;
    virtual void calculateJacobian(arguments_dict& args)=0;
    virtual std::tuple<double, double, double, double, double, double, double, double, double, double, double>
        calculateMetricsAtEOS(arguments_dict& args)=0;
    virtual std::tuple<double, double, double, double, double> calculateMetricsAtETS(arguments_dict& args)=0;
    virtual void normalReconstruction(arguments_dict& args)=0;
    virtual void calculateRhsL2Proj(arguments_dict& args)=0;
    virtual void calculateLumpedMassMatrix(arguments_dict& args)=0;
    virtual void assembleSpinUpSystem(arguments_dict& args)=0;
    virtual void FCTStep(arguments_dict& args)=0;
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

      void calculateResidual(arguments_dict& args)
      {
        double dt = args.scalar<double>("dt");
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<double>& mesh_velocity_dof = args.array<double>("mesh_velocity_dof");
        double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        xt::pyarray<double>& u_grad_test_ref = args.array<double>("u_grad_test_ref");
        xt::pyarray<double>& mesh_trial_trace_ref = args.array<double>("mesh_trial_trace_ref");
        xt::pyarray<double>& mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
        xt::pyarray<double>& dS_ref = args.array<double>("dS_ref");
        xt::pyarray<double>& u_trial_trace_ref = args.array<double>("u_trial_trace_ref");
        xt::pyarray<double>& u_grad_trial_trace_ref = args.array<double>("u_grad_trial_trace_ref");
        xt::pyarray<double>& u_test_trace_ref = args.array<double>("u_test_trace_ref");
        xt::pyarray<double>& u_grad_test_trace_ref = args.array<double>("u_grad_test_trace_ref");
        xt::pyarray<double>& normal_ref = args.array<double>("normal_ref");
        xt::pyarray<double>& boundaryJac_ref = args.array<double>("boundaryJac_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        int nElements_owned = args.scalar<int>("nElements_owned");
        double useMetrics = args.scalar<double>("useMetrics");
        const xt::pyarray<double>& q_vos = args.array<double>("q_vos");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        int degree_polynomial = args.scalar<int>("degree_polynomial");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& u_dof_old = args.array<double>("u_dof_old");
        xt::pyarray<double>& velocity = args.array<double>("velocity");
        xt::pyarray<double>& velocity_old = args.array<double>("velocity_old");
        xt::pyarray<double>& q_m = args.array<double>("q_m");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_mH = args.array<double>("q_mH");
        xt::pyarray<double>& q_dV = args.array<double>("q_dV");
        xt::pyarray<double>& q_dV_last = args.array<double>("q_dV_last");
        xt::pyarray<double>& cfl = args.array<double>("cfl");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
        const xt::pyarray<double>& ebqe_vos_ext = args.array<double>("ebqe_vos_ext");
        xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
        xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
        xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
        xt::pyarray<double>& ebqe_bc_flux_u_ext = args.array<double>("ebqe_bc_flux_u_ext");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& ebqe_H = args.array<double>("ebqe_H");
        xt::pyarray<double>& ebqe_flux = args.array<double>("ebqe_flux");
        int timeOrder = args.scalar<int>("timeOrder");
        int timeStage = args.scalar<int>("timeStage");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactRedist = args.scalar<double>("epsFactRedist");
        double lambdaFact = args.scalar<double>("lambdaFact");
        xt::pyarray<double>& min_distance = args.array<double>("min_distance");
        xt::pyarray<double>& max_distance = args.array<double>("max_distance");
        xt::pyarray<double>& mean_distance = args.array<double>("mean_distance");
        xt::pyarray<double>& volume_domain = args.array<double>("volume_domain");
        double norm_factor_lagged = args.scalar<double>("norm_factor_lagged");
        double VelMax = args.scalar<double>("VelMax");
        xt::pyarray<double>& projected_qx_tn = args.array<double>("projected_qx_tn");
        xt::pyarray<double>& projected_qy_tn = args.array<double>("projected_qy_tn");
        xt::pyarray<double>& projected_qz_tn = args.array<double>("projected_qz_tn");
        xt::pyarray<double>& projected_qx_tStar = args.array<double>("projected_qx_tStar");
        xt::pyarray<double>& projected_qy_tStar = args.array<double>("projected_qy_tStar");
        xt::pyarray<double>& projected_qz_tStar = args.array<double>("projected_qz_tStar");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
        xt::pyarray<double>& H_dof = args.array<double>("H_dof");
        int preRedistancingStage = args.scalar<int>("preRedistancingStage");
        xt::pyarray<double>& interface_locator = args.array<double>("interface_locator");
        double alpha = args.scalar<double>("alpha");
        min_distance.data()[0] = 1E10;
        max_distance.data()[0] = -1E10;
        mean_distance.data()[0] = 0.;
        volume_domain.data()[0] = 0.;

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
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],
                                    jacInv,
                                    u_grad_trial);
                // get the components of the normal reconstruction //
                ck.valFromDOF(projected_qx_tn.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              qxn);
                ck.valFromDOF(projected_qy_tn.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              qyn);
                ck.valFromDOF(projected_qz_tn.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              qzn);
                //ck.valFromDOF(projected_qx_tStar.data(),
                //            &u_l2g.data()[eN_nDOF_trial_element],
                //            &u_trial_ref.data()[k*nDOF_trial_element],
                //            qxnStar);
                //ck.valFromDOF(projected_qy_tStar.data(),
                //            &u_l2g.data()[eN_nDOF_trial_element],
                //            &u_trial_ref.data()[k*nDOF_trial_element],
                //            qynStar);
                //ck.valFromDOF(projected_qz_tStar.data(),
                //            &u_l2g.data()[eN_nDOF_trial_element],
                //            &u_trial_ref.data()[k*nDOF_trial_element],
                //            qznStar);
                // get the solution (of Newton's solver)
                ck.valFromDOF(u_dof.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              u);
                // get old solution
                ck.valFromDOF(u_dof_old.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              un);
                //get the solution gradients at quad points
                ck.gradFromDOF(u_dof.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial,
                               grad_u);
                ck.gradFromDOF(u_dof_old.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial,
                               grad_un);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
                  }
                double hK=(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN])/degree_polynomial;
                //VRANS
                porosity=1.0-q_vos.data()[eN_k];

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
                if (q_dV_last.data()[eN_k] <= -100)
                  q_dV_last.data()[eN_k] = dV;
                q_dV.data()[eN_k] = dV;

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
                        relative_velocity[I] = (velocity.data()[eN_k_nSpace+I]
                                                -MOVING_DOMAIN*mesh_velocity[I]);
                        relative_velocity_old[I] = (velocity_old.data()[eN_k_nSpace+I]
                                                    -MOVING_DOMAIN*mesh_velocity[I]);
                        // compute advection term //
                        //fnp1[I] = relative_velocity[I]*Snp1; //implicit advection via BDF1
                        //fnHalf[I] = 0.5*(relative_velocity[I]*Snp1
                        //               +relative_velocity_old[I]*Sn); //implicit advection via CN
                        fnHalf[I] = 0.5*relative_velocity[I]*porosity*(Snp1+Sn);
                      }

                    //////////////////////////////
                    // CALCULATE CELL BASED CFL //
                    //////////////////////////////
                    calculateCFL(elementDiameter.data()[eN]/degree_polynomial,relative_velocity,cfl.data()[eN_k]);
                    //calculateNonlinearCFL(elementDiameter.data()[eN]/degree_polynomial,
                    //                relative_velocity,
                    //                norm_factor_lagged,
                    //                epsFactHeaviside,
                    //                lambdaFact,
                    //                cfl.data()[eN_k]);

                    ///////////////
                    // CALCULATE min, max and mean distance //
                    ///////////////
                    if (eN<nElements_owned) // locally owned?
                      {
                        min_distance.data()[0] = fmin(min_distance.data()[0],fabs(u));
                        max_distance.data()[0] = fmax(max_distance.data()[0],fabs(u));
                        mean_distance.data()[0] += fabs(u)*dV;
                        //min_distance.data()[0] = fmin(min_distance.data()[0],u);
                        //max_distance.data()[0] = fmax(max_distance.data()[0],u);
                        //mean_distance.data()[0] += u*dV;
                        volume_domain.data()[0] += dV;
                      }
                    ///////////////////
                    // SAVE SOLUTION // for other models
                    ///////////////////
                    q_u.data()[eN_k] = u;
                    q_m.data()[eN_k] = porosity*Hnp1;
                    q_H.data()[eN_k] = Hnp1;
                    q_mH.data()[eN_k] = porosity*Hnp1; //porosity*H(\phi)=(1-q_vos.data())*H(\phi)
                    // gradient //
                    for (int I=0;I<nSpace;I++)
                      q_n.data()[eN_k_nSpace+I]  = grad_u[I];
                  }

                //////////////////////////////////////////////////
                // *************** LOOP ON DOFs *************** //
                //////////////////////////////////////////////////
                for(int i=0;i<nDOF_test_element;i++)
                  {
                    int eN_i=eN*nDOF_test_element+i;
                    int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                    register int i_nSpace=i*nSpace;
                    if (preRedistancingStage==1)
                      {
                        //////////////////////
                        // EIKONAL RESIDUAL //
                        //////////////////////
                        elementResidual_u[i] +=
                          alpha*(u_dof.data()[gi]-u_dof_old.data()[gi])*delta*u_test_dV[i] // BCs
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
                          sign = Sign(u_dof.data()[gi]);
                        else if (same_sign==1)
                          {
                            same_sign = sign == Sign(u_dof.data()[gi]) ? 1 : 0;
                            sign = Sign(u_dof.data()[gi]);
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
                        //                           grad_unHalf,
                        //                                   &u_grad_test_dV[i_nSpace])
                        // TARGET for PENALIZATION. This is EXPLICIT
                        //- lambda*ck.NumericalDiffusion(1.0,
                        //                           normalReconstruction,
                        //                                   &u_grad_test_dV[i_nSpace]);
                        //}
                      }
                  }//i
                if (preRedistancingStage==0)
                  {
                    if (same_sign == 0) // This cell contains the interface
                      {
                        for(int i=0;i<nDOF_test_element;i++)
                          {
                            int gi = offset_u+stride_u*u_l2g.data()[eN*nDOF_test_element+i];
                            interface_locator.data()[gi] = 1.0;
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
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                // distribute global residual for (lumped) mass matrix
                globalResidual.data()[gi] += elementResidual_u[i];
                if (preRedistancingStage==0)
                  H_dof.data()[gi] += element_rhs_L2proj_H[i]; // int(H*wi*dx)
              }//i
          }//elements
        if (preRedistancingStage==0)
          {
            // COMPUTE LUMPED L2 PROJECTION
            for (int i=0; i<numDOFs; i++)
              H_dof.data()[i] /= lumped_mass_matrix.data()[i];
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
            register int ebN = exteriorElementBoundariesArray.data()[ebNE],
              eN  = elementBoundaryElementsArray.data()[ebN*2+0],
              ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
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
                dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt
                      + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
                ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],
                                    jacInv_ext,
                                    u_grad_trial_trace);
                //solution at quad points
                ck.valFromDOF(u_dof.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],
                              u_ext);
                ck.valFromDOF(u_dof_old.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],
                              un_ext);
                ck.gradFromDOF(u_dof.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial_trace,
                               grad_u_ext);
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                //
                //load the boundary values
                //
                double hK = elementDiameter.data()[eN]/degree_polynomial;
                double epsHeaviside = epsFactHeaviside*hK;
                double Su_ext = smoothedSign(epsHeaviside,u_ext); //Sign(u_ext)
                double Sun_ext = smoothedSign(epsHeaviside,un_ext); //Sign(un_ext)
                // NOTE: ebqe_bc_u_ext is provided by the user as BCs for VOF (i.e., 0 or 1)
                // SuDBC = Sign(uBC)
                double SuBC = (USE_SIGN_FUNCTION == 0 ? ebqe_bc_u_ext.data()[ebNE_kb] :
                               2*ebqe_bc_u_ext.data()[ebNE_kb] - 1);
                if (IMPLICIT_BCs==1)
                  bc_u_ext = (isDOFBoundary_u.data()[ebNE_kb]*SuBC
                              +(1-isDOFBoundary_u.data()[ebNE_kb])*Su_ext);
                else
                  bc_u_ext = (isDOFBoundary_u.data()[ebNE_kb]*SuBC
                              +(1-isDOFBoundary_u.data()[ebNE_kb])*Sun_ext);
                //VRANS
                porosity_ext = 1.-ebqe_vos_ext.data()[ebNE_kb];

                //
                //moving mesh
                //
                double mesh_velocity[3];
                mesh_velocity[0] = xt_ext;
                mesh_velocity[1] = yt_ext;
                mesh_velocity[2] = zt_ext;

                for (int I=0;I<nSpace;I++)
                  df_ext[I] = porosity_ext*(ebqe_velocity_ext.data()[ebNE_kb_nSpace+I]
                                            - MOVING_DOMAIN*mesh_velocity[I]);
                //
                //calculate the numerical fluxes
                //
                exteriorNumericalAdvectiveFlux(isDOFBoundary_u.data()[ebNE_kb],
                                               isFluxBoundary_u.data()[ebNE_kb],
                                               normal,
                                               bc_u_ext, //{-1,1} or {0,1}
                                               ebqe_bc_flux_u_ext.data()[ebNE_kb],
                                               IMPLICIT_BCs == 1 ? Su_ext : Sun_ext,
                                               df_ext, //VRANS includes porosity
                                               flux_ext);
                ebqe_flux.data()[ebNE_kb] = flux_ext;

                ///////////////////
                // save solution // for other models? cek need to be consistent with numerical flux
                ///////////////////
                ebqe_u.data()[ebNE_kb] = u_ext;
                // TODO: do I need ebqe_m?
                ebqe_H.data()[ebNE_kb] = smoothedHeaviside(epsHeaviside,u_ext);
                // gradient //
                for (int I=0;I<nSpace;I++)
                  ebqe_n.data()[ebNE_kb_nSpace+I]  = grad_u_ext[I];

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
                globalResidual.data()[offset_u+stride_u*u_l2g.data()[eN_i]] += elementResidual_u[i];
              }//i
          }//ebNE
        // END OF BOUNDARY //
      }

      void calculateJacobian(arguments_dict& args)
      {
        double dt = args.scalar<double>("dt");
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<double>& mesh_velocity_dof = args.array<double>("mesh_velocity_dof");
        double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        xt::pyarray<double>& u_grad_test_ref = args.array<double>("u_grad_test_ref");
        xt::pyarray<double>& mesh_trial_trace_ref = args.array<double>("mesh_trial_trace_ref");
        xt::pyarray<double>& mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
        xt::pyarray<double>& dS_ref = args.array<double>("dS_ref");
        xt::pyarray<double>& u_trial_trace_ref = args.array<double>("u_trial_trace_ref");
        xt::pyarray<double>& u_grad_trial_trace_ref = args.array<double>("u_grad_trial_trace_ref");
        xt::pyarray<double>& u_test_trace_ref = args.array<double>("u_test_trace_ref");
        xt::pyarray<double>& u_grad_test_trace_ref = args.array<double>("u_grad_test_trace_ref");
        xt::pyarray<double>& normal_ref = args.array<double>("normal_ref");
        xt::pyarray<double>& boundaryJac_ref = args.array<double>("boundaryJac_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        double useMetrics = args.scalar<double>("useMetrics");
        const xt::pyarray<double>& q_vos = args.array<double>("q_vos");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        int degree_polynomial = args.scalar<int>("degree_polynomial");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& u_dof_old = args.array<double>("u_dof_old");
        xt::pyarray<double>& velocity = args.array<double>("velocity");
        xt::pyarray<double>& cfl = args.array<double>("cfl");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& globalJacobian = args.array<double>("globalJacobian");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
        const xt::pyarray<double>& ebqe_vos_ext = args.array<double>("ebqe_vos_ext");
        xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
        xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
        xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
        xt::pyarray<double>& ebqe_bc_flux_u_ext = args.array<double>("ebqe_bc_flux_u_ext");
        xt::pyarray<int>& csrColumnOffsets_eb_u_u = args.array<int>("csrColumnOffsets_eb_u_u");
        int timeOrder = args.scalar<int>("timeOrder");
        int timeStage = args.scalar<int>("timeStage");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactRedist = args.scalar<double>("epsFactRedist");
        double lambdaFact = args.scalar<double>("lambdaFact");
        int preRedistancingStage = args.scalar<int>("preRedistancingStage");
        double norm_factor_lagged = args.scalar<double>("norm_factor_lagged");
        double alpha = args.scalar<double>("alpha");
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
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref.data()[k];
                //get the trial function gradients
                ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],
                                    jacInv,
                                    u_grad_trial);
                //get the solution
                ck.valFromDOF(u_dof.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              u);
                ck.valFromDOF(u_dof_old.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              un);
                //get the solution gradients
                ck.gradFromDOF(u_dof.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial,
                               grad_u);
                ck.gradFromDOF(u_dof_old.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial,
                               grad_un);
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
                  }
                double hK=(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN])/degree_polynomial;
                //VRANS
                porosity=1.0-q_vos.data()[eN_k];

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
                        relative_velocity[I] = (velocity.data()[eN_k_nSpace+I]
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
                              time_derivative_jacobian*(u_trial_ref.data()[k*nDOF_trial_element+j]
                                                        *u_test_dV[i])
                              // IMPLICIT TERMS: ADVECTION, DIFFUSION
                              + 0.5*ck.AdvectionJacobian_weak(df,
                                                              u_trial_ref.data()[k*nDOF_trial_element+j],
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
                    globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] +=
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
            register int ebN = exteriorElementBoundariesArray.data()[ebNE];
            register int eN  = elementBoundaryElementsArray.data()[ebN*2+0],
              ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
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
                dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt
                      + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
                ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],
                                    jacInv_ext,
                                    u_grad_trial_trace);
                ck.valFromDOF(u_dof.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],
                              u_ext);
                ck.gradFromDOF(u_dof.data(),
                               &u_l2g.data()[eN_nDOF_trial_element],
                               u_grad_trial_trace,
                               grad_u_ext);
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                //
                //load the boundary values
                //
                double hK = elementDiameter.data()[eN]/degree_polynomial;
                double epsHeaviside = epsFactHeaviside*hK;
                double dSu_ext = smoothedDerivativeSign(epsHeaviside,u_ext);
                //VRANS
                porosity_ext = 1.-ebqe_vos_ext.data()[ebNE_kb];
                //
                //moving domain
                //
                double mesh_velocity[3];
                mesh_velocity[0] = xt_ext;
                mesh_velocity[1] = yt_ext;
                mesh_velocity[2] = zt_ext;
                //std::cout<<"ext J mesh_velocity"<<std::endl;
                for (int I=0;I<nSpace;I++)
                  df_ext[I] = porosity_ext*(ebqe_velocity_ext.data()[ebNE_kb_nSpace+I]
                                            - MOVING_DOMAIN*mesh_velocity[I]);
                //
                //calculate the numerical fluxes
                //
                exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u.data()[ebNE_kb],
                                                         isFluxBoundary_u.data()[ebNE_kb],
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
                                                                u_trial_trace_ref.data()[ebN_local_kb_j]);
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
                        globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_eb_u_u.data()[ebN_i_j]]+=
                          fluxJacobian_u_u[j]*u_test_dS[i];
                      }//j
                  }//i
              }//kb
          }//ebNE
      }//computeJacobian for MCorr with CLSVOF

      std::tuple<double, double, double, double, double, double, double, double, double, double, double>
          calculateMetricsAtEOS(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        int nElements_owned = args.scalar<int>("nElements_owned");
        int useMetrics = args.scalar<int>("useMetrics");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        double degree_polynomial = args.scalar<double>("degree_polynomial");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& u0_dof = args.array<double>("u0_dof");
        xt::pyarray<double>& u_exact = args.array<double>("u_exact");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        double global_I_err = 0.0;
        double global_sI_err = 0.0;
        double global_V = 0.0;
        double global_V0 = 0.0;
        double global_sV = 0.0;
        double global_sV0 = 0.0;
        double global_D_err = 0.0;
        double global_L2_err = 0.0;
        double global_L2Banded_err = 0.0;
        double global_area_band = 0.0;
        double global_sH_L2_err = 0.0;
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
                    // get functions at quad points
                    ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],uh);
                    ck.valFromDOF(u0_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],u0);
                    u = u_exact.data()[eN_k];
                    // get gradients
                    ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                    ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_uh);

                    double epsHeaviside = epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN])/degree_polynomial;
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
                global_V += cell_V;
                global_V0 += cell_V0;
                global_sV += cell_sV;
                global_sV0 += cell_sV0;
                // metrics //
                global_I_err    += cell_I_err;
                global_sI_err += cell_sI_err;
                global_D_err    += cell_D_err;
                global_L2_err += cell_L2_err;
                global_L2Banded_err += cell_L2Banded_err;
                global_area_band += cell_area_band;
                global_sH_L2_err += cell_sH_L2_err;
              }//elements
          }
        global_D_err *= 0.5;

        return std::tuple<double, double, double, double, double, double, double, double, double, double, double>(global_I_err, global_sI_err, global_V, global_V0, global_sV, global_sV0, global_D_err, global_L2_err, global_L2Banded_err, global_area_band, global_sH_L2_err);
      }

      std::tuple<double, double, double, double, double> calculateMetricsAtETS(arguments_dict& args)
      {
        double dt = args.scalar<double>("dt");
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        int nElements_owned = args.scalar<int>("nElements_owned");
        int useMetrics = args.scalar<int>("useMetrics");
        xt::pyarray<double>& q_vos = args.array<double>("q_vos");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        double degree_polynomial = args.scalar<double>("degree_polynomial");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& u_dof_old = args.array<double>("u_dof_old");
        xt::pyarray<double>& u0_dof = args.array<double>("u0_dof");
        xt::pyarray<double>& velocity = args.array<double>("velocity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& R_vector = args.array<double>("R_vector");
        xt::pyarray<double>& sR_vector = args.array<double>("sR_vector");
        double global_V = 0.0;
        double global_V0 = 0.0;
        double global_sV = 0.0;
        double global_sV0 = 0.0;
        double global_D_err = 0.0;
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
                // get functions at quad points
                ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],unp1);
                ck.valFromDOF(u_dof_old.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],un);
                ck.valFromDOF(u0_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],u0);
                // get gradients
                ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_unp1);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
                  }

                porosity = 1.0-q_vos.data()[eN_k];

                double epsHeaviside = epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN])/degree_polynomial;
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
                    Flux_np1[I] = velocity.data()[eN_k_nSpace+I]*Sunp1;
                    sFlux_np1[I] = velocity.data()[eN_k_nSpace+I]*sSunp1;
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
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                R_vector.data()[gi] += element_R[i];
                sR_vector.data()[gi] += element_sR[i];
              }
            if (eN<nElements_owned) // just consider the locally owned cells
              {
                global_V += cell_V;
                global_V0 += cell_V0;
                global_sV += cell_sV;
                global_sV0 += cell_sV0;
                global_D_err    += cell_D_err;
              }
          }//elements
        global_D_err *= 0.5;

        return std::tuple<double, double, double, double, double>(global_V, global_V0, global_sV, global_sV0, global_D_err);
      }

      void normalReconstruction(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& weighted_lumped_mass_matrix = args.array<double>("weighted_lumped_mass_matrix");
        xt::pyarray<double>& rhs_qx = args.array<double>("rhs_qx");
        xt::pyarray<double>& rhs_qy = args.array<double>("rhs_qy");
        xt::pyarray<double>& rhs_qz = args.array<double>("rhs_qz");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& weighted_mass_matrix = args.array<double>("weighted_mass_matrix");
        for (int i=0; i<numDOFs; i++)
          {
            weighted_lumped_mass_matrix.data()[i]=0.;
            rhs_qx.data()[i]=0.;
            rhs_qy.data()[i]=0.;
            rhs_qz.data()[i]=0.;
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
                                            mesh_dof.data(),
                                            mesh_l2g.data(),
                                            mesh_trial_ref.data(),
                                            mesh_grad_trial_ref.data(),
                                            jac,
                                            jacDet,
                                            jacInv,
                                            x,y,z);
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_u);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;

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
                        norm_grad_u*u_trial_ref.data()[k*nDOF_trial_element+j]*u_test_dV[i];
                  }
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index

                weighted_lumped_mass_matrix.data()[gi] += element_weighted_lumped_mass_matrix[i];
                // rhs for reconstruction via consistent mass matrix
                rhs_qx.data()[gi] += element_rhsx_normal_reconstruction[i];
                rhs_qy.data()[gi] += element_rhsy_normal_reconstruction[i];
                rhs_qz.data()[gi] += element_rhsz_normal_reconstruction[i];
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    int eN_i_j = eN_i*nDOF_trial_element+j;
                    weighted_mass_matrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]]
                      += element_weighted_mass_matrix[i][j];
                  }
              }//i
          }//elements
      }

      void calculateRhsL2Proj(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        double he_for_disc_ICs = args.scalar<double>("he_for_disc_ICs");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& rhs_l2_proj = args.array<double>("rhs_l2_proj");
        for (int i=0; i<numDOFs; i++)
          rhs_l2_proj.data()[i]=0.;
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
                                            mesh_dof.data(),
                                            mesh_l2g.data(),
                                            mesh_trial_ref.data(),
                                            mesh_grad_trial_ref.data(),
                                            jac,
                                            jacDet,
                                            jacInv,
                                            x,y,z);
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.valFromDOF(u_dof.data(),
                              &u_l2g.data()[eN_nDOF_trial_element],
                              &u_trial_ref.data()[k*nDOF_trial_element],
                              u);
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;

                for(int i=0;i<nDOF_test_element;i++)
                  element_rhs_l2_proj[i] += he_for_disc_ICs*u*u_test_dV[i];
                //element_rhs_l2_proj[i] += u*u_test_dV[i];
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                rhs_l2_proj.data()[gi] += element_rhs_l2_proj[i];
              }//i
          }//elements
      }

      void calculateLumpedMassMatrix(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
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
                                            mesh_dof.data(),
                                            mesh_l2g.data(),
                                            mesh_trial_ref.data(),
                                            mesh_grad_trial_ref.data(),
                                            jac,
                                            jacDet,
                                            jacInv,
                                            x,y,z);
                dV = fabs(jacDet)*dV_ref.data()[k];
                //precalculate test function products with integration weights for mass matrix terms
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;

                for(int i=0;i<nDOF_test_element;i++)
                  element_lumped_mass_matrix[i] += u_test_dV[i];
              } //k
            // DISTRIBUTE //
            for(int i=0;i<nDOF_test_element;i++)
              {
                int eN_i=eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                lumped_mass_matrix.data()[gi] += element_lumped_mass_matrix[i];
              }//i
          }//elements
      }

      void assembleSpinUpSystem(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
        xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
        xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
        int nElements_global = args.scalar<int>("nElements_global");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& uInitial = args.array<double>("uInitial");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& globalMassMatrix = args.array<double>("globalMassMatrix");
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
                                            mesh_dof.data(),
                                            mesh_l2g.data(),
                                            mesh_trial_ref.data(),
                                            mesh_grad_trial_ref.data(),
                                            jac,
                                            jacDet,
                                            jacInv,
                                            x,y,z);
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref.data()[k];
                //precalculate test function products with integration weights
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;

                //////////////////
                // LOOP ON DOFs //
                //////////////////
                for(int i=0;i<nDOF_test_element;i++)
                  {
                    elementResidual_u[i] += uInitial.data()[eN_k]*u_test_dV[i];
                    for(int j=0;j<nDOF_trial_element;j++)
                      {
                        elementMassMatrix_u_u[i][j] +=
                          u_trial_ref.data()[k*nDOF_trial_element+j]*u_test_dV[i];
                      }//j
                  }//i
              }//k
            //
            //load into element Jacobian into global Jacobian
            //
            for (int i=0;i<nDOF_test_element;i++)
              {
                int eN_i = eN*nDOF_test_element+i;
                int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                globalResidual.data()[gi] += elementResidual_u[i];
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    int eN_i_j = eN_i*nDOF_trial_element+j;
                    globalMassMatrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] +=
                      elementMassMatrix_u_u[i][j];
                  }//j
              }//i
          }//elements
      }

      void FCTStep(arguments_dict& args)
      {
        int NNZ = args.scalar<int>("NNZ");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
        xt::pyarray<double>& soln = args.array<double>("soln");
        xt::pyarray<double>& solH = args.array<double>("solH");
        xt::pyarray<double>& solL = args.array<double>("solL");
        xt::pyarray<double>& limited_solution = args.array<double>("limited_solution");
        xt::pyarray<int>& csrRowIndeces_DofLoops = args.array<int>("csrRowIndeces_DofLoops");
        xt::pyarray<int>& csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops");
        xt::pyarray<double>& MassMatrix = args.array<double>("matrix");
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
            double solHi = solH.data()[i];
            double solLi = solL.data()[i];
            double mi = lumped_mass_matrix.data()[i];

            double mini=-1.0, maxi=1.0; // global FCT
            //double mini=1.0E10, maxi=-1.0E10;
            double Pposi=0, Pnegi=0;
            // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
            for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
              {
                int j = csrColumnOffsets_DofLoops.data()[offset];
                // i-th row of flux correction matrix
                FluxCorrectionMatrix[ij] = ((i==j ? 1. : 0.)*mi - MassMatrix.data()[ij]) * (solH.data()[j]-solHi);

                //mini = fmin(mini,limited_solution.data()[j]);
                //maxi = fmax(maxi,limited_solution.data()[j]);

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
            for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
              {
                int j = csrColumnOffsets_DofLoops.data()[offset];
                ith_Limiter_times_FluxCorrectionMatrix +=
                  ((FluxCorrectionMatrix[ij]>0) ? std::min(Rposi,Rneg[j]) : std::min(Rnegi,Rpos[j]))
                  * FluxCorrectionMatrix[ij];
                //ith_Limiter_times_FluxCorrectionMatrix += FluxCorrectionMatrix[ij];
                //update ij
                ij+=1;
              }
            limited_solution.data()[i] = solL.data()[i] + 1./lumped_mass_matrix.data()[i]*ith_Limiter_times_FluxCorrectionMatrix;
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
