#ifndef Poisson_H1
#define Poisson_H1
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <set>
#include "CompKernel.h"
#include "ModelFactory.h"

#define GAMMA 1.4
#define Sign(z) (z >= 0.0 ? 1.0 : -1.0)

#define POWER_SMOOTHNESS_INDICATOR 2
#define IS_BETAij_ONE 0

/////////////////////
//ENTROPY FUNCTION //
/////////////////////
// Power entropy //
#define entropy_power 2. // phiL and phiR are dummy variables
#define ENTROPY(phi,dummyL,dummyR) 1./entropy_power*std::pow(fabs(phi),entropy_power)
#define DENTROPY(phi,dummyL,dummyR) entropy_power == 1. ? 0. : std::pow(fabs(phi),entropy_power-1.)*(phi >= 0. ? 1. : -1.)

#define ENTROPY_LOG(phi,phiL,phiR) std::log(fabs((phi-phiL)*(phiR-phi))+1E-14)
#define DENTROPY_LOG(phi,phiL,phiR) (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(fabs((phi-phiL)*(phiR-phi))+1E-14)

namespace proteus
{
  class Poisson_base
  {
    //The base class defining the interface
  public:
    virtual ~Poisson_base(){}
    virtual void calculateResidual(//element
                                           double dt,
                                           double * mesh_trial_ref,
                                           double * mesh_grad_trial_ref,
                                           double * mesh_dof,
                                           double * mesh_velocity_dof,
                                           int * mesh_l2g,
                                           double * dV_ref,
                                           double * u_trial_ref,
                                           double * u_grad_trial_ref,
                                           double * u_test_ref,
                                           double * u_grad_test_ref,
                                           double * mesh_trial_trace_ref,
                                           double * mesh_grad_trial_trace_ref,
                                           double * dS_ref,
                                           double * u_trial_trace_ref,
                                           double * u_grad_trial_trace_ref,
                                           double * u_test_trace_ref,
                                           double * u_grad_test_trace_ref,
                                           double * normal_ref,
                                           double * boundaryJac_ref,
                                           int nElements_global,
                                           double alphaBDF,
                                           int * u_l2g,
                                           double * elementDiameter,
                                           double * nodeDiametersArray,
                                           double * u_dof,
                                           double * r,
                                           int u_offset,
                                           int u_stride,
                                           int nExteriorElementBoundaries_global,
                                           int * exteriorElementBoundariesArray,
                                           int * elementBoundariesArray,
                                           int * elementBoundaryElementsArray,
                                           int * elementBoundaryLocalElementBoundariesArray,
                                           int u_ndofs,
                                           int NNZ,
                                           int * csrRowIndeces_DofLoops,
                                           int * csrColumnOffsets_DofLoops,
                                           int * csrRowIndeces_CellLoops_rho,
                                           int * csrColumnOffsets_CellLoops_rho,
                                           int * csrColumnOffsets_eb_CellLoops_rho,
                                           double * quantDOFs,
                                           double * ML,
                                           double* isActiveDOF,
                                           int USE_SBM)=0;
    virtual void calculateJacobian(//element
                                    double dt,
                                    double* mesh_trial_ref,
                                    double* mesh_grad_trial_ref,
                                    double* mesh_dof,
                                    double* mesh_velocity_dof,
                                    int* mesh_l2g,
                                    double* dV_ref,
                                    double* u_trial_ref,
                                    double* u_grad_trial_ref,
                                    double* u_test_ref,
                                    double* u_grad_test_ref,
                                    double* mesh_trial_trace_ref,
                                    double* mesh_grad_trial_trace_ref,
                                    double* dS_ref,
                                    double* u_trial_trace_ref,
                                    double* u_grad_trial_trace_ref,
                                    double* u_test_trace_ref,
                                    double* u_grad_test_trace_ref,
                                    double* normal_ref,
                                    double* boundaryJac_ref,
                                    int nElements_global,
                                    int* u_l2g,
                                    double* elementDiameter,
                                    double* u_dof,
                                    double* velocity,
                                    double* q_m_betaBDF,
                                    int* csrRowIndeces_u_u,
                                    int* csrColumnOffsets_u_u,
                                    int* csrColumnOffsets_eb_u_u,
                                    double* globalJacobian,
                                    int nExteriorElementBoundaries_global,
                                    int* exteriorElementBoundariesArray,
                                    int* elementBoundariesArray,
                                    int* elementBoundaryElementsArray,
                                    int* elementBoundaryLocalElementBoundariesArray,
                                    int USE_SBM
                                    )=0;
    virtual void calculateMassMatrix(//element
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
                                     int* u_l2g,
                                     double* elementDiameter,
                                     int degree_polynomial,
                                     double* u_dof,
                                     double* velocity,
                                     double* q_m_betaBDF,
                                     double* cfl,
                                     double* q_numDiff_u_last,
                                     int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                                     double* globalJacobian,
                                     int nExteriorElementBoundaries_global,
                                     int* exteriorElementBoundariesArray,
                                     int* elementBoundaryElementsArray,
                                     int* elementBoundaryLocalElementBoundariesArray,
                                     double* ebqe_velocity_ext,
                                     int* isDOFBoundary_u,
                                     double* ebqe_rd_u_ext,
                                     double* ebqe_bc_u_ext,
                                     int* csrColumnOffsets_eb_u_u,
                                     int PURE_BDF,
                                     int LUMPED_MASS_MATRIX
                                     )=0;
  };

  template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class Poisson : public Poisson_base
    {
    public:
      const int nDOF_test_X_trial_element;
      CompKernelType ck;

//      std::vector<int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;
      static const double C_sbm=1000;//penalty constant for sbm
      static const double beta_sbm=0.0;//tangent penalty constant for sbm
      static const int nParticles=1;
      double ball_center[3*nParticles];
      double ball_radius[nParticles];
      double ball_velocity[3*nParticles];
      double ball_angular_velocity[3*nParticles];
      static const int use_ball_as_particle=1;
      int has_surrogate_boundaries;

    Poisson():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck(),
      has_surrogate_boundaries(0)
      {
        ball_center[0]=0.0;
        ball_center[1]=0.0;
        ball_center[2]=0.0;
        ball_radius[0]=0.7;
        ball_velocity[0]=0.0;
        ball_velocity[1]=0.0;
        ball_velocity[2]=0.0;
        ball_angular_velocity[0]=0.0;
        ball_angular_velocity[1]=0.0;
        ball_angular_velocity[2]=0.0;
      }

      inline void evaluateCoefficients(const double v[nSpace],
                       const double& u,
                       const double grad_u[nSpace],
                       double& m,
                       double& dm,
                       double& H,
                       double dH[nSpace])
      {
    m = u;
    dm=1.0;
    H = 0.0;
    for (int I=0; I < nSpace; I++)
      {
        H += v[I]*grad_u[I];
        dH[I] = v[I];
      }
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
    void calculateCFL(const double& elementDiameter,
              const double df,
              double& cfl)
      {
    double h,nrm_v;
    h = elementDiameter;
    nrm_v=fabs(df);
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

      void exteriorNumericalFlux(const double n[nSpace],
                 const double& bc_u,
                 const double& u,
                 const double velocity[nSpace],
                 const double velocity_movingDomain[nSpace],
                 double& flux)
      {
    double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
    for (int I=0; I < nSpace; I++)
      {
        flow_fluid += n[I]*velocity[I];
        //flow_movingDomain -= n[I]*velocity_movingDomain[I];
      }
    flow_total = flow_fluid+flow_movingDomain;
    if (flow_total > 0.0)
      {
        flux = u*flow_movingDomain;
      }
    else
      {
        flux = bc_u*flow_movingDomain - flow_fluid*(u-bc_u);
      }
      }

      inline
    void exteriorNumericalFluxDerivative(const double n[nSpace],
                         const double velocity[nSpace],
                         const double velocity_movingDomain[nSpace],
                         double& dflux)
      {
    double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
    for (int I=0; I < nSpace; I++)
      {
        flow_fluid += n[I]*velocity[I];
        //flow_movingDomain -= n[I]*velocity_movingDomain[I];
      }
    flow_total=flow_fluid+flow_movingDomain;
    if (flow_total > 0.0)
      {
        dflux = flow_movingDomain;
      }
    else
      {
        dflux = -flow_fluid;
      }
      }

      inline
    double sign(const double phi,
            const double eps) //eps=epsFactRedistancing=0.33*he (for instance)
      {
    double H;
    if (phi > eps)
      H = 1.0;
    else if (phi < -eps)
      H = 0.0;
    else if (phi == 0.0)
      H = 0.5;
    else
      H = 0.5*(1.0 + phi/eps + std::sin(M_PI*phi/eps)/M_PI);
    return -1.0 + 2.0*H;
    //return (u > 0 ? 1. : -1.)*(u ==0 ? 0. : 1.);
    //double tol_sign 0.1;
    //return (u > tol_sign*epsFactRedistancing ? 1. : -1.)*((u > -tol_sign*epsFactRedistancing && u < tol_sign*epsFactRedistancing) ? 0. : 1.);
      }


      void calculateJacobian(//element
                 double dt,
                 double* mesh_trial_ref,
                 double* mesh_grad_trial_ref,
                 double* mesh_dof,
                 double* mesh_velocity_dof,
                 int* mesh_l2g,
                 double* dV_ref,
                 double* u_trial_ref,
                 double* u_grad_trial_ref,
                 double* u_test_ref,
                 double* u_grad_test_ref,
                 double* mesh_trial_trace_ref,
                 double* mesh_grad_trial_trace_ref,
                 double* dS_ref,
                 double* u_trial_trace_ref,
                 double* u_grad_trial_trace_ref,
                 double* u_test_trace_ref,
                 double* u_grad_test_trace_ref,
                 double* normal_ref,
                 double* boundaryJac_ref,
                 int nElements_global,
                 int* u_l2g,
                 double* elementDiameter,
                 double* u_dof,
                 double* velocity,
                 double* q_m_betaBDF,
                 int* csrRowIndeces_u_u,
                 int* csrColumnOffsets_u_u,
                 int* csrColumnOffsets_eb_u_u,
                 double* globalJacobian,
                 int nExteriorElementBoundaries_global,
                 int* exteriorElementBoundariesArray,
                 int* elementBoundariesArray,
                 int* elementBoundaryElementsArray,
                 int* elementBoundaryLocalElementBoundariesArray,
                 int USE_SBM)
      {

          double* node_coord=mesh_dof;
          std::vector<int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;
          ///////////////////////////////////////////////////////////////////////////////////////////////
          for(int eN=0;eN<nElements_global;eN++)
          {
          ///////////////////////////////////////////////////////////////////////////////////////////////
              register double  local_dd[nDOF_test_element][nDOF_trial_element];

              double element_active=1;//use 1 since by default it is assembled over all elements

              for (int i=0;i<nDOF_test_element;i++)
              {
                  for (int j=0;j<nDOF_trial_element;j++)
                  {
                          local_dd[i][j]=0.0;
                  }
              }
              if(USE_SBM>0)///////YYYYYYY: has to update every time since isActive is assigend to be 1 in getResidual
              {
                  ///////////////////////////////////////////////////////////////////////////////////////////////
                  /////YYYYYYYYY: this is a bug since maybe isActive is reset to be 0 if it is 1.
                  //                  for (int i=0;i<nDOF_test_element;i++)
                  //                  {
                  //                      isActiveDOF[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=0.0;//since it has 1 by default
                  //                      quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=0.0;
                  //                  }
                  ///////////////////////////////////////////////////////////////////////////////////////////////
                  double _distance[nDOF_mesh_trial_element]={0.0};
                  int pos_counter=0;
                  for (int I=0;I<nDOF_mesh_trial_element;I++)
                  {
                      if(use_ball_as_particle==1)
                      {
                          get_distance_to_ball(nParticles, ball_center, ball_radius,
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+0],
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+1],
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+2],
                                  _distance[I]);
                      }
                      else
                      {
                          //                          _distance[I] = phi_solid_nodes[mesh_l2g[eN*nDOF_mesh_trial_element+I]];
                      }
                      if ( _distance[I] >= 0)
                          pos_counter++;
                  }
                  if (pos_counter == 2)
                  {
                      element_active=0.0;
                      int opp_node=-1;
                      for (int I=0;I<nDOF_mesh_trial_element;I++)
                      {
                          if (_distance[I] < 0)
                          {
                              opp_node = I;
                          }
                      }
                      assert(opp_node >=0);
                      assert(opp_node <nDOF_mesh_trial_element);
                      int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                      surrogate_boundaries.push_back(ebN);
                      //now find which element neighbor this element is
                      if (eN == elementBoundaryElementsArray[eN*2+0])
                          surrogate_boundary_elements.push_back(1);
                      else
                          surrogate_boundary_elements.push_back(0);

                      //check which particle this surrogate edge is related to.
                      int j=-1;
                      if(use_ball_as_particle==1)
                      {
                          double middle_point_coord[3]={0.0};
                          double middle_point_distance;
                          if(opp_node == 0)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+1]);
                          }
                          else if(opp_node == 1)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+1]);
                          }
                          else if(opp_node == 2)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+1]);
                          }
                          j = get_distance_to_ball(nParticles, ball_center, ball_radius,
                                  middle_point_coord[0],middle_point_coord[1],middle_point_coord[2],
                                  middle_point_distance);

                      }
                      else
                      {
                          //                          //The method is to check one quadrature point inside of this element.
                          //                          //It works based on the assumption that the distance between any two particles
                          //                          //is larger than 2*h_min, otherwise it depends on the choice of the quadrature point
                          //                          //or one edge belongs to two particles .
                          //                          //But in any case, phi_s is well defined as the minimum.
                          //                          double distance=1e10, distance_to_ith_particle;
                          //                          for (int i=0;i<nParticles;++i)
                          //                          {
                          //                              distance_to_ith_particle=particle_signed_distances[i*nElements_global*nQuadraturePoints_element
                          //                                                                                 +eN*nQuadraturePoints_element
                          //                                                                                 +0];//0-th quadrature point
                          //                              if (distance_to_ith_particle<distance)
                          //                              {
                          //                                  distance = distance_to_ith_particle;
                          //                                  j = i;
                          //                              }
                          //                          }
                      }
                      surrogate_boundary_particle.push_back(j);
                  }
                  else if (pos_counter == 3)
                  {
                      element_active=1.0;
                  }
                  else
                  {
                      element_active=0.0;
                  }
              }
          ///////////////////////////////////////////////////////////////////////////////////////////////
              for (int k=0;k<nQuadraturePoints_element;k++)
              {
          ///////////////////////////////////////////////////////////////////////////////////////////////
                  //compute indeces and declare local storage
                  register int eN_k = eN*nQuadraturePoints_element+k,
                    eN_k_nSpace = eN_k*nSpace,
                    eN_nDOF_trial_element = eN*nDOF_trial_element;
                  register double
                    u_test_dV[nDOF_trial_element],
                    u_grad_trial[nDOF_trial_element*nSpace],
                    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                    dV,x,y,z,xt,yt,zt;
                  register double u_at_qp,grad_u_at_qp[nSpace], f_at_qp;
                  ck.calculateMapping_element(eN,
                                  k,
                                  node_coord,////////////////////////////use updated mesh
                                  mesh_l2g,
                                  mesh_trial_ref,
                                  mesh_grad_trial_ref,
                                  jac,
                                  jacDet,
                                  jacInv,x,y,z);
                  ck.calculateMappingVelocity_element(eN,
                                      k,
                                      mesh_velocity_dof,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      xt,yt,zt);
                  dV = fabs(jacDet)*dV_ref[k];

                  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);

                  for (int j=0;j<nDOF_trial_element;j++)
                      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
          ///////////////////////////////////////////////////////////////////////////////////////////////

                  for (int i=0;i<nDOF_test_element;i++)
                  {
                      for (int j=0;j<nDOF_trial_element;j++)
                      {
                              local_dd[i][j] += (u_grad_trial[i*nSpace+0]*u_grad_trial[j*nSpace+0]
                                               +u_grad_trial[i*nSpace+1]*u_grad_trial[j*nSpace+1])*dV;
                      }

                  }

          ///////////////////////////////////////////////////////////////////////////////////////////////
              }//end of eN_k
          ///////////////////////////////////////////////////////////////////////////////////////////////
              for(int i=0;i<nDOF_test_element;i++)
              {
                  int eN_i=eN*nDOF_test_element+i;
                  for (int j=0;j<nDOF_trial_element;j++)
                  {
                      int eN_i_j = eN_i*nDOF_trial_element+j;
                      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += element_active*local_dd[j][i];

                  }//j
              }
          }//end of eN

          if(USE_SBM>0)
          {
              //loop over the surrogate boundaries in SB method and assembly into jacobian
              //
              for (int ebN_s=0;ebN_s < surrogate_boundaries.size();ebN_s++)
              {
                  register int ebN = surrogate_boundaries[ebN_s],
                          eN = elementBoundaryElementsArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                          ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                          eN_nDOF_trial_element = eN*nDOF_trial_element;

//                  if (ebN >= nElementBoundaries_owned) continue;
                  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                  {
                      register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                              ebN_kb_nSpace = ebN_kb*nSpace,
                              ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                              ebN_local_kb_nSpace = ebN_local_kb*nSpace;

                      register double u_ext=0.0,
                              bc_u_ext=0.0,
                              bc_v_ext=0.0,
                              grad_u_ext[nSpace],
                              jac_ext[nSpace*nSpace],
                              jacDet_ext,
                              jacInv_ext[nSpace*nSpace],
                              boundaryJac[nSpace*(nSpace-1)],
                              metricTensor[(nSpace-1)*(nSpace-1)],
                              metricTensorDetSqrt,
                              u_grad_trial_trace[nDOF_trial_element*nSpace],
                              dS,
                              u_test_dS[nDOF_test_element],
                              normal[2],
                              x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                              u_grad_test_dS[nDOF_trial_element*nSpace],
                              G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty;
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
                      dS = metricTensorDetSqrt*dS_ref[kb];
                      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                      //compute shape and solution information
                      //shape
                      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                      //solution and gradients
                      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
                      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
                      //precalculate test function products with integration weights
                      for (int j=0;j<nDOF_trial_element;j++)
                      {
                          u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                          for (int I=0;I<nSpace;I++)
                              u_grad_test_dS[j*nSpace+I] = u_grad_trial_trace[j*nSpace+I]*dS;//cek hack, using trial
                      }
                      //
                      //load the boundary values
                      //
                      bc_u_ext = 0.0;
                      bc_v_ext = 0.0;
                      ck.calculateGScale(G,normal,h_penalty);
                      penalty = h_penalty;
                      //
                      //update the global Jacobian from the flux Jacobian
                      //

                      double dist = 0.0;
                      double distance[2], P_normal[2], P_tangent[2], normal_Omega[2]; // distance vector, normal and tangent of the physical boundary

                      if(use_ball_as_particle==1)
                      {
                          get_distance_to_ball(nParticles,ball_center,ball_radius,
                                  x_ext,y_ext,z_ext,
                                  dist);
                          get_normal_to_ith_ball(nParticles,ball_center,ball_radius,
                                  surrogate_boundary_particle[ebN_s],
                                  x_ext,y_ext,z_ext,
                                  P_normal[0],P_normal[1]);
                          get_velocity_to_ith_ball(nParticles,ball_center,ball_radius,
                                  ball_velocity, ball_angular_velocity,
                                  surrogate_boundary_particle[ebN_s],
                                  x_ext-dist*P_normal[0],
                                  y_ext-dist*P_normal[1],
                                  0.0,//z_ext,
                                  bc_u_ext,bc_v_ext);
                      }
                      else
                      {
//                          dist = ebq_global_phi_solid[ebN_kb];
//                          P_normal[0] = ebq_global_grad_phi_solid[ebN_kb*nSpace+0];
//                          P_normal[1] = ebq_global_grad_phi_solid[ebN_kb*nSpace+1];
//                          bc_u_ext = ebq_particle_velocity_solid [ebN_kb*nSpace+0];
//                          bc_v_ext = ebq_particle_velocity_solid [ebN_kb*nSpace+1];
                      }
                      distance[0] = -P_normal[0]*dist;
                      distance[1] = -P_normal[1]*dist;
                      P_tangent[0]= -P_normal[1];
                      P_tangent[1]= P_normal[0];

                      normal_Omega[0] = -P_normal[0];
                      normal_Omega[1] = -P_normal[1];

                      double C_adim = C_sbm/h_penalty;
                      double beta_adim = beta_sbm/h_penalty;

                      for (int i=0;i<nDOF_test_element;i++)
                      {
                          register int eN_i = eN*nDOF_test_element+i;
                          double phi_i = u_test_dS[i];
                          double* grad_phi_i = &u_grad_test_dS[i*nSpace+0];
                          const double grad_phi_i_dot_d = get_dot_product(grad_phi_i,distance);
                          const double grad_phi_i_dot_t = get_dot_product(P_tangent,grad_phi_i);

                          for (int j=0;j<nDOF_trial_element;j++)
                          {
                              register int ebN_i_j = (ebN*4
                                                      +surrogate_boundary_elements[ebN_s]*3/////////YYYYY
                                                      )*nDOF_test_X_trial_element
                                                      + i*nDOF_trial_element + j;

                              double phi_j = u_test_dS[j]/dS;
                              const double grad_phi_j[2]={u_grad_test_dS[j*nSpace+0]/dS,
                                                          u_grad_test_dS[j*nSpace+1]/dS};
                              const double grad_phi_j_dot_d = get_dot_product(distance, grad_phi_j);
                              const double grad_phi_j_dot_t = get_dot_product(P_tangent,grad_phi_j);

                              // (1)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] +=
                                      C_adim*(phi_i+grad_phi_i_dot_d)*(phi_j+grad_phi_j_dot_d);

                              // (2)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] -=
                                      phi_i * get_dot_product(grad_phi_j,normal);/////YY: not normal_Omega

                              // (3)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] -=
                                      get_dot_product(grad_phi_i,normal)*(phi_j+grad_phi_j_dot_d);/////YY: not normal_Omega
                          }//j
                      }//i
                  }//kb
              }//ebN_s
          }
          ///////////////////////////////////////////////////////////////////////////////////////////////
          return;
      }//computeJacobian

      double get_max_wave_speed(double ul, double ur, std::complex<double> w, std::complex<double> n)
      {
          return fmax(fabs(1.0-std::real(w*n)), fabs(-1.0-std::real(w*n)));
      }

      double get_sound_speed(double rho, double Mx, double My, double E)
      {
          double e = (E-(0.5*Mx*Mx+0.5*My*My)/(rho+1.0e-10))/(rho+1.0e-10);
          double p = (GAMMA-1.0)*rho*e;
          return std::sqrt(GAMMA*p/(rho+1.0e-10));
      }
      double get_rhs(double x, double y, double z)
      {
          return 2.0;//2*(1+x)*(1-x)+2*(1+y)*(1-y);
      }

      double get_cross_product(const double *u, const double *v)
      {
          return u[0]*v[1]-u[1]*v[0];
      }
      double get_dot_product(const double *u, const double *v)
      {
          return u[0]*v[0]+u[1]*v[1];
      }
      int get_distance_to_ball(int n_balls,double* ball_center, double* ball_radius, double x, double y, double z, double& distance)
      {
          int index = 0;
          int i = 0;
          distance = std::sqrt((ball_center[i*3+0]-x)*(ball_center[i*3+0]-x)
                                  +(ball_center[i*3+1]-y)*(ball_center[i*3+1]-y)
//                                  +(ball_center[i*3+2]-z)*(ball_center[i*3+2]-z)
                                  )-ball_radius[0];

          distance = -distance;
          return index;
      }

      void get_normal_to_ith_ball(int n_balls,double* ball_center, double* ball_radius,
                                  int I,
                                  double x, double y, double z,
                                  double& nx, double& ny)
      {
          double distance = std::sqrt((ball_center[I*3+0]-x)*(ball_center[I*3+0]-x)
                                    + (ball_center[I*3+1]-y)*(ball_center[I*3+1]-y)
//                                  + (ball_center[I*3+2]-z)*(ball_center[I*3+2]-z)
                            );
          nx = (x-ball_center[I*3+0])/(distance+1e-10);
          ny = (y-ball_center[I*3+1])/(distance+1e-10);

          nx = -nx;
          ny = -ny;
      }
      void get_velocity_to_ith_ball(int n_balls,double* ball_center, double* ball_radius,
                                    double* ball_velocity, double* ball_angular_velocity,
                                    int I,
                                    double x, double y, double z,
                                    double& vx, double& vy)
      {
          vx = ball_velocity[3*I + 0] - ball_angular_velocity[3*I + 2]*(y-ball_center[3*I + 1]);
          vy = ball_velocity[3*I + 1] + ball_angular_velocity[3*I + 2]*(x-ball_center[3*I + 0]);
      }
      void calculateResidual(//element
              double dt,
              double * mesh_trial_ref,
              double * mesh_grad_trial_ref,
              double * mesh_dof,
              double * mesh_velocity_dof,
              int * mesh_l2g,
              double * dV_ref,
              double * u_trial_ref,
              double * u_grad_trial_ref,
              double * u_test_ref,
              double * u_grad_test_ref,
              double * mesh_trial_trace_ref,
              double * mesh_grad_trial_trace_ref,
              double * dS_ref,
              double * u_trial_trace_ref,
              double * u_grad_trial_trace_ref,
              double * u_test_trace_ref,
              double * u_grad_test_trace_ref,
              double * normal_ref,
              double * boundaryJac_ref,
              int nElements_global,
              double alphaBDF,
              int * u_l2g,
              double * elementDiameter,
              double * nodeDiametersArray,
              double * u_dof,
              double * r,
              int u_offset,
              int u_stride,
              int nExteriorElementBoundaries_global,
              int * exteriorElementBoundariesArray,
              int * elementBoundariesArray,
              int * elementBoundaryElementsArray,
              int * elementBoundaryLocalElementBoundariesArray,
              int u_ndofs,
              int NNZ,
              int * csrRowIndeces_DofLoops,
              int * csrColumnOffsets_DofLoops,
              int * csrRowIndeces_CellLoops_rho,
              int * csrColumnOffsets_CellLoops_rho,
              int * csrColumnOffsets_eb_CellLoops_rho,
              double * quantDOFs,
              double * ML,
              double* isActiveDOF,
              int USE_SBM)
      {
          double* node_coord=mesh_dof;
          std::vector<int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;

          if(USE_SBM>0)
              for(int i=0;i<u_ndofs;i++)
              {
                  isActiveDOF[i]=0.0;//since it has 1 by default
                  quantDOFs[i]=0.0;
              }
          ///////////////////////////////////////////////////////////////////////////////////////////////
          for(int eN=0;eN<nElements_global;eN++)
          {
          ///////////////////////////////////////////////////////////////////////////////////////////////
              register double  local_dd[nDOF_test_element][nDOF_trial_element];
              register double  local_r[nDOF_test_element];

              double element_active=1;//use 1 since by default it is assembled over all elements

              for (int i=0;i<nDOF_test_element;i++)
              {
                  for (int j=0;j<nDOF_trial_element;j++)
                  {
                          local_dd[i][j]=0.0;
                  }
                  local_r[i]=0.0;
              }
              ///////////////////////////////////////////////////////////////////////////////////////////////check surrogate boundary
              //1. define offset_u, stride_u
              //2. delete offset_v and stride_v since Poisson equation has 1 variable
              //3. use u_l2g
              if(USE_SBM>0)///////YYYYYYY: has to update every time since isActive is assigend to be 1 in getResidual
              {
                  const int stride_u=u_stride,offset_u=u_offset;
                  ///////////////////////////////////////////////////////////////////////////////////////////////
                  /////YYYYYYYYY: this is a bug since maybe isActive is reset to be 0 if it is 1.
//                  for (int i=0;i<nDOF_test_element;i++)
//                  {
//                      isActiveDOF[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=0.0;//since it has 1 by default
//                      quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=0.0;
//                  }
                  ///////////////////////////////////////////////////////////////////////////////////////////////
                  double _distance[nDOF_mesh_trial_element]={0.0};
                  int pos_counter=0;
                  for (int I=0;I<nDOF_mesh_trial_element;I++)
                  {
                      if(use_ball_as_particle==1)
                      {
                          get_distance_to_ball(nParticles, ball_center, ball_radius,
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+0],
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+1],
                                  mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+I]+2],
                                  _distance[I]);
                      }
                      else
                      {
//                          _distance[I] = phi_solid_nodes[mesh_l2g[eN*nDOF_mesh_trial_element+I]];
                      }
                      if ( _distance[I] >= 0)
                          pos_counter++;
                  }
                  if (pos_counter == 2)
                  {
                      element_active=0.0;
                      int opp_node=-1;
                      for (int I=0;I<nDOF_mesh_trial_element;I++)
                      {
                          quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + I]] = 2.0;
                          if (_distance[I] < 0)
                          {
                              opp_node = I;
                              quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + I]] = 1.0;
                          }
                      }
                      assert(opp_node >=0);
                      assert(opp_node <nDOF_mesh_trial_element);
                      int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                      surrogate_boundaries.push_back(ebN);
                      //now find which element neighbor this element is
                      if (eN == elementBoundaryElementsArray[eN*2+0])
                          surrogate_boundary_elements.push_back(1);
                      else
                          surrogate_boundary_elements.push_back(0);

                      //check which particle this surrogate edge is related to.
                      int j=-1;
                      if(use_ball_as_particle==1)
                      {
                          double middle_point_coord[3]={0.0};
                          double middle_point_distance;
                          if(opp_node == 0)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+1]);
                          }
                          else if(opp_node == 1)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+2]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+1]);
                          }
                          else if(opp_node == 2)
                          {
                              middle_point_coord[0] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+0]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+0]);
                              middle_point_coord[1] = 0.5*(mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+0]+1]+mesh_dof[3*mesh_l2g[eN*nDOF_mesh_trial_element+1]+1]);
                          }
                          j = get_distance_to_ball(nParticles, ball_center, ball_radius,
                                  middle_point_coord[0],middle_point_coord[1],middle_point_coord[2],
                                  middle_point_distance);

                      }
                      else
                      {
//                          //The method is to check one quadrature point inside of this element.
//                          //It works based on the assumption that the distance between any two particles
//                          //is larger than 2*h_min, otherwise it depends on the choice of the quadrature point
//                          //or one edge belongs to two particles .
//                          //But in any case, phi_s is well defined as the minimum.
//                          double distance=1e10, distance_to_ith_particle;
//                          for (int i=0;i<nParticles;++i)
//                          {
//                              distance_to_ith_particle=particle_signed_distances[i*nElements_global*nQuadraturePoints_element
//                                                                                 +eN*nQuadraturePoints_element
//                                                                                 +0];//0-th quadrature point
//                              if (distance_to_ith_particle<distance)
//                              {
//                                  distance = distance_to_ith_particle;
//                                  j = i;
//                              }
//                          }
                      }
                      surrogate_boundary_particle.push_back(j);
                  }
                  else if (pos_counter == 3)
                  {
                      element_active=1.0;
                      for (int i=0;i<nDOF_test_element;i++)
                      {
                          isActiveDOF[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=1.0;
                      }
                  }
                  else
                  {
                      element_active=0.0;
                  }
              }
          ///////////////////////////////////////////////////////////////////////////////////////////////
              for (int k=0;k<nQuadraturePoints_element;k++)
              {
          ///////////////////////////////////////////////////////////////////////////////////////////////
                  //compute indeces and declare local storage
                  register int eN_k = eN*nQuadraturePoints_element+k,
                    eN_k_nSpace = eN_k*nSpace,
                    eN_nDOF_trial_element = eN*nDOF_trial_element;
                  register double
                    u_test_dV[nDOF_trial_element],
                    u_grad_trial[nDOF_trial_element*nSpace],
                    jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace],
                    dV,x,y,z,xt,yt,zt;
                  register double u_at_qp,grad_u_at_qp[nSpace], f_at_qp;
                  ck.calculateMapping_element(eN,
                                  k,
                                  node_coord,////////////////////////////use updated mesh
                                  mesh_l2g,
                                  mesh_trial_ref,
                                  mesh_grad_trial_ref,
                                  jac,
                                  jacDet,
                                  jacInv,x,y,z);
                  ck.calculateMappingVelocity_element(eN,
                                      k,
                                      mesh_velocity_dof,
                                      mesh_l2g,
                                      mesh_trial_ref,
                                      xt,yt,zt);
                  dV = fabs(jacDet)*dV_ref[k];

                  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);

                  //get the solution
                  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u_at_qp);
                  //get the solution gradients
                  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_at_qp);
                  for (int j=0;j<nDOF_trial_element;j++)
                      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
          ///////////////////////////////////////////////////////////////////////////////////////////////

                  for (int i=0;i<nDOF_test_element;i++)
                  {

                      f_at_qp = get_rhs(x,y,z);
                      local_r[i] += (u_grad_trial[i*nSpace+0]*grad_u_at_qp[0]
                                    +u_grad_trial[i*nSpace+1]*grad_u_at_qp[1])*dV
                                    -f_at_qp*u_test_dV[i];
                  }

          ///////////////////////////////////////////////////////////////////////////////////////////////
              }//end of eN_k

          ///////////////////////////////////////////////////////////////////////////////////////////////
              for(int i=0;i<nDOF_test_element;i++)
              {
                  int eN_i=eN*nDOF_test_element+i;
                  r[u_offset+u_stride*u_l2g[eN_i]] += element_active*local_r[i];
//                  for (int j=0;j<nDOF_trial_element;j++)
//                  {
//                      int eN_i_j = eN_i*nDOF_trial_element+j;
////                      Cx[csrRowIndeces_CellLoops_rho[eN_i] + csrColumnOffsets_CellLoops_rho[eN_i_j]]+= local_c_x[i][j];//This is why to have NNZ nonzeros.
////                      Cy[csrRowIndeces_CellLoops_rho[eN_i] + csrColumnOffsets_CellLoops_rho[eN_i_j]]+= local_c_y[i][j];
////                      Cx_T[csrRowIndeces_CellLoops_rho[eN_i] + csrColumnOffsets_CellLoops_rho[eN_i_j]]+= local_c_x[j][i];
////                      Cy_T[csrRowIndeces_CellLoops_rho[eN_i] + csrColumnOffsets_CellLoops_rho[eN_i_j]]+= local_c_y[j][i];
//
////                      DD[csrRowIndeces_CellLoops_rho[eN_i] + csrColumnOffsets_CellLoops_rho[eN_i_j]]+= local_dd[j][i];
//
//                  }//j
              }
          }//end of eN
          if(USE_SBM>0)
          {
              for (int ebN_s=0;ebN_s < surrogate_boundaries.size();ebN_s++)
              {
                  register int ebN = surrogate_boundaries[ebN_s],
                    eN = elementBoundaryElementsArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                    eN_nDOF_trial_element = eN*nDOF_trial_element;
                  register double elementResidual_u[nDOF_test_element];
//                  if (ebN >= nElementBoundaries_owned) continue;/////for parallel; only loop one time across all processors
                  for (int i=0;i<nDOF_test_element;i++)
                  {
                      elementResidual_u[i]=0.0;
                  }
                  for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                  {
                      register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                        ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                        ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                      register double
                        u_ext=0.0,
                        bc_u_ext=0.0,
                        bc_v_ext=0.0,
                        grad_u_ext[nSpace],
                        jac_ext[nSpace*nSpace],
                        jacDet_ext,
                        jacInv_ext[nSpace*nSpace],
                        boundaryJac[nSpace*(nSpace-1)],
                        metricTensor[(nSpace-1)*(nSpace-1)],
                        metricTensorDetSqrt,
                        dS,
                        u_test_dS[nDOF_test_element],
                        u_grad_test_dS[nDOF_trial_element*nSpace],
                        u_grad_trial_trace[nDOF_trial_element*nSpace],
                        normal[2],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                        G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty,penalty;
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
                      dS = metricTensorDetSqrt*dS_ref[kb];
                      //get the metric tensor
                      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                      //compute shape and solution information
                      //shape
                      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                      //solution and gradients
                      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);

                      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
                      //precalculate test function products with integration weights
                      for (int j=0;j<nDOF_trial_element;j++)
                      {
                          u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                          for (int I=0;I<nSpace;I++)
                            u_grad_test_dS[j*nSpace+I] = u_grad_trial_trace[j*nSpace+I]*dS;//cek hack, using trial
                      }

                      double dist = 0.0;
                      double distance[2], P_normal[2], P_tangent[2], normal_Omega[2]; // distance vector, normal and tangent of the physical boundary

                      if(use_ball_as_particle==1)
                      {
                          get_distance_to_ball(nParticles,ball_center,ball_radius,
                                               x_ext,y_ext,z_ext,
                                               dist);
                          get_normal_to_ith_ball(nParticles,ball_center,ball_radius,
                                                 surrogate_boundary_particle[ebN_s],
                                                 x_ext,y_ext,z_ext,
                                                 P_normal[0],P_normal[1]);
                          get_velocity_to_ith_ball(nParticles,ball_center,ball_radius,
                                                   ball_velocity,ball_angular_velocity,
                                                   surrogate_boundary_particle[ebN_s],
                                                   x_ext-dist*P_normal[0],//corresponding point on the boundary of the particle
                                                   y_ext-dist*P_normal[1],
                                                   0.0,//z_ext,
                                                   bc_u_ext,bc_v_ext);
                      }
                      else
                      {
//                          dist = ebq_global_phi_solid[ebN_kb];
//                          P_normal[0] = ebq_global_grad_phi_solid[ebN_kb*nSpace+0];
//                          P_normal[1] = ebq_global_grad_phi_solid[ebN_kb*nSpace+1];
//                          bc_u_ext = ebq_particle_velocity_solid [ebN_kb*nSpace+0];
//                          bc_v_ext = ebq_particle_velocity_solid [ebN_kb*nSpace+1];

                      }
                      ck.calculateGScale(G,normal,h_penalty);
                      //
                      //update the element and global residual storage
                      //

                      assert(dist>0.0);
                      assert(h_penalty>0.0);

                      distance[0] = -P_normal[0]*dist;
                      distance[1] = -P_normal[1]*dist;
                      P_tangent[0] = -P_normal[1];
                      P_tangent[1] = P_normal[0];
                      normal_Omega[0] = -P_normal[0];
                      normal_Omega[1] = -P_normal[1];

                      const double C_adim = C_sbm/h_penalty;
                      double beta_adim = beta_sbm/h_penalty;

                      const double grad_u_d = get_dot_product(distance,grad_u_ext);

                      double res;
                      const double u_m_uD = u_ext - bc_u_ext;
                      const double grad_u_t = get_dot_product(P_tangent,grad_u_ext);

                      for (int i=0;i<nDOF_test_element;i++)
                        {
                          int eN_i = eN*nDOF_test_element+i;

                          int GlobPos_u = u_offset+u_stride*u_l2g[eN_i];

                          double phi_i = u_test_dS[i];
                          double Gxphi_i = u_grad_test_dS[i*nSpace+0];
                          double Gyphi_i = u_grad_test_dS[i*nSpace+1];
                          double *grad_phi_i = &u_grad_test_dS[i*nSpace+0];
                          const double grad_phi_i_dot_d =  get_dot_product(distance,grad_phi_i);
                          const double grad_phi_i_dot_t =  get_dot_product(P_tangent,grad_phi_i);

                          // (1)
                          r[GlobPos_u] += C_adim*(phi_i+grad_phi_i_dot_d)*(u_m_uD+grad_u_d);

                          // (2)
                          r[GlobPos_u] -= phi_i* get_dot_product(grad_u_ext,normal);////YY: not normal_Omega

                          // (3)
                          r[GlobPos_u] -= get_dot_product(grad_phi_i,normal)*(u_m_uD+grad_u_d);////YY: not normal_Omega
                        }//i

                    }//kb

              }//ebN_s
              //
            }
          ///////////////////////////////////////////////////////////////////////////////////////////////
          return;
}

      void calculateMassMatrix(//element
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
                   int* u_l2g,
                   double* elementDiameter,
                   int degree_polynomial,
                   double* u_dof,
                   double* velocity,
                   double* q_m_betaBDF,
                   double* cfl,
                   double* q_numDiff_u_last,
                   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                   double* globalJacobian,
                   int nExteriorElementBoundaries_global,
                   int* exteriorElementBoundariesArray,
                   int* elementBoundaryElementsArray,
                   int* elementBoundaryLocalElementBoundariesArray,
                   double* ebqe_velocity_ext,
                   int* isDOFBoundary_u,
                   double* ebqe_rd_u_ext,
                   double* ebqe_bc_u_ext,
                   int* csrColumnOffsets_eb_u_u,
                   int PURE_BDF,
                   int LUMPED_MASS_MATRIX)
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
        register double u=0.0,
          grad_u[nSpace],
          m=0.0,dm=0.0,
          H=0.0,dH[nSpace],
          f[nSpace],df[nSpace],//MOVING_MESH
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
        //get the solution gradients
        ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
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
                     u,
                     grad_u,
                     m,
                     dm,
                     H,
                     dH);
        //
        //moving mesh
        //
        double mesh_velocity[3];
        mesh_velocity[0] = xt;
        mesh_velocity[1] = yt;
        mesh_velocity[2] = zt;
        for (int I=0;I<nSpace;I++)
          {
            f[I] = -MOVING_DOMAIN*m*mesh_velocity[I];
            df[I] = -MOVING_DOMAIN*dm*mesh_velocity[I];
          }
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
            //int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
            //Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[eN_k_i_nSpace]);
            register int i_nSpace = i*nSpace;
            Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[i_nSpace]) + MOVING_DOMAIN*ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);

          }
        //calculate the Jacobian of strong residual
        for (int j=0;j<nDOF_trial_element;j++)
          {
            //int eN_k_j=eN_k*nDOF_trial_element+j;
            //int eN_k_j_nSpace = eN_k_j*nSpace;
            int j_nSpace = j*nSpace;
            dpdeResidual_u_u[j]=ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
              ck.HamiltonianJacobian_strong(dH,&u_grad_trial[j_nSpace]) +
              MOVING_DOMAIN*ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);

          }
        //tau and tau*Res
        double subgridErrorVelocity[nSpace];
        for (int I=0;I<nSpace;I++)
          subgridErrorVelocity[I] = dH[I] - MOVING_DOMAIN*df[I];

        calculateSubgridError_tau(elementDiameter[eN],
                      dm_t,
                      subgridErrorVelocity,//dH,
                      cfl[eN_k],
                      tau0);

        calculateSubgridError_tau(Ct_sge,
                      G,
                      dm_t,
                      subgridErrorVelocity,//dH,
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
            if (LUMPED_MASS_MATRIX==1)
              {
                if (i==j)
                  elementJacobian_u_u[i][j] += u_test_dV[i];
              }
            else
              {
                elementJacobian_u_u[i][j] +=
                  dt*ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]);
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
      }//computeMassMatrix

    };//Poisson

  inline Poisson_base* newPoisson(int nSpaceIn,
                            int nQuadraturePoints_elementIn,
                            int nDOF_mesh_trial_elementIn,
                            int nDOF_trial_elementIn,
                            int nDOF_test_elementIn,
                            int nQuadraturePoints_elementBoundaryIn,
                            int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<Poisson_base,Poisson,CompKernel>(nSpaceIn,
                                                                                   nQuadraturePoints_elementIn,
                                                                                   nDOF_mesh_trial_elementIn,
                                                                                   nDOF_trial_elementIn,
                                                                                   nDOF_test_elementIn,
                                                                                   nQuadraturePoints_elementBoundaryIn,
                                                                                   CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<Poisson_base,Poisson,CompKernel>(nSpaceIn,
                                                                                 nQuadraturePoints_elementIn,
                                                                                 nDOF_mesh_trial_elementIn,
                                                                                 nDOF_trial_elementIn,
                                                                                 nDOF_test_elementIn,
                                                                                 nQuadraturePoints_elementBoundaryIn,
                                                                                 CompKernelFlag);
    /* return proteus::chooseAndAllocateDiscretization<Poisson_base,Poisson>(nSpaceIn, */
    /*                                                                         nQuadraturePoints_elementIn, */
    /*                                                                         nDOF_mesh_trial_elementIn, */
    /*                                                                         nDOF_trial_elementIn, */
    /*                                                                         nDOF_test_elementIn, */
    /*                                                                         nQuadraturePoints_elementBoundaryIn, */
    /*                                                                         CompKernelFlag); */
  }
}//proteus
#endif
