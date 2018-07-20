#ifndef Poisson_H1
#define Poisson_H1
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <set>
#include <cstring>
#include "CompKernel.h"
#include "ModelFactory.h"


#define ASSERT_WITH_MSG(cond, msg) do \
        { if (!(cond)) { std::ostringstream str; str << msg; std::cerr << str.str();        std::abort(); } \
        } while(0)

#define USE_H_PERP 0


namespace proteus
{
  class Poisson_base
  {
    //The base class defining the interface
  public:
    virtual ~Poisson_base(){}
    virtual void calculateResidual(//element
                                           int rank,
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
                                           int nNodes_owned,
                                           int nElementBoundaries_owned,
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
                                    int u_ndofs,
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
                                    int nNodes_owned,
                                    int nElementBoundaries_owned,
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
        ball_radius[0]=0.2;
        ball_velocity[0]=0.0;
        ball_velocity[1]=0.0;
        ball_velocity[2]=0.0;
        ball_angular_velocity[0]=0.0;
        ball_angular_velocity[1]=0.0;
        ball_angular_velocity[2]=0.0;
      }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
                 int u_ndofs,
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
                 int nNodes_owned,
                 int nElementBoundaries_owned,
                 int USE_SBM)
      {

          double* node_coord=mesh_dof;
          std::vector<int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;
          ///////////////////////////////////////////////////////////////////////////////////////////////
          for(int eN=0;eN<nElements_global;eN++)
          {
          ///////////////////////////////////////////////////////////////////////////////////////////////
              register double  local_dd[nDOF_test_element][nDOF_trial_element];

              int element_active=1;//use 1 since by default it is assembled over all elements

              for (int i=0;i<nDOF_test_element;i++)
              {
                  for (int j=0;j<nDOF_trial_element;j++)
                  {
                          local_dd[i][j]=0.0;
                  }
              }
              if(USE_SBM>0)///////YYYYYYY: has to update every time since isActive is assin to be 1 in getResidual
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
                      element_active=0;
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

                      ///////For parallel. If no nodes of this edge are owned by this processor,
                      ///////there is no contribution of this edge integral to Jacobian and residual.
                      ///////As the surrogate boundary integral, maybe this edge is owned by this process. So add assert in else.

                      if(mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+1)%3]<nNodes_owned || mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+2)%3]<nNodes_owned)
                      {
                          int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                          surrogate_boundaries.push_back(ebN);
                          //now find which element neighbor this element is
                          if (eN == elementBoundaryElementsArray[ebN*2+0])/////////YY: should be ebN
                              surrogate_boundary_elements.push_back(1);
                          else if(eN == elementBoundaryElementsArray[ebN*2+1])
                              surrogate_boundary_elements.push_back(0);
                          else
                              assert(0);

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
                      }else{
                          int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                          assert(ebN>=nElementBoundaries_owned);
                      }
                  }
                  else if (pos_counter == 3)
                  {
                      element_active=1;
                  }
                  else
                  {
                      element_active=0;
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
                    dV,x,y,z;
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
//                  ck.calculateMappingVelocity_element(eN,
//                                      k,
//                                      mesh_velocity_dof,
//                                      mesh_l2g,
//                                      mesh_trial_ref,
//                                      xt,yt,zt);
                  dV = fabs(jacDet)*dV_ref[k];

                  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
          ///////////////////////////////////////////////////////////////////////////////////////////////

                  for (int i=0;i<nDOF_test_element;i++)
                  {
                      const double grad_phi_i[2] = {u_grad_trial[i*nSpace+0],
                                                    u_grad_trial[i*nSpace+1]};
                      for (int j=0;j<nDOF_trial_element;j++)
                      {
                          const double grad_phi_j[2] = {u_grad_trial[j*nSpace+0],
                                                        u_grad_trial[j*nSpace+1]};
                          local_dd[i][j] += (grad_phi_i[0]*grad_phi_j[0]
                                            +grad_phi_i[1]*grad_phi_j[1])*dV;
                      }

                  }

          ///////////////////////////////////////////////////////////////////////////////////////////////
              }//end of eN_k
          ///////////////////////////////////////////////////////////////////////////////////////////////
              if(element_active)
                  for(int i=0;i<nDOF_test_element;i++)
                  {
                      int eN_i=eN*nDOF_test_element+i;
                      for (int j=0;j<nDOF_trial_element;j++)
                      {
                          int eN_i_j = eN_i*nDOF_trial_element+j;
                          globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += local_dd[j][i];

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
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////get h_perp
                  double h_perp = 1.0e10;
                  if(USE_H_PERP)
                  {
                      for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                      {
                          register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                                  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                                  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                          register double
                          jac_ext[nSpace*nSpace]={0.0},
                          jacDet_ext,
                          jacInv_ext[nSpace*nSpace]={0.0},
                          boundaryJac[nSpace*(nSpace-1)]={0.0},
                          metricTensor[(nSpace-1)*(nSpace-1)]={0.0},
                          metricTensorDetSqrt,
                          dS,
                          normal[2],
                          x_ext,y_ext,z_ext,
                          G[nSpace*nSpace],G_dd_G,tr_G,h_penalty;

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
                          ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                          ck.calculateGScale(G,normal,h_penalty);////////YY: should use elementDiameter????
                          //
                          //update h_perp
                          //
                          double dist=1e10;
                          if(use_ball_as_particle==1)
                          {
                              get_distance_to_ball(nParticles,ball_center,ball_radius,
                                      x_ext,y_ext,z_ext,
                                      dist);
                          }else{
                              ASSERT_WITH_MSG(0, "use ball as particle");
                          }
                          assert(dist>0.0);
                          assert(h_penalty>0.0);

                          if (h_penalty < dist)
                          {
                              h_penalty = dist;
                          }
                          if(h_penalty<h_perp)
                          {
                              h_perp = h_penalty;
                          }
                      }//end-kb
                  }//end of h_perp
                  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                  {
                      register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                              ebN_kb_nSpace = ebN_kb*nSpace,
                              ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                              ebN_local_kb_nSpace = ebN_local_kb*nSpace;

                      register double u_ext=0.0,
                              bc_u_ext=0.0,
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
                              x_ext,y_ext,z_ext,
                              u_grad_test_dS[nDOF_trial_element*nSpace],
                              G[nSpace*nSpace],G_dd_G,tr_G,h_penalty;
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
//                      ck.calculateMappingVelocity_elementBoundary(eN,
//                              ebN_local,
//                              kb,
//                              ebN_local_kb,
//                              mesh_velocity_dof,
//                              mesh_l2g,
//                              mesh_trial_trace_ref,
//                              xt_ext,yt_ext,zt_ext,
//                              normal,
//                              boundaryJac,
//                              metricTensor,
//                              integralScaling);
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
                      ck.calculateGScale(G,normal,h_penalty);
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
                          get_u_on_ith_ball(nParticles,ball_center,ball_radius,
                                  ball_velocity, ball_angular_velocity,
                                  surrogate_boundary_particle[ebN_s],
                                  x_ext-dist*P_normal[0],
                                  y_ext-dist*P_normal[1],
                                  0.0,//z_ext,
                                  bc_u_ext);
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
                      assert(dist>0.0);
                      assert(h_penalty>0.0);
                      if (h_penalty < dist)
                          h_penalty = dist;
                      double C_adim = C_sbm/h_penalty;
                      double beta_adim = beta_sbm/h_penalty;
                      if(USE_H_PERP)
                      {
                          C_adim = C_sbm/h_perp;
                          beta_adim = beta_sbm/h_perp;
                      }
                      for (int i=0;i<nDOF_test_element;i++)
                      {
                          const int eN_i = eN*nDOF_test_element+i;/////use eN to get index of dofs

                          double phi_i = u_test_trace_ref[ebN_local_kb*nDOF_test_element+i];
                          double grad_phi_i[2] = {u_grad_trial_trace[i*nSpace+0],
                                                  u_grad_trial_trace[i*nSpace+1]};

                          const double grad_phi_i_dot_d =  get_dot_product(distance,grad_phi_i);
                          const double grad_phi_i_dot_t =  get_dot_product(P_tangent,grad_phi_i);

                          for (int j=0;j<nDOF_trial_element;j++)
                          {
//                              assert(surrogate_boundary_elements[ebN_s]==0);//////YY: some values are 1
                              register int ebN_i_j = (ebN*4
                                                      +surrogate_boundary_elements[ebN_s]*3/////////YY: bug
                                                      )*nDOF_test_X_trial_element
                                                      + i*nDOF_trial_element + j;

                              double phi_j = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j];
                              double grad_phi_j[2] = {u_grad_trial_trace[j*nSpace+0],
                                                      u_grad_trial_trace[j*nSpace+1]};

                              const double grad_phi_j_dot_d =  get_dot_product(distance,grad_phi_j);
                              const double grad_phi_j_dot_t =  get_dot_product(P_tangent,grad_phi_j);

                              // (1)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] +=
                                      C_adim*(phi_i+grad_phi_i_dot_d)*(phi_j+grad_phi_j_dot_d)*dS;

                              // (2)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] -=
                                      phi_i * get_dot_product(grad_phi_j,normal)*dS;/////YY: not normal_Omega

                              // (3)
                              globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] -=
                                      get_dot_product(grad_phi_i,normal)*(phi_j+grad_phi_j_dot_d)*dS;/////YY: not normal_Omega
                          }//j
                      }//i
                  }//kb
              }//ebN_s
          }
          ///////////////////////////////////////////////////////////////////////////////////////////////
          return;
      }//computeJacobian

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
      void get_u_on_ith_ball(int n_balls,double* ball_center, double* ball_radius,
                                    double* ball_velocity, double* ball_angular_velocity,
                                    int I,
                                    double x, double y, double z,
                                    double& u)
      {
          u = 0.0;
      }
      double get_u(const double x, const double y)
      {
          const double R2 = 0.2*0.2;
          const double r2 = x*x+y*y;
          if (r2<= R2)
              return 0.5*(R2-r2);
          else
              return 0.0;
      }
      void get_gradient_u(const double x, const double y, double& ux, double& uy)
      {
          const double R2 = 0.2*0.2;
          const double r2 = x*x+y*y;
          if (r2<= R2)
          {
              ux = -x;
              uy = -y;
          }
          else
          {
              ux = 0.0;
              uy = 0.0;
          }
      }
      void test_element(int eN, int ebN, double* node_xyz, int* ele_nodes)
      {
          const int i0=ele_nodes[3*eN+0], i1=ele_nodes[3*eN+1], i2=ele_nodes[3*eN+2];
          const double *xyz0=&node_xyz[3*i0];
          const double *xyz1=&node_xyz[3*i1];
          const double *xyz2=&node_xyz[3*i2];
          std::cout<<"element= "<<eN<<"\t"
                  <<"("<<xyz0[0]<<","<<xyz0[1]<<"),"
                  <<"("<<xyz1[0]<<","<<xyz1[1]<<"),"
                  <<"("<<xyz2[0]<<","<<xyz2[1]<<")"
                  <<"edges="<<ebN<<std::endl;
      }
      void calculateResidual(//element
              int rank,
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
              int nNodes_owned,
              int nElementBoundaries_owned,
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
          std::vector<long int> surrogate_boundaries, surrogate_boundary_elements, surrogate_boundary_particle;
          double error_u_L2=0.0,error_u_H1=0.0;

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
              register double  local_r[nDOF_test_element];
              double local_error_u_L2=0.0,local_error_u_H1=0.0;

              int element_active=1;//use 1 since by default it is assembled over all elements

              for (int i=0;i<nDOF_test_element;i++)
              {
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
                      element_active=0;
                      int opp_node=-1;
                      for (int I=0;I<nDOF_mesh_trial_element;I++)
                      {
                          quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + I]] = 2.0;////////YY: for test
                          if (_distance[I] < 0)
                          {
                              opp_node = I;
                              quantDOFs[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + I]] = 1.0;////////YY: for test
                          }
                      }
                      assert(opp_node >=0);
                      assert(opp_node <nDOF_mesh_trial_element);
                      //For parallel. Two reasons:
                      //if none of nodes of this edge is owned by this processor,
                      //1. The surrogate_boundary_elements corresponding to this edge is -1, which gives 0 JacDet and infty h_penalty.
                      //2. there is no contribution of the integral over this edge to Jacobian and residual.
                      if(mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+1)%3]<nNodes_owned || mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+2)%3]<nNodes_owned)
                      {
                          int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                          surrogate_boundaries.push_back(ebN);
                          //now find which element neighbor this element is
                          if (eN == elementBoundaryElementsArray[ebN*2+0])/////////YY: sould be ebN
                              surrogate_boundary_elements.push_back(1);
                          else if(eN == elementBoundaryElementsArray[ebN*2+1])/////////YY: sould be ebN
                              surrogate_boundary_elements.push_back(0);
                          else
                              assert(0);
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
                      }else{
                          //If the integral over the surrogate boundary is needed, we have to make sure all edges are in surrogate_boundaries,
                          //which is based on the assumption that if none of its nodes is owned by the processor, then the edge is not owned
                          //by this processor. This assert is used to make sure this.
                          int ebN = elementBoundariesArray[eN*nDOF_mesh_trial_element+opp_node];//only works for simplices
                          ASSERT_WITH_MSG(ebN>=nElementBoundaries_owned,
                                  ebN<<","<<mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+1)%3]<<","<<mesh_l2g[eN*nDOF_mesh_trial_element+(opp_node+2)%3]
                                  <<nElementBoundaries_owned<<","<<nNodes_owned);
                      }
                  }
                  else if (pos_counter == 3)
                  {
                      element_active=1;
                      for (int i=0;i<nDOF_test_element;i++)
                      {
                          isActiveDOF[offset_u+stride_u*u_l2g[eN*nDOF_trial_element + i]]=1.0;
                      }
                  }
                  else
                  {
                      element_active=0;
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
                      const double phi_i = u_test_ref[k*nDOF_trial_element+i];
                      const double grad_phi_i[2] = {u_grad_trial[i*nSpace+0],
                                                    u_grad_trial[i*nSpace+1]};

                      f_at_qp = get_rhs(x,y,z);
                      local_r[i] += (grad_phi_i[0]*grad_u_at_qp[0]
                                    +grad_phi_i[1]*grad_u_at_qp[1])*dV
                                    -f_at_qp*phi_i*dV;
                  }
                  ///////////////////////////////////////////////////////////////////////////////////////////////
                  local_error_u_L2 += std::pow(get_u(x,y)-u_at_qp, 2.0)*dV;
                  double ux,uy;
                  get_gradient_u(x,y,ux,uy);

                  local_error_u_H1 += std::pow(ux-grad_u_at_qp[0], 2.0)*dV+std::pow(uy-grad_u_at_qp[1], 2.0)*dV;

          ///////////////////////////////////////////////////////////////////////////////////////////////
              }//end of eN_k

              ///////////////////////////////////////////////////////////////////////////////////////////////
              if(element_active)
              {
                  for(int i=0;i<nDOF_test_element;i++)
                  {
                      int eN_i=eN*nDOF_test_element+i;
                      r[u_offset+u_stride*u_l2g[eN_i]] += local_r[i];
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
                      quantDOFs[u_offset+u_stride*u_l2g[eN_i]]=rank;
                  }
                  //////////////////////////////////////////////////////////////YY: In theory, the error is over active elements
                  error_u_L2 += local_error_u_L2;
                  error_u_H1 += local_error_u_H1;
              }
              for(int i=0;i<nDOF_test_element;i++)
              {
                  int eN_i=eN*nDOF_test_element+i;
                  quantDOFs[u_offset+u_stride*u_l2g[eN_i]]=rank;
              }
//              ///////////////////////////////////////////////////////////////////////////////////////////////
//              error_u_L2 += local_error_u_L2;
//              error_u_H1 += local_error_u_H1;

          }//end of eN
          std::cout<<"YY-pdb L^2 error of u: "
                  <<std::sqrt(error_u_L2)<<" "
                  <<std::endl;
          std::cout<<"YY-pdb H^1-seminorm error of u: "
                  <<std::sqrt(error_u_H1)<<" "
                  <<std::endl;
          if(USE_SBM>0)
          {
              if (rank==0)
              {

//                  test_element(57,0,mesh_dof,mesh_l2g);
                  for (int ebN_s=0;ebN_s < surrogate_boundaries.size();ebN_s++)
                  {
                      register int ebN = surrogate_boundaries[ebN_s],
                              eN = elementBoundaryElementsArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                              ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+surrogate_boundary_elements[ebN_s]];
                      if (eN<0)
                      {
                          std::cout<<rank<<","
                                  <<ebN_s<<","
                                  <<ebN<<","
                                  <<eN<<","<<std::endl;
                          test_element(eN,ebN_local,mesh_dof,mesh_l2g);
                      }
                  }
              }
              for (int ebN_s=0;ebN_s < surrogate_boundaries.size();ebN_s++)
              {
                  register int ebN = surrogate_boundaries[ebN_s],
                          eN = elementBoundaryElementsArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                          ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+surrogate_boundary_elements[ebN_s]],
                          eN_nDOF_trial_element = eN*nDOF_trial_element;
                  register double elementResidual_u[nDOF_test_element];
                  // YY: This is wrong for assembling since the associated dofs may by owned by this processor
                  // and the edge has contribution to the dof even it is not owned.
                  //                  if (ebN >= nElementBoundaries_owned) continue;/////for parallel;

                  for (int i=0;i<nDOF_test_element;i++)
                  {
                      elementResidual_u[i]=0.0;
                  }

                  ///////////////////////////////////////////////////////////////////////////////////////////////////////get h_perp
                  double h_perp = 1.0e10;
                  if(USE_H_PERP)
                  {
                      for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                      {
                          register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                                  ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                                  ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                          register double
                          jac_ext[nSpace*nSpace]={0.0},
                          jacDet_ext,
                          jacInv_ext[nSpace*nSpace]={0.0},
                          boundaryJac[nSpace*(nSpace-1)]={0.0},
                          metricTensor[(nSpace-1)*(nSpace-1)]={0.0},
                          metricTensorDetSqrt,
                          dS,
                          normal[2],
                          x_ext,y_ext,z_ext,
                          G[nSpace*nSpace],G_dd_G,tr_G,h_penalty;

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
                          ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                          ck.calculateGScale(G,normal,h_penalty);////////YY: should use elementDiameter????
                                  //
                                  //update h_perp
                                  //
                          double dist=1e10;
                          if(use_ball_as_particle==1)
                          {
                              get_distance_to_ball(nParticles,ball_center,ball_radius,
                                      x_ext,y_ext,z_ext,
                                      dist);
                          }else{
                              ASSERT_WITH_MSG(0, "use ball as particle");
                          }
                          assert(dist>0.0);
                          assert(h_penalty>0.0);

                          if (h_penalty < dist)
                          {
                              h_penalty = dist;
                          }
                          if(h_penalty<h_perp)
                          {
                              h_perp = h_penalty;
                          }
                      }//end-kb
                  }//end of h_perp
                  ///////////////////////////////////////////////////////////////////////////////////////////////////////
                  for (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                  {
                      register int ebN_kb = ebN*nQuadraturePoints_elementBoundary+kb,
                        ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                        ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                      register double
                        u_ext=0.0,
                        bc_u_ext=0.0,
                        grad_u_ext[nSpace],
                        jac_ext[nSpace*nSpace]={0.0},
                        jacDet_ext=M_PI,
                        jacInv_ext[nSpace*nSpace]={0.0},
                        boundaryJac[nSpace*(nSpace-1)]={0.0},
                        metricTensor[(nSpace-1)*(nSpace-1)]={0.0},
                        metricTensorDetSqrt,
                        dS,
                        u_test_dS[nDOF_test_element]={0.0},
                        u_grad_test_dS[nDOF_trial_element*nSpace]={0.0},
                        u_grad_trial_trace[nDOF_trial_element*nSpace]={0.0},
                        normal[2],x_ext,y_ext,z_ext,
                        G[nSpace*nSpace],G_dd_G,tr_G,h_penalty;
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
//                      ck.calculateMappingVelocity_elementBoundary(eN,
//                                                                  ebN_local,
//                                                                  kb,
//                                                                  ebN_local_kb,
//                                                                  mesh_velocity_dof,
//                                                                  mesh_l2g,
//                                                                  mesh_trial_trace_ref,
//                                                                  xt_ext,yt_ext,zt_ext,
//                                                                  normal,
//                                                                  boundaryJac,
//                                                                  metricTensor,
//                                                                  integralScaling);
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
                          get_u_on_ith_ball(nParticles,ball_center,ball_radius,
                                                   ball_velocity,ball_angular_velocity,
                                                   surrogate_boundary_particle[ebN_s],
                                                   x_ext-dist*P_normal[0],//corresponding point on the boundary of the particle
                                                   y_ext-dist*P_normal[1],
                                                   0.0,//z_ext,
                                                   bc_u_ext);
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

                      ASSERT_WITH_MSG(dist>0.0, dist);
                      ASSERT_WITH_MSG(h_penalty>0.0, rank<<","<<h_penalty<<","<<metricTensorDetSqrt<<","<<jacDet_ext<<","<<ebN_s;
                      /*test_element(eN,ebN_local,mesh_dof,mesh_l2g)*/
                      );
                      if (h_penalty < dist)
                      {
                          h_penalty = dist;
                      }
                      distance[0] = -P_normal[0]*dist;
                      distance[1] = -P_normal[1]*dist;
                      P_tangent[0] = -P_normal[1];
                      P_tangent[1] = P_normal[0];
                      normal_Omega[0] = -P_normal[0];
                      normal_Omega[1] = -P_normal[1];

                      double C_adim = C_sbm/h_penalty;
                      double beta_adim = beta_sbm/h_penalty;
                      if(USE_H_PERP)////////use h_perp instead
                      {
                          C_adim = C_sbm/h_perp;
                          beta_adim = beta_sbm/h_perp;
                      }

                      const double grad_u_d = get_dot_product(distance,grad_u_ext);

                      double res;
                      const double u_m_uD = u_ext - bc_u_ext;
                      const double grad_u_t = get_dot_product(P_tangent,grad_u_ext);
                      if(0)
                      {
                          double x1 = mesh_dof[3*mesh_l2g[eN*3+0]+0], y1 = mesh_dof[3*mesh_l2g[eN*3+0]+1];
                          double x2 = mesh_dof[3*mesh_l2g[eN*3+1]+0], y2 = mesh_dof[3*mesh_l2g[eN*3+1]+1];
                          double x3 = mesh_dof[3*mesh_l2g[eN*3+2]+0], y3 = mesh_dof[3*mesh_l2g[eN*3+2]+1];
                          std::cout<<"yyPDB-Surrogate bc: ";
                          if(ebN_local==0)
                          {
                              std::cout<<x2<<"\t"
                                       <<y2<<"\t"
                                       <<x3<<"\t"
                                       <<y3<<"\t";
                          }else if(ebN_local==1){

                              std::cout<<x3<<"\t"
                                       <<y3<<"\t"
                                       <<x1<<"\t"
                                       <<y1<<"\t";
                          }else if(ebN_local==2){

                              std::cout<<x1<<"\t"
                                       <<y1<<"\t"
                                       <<x2<<"\t"
                                       <<y2<<"\t";
                          }
                          std::cout<<distance[0]<<"\t"
                                  <<distance[1]<<"\t"
                                  <<P_normal[0]<<"\t"
                                  <<P_normal[1]<<"\t"
                                  <<normal[0]<<"\t"
                                  <<normal[1]<<"\t"
                                  <<x_ext<<"\t"
                                  <<y_ext<<"\t"
                                  <<"\n";
                      }
                      for (int i=0;i<nDOF_test_element;i++)
                      {
                          const int eN_i = eN*nDOF_test_element+i;/////use eN to get index of dofs
                          const int GlobPos_u = u_offset+u_stride*u_l2g[eN_i];

                          double phi_i = u_test_trace_ref[ebN_local_kb*nDOF_test_element+i];
                          double grad_phi_i[2] = {u_grad_trial_trace[i*nSpace+0],
                                                  u_grad_trial_trace[i*nSpace+1]};

                          const double grad_phi_i_dot_d =  get_dot_product(distance,grad_phi_i);
                          const double grad_phi_i_dot_t =  get_dot_product(P_tangent,grad_phi_i);

                          // (1)
                          r[GlobPos_u] += C_adim*(phi_i+grad_phi_i_dot_d)*(u_m_uD+grad_u_d)*dS;

                          // (2)
                          r[GlobPos_u] -= phi_i* get_dot_product(grad_u_ext,normal)*dS;////YY: for consistency; not normal_Omega

                          // (3)
                          r[GlobPos_u] -= get_dot_product(grad_phi_i,normal)*(u_m_uD+grad_u_d)*dS;////YY: for consistency; not normal_Omega
                      }//i

                    }//kb

              }//ebN_s
            }//surrogate bc
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
          assert(0);
      }//computeMassMatrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
