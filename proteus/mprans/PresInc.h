#ifndef PRESINC_H
#define PRESINC_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
  class cppPresInc_base
  {
  public:
    virtual ~cppPresInc_base(){}
    virtual void calculateResidual(//element
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
                                   int* isDOFBoundary,
                                   int* isFluxBoundary,
                                   int* u_l2g,
                                   double* u_dof,
                                   double alphaBDF,
                                   double* q_vf,
                                   double* q_divU,
                                   double* q_vs,
                                   double* q_vos,
                                   double rho_s,
                                   double* q_rho_f,
                                   double rho_s_min,
                                   double rho_f_min,
                                   double* ebqe_vf,
                                   double* ebqe_vs,
                                   double* ebqe_vos,
                                   double* ebqe_rho_f,
                                   double* q_u,
                                   double* q_grad_u,
                                   double* ebqe_u,
                                   double* ebqe_grad_u,
                                   double* ebqe_bc_u_ext,
                                   double* ebqe_adv_flux,
                                   double* ebqe_diff_flux,
                                   double* bc_adv_flux,
                                   double* bc_diff_flux,
                                   int offset_u,
                                   int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   int INTEGRATE_BY_PARTS_DIV_U)=0;
    virtual void calculateJacobian(//element
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
                                   int* isDOFBoundary,
                                   int* isFluxBoundary,
                                   int* u_l2g,
                                   double* u_dof,
                                   double alphaBDF,
                                   double* q_vf,
                                   double* q_vs,
                                   double* q_vos,
                                   double rho_s,
                                   double* q_rho_f,
                                   double rho_s_min,
                                   double rho_f_min,
                                   double* ebqe_vf,
                                   double* ebqe_vs,
                                   double* ebqe_vos,
                                   double* ebqe_rho_f,
                                   int* csrRowIndeces_u_u,
                                   int* csrColumnOffsets_u_u,
                                   double* globalJacobian,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   int* csrColumnOffsets_eb_u_u)=0;
  };

  template<class CompKernelType,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nQuadraturePoints_elementBoundary>
  class cppPresInc : public cppPresInc_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    cppPresInc():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}
    inline
      void evaluateCoefficients(const double& alphaBDF,
                                const double vf[nSpace],
                                const double vs[nSpace],
                                const double& vos,
                                const double& rhos_min,
                                const double& rhof_min,
                                double f[nSpace],
                                double& a)
    {
      for (int I=0;I<nSpace;I++)
        f[I] = (1.0-vos)*vf[I] + vos*vs[I];
      a = (1.0-vos)/(rhof_min*alphaBDF) + vos/(rhos_min*alphaBDF);
    }

    inline
      void exteriorNumericalAdvectiveFlux(const int& isFluxBoundary,
                                          const double& bc_flux,
                                          const double n[nSpace],
                                          const double f[nSpace],
                                          double& flux)
    {
      if (isFluxBoundary == 1)
        flux = bc_flux;
      else
        {
          flux = 0.0;
          for (int I=0; I < nSpace; I++)
            flux += n[I]*f[I];
        }
    }

    inline
    void exteriorNumericalDiffusiveFlux(const int& isDOFBoundary,
                                        const int& isFluxBoundary,
                                        const double n[nSpace],
                                        const double& a,
                                        const double grad_potential[nSpace],
                                        const double& u,
                                        const double& bc_u,
                                        const double& bc_flux,
                                        const double& penalty,
                                        double& flux)
    {
      if(isFluxBoundary == 1)
        {
          flux = bc_flux;
        }
      else if(isDOFBoundary == 1)
        {
          flux = 0.0;
          for(int I=0;I<nSpace;I++)
            flux-= a*grad_potential[I]*n[I];
          flux += a*penalty*(u-bc_u);
        }
      else
        {
          std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
          flux = 0.0;
        }
    }

    inline
    double ExteriorNumericalDiffusiveFluxJacobian(const int& isDOFBoundary,
                                                  const int& isFluxBoundary,
                                                  const double n[nSpace],
                                                  const double& a,
                                                  const double& v,
                                                  const double grad_v[nSpace],
                                                  const double& penalty)
    {
      double tmp=0.0;
      if(isFluxBoundary==0 && isDOFBoundary==1)
        {
          for(int I=0;I<nSpace;I++)
            tmp -= a*grad_v[I]*n[I];
          tmp +=a*penalty*v;
        }
      return tmp;
    }

    inline void calculateElementResidual(//element
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
                                         double* u_dof,
                                         double alphaBDF,
                                         double* q_vf,
                                         double* q_divU,
                                         double* q_vs,
                                         double* q_vos,
                                         double rho_s,
                                         double* q_rho_f,
                                         double rho_s_min,
                                         double rho_f_min,
                                         double* ebqe_vf,
                                         double* ebqe_vs,
                                         double* ebqe_vos,
                                         double* ebqe_rho_f,
                                         double* q_u,
                                         double* q_grad_u,
                                         double* ebqe_u,
                                         double* ebqe_grad_u,
                                         int offset_u,
                                         int stride_u,
                                         double* elementResidual_u,
                                         int nExteriorElementBoundaries_global,
                                         int* exteriorElementBoundariesArray,
                                         int* elementBoundaryElementsArray,
                                         int* elementBoundaryLocalElementBoundariesArray,
                                         double* element_u,
                                         int eN,
                                         double compatibility_condition,
                                         int INTEGRATE_BY_PARTS_DIV_U)
    {
      for (int i=0;i<nDOF_test_element;i++)
        {
          elementResidual_u[i]=0.0;
        }//i
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
          //compute indeces and declare local storage
          register int eN_k = eN*nQuadraturePoints_element+k,
            eN_k_nSpace = eN_k*nSpace;
            //eN_nDOF_trial_element = eN*nDOF_trial_element;
          register double u=0.0,grad_u[nSpace],
            a=0.0,
            f[nSpace],
            jac[nSpace*nSpace],
            jacDet,
            jacInv[nSpace*nSpace],
            u_grad_trial[nDOF_trial_element*nSpace],
            u_test_dV[nDOF_trial_element],
            u_grad_test_dV[nDOF_test_element*nSpace],
            dV,x,y,z,
            G[nSpace*nSpace],G_dd_G,tr_G;
          //
          //compute solution and gradients at quadrature points
          //
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
          ck.calculateG(jacInv,G,G_dd_G,tr_G);
          //get the trial function gradients
          ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
          //get the solution
          ck.valFromElementDOF(element_u,&u_trial_ref[k*nDOF_trial_element],u);
          //get the solution gradients
          ck.gradFromElementDOF(element_u,u_grad_trial,grad_u);
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
          evaluateCoefficients(alphaBDF,
                               &q_vf[eN_k_nSpace],
                               &q_vs[eN_k_nSpace],
                               q_vos[eN_k],
                               rho_s_min,
                               rho_f_min,
                               f,
                               a);
          //
          //update element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              //register int eN_k_i=eN_k*nDOF_test_element+i;
              //register int eN_k_i_nSpace = eN_k_i*nSpace;
              register int  i_nSpace=i*nSpace;
              elementResidual_u[i] +=
                (INTEGRATE_BY_PARTS_DIV_U == 1 ? ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) : q_divU[eN_k]*u_test_dV[i])
                //                + compatibility_condition*u_test_dV[i] // mql: to make the system solvable if int(div(u))!=0
                + ck.NumericalDiffusion(a,grad_u,&u_grad_test_dV[i_nSpace]);
            }//i
          //
          //save momentum for time history and velocity for subgrid error
          //save solution for other models
          //

          q_u[eN_k] = u;
          for (int I=0;I<nSpace;I++)
            q_grad_u[eN_k_nSpace+I] = grad_u[I];
        }
    }

    void calculateResidual(//element
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
                           int* isDOFBoundary,
                           int* isFluxBoundary,
                           int* u_l2g,
                           double* u_dof,
                           double alphaBDF,
                           double* q_vf,
                           double* q_divU,
                           double* q_vs,
                           double* q_vos,
                           double rho_s,
                           double* q_rho_f,
                           double rho_s_min,
                           double rho_f_min,
                           double* ebqe_vf,
                           double* ebqe_vs,
                           double* ebqe_vos,
                           double* ebqe_rho_f,
                           double* q_u,
                           double* q_grad_u,
                           double* ebqe_u,
                           double* ebqe_grad_u,
                           double* ebqe_bc_u_ext,
                           double* ebqe_adv_flux,
                           double* ebqe_diff_flux,
                           double* bc_adv_flux,
                           double* bc_diff_flux,
                           int offset_u,
                           int stride_u,
                           double* globalResidual,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           int INTEGRATE_BY_PARTS_DIV_U)
    {
      double compatibility_condition=0.;
      // COMPUTE COMPATIBILITY CONSTANT
      // mql: Modify the rhs (by adding a constant) so that the Poission system is solvable (assume diffusive flux = 0).
      // Note that this is equivalent to consider a (not known) diffusive flux != 0 s.t. the system is solvable.
      for(int eN=0;eN<nElements_global;eN++)
        {
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              register int eN_k = eN*nQuadraturePoints_element+k;
              register double
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
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
              compatibility_condition -= q_divU[eN_k]*dV;
            }
        }
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
      for(int eN=0;eN<nElements_global;eN++)
        {
          //declare local storage for element residual and initialize
          register double elementResidual_u[nDOF_test_element],element_u[nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              element_u[i] = u_dof[u_l2g[eN_i]];
            }//i
          calculateElementResidual(mesh_trial_ref,
                                   mesh_grad_trial_ref,
                                   mesh_dof,
                                   mesh_l2g,
                                   dV_ref,
                                   u_trial_ref,
                                   u_grad_trial_ref,
                                   u_test_ref,
                                   u_grad_test_ref,
                                   mesh_trial_trace_ref,
                                   mesh_grad_trial_trace_ref,
                                   dS_ref,
                                   u_trial_trace_ref,
                                   u_grad_trial_trace_ref,
                                   u_test_trace_ref,
                                   u_grad_test_trace_ref,
                                   normal_ref,
                                   boundaryJac_ref,
                                   nElements_global,
                                   u_l2g,
                                   u_dof,
                                   alphaBDF,
                                   q_vf,
                                   q_divU,
                                   q_vs,
                                   q_vos,
                                   rho_s,
                                   q_rho_f,
                                   rho_s_min,
                                   rho_f_min,
                                   ebqe_vf,
                                   ebqe_vs,
                                   ebqe_vos,
                                   ebqe_rho_f,
                                   q_u,
                                   q_grad_u,
                                   ebqe_u,
                                   ebqe_grad_u,
                                   offset_u,
                                   stride_u,
                                   elementResidual_u,
                                   nExteriorElementBoundaries_global,
                                   exteriorElementBoundariesArray,
                                   elementBoundaryElementsArray,
                                   elementBoundaryLocalElementBoundariesArray,
                                   element_u,
                                   eN,
                                   compatibility_condition,
                                   INTEGRATE_BY_PARTS_DIV_U);
          //
          //load element into global residual and save element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              globalResidual[offset_u+stride_u*u_l2g[eN_i]]+=elementResidual_u[i];
            }//i
        }//elements
      //
      //loop over exterior element boundaries to calculate levelset gradient
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
        {
          register int ebN = exteriorElementBoundariesArray[ebNE],
            eN  = elementBoundaryElementsArray[ebN*2+0],
            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
            //eN_nDOF_trial_element = eN*nDOF_trial_element;
          register double elementResidual_u[nDOF_test_element];
          double element_u[nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              element_u[i] = u_dof[u_l2g[eN_i]];
              elementResidual_u[i] = 0.0;
            }//i
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                ebNE_kb_nSpace = ebNE_kb*nSpace,
                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                ebN_local_kb_nSpace = ebN_local_kb*nSpace;
              register double h_b=0.0,
                penalty=0.0,
                u_ext=0.0,
                bc_u_ext=0.0,
                adv_flux_ext=0.0,
                diff_flux_ext=0.0,
                a_ext,
                f_ext[nSpace],
                grad_u_ext[nSpace],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,
                u_test_dS[nDOF_test_element],
                u_grad_trial_trace[nDOF_trial_element*nSpace],
                u_grad_test_dS[nDOF_test_element*nSpace],
                normal[nSpace],x_ext,y_ext,z_ext,
                G[nSpace*nSpace],G_dd_G,tr_G;
              //
              //calculate the solution and gradients at quadrature points
              //
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
              dS = metricTensorDetSqrt*dS_ref[kb];
              //get the metric tensor
              //cek todo use symmetry
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              ck.calculateGScale(G,normal,h_b);
              penalty = 10.0/h_b;
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
              //solution and gradients
              ck.valFromElementDOF(element_u,&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
              ck.gradFromElementDOF(element_u,u_grad_trial_trace,grad_u_ext);
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
              bc_u_ext = isDOFBoundary[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary[ebNE_kb])*u_ext;
              //
              //calculate the pde coefficients using the solution and the boundary values for the solution
              //
              evaluateCoefficients(alphaBDF,
                                   &ebqe_vf[ebNE_kb_nSpace],
                                   &ebqe_vs[ebNE_kb_nSpace],
                                   q_vos[ebNE_kb],
                                   rho_s_min,
                                   rho_f_min,
                                   f_ext,
                                   a_ext);
              ebqe_u[ebNE_kb] = u_ext;
              for (int I=0;I<nSpace;I++)
                ebqe_grad_u[ebNE_kb_nSpace+I] = grad_u_ext[I];
              //
              //calculate the numerical fluxes
              //
              exteriorNumericalAdvectiveFlux(isFluxBoundary[ebNE_kb],
                                             bc_adv_flux[ebNE_kb],
                                             normal,
                                             f_ext,
                                             adv_flux_ext); //=f.normal = [(1-vos)*vf + vos*vs].normal
              exteriorNumericalDiffusiveFlux(isDOFBoundary[ebNE_kb],
                                             isFluxBoundary[ebNE_kb],
                                             normal,
                                             a_ext,
                                             grad_u_ext,
                                             u_ext,
                                             bc_u_ext,
                                             bc_diff_flux[ebNE_kb],
                                             penalty,
                                             diff_flux_ext);
              if (isDOFBoundary[ebNE_kb] != 1)
                diff_flux_ext = 0.0; // mql: don't consider diffusive flux unless Dirichlet BC
              //if(isFluxBoundary[ebNE_kb] == 1)
              //{
              //  adv_flux_ext = 0.0;
              //  diff_flux_ext = 0.0;
              //}
              ebqe_adv_flux[ebNE_kb] = adv_flux_ext;
              ebqe_diff_flux[ebNE_kb] = diff_flux_ext;
              //
              //update residuals
              //
              for (int i=0;i<nDOF_test_element;i++)
                {
                  elementResidual_u[i] +=
                    (INTEGRATE_BY_PARTS_DIV_U == 1 ? ck.ExteriorElementBoundaryFlux(adv_flux_ext,u_test_dS[i]) : 0.)
                    + ck.ExteriorElementBoundaryFlux(diff_flux_ext,u_test_dS[i]) // mql: just != 0 if Dirichlet BC
                    + ck.ExteriorElementBoundaryScalarDiffusionAdjoint(isDOFBoundary[ebNE_kb],
                                                                       0, // mql: if Dirichlet BCs, then the flux BCs don't matter
                                                                       //isFluxBoundary[ebNE_kb],
                                                                       1.0,
                                                                       u_ext,
                                                                       bc_u_ext,
                                                                       normal,
                                                                       a_ext,
                                                                       &u_grad_test_dS[i*nSpace]);
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

    inline void calculateElementJacobian(//element
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
                                         double* u_dof,
                                         double alphaBDF,
                                         double* q_vf,
                                         double* q_vs,
                                         double* q_vos,
                                         double rho_s,
                                         double* q_rho_f,
                                         double rho_s_min,
                                         double rho_f_min,
                                         double* elementJacobian_u_u,
                                         double* element_u,
                                         int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
        for (int j=0;j<nDOF_trial_element;j++)
          {
            elementJacobian_u_u[i*nDOF_trial_element+j]=0.0;
          }
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
          int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
            eN_k_nSpace = eN_k*nSpace;
            //eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

          //declare local storage
          register double u=0.0,
            grad_u[nSpace],
            f[nSpace],
            a=0.0,
            jac[nSpace*nSpace],
            jacDet,
            jacInv[nSpace*nSpace],
            u_grad_trial[nDOF_trial_element*nSpace],
            dV,
            u_test_dV[nDOF_test_element],
            u_grad_test_dV[nDOF_test_element*nSpace],
            x,y,z,
            G[nSpace*nSpace],G_dd_G,tr_G;
          //
          //calculate solution and gradients at quadrature points
          //
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
          ck.calculateG(jacInv,G,G_dd_G,tr_G);
          //get the trial function gradients
          ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
          //get the solution
          ck.valFromElementDOF(element_u,&u_trial_ref[k*nDOF_trial_element],u);
          //get the solution gradients
          ck.gradFromElementDOF(element_u,u_grad_trial,grad_u);
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
          evaluateCoefficients(alphaBDF,
                               &q_vf[eN_k_nSpace],
                               &q_vs[eN_k_nSpace],
                               q_vos[eN_k],
                               rho_s_min,
                               rho_f_min,
                               f,
                               a);
          for(int i=0;i<nDOF_test_element;i++)
            {
              //int eN_k_i=eN_k*nDOF_test_element+i;
              //int eN_k_i_nSpace=eN_k_i*nSpace;
              int i_nSpace=i*nSpace;
              for(int j=0;j<nDOF_trial_element;j++)
                {
                  //int eN_k_j=eN_k*nDOF_trial_element+j;
                  //int eN_k_j_nSpace = eN_k_j*nSpace;
                  int j_nSpace = j*nSpace;
                  elementJacobian_u_u[i*nDOF_trial_element+j] +=
                    ck.NumericalDiffusionJacobian(a,&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]);
                }//j
            }//i
        }//k
    }

    void calculateJacobian(//element
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
                           int* isDOFBoundary,
                           int* isFluxBoundary,
                           int* u_l2g,
                           double* u_dof,
                           double alphaBDF,
                           double* q_vf,
                           double* q_vs,
                           double* q_vos,
                           double rho_s,
                           double* q_rho_f,
                           double rho_s_min,
                           double rho_f_min,
                           double* ebqe_vf,
                           double* ebqe_vs,
                           double* ebqe_vos,
                           double* ebqe_rho_f,
                           int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                           double* globalJacobian,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           int* csrColumnOffsets_eb_u_u)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
        {
          register double  elementJacobian_u_u[nDOF_test_element*nDOF_trial_element],element_u[nDOF_trial_element];
          for (int j=0;j<nDOF_trial_element;j++)
            {
              register int eN_j = eN*nDOF_trial_element+j;
              element_u[j] = u_dof[u_l2g[eN_j]];
            }
          calculateElementJacobian(mesh_trial_ref,
                                   mesh_grad_trial_ref,
                                   mesh_dof,
                                   mesh_l2g,
                                   dV_ref,
                                   u_trial_ref,
                                   u_grad_trial_ref,
                                   u_test_ref,
                                   u_grad_test_ref,
                                   mesh_trial_trace_ref,
                                   mesh_grad_trial_trace_ref,
                                   dS_ref,
                                   u_trial_trace_ref,
                                   u_grad_trial_trace_ref,
                                   u_test_trace_ref,
                                   u_grad_test_trace_ref,
                                   normal_ref,
                                   boundaryJac_ref,
                                   nElements_global,
                                   u_l2g,
                                   u_dof,
                                   alphaBDF,
                                   q_vf,
                                   q_vs,
                                   q_vos,
                                   rho_s,
                                   q_rho_f,
                                   rho_s_min,
                                   rho_f_min,
                                   elementJacobian_u_u,
                                   element_u,
                                   eN);
          //
          //load into element Jacobian into global Jacobian
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i*nDOF_trial_element+j];
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

              register double h_b=0.0,
                u_ext=0.0,
                grad_u_ext[nSpace],
                m_ext=0.0,
                dm_ext=0.0,
                a_ext=0.0,
                f_ext[nSpace],
                df_ext[nSpace],
                dflux_u_u_ext=0.0,
                bc_u_ext=0.0,
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
                u_grad_test_dS[nDOF_test_element*nSpace],
                normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                penalty=0.0,
                //
                G[nSpace*nSpace],G_dd_G,tr_G;
              //
              //calculate the solution and gradients at quadrature points
              //
              u_ext=0.0;
              for (int I=0;I<nSpace;I++)
                {
                  grad_u_ext[I] = 0.0;
                }
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
              dS = metricTensorDetSqrt*dS_ref[kb];
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              ck.calculateGScale(G,normal,h_b);
              penalty=10.0/h_b;
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
              //calculate the internal and external trace of the pde coefficients
              //
              evaluateCoefficients(alphaBDF,
                                   &ebqe_vf[ebNE_kb_nSpace],
                                   &ebqe_vs[ebNE_kb_nSpace],
                                   q_vos[ebNE_kb],
                                   rho_s_min,
                                   rho_f_min,
                                   f_ext,
                                   a_ext);
              //
              //update the global Jacobian from the flux Jacobian
              //
              for (int i=0;i<nDOF_test_element;i++)
                {
                  register int eN_i = eN*nDOF_test_element+i;
                  //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j,
                        ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j,
                        j_nSpace = j*nSpace;

                      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] +=
                        ExteriorNumericalDiffusiveFluxJacobian(isDOFBoundary[ebNE_kb],
                                                               isFluxBoundary[ebNE_kb],
                                                               normal,
                                                               a_ext,
                                                               u_trial_trace_ref[ebN_local_kb_j],
                                                               &u_grad_trial_trace[j_nSpace],
                                                               penalty)*u_test_dS[i]
                        +
                        ck.ExteriorElementBoundaryScalarDiffusionAdjointJacobian
                        (isDOFBoundary[ebNE_kb],
                         //isFluxBoundary[ebNE_kb],
                         0, // mql: if Dirichlet BCs, then the flux BCs don't matte
                         1.0,
                         u_trial_trace_ref[ebN_local_kb_j],
                         normal,
                         a_ext,
                         &u_grad_test_dS[i*nSpace]);
                    }//j
                }//i
            }//kb
        }//ebNE
    }//computeJacobian
  };//cppPresInc

  inline cppPresInc_base* newPresInc(int nSpaceIn,
                                 int nQuadraturePoints_elementIn,
                                 int nDOF_mesh_trial_elementIn,
                                 int nDOF_trial_elementIn,
                                 int nDOF_test_elementIn,
                                 int nQuadraturePoints_elementBoundaryIn,
                                 int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppPresInc_base,cppPresInc,CompKernel>(nSpaceIn,
                                                                                               nQuadraturePoints_elementIn,
                                                                                               nDOF_mesh_trial_elementIn,
                                                                                               nDOF_trial_elementIn,
                                                                                               nDOF_test_elementIn,
                                                                                               nQuadraturePoints_elementBoundaryIn,
                                                                                               CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppPresInc_base,cppPresInc,CompKernel>(nSpaceIn,
                                                                                             nQuadraturePoints_elementIn,
                                                                                             nDOF_mesh_trial_elementIn,
                                                                                             nDOF_trial_elementIn,
                                                                                             nDOF_test_elementIn,
                                                                                             nQuadraturePoints_elementBoundaryIn,
                                                                                             CompKernelFlag);
  }
}//proteus
#endif
