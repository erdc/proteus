#include "cRANS2P.h"
#define nComponents 1
#define nSpace 3
#define nQuadraturePoints_element 5
#define nDOF_trial_element 4
#define nDOF_test_element 4

void calculateElementResidual(int nElements_global,
                              double* p, double* u, double *v, double* w,
                              double* p_dof, double* u_dof, double* v_dof, double* w_dof,
                              int* p_l2g, int* u_l2g, int* v_l2g, int* w_l2g,
                              double* p_trial, double* u_trial, double* v_trial, double* w_trial,
                              double* grad_p, double* grad_u, double* grad_v, double* grad_w,
                              double* p_grad_trial, double* u_grad_trial, double* v_grad_trial, double* w_grad_trial,
                              double eps_rho,
                              double eps_mu,
                              double sigma,
                              double rho_0,
                              double nu_0,
                              double rho_1,
                              double nu_1,
                              double* g,
                              double* phi,
                              double* n,
                              double* kappa,
                              double* p,
                              double* grad_p,
                              double* u,
                              double* v,
                              double* w,
                              double*  mom_u_acc,
                              double*  dmom_u_acc_u,
                              double*  mom_v_acc,
                              double*  dmom_v_acc_v,
                              double*  mom_w_acc,
                              double*  dmom_w_acc_w,
                              double*  mass_adv,
                              double*  dmass_adv_u,
                              double*  dmass_adv_v,
                              double*  dmass_adv_w,
                              double*  mom_u_adv,
                              double*  dmom_u_adv_u,
                              double*  dmom_u_adv_v,
                              double*  dmom_u_adv_w,
                              double*  mom_v_adv,
                              double*  dmom_v_adv_u,
                              double*  dmom_v_adv_v,
                              double*  dmom_v_adv_w,
                              double*  mom_w_adv,
                              double*  dmom_w_adv_u,
                              double*  dmom_w_adv_v,
                              double*  dmom_w_adv_w,
                              double*  mom_u_diff_ten,
                              double*  mom_v_diff_ten,
                              double*  mom_w_diff_ten,
                              double*  mom_uv_diff_ten,
                              double*  mom_uw_diff_ten,
                              double*  mom_vu_diff_ten,
                              double*  mom_vw_diff_ten,
                              double*  mom_wu_diff_ten,
                              double*  mom_wv_diff_ten,
                              double*  mom_u_source,
                              double*  mom_v_source,
                              double*  mom_w_source,
                              double*  mom_u_ham,
                              double*  dmom_u_ham_grad_p,
                              double*  mom_v_ham,
                              double*  dmom_v_ham_grad_p,
                              double*  mom_w_ham,
                              double*  dmom_w_ham_grad_p,
                              double dt,
                              double* mom_u_acc_t, double* dmom_u_acc_u_t,
                              double* mom_v_acc_t, double* dmom_v_acc_v_t,
                              double* mom_w_acc_t, double* dmom_w_acc_w_t,
                              double* pdeResidual_p, double* dpdeResidual_p_u, double* dpdeResidual_p_v, double* dpdeResidual_p_w,
                              double* pdeResidual_u, double* dpdeResidual_u_p, double* dpdeResidual_u_u, double* dpdeResidual_u_v, double* dpdeResidual_u_w,
                              double* pdeResidual_v, double* dpdeResidual_v_p, double* dpdeResidual_v_u, double* dpdeResidual_v_v, double* dpdeResidual_v_w,
                              double* pdeResidual_w, double* dpdeResidual_w_p, double* dpdeResidual_w_u, double* dpdeResidual_w_v, double* dpdeResidual_w_w,
                              double* Lstar_u_p, double* Lstar_v_p, double* Lstar_w_p,
                              double* Lstar_u_u, double* Lstar_v_v, double* Lstar_w_w,
                              double* Lstar_p_u, double* Lstar_p_v, double* Lstar_p_w,
                              double* velocity_old,
                              double* cfl,
                              double* subgridError_p, double* dsubgridError_p_u, double* dsubgridError_p_v, double* dsubgridError_p_w,
                              double* subgridError_u, double* dsubgridError_u_p, double* dsubgridError_u_u, double* dsubgridError_u_v, double* dsubgridError_u_w,
                              double* subgridError_v, double* dsubgridError_v_p, double* dsubgridError_v_u, double* dsubgridError_v_v, double* dsubgridError_v_w,
                              double* subgridError_w, double* dsubgridError_w_p, double* dsubgridError_w_u, double* dsubgridError_w_v, double* dsubgridError_w_w,
                              double* numDiff_p, double* numDiff_u, double* numDiff_v, double* numDiff_w,
                              double* elementResidual_p, double* elementResidual_u, double* elementResidual_v, double* elementResidual_w)
{
  //loop over elements
  for(eN=0;eN<nElements_global;eN++)
    {
      for  (k=0;k<nQuadraturePoints_element;k++)
        {
          //
          //solution and gradients at quadrature points
          //
          for (j=0;j<nDOF_trial_element;j++)
            {
              for (t=0;t<nComponents;t++)
                {
                  valFromDOF(eN,k,j,t,p,p_dof,p_l2g,p_trial);
                  valFromDOF(eN,k,j,t,u,u_dof,u_l2g,u_trial);
                  valFromDOF(eN,k,j,t,v,v_dof,v_l2g,v_trial);
                  valFromDOF(eN,k,j,t,w,w_dof,w_l2g,w_trial);
                  for  (I=0;I<nSpace;I++)
                    {
                      gradFromDOF(eN,k,j,t,I,grad_p,p_dof,p_l2g,p_grad_trial);
                      gradFromDOF(eN,k,j,t,I,grad_u,u_dof,u_l2g,u_grad_trial);
                      gradFromDOF(eN,k,j,t,I,grad_v,v_dof,v_l2g,v_grad_trial);
                      gradFromDOF(eN,k,j,t,I,grad_w,w_dof,w_l2g,w_grad_trial);
                    }
                }
            }
          //
          //pde coefficients at quadrature points
          //
          K = eN*nQuadraturePoints_element+k;
          evaluateCoefficients(K,
                               eps_rho,
                               eps_mu,
                               sigma,
                               rho_0,
                               nu_0,
                               rho_1,
                               nu_1,
                               g,
                               phi,
                               n,
                               kappa,
                               p,
                               grad_p,
                               u,
                               v,
                               w,
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
                               mom_u_diff_ten,
                               mom_v_diff_ten,
                               mom_w_diff_ten,
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
                               mom_v_ham,
                               dmom_v_ham_grad_p,
                               mom_w_ham,
                               dmom_w_ham_grad_p);          
          //
          //moving mesh
          //
          //omit for now
          //
          //time integration
          //
          backwardEuler(eN,k,dt,mom_u_acc_old,mom_u_acc,dmom_u_acc_u,mom_u_acc_t,dmom_u_acc_u_t);
          backwardEuler(eN,k,dt,mom_v_acc_old,mom_v_acc,dmom_v_acc_v,mom_v_acc_t,dmom_v_acc_v_t);
          backwardEuler(eN,k,dt,mom_w_acc_old,mom_w_acc,dmom_w_acc_v,mom_w_acc_t,dmom_w_acc_w_t);
          //
          //subgrid error (strong residual, adjoint, jacobian of strong residual)
          //
          //strong residual
          for  (I=0;I<nSpace;I++)
            {
              updateAdvection_strong(eN,k,I,dmass_adv_u_sge,grad_u,pdeResidual_p);
              updateAdvection_strong(eN,k,I,dmass_adv_v_sge,grad_v,pdeResidual_p);
              updateAdvection_strong(eN,k,I,dmass_adv_w_sge,grad_w,pdeResidual_p);
              
              updateAdvection_strong(eN,k,I,dmom_adv_u_sge,grad_u,pdeResidual_u);
              updateAdvection_strong(eN,k,I,dmom_adv_v_sge,grad_v,pdeResidual_v);
              updateAdvection_strong(eN,k,I,dmom_adv_w_sge,grad_w,pdeResidual_w);
              
              updateHamiltonian_strong(eN,k,I,dmom_u_ham_grad_p_sge,grad_p,pdeResidual_u);
              updateHamiltonian_strong(eN,k,I,dmom_v_ham_grad_p_sge,grad_p,pdeResidual_v);
              updateHamiltonian_strong(eN,k,I,dmom_w_ham_grad_p_sge,grad_p,pdeResidual_w);
            }
          updateReaction_strong(eN,k,mom_u_source,pdeResidual_u);
          updateReaction_strong(eN,k,mom_v_source,pdeResidual_v);
          updateReaction_strong(eN,k,mom_w_source,pdeResidual_w);
          
          updateMass_strong(eN,k,mom_u_acc_t,pdeResidual_u);
          updateMass_strong(eN,k,mom_v_acc_t,pdeResidual_v);
          updateMass_strong(eN,k,mom_w_acc_t,pdeResidual_w); 
          //adjoint
          for (i=0;i<nDOF_test_element;i++)
            {
              for  (I=0;I<nSpace;I++)
                {
                  updateAdvection_adjoint(eN,k,i,I,dmass_adv_u_sge,grad_test_p_dV,Lstar_u_p);
                  updateAdvection_adjoint(eN,k,i,I,dmass_adv_v_sge,grad_test_p_dV,Lstar_v_p);
                  updateAdvection_adjoint(eN,k,i,I,dmass_adv_w_sge,grad_test_p_dV,Lstar_w_p);
                  
                  updateAdvection_adjoint(eN,k,i,I,dmom_adv_u_sge,grad_test_u_dV,Lstar_u_u);
                  updateAdvection_adjoint(eN,k,i,I,dmom_adv_v_sge,grad_test_v_dV,Lstar_v_v);
                  updateAdvection_adjoint(eN,k,i,I,dmom_adv_w_sge,grad_test_w_dV,Lstar_w_w);

                  updateHamiltonian_adjoint(eN,k,i,I,dmom_u_ham_grad_p_sge,grad_test_p_dV,Lstar_p_u);
                  updateHamiltonian_adjoint(eN,k,i,I,dmom_v_ham_grad_p_sge,grad_test_p_dV,Lstar_p_v);
                  updateHamiltonian_adjoint(eN,k,i,I,dmom_w_ham_grad_p_sge,grad_test_p_dV,Lstar_p_w);
                }
            }
          //Jacobian of strong residual
          for (j=0;j<nDOF_trial_element;j++)
            {
              for  (I=0;I<nSpace;I++)
                {
                  updateAdvectionJacobina_strong(eN,k,j,I,dmass_adv_u_sge,grad_trial_u,dpdeResidual_p_u);
                  updateAdvectionJacobina_strong(eN,k,j,I,dmass_adv_v_sge,grad_trial_v,dpdeResidual_p_v);
                  updateAdvectionJacobina_strong(eN,k,j,I,dmass_adv_w_sge,grad_trial_w,dpdeResidual_p_w);

                  updateAdvectionJacobina_strong(eN,k,j,I,dmom_adv_u_sge,grad_trial_u,dpdeResidual_u_u);
                  updateAdvectionJacobina_strong(eN,k,j,I,dmom_adv_v_sge,grad_trial_v,dpdeResidual_v_v);
                  updateAdvectionJacobina_strong(eN,k,j,I,dmom_adv_w_sge,grad_trial_w,dpdeResidual_w_w);

                  updateHamiltonianJacobina_strong(eN,k,j,I,dmom_u_ham_grad_p_sge,grad_trial_p,dpdeResidual_u_p);
                  updateHamiltonianJacobina_strong(eN,k,j,I,dmom_v_ham_grad_p_sge,grad_trial_p,dpdeResidual_v_p);
                  updateHamiltonianJacobina_strong(eN,k,j,I,dmom_w_ham_grad_p_sge,grad_trial_p,dpdeResidual_w_p);
                }
              updateMassJacobina_strong(eN,k,j,mom_u_acc_t,dpdeResidual_u_u);
              updateMassJacobina_strong(eN,k,j,mom_v_acc_t,dpdeResidual_v_v);
              updateMassJacobina_strong(eN,k,j,mom_w_acc_t,dpdeResidual_w_w); 
            }
          //tau and tau*Res
          calculateSubgridError_tau(eN,k,nSpace,hFactor,elementDiameter,dmom_u_acc_u_t,dmom_u_acc_u,velocity_old,tau_0,tau_1,cfl);
          calculateSubgridError_tauRes(eN,k,
                                       tau_0,
                                       tau_1,
                                       pdeResidual_p,
                                       dpdeResidual_p_u,
                                       dpdeResidual_p_v,
                                       dpdeResidual_p_w,
                                       pdeResidual_u,
                                       dpdeResidual_u_p,
                                       dpdeResidual_u_u,
                                       dpdeResidual_u_v,
                                       dpdeResidual_u_w,
                                       pdeResidual_v,
                                       dpdeResidual_v_p,
                                       dpdeResidual_v_u,
                                       dpdeResidual_v_v,
                                       dpdeResidual_v_w,
                                       pdeResidual_w,
                                       dpdeResidual_w_p,
                                       dpdeResidual_w_u,
                                       dpdeResidual_w_v,
                                       dpdeResidual_w_w,
                                       subgridError_p,
                                       dsubgridError_p_u,
                                       dsubgridError_p_v,
                                       dsubgridError_p_w,
                                       subgridError_u,
                                       dsubgridError_u_p,
                                       dsubgridError_u_u,
                                       dsubgridError_u_v,
                                       dsubgridError_u_w,
                                       subgridError_v,
                                       dsubgridError_v_p,
                                       dsubgridError_v_u,
                                       dsubgridError_v_v,
                                       dsubgridError_v_w,
                                       subgridError_w,
                                       dsubgridError_w_p,
                                       dsubgridError_w_u,
                                       dsubgridError_w_v,
                                       dsubgridError_w_w);
          //
          //shock capturing diffusion
          //
          calculateNumericalDiffusion(eN,k,nSpace,shockCapturingDiffusion,elementDiameter,pdeResidual_p,grad_p,numDiff_p);
          calculateNumericalDiffusion(eN,k,nSpace,shockCapturingDiffusion,elementDiameter,pdeResidual_u,grad_u,numDiff_u);
          calculateNumericalDiffusion(eN,k,nSpace,shockCapturingDiffusion,elementDiameter,pdeResidual_v,grad_v,numDiff_v);
          calculateNumericalDiffusion(eN,k,nSpace,shockCapturingDiffusion,elementDiameter,pdeResidual_w,grad_w,numDiff_w);
/*           // */
/*           //element residual */
/*           // */
/*           for(i=0;i<nDOF_test_element;i++) */
/*             { */
/*               for (I=0;I<nSpace;I++) */
/*                 { */
/*                   updateAdvection_weak(eN,k,i,I,mass_adv,grad_test_dV,elementResidual_p); */
/*                   updateAdvection_weak(eN,k,i,I,mom_u_adv,grad_test_dV,elementResidual_u); */
/*                   updateAdvection_weak(eN,k,i,I,mom_v_adv,grad_test_dV,elementResidual_v); */
/*                   updateAdvection_weak(eN,k,i,I,mom_w_adv,grad_test_dV,elementResidual_w); */
                  
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_u_u_rowptr,sdInfo_u_u_colind,mom_u_diff_ten,grad_u,grad_test_dV,elementResidual_u); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_u_v_rowptr,sdInfo_u_v_colind,mom_uv_diff_ten,grad_v,grad_test_dV,elementResidual_u); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_u_w_rowptr,sdInfo_u_w_colind,mom_uw_diff_ten,grad_w,grad_test_dV,elementResidual_u); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_v_v_rowptr,sdInfo_v_v_colind,mom_v_diff_ten,grad_v,grad_test_dV,elementResidual_v); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_v_u_rowptr,sdInfo_v_u_colind,mom_vu_diff_ten,grad_u,grad_test_dV,elementResidual_v); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_v_w_rowptr,sdInfo_v_w_colind,mom_vw_diff_ten,grad_w,grad_test_dV,elementResidual_v); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_w_w_rowptr,sdInfo_w_w_colind,mom_w_diff_ten,grad_w,grad_test_dV,elementResidual_w); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_w_u_rowptr,sdInfo_w_u_colind,mom_wu_diff_ten,grad_u,grad_test_dV,elementResidual_w); */
/*                   updateDiffusion_weak(eN,k,i,I,sdInfo_w_v_rowptr,sdInfo_w_v_colind,mom_wv_diff_ten,grad_v,grad_test_dV,elementResidual_w); */
/*                 } */
/*               updateReaction_weak(eN,k,i,mom_u_source,test_dV,elementResidual_u); */
/*               updateReaction_weak(eN,k,i,mom_v_source,test_dV,elementResidual_v); */
/*               updateReaction_weak(eN,k,i,mom_w_source,test_dV,elementResidual_w); */
              
/*               updateHamiltonian_weak(eN,k,i,mom_u_ham,test_dV,elementResidual_u); */
/*               updateHamiltonian_weak(eN,k,i,mom_v_ham,test_dV,elementResidual_v); */
/*               updateHamiltonian_weak(eN,k,i,mom_w_ham,test_dV,elementResidual_w); */
              
/*               updateSubgridError(eN,k,i,subgridError_u,Lstar_u_p,elementResidual_p); */
/*               updateSubgridError(eN,k,i,subgridError_v,Lstar_v_p,elementResidual_p); */
/*               updateSubgridError(eN,k,i,subgridError_w,Lstar_w_p,elementResidual_p); */
              
/*               updateSubgridError(eN,k,i,subgridError_p,Lstar_p_u,elementResidual_u); */
/*               updateSubgridError(eN,k,i,subgridError_u,Lstar_u_u,elementResidual_u); */
              
/*               updateSubgridError(eN,k,i,subgridError_p,Lstar_p_v,elementResidual_v); */
/*               updateSubgridError(eN,k,i,subgridError_v,Lstar_v_v,elementResidual_v); */
              
/*               updateSubgridError(eN,k,i,subgridError_p,Lstar_p_w,elementResidual_w); */
/*               updateSubgridError(eN,k,i,subgridError_w,Lstar_w_w,elementResidual_w); */
              
/*               updateNumericalDiffusion(eN,k,i,numDiff_p,grad_p,grad_test_dV,elementResidual_p); */
/*               updateNumericalDiffusion(eN,k,i,numDiff_u,grad_u,grad_test_dV,elementResidual_u); */
/*               updateNumericalDiffusion(eN,k,i,numDiff_v,grad_v,grad_test_dV,elementResidual_v); */
/*               updateNumericalDiffusion(eN,k,i,numDiff_w,grad_w,grad_test_dV,elementResidual_w); */
/*             } */
/*         } */
/*       //loop over element boundaries */
/*       for(ebN=0;ebN<nElementBoundaries_element;ebN++) */
/*         { */
/*           for(kb=0;kb<nQuadraturePoints_elementBoundary;kb++) */
/*             { */
/*               // */
/*               //solution and gradients at quadrature points */
/*               // */
/*               for (j=0;j<nDOF_trial_element;j++) */
/*                 { */
/*                   for (t=0;t<nComponents;t++) */
/*                     { */
/*                       valFromDOF_trace(eN,ebN,kb,j,t,p_trace,p_dof,p_l2g,p_trial_trace); */
/*                       valFromDOF_trace(eN,ebN,kb,j,t,u_trace,u_dof,u_l2g,u_trial_trace); */
/*                       valFromDOF_trace(eN,ebN,kb,j,t,v_trace,v_dof,v_l2g,v_trial_trace); */
/*                       valFromDOF_trace(eN,ebN,kb,j,t,w_trace,w_dof,w_l2g,w_trial_trace); */
/*                       for  (I=0;I<nSpace;I++) */
/*                         { */
/*                           gradFromDOF_trace(eN,ebN,kb,j,t,I,grad_p_trace,p_dof,p_l2g,p_grad_trial_trace); */
/*                           gradFromDOF_trace(eN,ebN,kb,j,t,I,grad_u_trace,u_dof,u_l2g,u_grad_trial_trace); */
/*                           gradFromDOF_trace(eN,ebN,kb,j,t,I,grad_v_trace,v_dof,v_l2g,v_grad_trial_trace); */
/*                           gradFromDOF_trace(eN,ebN,kb,j,t,I,grad_w_trace,w_dof,w_l2g,w_grad_trial_trace); */
/*                         } */
/*                     } */
/*                 } */
/*               // */
/*               //pde coefficients at quadrature points */
/*               // */
/*               K = eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+ebN*nQuadraturePoints_elementBoundary+kb; */
/*               evaluateCoefficients(K, */
/*                                    eps_rho, */
/*                                    eps_mu, */
/*                                    sigma, */
/*                                    rho_0, */
/*                                    nu_0, */
/*                                    rho_1, */
/*                                    nu_1, */
/*                                    g, */
/*                                    phi, */
/*                                    n, */
/*                                    kappa, */
/*                                    p, */
/*                                    grad_p, */
/*                                    u, */
/*                                    v, */
/*                                    w, */
/*                                    mom_u_acc, */
/*                                    dmom_u_acc_u, */
/*                                    mom_v_acc, */
/*                                    dmom_v_acc_v, */
/*                                    mom_w_acc, */
/*                                    dmom_w_acc_w, */
/*                                    mass_adv, */
/*                                    dmass_adv_u, */
/*                                    dmass_adv_v, */
/*                                    dmass_adv_w, */
/*                                    mom_u_adv, */
/*                                    dmom_u_adv_u, */
/*                                    dmom_u_adv_v, */
/*                                    dmom_u_adv_w, */
/*                                    mom_v_adv, */
/*                                    dmom_v_adv_u, */
/*                                    dmom_v_adv_v, */
/*                                    dmom_v_adv_w, */
/*                                    mom_w_adv, */
/*                                    dmom_w_adv_u, */
/*                                    dmom_w_adv_v, */
/*                                    dmom_w_adv_w, */
/*                                    mom_u_diff_ten, */
/*                                    mom_v_diff_ten, */
/*                                    mom_w_diff_ten, */
/*                                    mom_uv_diff_ten, */
/*                                    mom_uw_diff_ten, */
/*                                    mom_vu_diff_ten, */
/*                                    mom_vw_diff_ten, */
/*                                    mom_wu_diff_ten, */
/*                                    mom_wv_diff_ten, */
/*                                    mom_u_source, */
/*                                    mom_v_source, */
/*                                    mom_w_source, */
/*                                    mom_u_ham, */
/*                                    dmom_u_ham_grad_p, */
/*                                    mom_v_ham, */
/*                                    dmom_v_ham_grad_p, */
/*                                    mom_w_ham, */
/*                                    dmom_w_ham_grad_p);      */
/*               // */
/*               //moving domain */
/*               // */
/*               //skip for now */
/*               // */
/*               //time integration */
/*               // */
/*               //skip for now */
/*               // */
/*               //average velocities */
/*               // */
/*               //skip for now */
/*             }//quad */
/*         }//element boundaries per element */
/*     }//elements */
/*   for (ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) */
/*     { */
/*       ebN = exteriorElementBoundariesArray[ebNE]; */
/*       eN  = elementBoundaryElementsArray[ebN*2+0]; */
/*       ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0]; */
/*       for  (kb=0;k<nQuadraturePoints_elementBoundary;kb++) */
/*         { */
/*           // */
/*           //solution and gradients at quadrature points */
/*           // */
/*           for (j=0;j<nDOF_trial_element;j++) */
/*             { */
/*               for (t=0;t<nComponents;t++) */
/*                 { */
/*                   uFromDOF_trace_ext(ebNE,kb,j,t,p_trace_ext,p_dof,p_l2g,p_trial_trace); */
/*                   uFromDOF_trace_ext(ebNE,kb,j,t,u_trace_ext,u_dof,u_l2g,u_trial_trace); */
/*                   uFromDOF_trace_ext(ebNE,kb,j,t,v_trace_ext,v_dof,v_l2g,v_trial_trace); */
/*                   uFromDOF_trace_ext(ebNE,kb,j,t,w_trace_ext,w_dof,w_l2g,w_trial_trace); */
/*                   for  (I=0;I<nSpace;I++) */
/*                     { */
/*                       grad_uFromDOF_trace_ext(ebNE,kb,j,t,I,grad_p_trace_ext,p_dof,p_l2g,p_grad_trial_trace); */
/*                       grad_uFromDOF_trace_ext(ebNE,kb,j,t,I,grad_u_trace_ext,u_dof,u_l2g,u_grad_trial_trace); */
/*                       grad_uFromDOF_trace_ext(ebNE,kb,j,t,I,grad_v_trace_ext,v_dof,v_l2g,v_grad_trial_trace); */
/*                       grad_uFromDOF_trace_ext(ebNE,kb,j,t,I,grad_w_trace_ext,w_dof,w_l2g,w_grad_trial_trace); */
/*                     } */
/*                 } */
/*             } */
/*           // */
/*           //pde coefficients */
/*           // */
/*           // */
/*           //numerical flux */
/*           // */
/*           cfemIntegrals.updateExteriorElementBoundaryFlux(self.mesh.exteriorElementBoundariesArray, */
/*                                                           self.mesh.elementBoundaryElementsArray, */
/*                                                           self.mesh.elementBoundaryLocalElementBoundariesArray, */
/*                                                           self.ebqe[('advectiveFlux',ci)], */
/*                                                           self.ebqe[('w*dS_f',ci)], */
/*                                                           self.elementResidual[ci]) */
/*             cfemIntegrals.updateExteriorElementBoundaryFlux(self.mesh.exteriorElementBoundariesArray, */
/*                                                             self.mesh.elementBoundaryElementsArray, */
/*                                                             self.mesh.elementBoundaryLocalElementBoundariesArray, */
/*                                                             self.ebqe[('diffusiveFlux',ck,ci)], */
/*                                                             self.ebqe[('w*dS_a',ck,ci)], */
/*                                                             self.elementResidual[ci]) */
/*             cfemIntegrals.updateExteriorElementBoundaryDiffusionAdjoint_sd(self.coefficients.sdInfo[(ci,ck)][0],self.coefficients.sdInfo[(ci,ck)][1], */
/*                                                                            self.numericalFlux.isDOFBoundary[ck], */
/*                                                                            self.mesh.exteriorElementBoundariesArray, */
/*                                                                            self.mesh.elementBoundaryElementsArray, */
/*                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray, */
/*                                                                            self.numericalFlux.boundaryAdjoint_sigma, */
/*                                                                            self.ebqe[('u',ck)], */
/*                                                                            self.numericalFlux.ebqe[('u',ck)], */
/*                                                                            self.ebqe['n'], */
/*                                                                            self.ebqe[('a',ci,ck)], */
/*                                                                            self.ebqe[('grad(v)',ci)],#cek grad w */
/*                                                                            self.ebqe[('dS_u',ci)], */
/*                                                                            self.elementResidual[ci]) */
/*             }//qu */
/*     }//exterior element boundaires */
}

