#ifndef PRESINIT_H
#define PRESINIT_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
  class cppPresInit_base
  {
  public:
    virtual ~cppPresInit_base(){}
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
                                   double useMetrics,
                                   double epsFactHeaviside,
                                   double epsFactDirac,
                                   double epsFactDiffusion,
                                   int* u_l2g,
                                   double* elementDiameter,
                                   double* nodeDiametersArray,
                                   double* u_dof,
                                   double* q_phi,
                                   double* q_normal_phi,
                                   double* ebqe_phi,
                                   double* ebqe_normal_phi,
                                   double* q_H,
                                   double* q_u,
                                   double* q_n,
                                   double* ebqe_u,
                                   double* ebqe_n,
                                   double* q_r,
                                   double* q_vos,
                                   int offset_u, int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray)=0;
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
                                   double useMetrics,
                                   double epsFactHeaviside,
                                   double epsFactDirac,
                                   double epsFactDiffusion,
                                   int* u_l2g,
                                   double* elementDiameter,
                                   double* nodeDiametersArray,
                                   double* u_dof,
                                   double* q_phi,
                                   double* q_normal_phi,
                                   double* q_H,
                                   double* q_vos,
                                   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                                   double* globalJacobian)=0;
  };

  template<class CompKernelType,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nQuadraturePoints_elementBoundary>
  class cppPresInit : public cppPresInit_base
  {
  public:
    CompKernelType ck;
    cppPresInit():ck()
    {}
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
          HI= phi - eps +       0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
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
    void evaluateCoefficients(const double& epsHeaviside,
                              const double& epsDirac,
                              const double& phi,
                              const double& H,
                              const double& u,
                              const double& porosity,
                              double& r,
                              double& dr)
    {
      r = porosity*smoothedHeaviside(epsHeaviside,phi+u) - H;
      dr = porosity*smoothedDirac(epsDirac,phi+u);
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
                                         double useMetrics,
                                         double epsFactHeaviside,
                                         double epsFactDirac,
                                         double epsFactDiffusion,
                                         int* u_l2g,
                                         double* elementDiameter,
                                         double* nodeDiametersArray,
                                         double* u_dof,
                                         double* q_phi,
                                         double* q_normal_phi,
                                         double* ebqe_phi,
                                         double* ebqe_normal_phi,
                                         double* q_H,
                                         double* q_u,
                                         double* q_n,
                                         double* ebqe_u,
                                         double* ebqe_n,
                                         double* q_r,
                                         double* q_vos,
                                         int offset_u, int stride_u,
                                         double* elementResidual_u,
                                         int nExteriorElementBoundaries_global,
                                         int* exteriorElementBoundariesArray,
                                         int* elementBoundaryElementsArray,
                                         int* elementBoundaryLocalElementBoundariesArray,
                                         double* element_u,
                                         int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
        {
          elementResidual_u[i]=0.0;
        }//i
      double epsHeaviside,epsDirac,epsDiffusion,norm;
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
          //compute indeces and declare local storage
          register int eN_k = eN*nQuadraturePoints_element+k,
            eN_k_nSpace = eN_k*nSpace;
            //eN_nDOF_trial_element = eN*nDOF_trial_element;
          register double u=0.0,grad_u[nSpace],
            r=0.0,dr=0.0,
            jac[nSpace*nSpace],
            jacDet,
            jacInv[nSpace*nSpace],
            u_grad_trial[nDOF_trial_element*nSpace],
            u_test_dV[nDOF_trial_element],
            u_grad_test_dV[nDOF_test_element*nSpace],
            dV,x,y,z,
            G[nSpace*nSpace],G_dd_G,tr_G,h_phi;
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
          ck.calculateH_element(eN,
                                k,
                                nodeDiametersArray,
                                mesh_l2g,
                                mesh_trial_ref,
                                h_phi);
          //get the physical integration weight
          dV = fabs(jacDet)*dV_ref[k];
          ck.calculateG(jacInv,G,G_dd_G,tr_G);

          /* double dir[nSpace]; */
          /* double norm = 1.0e-8; */
          /* for (int I=0;I<nSpace;I++) */
          /*   norm += q_normal_phi[eN_k_nSpace+I]*q_normal_phi[eN_k_nSpace+I]; */
          /* norm = sqrt(norm);    */
          /* for (int I=0;I<nSpace;I++) */
          /*   dir[I] = q_normal_phi[eN_k_nSpace+I]/norm; */
          /* ck.calculateGScale(G,dir,h_phi); */

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
          epsHeaviside = epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          epsDirac     = epsFactDirac*    (useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          epsDiffusion = epsFactDiffusion*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          // *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          evaluateCoefficients(epsHeaviside,
                               epsDirac,
                               q_phi[eN_k],
                               q_H[eN_k],
                               u,
                               1.0 - q_vos[eN_k],
                               r,
                               dr);
          //
          //update element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              //register int eN_k_i=eN_k*nDOF_test_element+i;
              //register int eN_k_i_nSpace = eN_k_i*nSpace;
              register int  i_nSpace=i*nSpace;

              elementResidual_u[i] += ck.Reaction_weak(r,u_test_dV[i]) +
                ck.NumericalDiffusion(epsDiffusion,grad_u,&u_grad_test_dV[i_nSpace]);
            }//i
          //
          //save momentum for time history and velocity for subgrid error
          //save solution for other models
          //

          q_r[eN_k] = r;
          q_u[eN_k] = u;


          norm = 1.0e-8;
          for (int I=0;I<nSpace;I++)
            norm += grad_u[I]*grad_u[I];
          norm = sqrt(norm);
          for(int I=0;I<nSpace;I++)
            q_n[eN_k_nSpace+I] = grad_u[I]/norm;
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
                           double useMetrics,
                           double epsFactHeaviside,
                           double epsFactDirac,
                           double epsFactDiffusion,
                           int* u_l2g,
                           double* elementDiameter,
                           double* nodeDiametersArray,
                           double* u_dof,
                           double* q_phi,
                           double* q_normal_phi,
                           double* ebqe_phi,
                           double* ebqe_normal_phi,
                           double* q_H,
                           double* q_u,
                           double* q_n,
                           double* ebqe_u,
                           double* ebqe_n,
                           double* q_r,
                           double* q_vos,
                           int offset_u, int stride_u,
                           double* globalResidual,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray)
    {
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
                                   useMetrics,
                                   epsFactHeaviside,
                                   epsFactDirac,
                                   epsFactDiffusion,
                                   u_l2g,
                                   elementDiameter,
                                   nodeDiametersArray,
                                   u_dof,
                                   q_phi,
                                   q_normal_phi,
                                   ebqe_phi,
                                   ebqe_normal_phi,
                                   q_H,
                                   q_u,
                                   q_n,
                                   ebqe_u,
                                   ebqe_n,
                                   q_r,
                                   q_vos,
                                   offset_u,stride_u,
                                   elementResidual_u,
                                   nExteriorElementBoundaries_global,
                                   exteriorElementBoundariesArray,
                                   elementBoundaryElementsArray,
                                   elementBoundaryLocalElementBoundariesArray,
                                   element_u,
                                   eN);
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
          //register double elementResidual_u[nDOF_test_element];
          double element_u[nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            {
              register int eN_i=eN*nDOF_test_element+i;
              element_u[i] = u_dof[u_l2g[eN_i]];
            }//i
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                ebNE_kb_nSpace = ebNE_kb*nSpace,
                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                ebN_local_kb_nSpace = ebN_local_kb*nSpace;
              register double u_ext=0.0,
                grad_u_ext[nSpace],
                //m_ext=0.0,
                //dm_ext=0.0,
                //H_ext=0.0,
                //dH_ext[nSpace],
                //flux_ext=0.0,
                //bc_u_ext=0.0,
                //bc_grad_u_ext[nSpace],
                //bc_m_ext=0.0,
                //bc_dm_ext=0.0,
                //bc_H_ext=0.0,
                //bc_dH_ext[nSpace],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,
                //u_test_dS[nDOF_test_element],
                u_grad_trial_trace[nDOF_trial_element*nSpace],
                normal[nSpace],x_ext,y_ext,z_ext,
                G[nSpace*nSpace],G_dd_G,tr_G,norm;
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
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
              //solution and gradients
              ck.valFromElementDOF(element_u,&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
              ck.gradFromElementDOF(element_u,u_grad_trial_trace,grad_u_ext);

              ebqe_u[ebNE_kb] = u_ext;
              norm = 1.0e-8;
              for (int I=0;I<nSpace;I++)
                norm += grad_u_ext[I]*grad_u_ext[I];
              norm = sqrt(norm);
              for (int I=0;I<nSpace;I++)
                ebqe_n[ebNE_kb_nSpace+I] = grad_u_ext[I]/norm;
            }//kb
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
                                         double useMetrics,
                                         double epsFactHeaviside,
                                         double epsFactDirac,
                                         double epsFactDiffusion,
                                         int* u_l2g,
                                         double* elementDiameter,
                                         double* nodeDiametersArray,
                                         double* u_dof,
                                         // double* u_trial,
                                         // double* u_grad_trial,
                                         // double* u_test_dV,
                                         // double* u_grad_test_dV,
                                         double* q_phi,
                                         double* q_normal_phi,
                                         double* q_H,
                                         double* q_vos,
                                         double* elementJacobian_u_u,
                                         double* element_u,
                                         int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
        for (int j=0;j<nDOF_trial_element;j++)
          {
            elementJacobian_u_u[i*nDOF_trial_element+j]=0.0;
          }
      double epsHeaviside,epsDirac,epsDiffusion;
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
          int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
            eN_k_nSpace = eN_k*nSpace;
            //eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

          //declare local storage
          register double u=0.0,
            grad_u[nSpace],
            r=0.0,dr=0.0,
            jac[nSpace*nSpace],
            jacDet,
            jacInv[nSpace*nSpace],
            u_grad_trial[nDOF_trial_element*nSpace],
            dV,
            u_test_dV[nDOF_test_element],
            u_grad_test_dV[nDOF_test_element*nSpace],
            x,y,z,
            G[nSpace*nSpace],G_dd_G,tr_G,h_phi;
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
          ck.calculateH_element(eN,
                                k,
                                nodeDiametersArray,
                                mesh_l2g,
                                mesh_trial_ref,
                                h_phi);
          //get the physical integration weight
          dV = fabs(jacDet)*dV_ref[k];
          ck.calculateG(jacInv,G,G_dd_G,tr_G);

          /* double dir[nSpace]; */
          /* double norm = 1.0e-8; */
          /* for (int I=0;I<nSpace;I++) */
          /*   norm += q_normal_phi[eN_k_nSpace+I]*q_normal_phi[eN_k_nSpace+I]; */
          /* norm = sqrt(norm); */
          /* for (int I=0;I<nSpace;I++) */
          /*   dir[I] = q_normal_phi[eN_k_nSpace+I]/norm; */
          /* ck.calculateGScale(G,dir,h_phi); */


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
          epsHeaviside=epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          epsDirac    =epsFactDirac*    (useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          epsDiffusion=epsFactDiffusion*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          //    *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
          evaluateCoefficients(epsHeaviside,
                               epsDirac,
                               q_phi[eN_k],
                               q_H[eN_k],
                               u,
                               1.0 - q_vos[eN_k],
                               r,
                               dr);
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

                  elementJacobian_u_u[i*nDOF_trial_element+j] += ck.ReactionJacobian_weak(dr,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) +
                    ck.NumericalDiffusionJacobian(epsDiffusion,&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]);
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
                           double useMetrics,
                           double epsFactHeaviside,
                           double epsFactDirac,
                           double epsFactDiffusion,
                           int* u_l2g,
                           double* elementDiameter,
                           double* nodeDiametersArray,
                           double* u_dof,
                           // double* u_trial,
                           // double* u_grad_trial,
                           // double* u_test_dV,
                           // double* u_grad_test_dV,
                           double* q_phi,
                           double* q_normal_phi,
                           double* q_H,
                           double* q_vos,
                           int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                           double* globalJacobian)
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
                                   useMetrics,
                                   epsFactHeaviside,
                                   epsFactDirac,
                                   epsFactDiffusion,
                                   u_l2g,
                                   elementDiameter,
                                   nodeDiametersArray,
                                   u_dof,
                                   q_phi,
                                   q_normal_phi,
                                   q_H,
                                   q_vos,
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
    }//computeJacobian
  };//cppPresInit

  inline cppPresInit_base* newPresInit(int nSpaceIn,
                                 int nQuadraturePoints_elementIn,
                                 int nDOF_mesh_trial_elementIn,
                                 int nDOF_trial_elementIn,
                                 int nDOF_test_elementIn,
                                 int nQuadraturePoints_elementBoundaryIn,
                                 int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppPresInit_base,cppPresInit,CompKernel>(nSpaceIn,
                                                                                               nQuadraturePoints_elementIn,
                                                                                               nDOF_mesh_trial_elementIn,
                                                                                               nDOF_trial_elementIn,
                                                                                               nDOF_test_elementIn,
                                                                                               nQuadraturePoints_elementBoundaryIn,
                                                                                               CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppPresInit_base,cppPresInit,CompKernel>(nSpaceIn,
                                                                                             nQuadraturePoints_elementIn,
                                                                                             nDOF_mesh_trial_elementIn,
                                                                                             nDOF_trial_elementIn,
                                                                                             nDOF_test_elementIn,
                                                                                             nQuadraturePoints_elementBoundaryIn,
                                                                                             CompKernelFlag);
  }
}//proteus
#endif
