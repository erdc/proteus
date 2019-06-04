#ifndef ADDEDMASS_H
#define ADDEDMASS_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

using namespace std;

namespace proteus
{
  class cppAddedMass_base
  {
  public:
    virtual ~cppAddedMass_base(){}
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
                                   int nElementBoundaries_owned,
                                   int* u_l2g,
                                   double* u_dof,
                                   double* q_rho,
                                   int offset_u,
                                   int stride_u,
                                   double* globalResidual,
                                   int nExteriorElementBoundaries_global,
                                   int* exteriorElementBoundariesArray,
                                   int* elementBoundaryElementsArray,
                                   int* elementBoundaryLocalElementBoundariesArray,
                                   int* elementBoundaryMaterialTypesArray,
                                   double* Aij,
                                   int added_mass_i,
                                   double* barycenters,
                                   int* flags_rigidbody)=0;
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
                                   int* u_l2g,
                                   double* u_dof,
                                   double* q_rho,
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
  class cppAddedMass : public cppAddedMass_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    cppAddedMass():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}
    inline
      void evaluateCoefficients(const double& rho,
                                double& a)
    {
      a = 1.0/rho;
    }

    inline
      void exteriorNumericalDiffusiveFlux(const double n[nSpace],
                                          const double a[nSpace],
                                          int isBodyBoundary,
                                          double& flux)
    {
      flux=0.0;
      if (isBodyBoundary == 1) {
        for (int I=0;I<nSpace;I++) {
          flux -= a[I]*n[I];
        }
      }
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
                                         double* q_rho,
                                         int offset_u,
                                         int stride_u,
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
          evaluateCoefficients(q_rho[eN_k], a);
          //
          //update element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              //register int eN_k_i=eN_k*nDOF_test_element+i;
              //register int eN_k_i_nSpace = eN_k_i*nSpace;
              register int  i_nSpace=i*nSpace;
              elementResidual_u[i] += ck.NumericalDiffusion(a,grad_u,&u_grad_test_dV[i_nSpace]);
            }//i
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
                           int nElementBoundaries_owned,
                           int* u_l2g,
                           double* u_dof,
                           double* q_rho,
                           int offset_u,
                           int stride_u,
                           double* globalResidual,
                           int nExteriorElementBoundaries_global,
                           int* exteriorElementBoundariesArray,
                           int* elementBoundaryElementsArray,
                           int* elementBoundaryLocalElementBoundariesArray,
                           int* elementBoundaryMaterialTypesArray,
                           double* Aij,
                           int added_mass_i,
                           double* barycenters,
                           int* flags_rigidbody)
    {
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
                                   q_rho,
                                   offset_u,
                                   stride_u,
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
              register double penalty=0.0,
                u_ext=0.0,
                bc_u_ext=0.0,
                adv_flux_ext=0.0,
                diff_flux_ext=0.0,
                a_ext,
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
                normal[nSpace],x_ext,y_ext,z_ext=0.0,
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
              ck.calculateGScale(G,normal,penalty);
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
              //solution and gradients
              ck.valFromElementDOF(element_u,&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
              //ck.gradFromElementDOF(element_u,u_grad_trial_trace,grad_u_ext);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
                  for (int I=0;I<nSpace;I++)
                    u_grad_test_dS[j*nSpace+I] = u_grad_trial_trace[j*nSpace+I]*dS;//cek hack, using trial
                }
              //
              //calculate the numerical fluxes
              //
              int eBMT = elementBoundaryMaterialTypesArray[ebN];
              double rx, ry, rz;
              rx = x_ext-barycenters[3*eBMT+0];
              ry = y_ext-barycenters[3*eBMT+1];
              rz = z_ext-barycenters[3*eBMT+2];
              double added_mass_a[3] = {0.0, 0.0, 0.0};
	      if (eBMT > 0)
		{
		  switch (added_mass_i)
		    {
		    case 0:
		      added_mass_a[0] = 1.0;
		      break;
		    case 1:
		      added_mass_a[1] = 1.0;
		      break;
		    case 2:
		      added_mass_a[2] = 1.0;
		      break;
		    case 3:
		      added_mass_a[1] = -rz;
		      added_mass_a[2] =  ry;
		      break;
		    case 4:
		      added_mass_a[0] =  rz;
		      added_mass_a[2] = -rx;
		      break;
		    case 5:
		      added_mass_a[0] = -ry;
		      added_mass_a[1] =  rx;
		      break;
		    default:
		      assert(0);
		    }
		}
              // normalise unit accelerations (necessary for angular ones)
	      //I think we want the angular acceleration to be 1
	      //but the flux uses whatever the linear acceleration works
	      //out to be, so I'm commenting this out for now
              /* double added_mass_a_tot = sqrt(added_mass_a[0]*added_mass_a[0]+ */
              /*                                added_mass_a[1]*added_mass_a[1]+ */
              /*                                added_mass_a[2]*added_mass_a[2]); */
              /* added_mass_a[0] = added_mass_a[0]/added_mass_a_tot; */
              /* added_mass_a[1] = added_mass_a[1]/added_mass_a_tot; */
              /* added_mass_a[2] = added_mass_a[2]/added_mass_a_tot; */
	      
              exteriorNumericalDiffusiveFlux(normal,
                                             added_mass_a,
                                             flags_rigidbody[eBMT],
                                             diff_flux_ext);
              //
              //update residuals
              //
              for (int i=0;i<nDOF_test_element;i++)
                {
                  elementResidual_u[i] +=
                    + ck.ExteriorElementBoundaryFlux(diff_flux_ext,u_test_dS[i]);
                }//i
              //calculate Aij
              if (ebN < nElementBoundaries_owned)
                {
                  double px, py, pz;
                  px = u_ext*normal[0];
                  py = u_ext*normal[1];
                  if (nSpace==3)
                    pz = u_ext*normal[2];
                  else
                    pz=0.0;
                  Aij[36*eBMT+added_mass_i+6*0] += px*dS;
                  Aij[36*eBMT+added_mass_i+6*1] += py*dS;
                  Aij[36*eBMT+added_mass_i+6*2] += pz*dS;
                  Aij[36*eBMT+added_mass_i+6*3] += (ry*pz-rz*py)*dS;
                  Aij[36*eBMT+added_mass_i+6*4] += (rz*px-rx*pz)*dS;
                  Aij[36*eBMT+added_mass_i+6*5] += (rx*py-ry*px)*dS;
                }
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
                                         double* q_rho,
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
          evaluateCoefficients(q_rho[eN_k], a);
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
                           int* u_l2g,
                           double* u_dof,
                           double* q_rho,
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
                                   q_rho,
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
  };//cppAddedMass

  inline cppAddedMass_base* newAddedMass(int nSpaceIn,
					 int nQuadraturePoints_elementIn,
					 int nDOF_mesh_trial_elementIn,
					 int nDOF_trial_elementIn,
					 int nDOF_test_elementIn,
					 int nQuadraturePoints_elementBoundaryIn,
					 int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppAddedMass_base,cppAddedMass,CompKernel>(nSpaceIn,
												   nQuadraturePoints_elementIn,
												   nDOF_mesh_trial_elementIn,
												   nDOF_trial_elementIn,
												   nDOF_test_elementIn,
												   nQuadraturePoints_elementBoundaryIn,
												   CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppAddedMass_base,cppAddedMass,CompKernel>(nSpaceIn,
												 nQuadraturePoints_elementIn,
												 nDOF_mesh_trial_elementIn,
												 nDOF_trial_elementIn,
												 nDOF_test_elementIn,
												 nQuadraturePoints_elementBoundaryIn,
												 CompKernelFlag);
  }
}//proteus
#endif
