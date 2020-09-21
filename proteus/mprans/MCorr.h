#ifndef MCorr_H
#define MCorr_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "equivalent_polynomials.h"
#include PROTEUS_LAPACK_H
#include "ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"

namespace py = pybind11;

namespace proteus
{

  template<int nSpace, int nP, int nQ, int nEBQ>
  using GeneralizedFunctions = equivalent_polynomials::GeneralizedFunctions_mix<nSpace, nP, nQ, nEBQ>;
  //using GeneralizedFunctions = equivalent_polynomials::Regularized<nSpace, nP, nQ>;
  //using GeneralizedFunctions = equivalent_polynomials::EquivalentPolynomials<nSpace, nP, nQ>;

  class MCorr_base
  {
  public:
    std::valarray<double> Rpos, Rneg;
    std::valarray<double> FluxCorrectionMatrix;
    virtual ~MCorr_base(){}
    virtual void calculateResidual(arguments_dict& args, bool useExact)=0;
    virtual void calculateJacobian(arguments_dict& args, bool useExact)=0;
    virtual void elementSolve(arguments_dict& args)=0;
    virtual void elementConstantSolve(arguments_dict& args)=0;
    virtual std::tuple<double, double> globalConstantRJ(arguments_dict& args)=0;
    virtual double calculateMass(arguments_dict& args, bool useExact)=0;
    virtual void setMassQuadrature(arguments_dict& args, bool useExact)=0;
    virtual void FCTStep(arguments_dict& args)=0;
    virtual void calculateMassMatrix(arguments_dict& args)=0;
    virtual void setMassQuadratureEdgeBasedStabilizationMethods(arguments_dict& args,
                                                                bool useExact)=0;
  };

  template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class MCorr : public MCorr_base
    {
    public:
      CompKernelType ck;
      GeneralizedFunctions<nSpace,2,nQuadraturePoints_element,nQuadraturePoints_elementBoundary> gf;
      GeneralizedFunctions<nSpace,2,nDOF_trial_element,nQuadraturePoints_elementBoundary> gf_nodes;
    MCorr():ck()
        {}

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
        r = porosity*(gf.H(epsHeaviside,phi+u) - H);
        dr = porosity*gf.D(epsDirac,phi+u);
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
                                           double* q_porosity,
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
            gf.set_quad(k);
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
                                 q_porosity[eN_k],
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
      void calculateResidual(arguments_dict& args,
                             bool useExact)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& x_ref = args.array<double>("x_ref");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<int>& r_l2g = args.array<int>("r_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& phi_dof = args.array<double>("phi_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
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
        gf.useExact = useExact;
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for element residual and initialize
            register double elementResidual_u[nDOF_test_element],element_u[nDOF_trial_element],element_phi[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i=eN*nDOF_test_element+i;
                element_u[i] = u_dof.data()[u_l2g.data()[eN_i]];
                element_phi[i] = phi_dof.data()[u_l2g.data()[eN_i]] + element_u[i];
              }//i
            double element_nodes[nDOF_mesh_trial_element*3];
            for (int i=0;i<nDOF_mesh_trial_element;i++)
              {
                register int eN_i=eN*nDOF_mesh_trial_element+i;
                for(int I=0;I<3;I++)
                  element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
	      }//i
            gf.calculate(element_phi, element_nodes, x_ref.data(),false);
	    calculateElementResidual(mesh_trial_ref.data(),
				     mesh_grad_trial_ref.data(),
				     mesh_dof.data(),
				     mesh_l2g.data(),
				     dV_ref.data(),
				     u_trial_ref.data(),
				     u_grad_trial_ref.data(),
				     u_test_ref.data(),
				     u_grad_test_ref.data(),
				     mesh_trial_trace_ref.data(),
				     mesh_grad_trial_trace_ref.data(),
				     dS_ref.data(),
				     u_trial_trace_ref.data(),
				     u_grad_trial_trace_ref.data(),
				     u_test_trace_ref.data(),
				     u_grad_test_trace_ref.data(),
				     normal_ref.data(),
				     boundaryJac_ref.data(),
				     nElements_global,
				     useMetrics,
				     epsFactHeaviside,
				     epsFactDirac,
				     epsFactDiffusion,
				     u_l2g.data(),
				     elementDiameter.data(),
				     nodeDiametersArray.data(),
				     u_dof.data(),
				     q_phi.data(),
				     q_normal_phi.data(),
				     ebqe_phi.data(),
				     ebqe_normal_phi.data(),
				     q_H.data(),
				     q_u.data(),
				     q_n.data(),
				     ebqe_u.data(),
				     ebqe_n.data(),
				     q_r.data(),
				     q_porosity.data(),
				     offset_u,stride_u,
				     elementResidual_u,
				     nExteriorElementBoundaries_global,
				     exteriorElementBoundariesArray.data(),
				     elementBoundaryElementsArray.data(),
				     elementBoundaryLocalElementBoundariesArray.data(),
				     element_u,
				     eN);
	    //
	    //load element into global residual and save element residual
	    //
	    for(int i=0;i<nDOF_test_element;i++)
	      {
		register int eN_i=eN*nDOF_test_element+i;

                globalResidual.data()[offset_u+stride_u*r_l2g.data()[eN_i]]+=elementResidual_u[i];
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
            register int ebN = exteriorElementBoundariesArray.data()[ebNE],
              eN  = elementBoundaryElementsArray.data()[ebN*2+0],
              ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0];
            //eN_nDOF_trial_element = eN*nDOF_trial_element;
            //register double elementResidual_u[nDOF_test_element];
            double element_u[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i=eN*nDOF_test_element+i;
                element_u[i] = u_dof.data()[u_l2g.data()[eN_i]];
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
                dS = metricTensorDetSqrt*dS_ref.data()[kb];
                //get the metric tensor
                //cek todo use symmetry
                ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                //compute shape and solution information
                //shape
                ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                //solution and gradients
                ck.valFromElementDOF(element_u,&u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],u_ext);
                ck.gradFromElementDOF(element_u,u_grad_trial_trace,grad_u_ext);

                ebqe_u.data()[ebNE_kb] = u_ext;
                norm = 1.0e-8;
                for (int I=0;I<nSpace;I++)
                  norm += grad_u_ext[I]*grad_u_ext[I];
                norm = sqrt(norm);
                for (int I=0;I<nSpace;I++)
                  ebqe_n.data()[ebNE_kb_nSpace+I] = grad_u_ext[I]/norm;
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
                                           double* q_porosity,
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
            gf.set_quad(k);
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
                                 q_porosity[eN_k],
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
                    elementJacobian_u_u[i*nDOF_trial_element+j] +=
                      ck.ReactionJacobian_weak(dr,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) +
                      ck.NumericalDiffusionJacobian(epsDiffusion,&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]);
                  }//j
              }//i
          }//k
      }

      void calculateJacobian(arguments_dict& args,
                             bool useExact)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& x_ref = args.array<double>("x_ref");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& phi_dof = args.array<double>("phi_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& globalJacobian = args.array<double>("globalJacobian");
        //
        //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
        //
        gf.useExact = useExact;
        for(int eN=0;eN<nElements_global;eN++)
          {
            register double  elementJacobian_u_u[nDOF_test_element*nDOF_trial_element],element_u[nDOF_trial_element],element_phi[nDOF_trial_element];
            for (int j=0;j<nDOF_trial_element;j++)
              {
                register int eN_j = eN*nDOF_trial_element+j;
                element_u[j] = u_dof.data()[u_l2g.data()[eN_j]];
                element_phi[j] = phi_dof.data()[u_l2g.data()[eN_j]] + element_u[j];
              }
            double element_nodes[nDOF_mesh_trial_element*3];
            for (int i=0;i<nDOF_mesh_trial_element;i++)
              {
                register int eN_i=eN*nDOF_mesh_trial_element+i;
                for(int I=0;I<3;I++)
                  element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
	      }//i
          gf.calculate(element_phi, element_nodes, x_ref.data(),false);
	    calculateElementJacobian(mesh_trial_ref.data(),
				     mesh_grad_trial_ref.data(),
				     mesh_dof.data(),
				     mesh_l2g.data(),
				     dV_ref.data(),
				     u_trial_ref.data(),
				     u_grad_trial_ref.data(),
				     u_test_ref.data(),
				     u_grad_test_ref.data(),
				     mesh_trial_trace_ref.data(),
				     mesh_grad_trial_trace_ref.data(),
				     dS_ref.data(),
				     u_trial_trace_ref.data(),
				     u_grad_trial_trace_ref.data(),
				     u_test_trace_ref.data(),
				     u_grad_test_trace_ref.data(),
				     normal_ref.data(),
				     boundaryJac_ref.data(),
				     nElements_global,
				     useMetrics,
				     epsFactHeaviside,
				     epsFactDirac,
				     epsFactDiffusion,
				     u_l2g.data(),
				     elementDiameter.data(),
				     nodeDiametersArray.data(),
				     u_dof.data(),
				     q_phi.data(),
				     q_normal_phi.data(),
				     q_H.data(),
				     q_porosity.data(),
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

                    globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i*nDOF_trial_element+j];
                  }//j
              }//i
          }//elements
      }//computeJacobian
      void elementSolve(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        int maxIts = args.scalar<int>("maxIts");
        double atol = args.scalar<double>("atol");
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
            register double element_u[nDOF_test_element],
              element_du[nDOF_test_element],
              elementResidual_u[nDOF_test_element],
              elementJacobian_u_u[nDOF_test_element*nDOF_trial_element],scale=1.0;
            register PROTEUS_LAPACK_INTEGER elementPivots[nDOF_test_element],
              elementColPivots[nDOF_test_element];
            //double epsHeaviside,epsDirac,epsDiffusion;
            for (int i=0;i<nDOF_test_element;i++)
              {
                element_u[i]=0.0;
              }//i
            calculateElementResidual(mesh_trial_ref.data(),
                                     mesh_grad_trial_ref.data(),
                                     mesh_dof.data(),
                                     mesh_l2g.data(),
                                     dV_ref.data(),
                                     u_trial_ref.data(),
                                     u_grad_trial_ref.data(),
                                     u_test_ref.data(),
                                     u_grad_test_ref.data(),
                                     mesh_trial_trace_ref.data(),
                                     mesh_grad_trial_trace_ref.data(),
                                     dS_ref.data(),
                                     u_trial_trace_ref.data(),
                                     u_grad_trial_trace_ref.data(),
                                     u_test_trace_ref.data(),
                                     u_grad_test_trace_ref.data(),
                                     normal_ref.data(),
                                     boundaryJac_ref.data(),
                                     nElements_global,
                                     useMetrics,
                                     epsFactHeaviside,
                                     epsFactDirac,
                                     epsFactDiffusion,
                                     u_l2g.data(),
                                     elementDiameter.data(),
                                     nodeDiametersArray.data(),
                                     u_dof.data(),
                                     q_phi.data(),
                                     q_normal_phi.data(),
                                     ebqe_phi.data(),
                                     ebqe_normal_phi.data(),
                                     q_H.data(),
                                     q_u.data(),
                                     q_n.data(),
                                     ebqe_u.data(),
                                     ebqe_n.data(),
                                     q_r.data(),
                                     q_porosity.data(),
                                     offset_u,stride_u,
                                     elementResidual_u,
                                     nExteriorElementBoundaries_global,
                                     exteriorElementBoundariesArray.data(),
                                     elementBoundaryElementsArray.data(),
                                     elementBoundaryLocalElementBoundariesArray.data(),
                                     element_u,
                                     eN);
            //compute l2 norm
            double resNorm=0.0;
            for (int i=0;i<nDOF_test_element;i++)
              {
                resNorm += elementResidual_u[i];
              }//i
            resNorm = fabs(resNorm);
            //now do Newton
            int its=0;
            //std::cout<<"element "<<eN<<std::endl;
            //std::cout<<"resNorm0 "<<resNorm<<std::endl;
            while (resNorm  >= atol && its < maxIts)
              {
                its+=1;
                calculateElementJacobian(mesh_trial_ref.data(),
                                         mesh_grad_trial_ref.data(),
                                         mesh_dof.data(),
                                         mesh_l2g.data(),
                                         dV_ref.data(),
                                         u_trial_ref.data(),
                                         u_grad_trial_ref.data(),
                                         u_test_ref.data(),
                                         u_grad_test_ref.data(),
                                         mesh_trial_trace_ref.data(),
                                         mesh_grad_trial_trace_ref.data(),
                                         dS_ref.data(),
                                         u_trial_trace_ref.data(),
                                         u_grad_trial_trace_ref.data(),
                                         u_test_trace_ref.data(),
                                         u_grad_test_trace_ref.data(),
                                         normal_ref.data(),
                                         boundaryJac_ref.data(),
                                         nElements_global,
                                         useMetrics,
                                         epsFactHeaviside,
                                         epsFactDirac,
                                         epsFactDiffusion,
                                         u_l2g.data(),
                                         elementDiameter.data(),
                                         nodeDiametersArray.data(),
                                         u_dof.data(),
                                         q_phi.data(),
                                         q_normal_phi.data(),
                                         q_H.data(),
                                         q_porosity.data(),
                                         elementJacobian_u_u,
                                         element_u,
                                         eN);
                for (int i=0;i<nDOF_test_element;i++)
                  {
                    element_du[i] = -elementResidual_u[i];
                    elementPivots[i] = ((PROTEUS_LAPACK_INTEGER)0);
                    elementColPivots[i]=((PROTEUS_LAPACK_INTEGER)0);
                    /* std::cout<<"element jacobian"<<std::endl; */
                    /* for (int j=0;j<nDOF_test_element;j++) */
                    /*   { */
                    /*     std::cout<<elementJacobian_u_u[i*nDOF_trial_element+j]<<'\t'; */
                    /*   } */
                    /* std::cout<<std::endl; */
                  }//i
                //factor
                PROTEUS_LAPACK_INTEGER La_N=((PROTEUS_LAPACK_INTEGER)nDOF_test_element),
                  INFO=0;
                dgetc2_(&La_N,
                        elementJacobian_u_u,
                        &La_N,
                        elementPivots,
                        elementColPivots,
                        &INFO);
                //solve
                dgesc2_(&La_N,
                        elementJacobian_u_u,
                        &La_N,
                        element_du,
                        elementPivots,
                        elementColPivots,
                        &scale);
                double resNormNew = resNorm,lambda=1.0;
                int lsIts=0;
                while (resNormNew > 0.99*resNorm && lsIts < 100)
                  {
                    //apply correction
                    for (int i=0;i<nDOF_test_element;i++)
                      {
                        element_u[i] += lambda*element_du[i];
                      }//i
                    lambda /= 2.0;
                    //compute new residual
                    calculateElementResidual(mesh_trial_ref.data(),
                                             mesh_grad_trial_ref.data(),
                                             mesh_dof.data(),
                                             mesh_l2g.data(),
                                             dV_ref.data(),
                                             u_trial_ref.data(),
                                             u_grad_trial_ref.data(),
                                             u_test_ref.data(),
                                             u_grad_test_ref.data(),
                                             mesh_trial_trace_ref.data(),
                                             mesh_grad_trial_trace_ref.data(),
                                             dS_ref.data(),
                                             u_trial_trace_ref.data(),
                                             u_grad_trial_trace_ref.data(),
                                             u_test_trace_ref.data(),
                                             u_grad_test_trace_ref.data(),
                                             normal_ref.data(),
                                             boundaryJac_ref.data(),
                                             nElements_global,
                                             useMetrics,
                                             epsFactHeaviside,
                                             epsFactDirac,
                                             epsFactDiffusion,
                                             u_l2g.data(),
                                             elementDiameter.data(),
                                             nodeDiametersArray.data(),
                                             u_dof.data(),
                                             q_phi.data(),
                                             q_normal_phi.data(),
                                             ebqe_phi.data(),
                                             ebqe_normal_phi.data(),
                                             q_H.data(),
                                             q_u.data(),
                                             q_n.data(),
                                             ebqe_u.data(),
                                             ebqe_n.data(),
                                             q_r.data(),
                                             q_porosity.data(),
                                             offset_u,stride_u,
                                             elementResidual_u,
                                             nExteriorElementBoundaries_global,
                                             exteriorElementBoundariesArray.data(),
                                             elementBoundaryElementsArray.data(),
                                             elementBoundaryLocalElementBoundariesArray.data(),
                                             element_u,
                                             eN);
                    lsIts +=1;
                    //compute l2 norm
                    resNormNew=0.0;
                    for (int i=0;i<nDOF_test_element;i++)
                      {
                        resNormNew += elementResidual_u[i];
                        std::cout<<"element_u["<<i<<"] "<<element_u[i]<<std::endl;
                        std::cout<<"elementResidual_u["<<i<<"] "<<elementResidual_u[i]<<std::endl;
                      }//i
                    resNormNew = fabs(resNormNew);
                  }
                resNorm = resNormNew;
                std::cout<<"INFO "<<INFO<<std::endl;
                std::cout<<"resNorm["<<its<<"] "<<resNorm<<std::endl;
              }
          }//elements
      }
      void elementConstantSolve(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        int maxIts = args.scalar<int>("maxIts");
        double atol = args.scalar<double>("atol");
        for(int eN=0;eN<nElements_global;eN++)
          {
            //declare local storage for element residual and initialize
            register double element_u[nDOF_test_element],elementConstant_u,
              elementResidual_u[nDOF_test_element],elementConstantResidual,
              elementJacobian_u_u[nDOF_test_element*nDOF_trial_element],elementConstantJacobian,resNorm;
            elementConstant_u=0.0;
            for (int i=0;i<nDOF_test_element;i++)
              {
                element_u[i]=elementConstant_u;
              }//i
            calculateElementResidual(mesh_trial_ref.data(),
                                     mesh_grad_trial_ref.data(),
                                     mesh_dof.data(),
                                     mesh_l2g.data(),
                                     dV_ref.data(),
                                     u_trial_ref.data(),
                                     u_grad_trial_ref.data(),
                                     u_test_ref.data(),
                                     u_grad_test_ref.data(),
                                     mesh_trial_trace_ref.data(),
                                     mesh_grad_trial_trace_ref.data(),
                                     dS_ref.data(),
                                     u_trial_trace_ref.data(),
                                     u_grad_trial_trace_ref.data(),
                                     u_test_trace_ref.data(),
                                     u_grad_test_trace_ref.data(),
                                     normal_ref.data(),
                                     boundaryJac_ref.data(),
                                     nElements_global,
                                     useMetrics,
                                     epsFactHeaviside,
                                     epsFactDirac,
                                     epsFactDiffusion,
                                     u_l2g.data(),
                                     elementDiameter.data(),
                                     nodeDiametersArray.data(),
                                     u_dof.data(),
                                     q_phi.data(),
                                     q_normal_phi.data(),
                                     ebqe_phi.data(),
                                     ebqe_normal_phi.data(),
                                     q_H.data(),
                                     q_u.data(),
                                     q_n.data(),
                                     ebqe_u.data(),
                                     ebqe_n.data(),
                                     q_r.data(),
                                     q_porosity.data(),
                                     offset_u,stride_u,
                                     elementResidual_u,
                                     nExteriorElementBoundaries_global,
                                     exteriorElementBoundariesArray.data(),
                                     elementBoundaryElementsArray.data(),
                                     elementBoundaryLocalElementBoundariesArray.data(),
                                     element_u,
                                     eN);
            //compute l2 norm
            elementConstantResidual=0.0;
            for (int i=0;i<nDOF_test_element;i++)
              {
                elementConstantResidual += elementResidual_u[i];
              }//i
            resNorm = fabs(elementConstantResidual);
            //now do Newton
            int its=0;
            //std::cout<<"element "<<eN<<std::endl;
            //std::cout<<"resNorm0 "<<resNorm<<std::endl;
            while (resNorm >= atol && its < maxIts)
              {
                its+=1;
                calculateElementJacobian(mesh_trial_ref.data(),
                                         mesh_grad_trial_ref.data(),
                                         mesh_dof.data(),
                                         mesh_l2g.data(),
                                         dV_ref.data(),
                                         u_trial_ref.data(),
                                         u_grad_trial_ref.data(),
                                         u_test_ref.data(),
                                         u_grad_test_ref.data(),
                                         mesh_trial_trace_ref.data(),
                                         mesh_grad_trial_trace_ref.data(),
                                         dS_ref.data(),
                                         u_trial_trace_ref.data(),
                                         u_grad_trial_trace_ref.data(),
                                         u_test_trace_ref.data(),
                                         u_grad_test_trace_ref.data(),
                                         normal_ref.data(),
                                         boundaryJac_ref.data(),
                                         nElements_global,
                                         useMetrics,
                                         epsFactHeaviside,
                                         epsFactDirac,
                                         epsFactDiffusion,
                                         u_l2g.data(),
                                         elementDiameter.data(),
                                         nodeDiametersArray.data(),
                                         u_dof.data(),
                                         q_phi.data(),
                                         q_normal_phi.data(),
                                         q_H.data(),
                                         q_porosity.data(),
                                         elementJacobian_u_u,
                                         element_u,
                                         eN);
                elementConstantJacobian=0.0;
                for (int i=0;i<nDOF_test_element;i++)
                  {
                    for (int j=0;j<nDOF_test_element;j++)
                      {
                        elementConstantJacobian += elementJacobian_u_u[i*nDOF_trial_element+j];
                      }
                  }//i
                std::cout<<"elementConstantJacobian "<<elementConstantJacobian<<std::endl;
                //apply correction
                elementConstant_u -= elementConstantResidual/(elementConstantJacobian+1.0e-8);
                for (int i=0;i<nDOF_test_element;i++)
                  {
                    element_u[i] = elementConstant_u;
                  }//i
                //compute new residual
                calculateElementResidual(mesh_trial_ref.data(),
                                         mesh_grad_trial_ref.data(),
                                         mesh_dof.data(),
                                         mesh_l2g.data(),
                                         dV_ref.data(),
                                         u_trial_ref.data(),
                                         u_grad_trial_ref.data(),
                                         u_test_ref.data(),
                                         u_grad_test_ref.data(),
                                         mesh_trial_trace_ref.data(),
                                         mesh_grad_trial_trace_ref.data(),
                                         dS_ref.data(),
                                         u_trial_trace_ref.data(),
                                         u_grad_trial_trace_ref.data(),
                                         u_test_trace_ref.data(),
                                         u_grad_test_trace_ref.data(),
                                         normal_ref.data(),
                                         boundaryJac_ref.data(),
                                         nElements_global,
                                         useMetrics,
                                         epsFactHeaviside,
                                         epsFactDirac,
                                         epsFactDiffusion,
                                         u_l2g.data(),
                                         elementDiameter.data(),
                                         nodeDiametersArray.data(),
                                         u_dof.data(),
                                         q_phi.data(),
                                         q_normal_phi.data(),
                                         ebqe_phi.data(),
                                         ebqe_normal_phi.data(),
                                         q_H.data(),
                                         q_u.data(),
                                         q_n.data(),
                                         ebqe_u.data(),
                                         ebqe_n.data(),
                                         q_r.data(),
                                         q_porosity.data(),
                                         offset_u,stride_u,
                                         elementResidual_u,
                                         nExteriorElementBoundaries_global,
                                         exteriorElementBoundariesArray.data(),
                                         elementBoundaryElementsArray.data(),
                                         elementBoundaryLocalElementBoundariesArray.data(),
                                         element_u,
                                         eN);
                //compute l2 norm
                elementConstantResidual=0.0;
                for (int i=0;i<nDOF_test_element;i++)
                  {
                    elementConstantResidual += elementResidual_u[i];
                  }//i
                resNorm = fabs(elementConstantResidual);
                std::cout<<"resNorm["<<its<<"] "<<resNorm<<std::endl;
              }
          }//elements
      }

      std::tuple<double, double> globalConstantRJ(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
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
        int nElements_owned = args.scalar<int>("nElements_owned");
        double useMetrics = args.scalar<double>("useMetrics");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        int maxIts = args.scalar<int>("maxIts");
        double atol = args.scalar<double>("atol");
        double constant_u = args.scalar<double>("constant_u");
        register double element_u[nDOF_test_element],
          elementResidual_u[nDOF_test_element],
          elementJacobian_u_u[nDOF_test_element*nDOF_trial_element];
        double constantResidual = 0.0;
        double constantJacobian = 0.0;
        for (int i=0;i<nDOF_trial_element;i++)
          {
            element_u[i]=constant_u;
          }//i
        //compute residual and Jacobian
        for(int eN=0;eN<nElements_owned;eN++)
          {
            calculateElementResidual(mesh_trial_ref.data(),
                                     mesh_grad_trial_ref.data(),
                                     mesh_dof.data(),
                                     mesh_l2g.data(),
                                     dV_ref.data(),
                                     u_trial_ref.data(),
                                     u_grad_trial_ref.data(),
                                     u_test_ref.data(),
                                     u_grad_test_ref.data(),
                                     mesh_trial_trace_ref.data(),
                                     mesh_grad_trial_trace_ref.data(),
                                     dS_ref.data(),
                                     u_trial_trace_ref.data(),
                                     u_grad_trial_trace_ref.data(),
                                     u_test_trace_ref.data(),
                                     u_grad_test_trace_ref.data(),
                                     normal_ref.data(),
                                     boundaryJac_ref.data(),
                                     nElements_owned,
                                     useMetrics,
                                     epsFactHeaviside,
                                     epsFactDirac,
                                     epsFactDiffusion,
                                     u_l2g.data(),
                                     elementDiameter.data(),
                                     nodeDiametersArray.data(),
                                     u_dof.data(),
                                     q_phi.data(),
                                     q_normal_phi.data(),
                                     ebqe_phi.data(),
                                     ebqe_normal_phi.data(),
                                     q_H.data(),
                                     q_u.data(),
                                     q_n.data(),
                                     ebqe_u.data(),
                                     ebqe_n.data(),
                                     q_r.data(),
                                     q_porosity.data(),
                                     offset_u,stride_u,
                                     elementResidual_u,
                                     nExteriorElementBoundaries_global,
                                     exteriorElementBoundariesArray.data(),
                                     elementBoundaryElementsArray.data(),
                                     elementBoundaryLocalElementBoundariesArray.data(),
                                     element_u,
                                     eN);
            //compute l2 norm
            for (int i=0;i<nDOF_test_element;i++)
              {
                constantResidual += elementResidual_u[i];
              }//i
            calculateElementJacobian(mesh_trial_ref.data(),
                                     mesh_grad_trial_ref.data(),
                                     mesh_dof.data(),
                                     mesh_l2g.data(),
                                     dV_ref.data(),
                                     u_trial_ref.data(),
                                     u_grad_trial_ref.data(),
                                     u_test_ref.data(),
                                     u_grad_test_ref.data(),
                                     mesh_trial_trace_ref.data(),
                                     mesh_grad_trial_trace_ref.data(),
                                     dS_ref.data(),
                                     u_trial_trace_ref.data(),
                                     u_grad_trial_trace_ref.data(),
                                     u_test_trace_ref.data(),
                                     u_grad_test_trace_ref.data(),
                                     normal_ref.data(),
                                     boundaryJac_ref.data(),
                                     nElements_owned,
                                     useMetrics,
                                     epsFactHeaviside,
                                     epsFactDirac,
                                     epsFactDiffusion,
                                     u_l2g.data(),
                                     elementDiameter.data(),
                                     nodeDiametersArray.data(),
                                     u_dof.data(),
                                     q_phi.data(),
                                     q_normal_phi.data(),
                                     q_H.data(),
                                     q_porosity.data(),
                                     elementJacobian_u_u,
                                     element_u,
                                     eN);
            for (int i=0;i<nDOF_test_element;i++)
              {
                for (int j=0;j<nDOF_test_element;j++)
                  {
                    constantJacobian += elementJacobian_u_u[i*nDOF_trial_element+j];
                  }
              }//i
          }
          return std::tuple<double, double>(constantResidual, constantJacobian);
      }

      double calculateMass(arguments_dict& args,
                           bool useExact)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& x_ref = args.array<double>("x_ref");
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
        int nElements_owned = args.scalar<int>("nElements_owned");
        double useMetrics = args.scalar<double>("useMetrics");
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& phi_dof = args.array<double>("phi_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        double globalMass = 0.0;
        gf.useExact=useExact;
        for(int eN=0;eN<nElements_owned;eN++)
          {
            double epsHeaviside;
            //loop over quadrature points and compute integrands
            //declare local storage for element residual and initialize
            register double element_phi[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i=eN*nDOF_test_element+i;
                element_phi[i] = phi_dof.data()[u_l2g.data()[eN_i]];
              }//i
            double element_nodes[nDOF_mesh_trial_element*3];
            for (int i=0;i<nDOF_mesh_trial_element;i++)
              {
                register int eN_i=eN*nDOF_mesh_trial_element+i;
                for(int I=0;I<3;I++)
                  element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
	      }//i
            gf.calculate(element_phi, element_nodes, x_ref.data(),false);
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace;
		//eN_nDOF_trial_element = eN*nDOF_trial_element;
		//double u=0.0,grad_u[nSpace],r=0.0,dr=0.0;
		double jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  //u_grad_trial[nDOF_trial_element*nSpace],
		  //u_test_dV[nDOF_trial_element],
		  //u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z,
		  G[nSpace*nSpace],G_dd_G,tr_G,h_phi;
                gf.set_quad(k);
                //
                //compute solution and gradients at quadrature points
                //
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
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.calculateG(jacInv,G,G_dd_G,tr_G);
                /* double dir[nSpace]; */
                /* double norm = 1.0e-8; */
                /* for (int I=0;I<nSpace;I++) */
                /*        norm += q_normal_phi.data()[eN_k_nSpace+I]*q_normal_phi.data()[eN_k_nSpace+I]; */
                /* norm = sqrt(norm); */
                /* for (int I=0;I<nSpace;I++) */
                /*        dir[I] = q_normal_phi.data()[eN_k_nSpace+I]/norm; */
                /* ck.calculateGScale(G,dir,h_phi); */
                epsHeaviside=epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                globalMass += q_porosity[eN_k]*gf.H(epsHeaviside,q_phi.data()[eN_k])*dV;
              }//k
          }//elements
          return globalMass;
      }

      void setMassQuadrature(arguments_dict& args,
                             bool useExact)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& x_ref = args.array<double>("x_ref");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& phi_l2g = args.array<int>("phi_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& phi_dof = args.array<double>("phi_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& H_dof = args.array<double>("H_dof");
        gf.useExact=useExact;
        gf_nodes.useExact=useExact;
        for(int eN=0;eN<nElements_global;eN++)
          {
            double epsHeaviside;
            //loop over quadrature points and compute integrands
            //declare local storage for element residual and initialize
            register double element_phi[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i=eN*nDOF_test_element+i;
                element_phi[i] = phi_dof.data()[phi_l2g.data()[eN_i]];
              }//i
            double element_nodes[nDOF_mesh_trial_element*3];
            for (int i=0;i<nDOF_mesh_trial_element;i++)
              {
                register int eN_i=eN*nDOF_mesh_trial_element+i;
                for(int I=0;I<3;I++)
                  element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
	      }//i
            gf.calculate(element_phi, element_nodes, x_ref.data(),false);
            gf_nodes.calculate(element_phi, element_nodes, element_nodes,false);
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace;
		//eN_nDOF_trial_element = eN*nDOF_trial_element;
		//register double u=0.0,grad_u[nSpace],r=0.0,dr=0.0;
		register double jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  //u_grad_trial[nDOF_trial_element*nSpace],
		  //u_test_dV[nDOF_trial_element],
		  //u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z,
		  G[nSpace*nSpace],G_dd_G,tr_G,h_phi;
		//
		//compute solution and gradients at quadrature points
		//
                gf.set_quad(k);
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
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.calculateG(jacInv,G,G_dd_G,tr_G);

                /* double dir[nSpace]; */
                /* double norm = 1.0e-8; */
                /* for (int I=0;I<nSpace;I++) */
                /*        norm += q_normal_phi.data()[eN_k_nSpace+I]*q_normal_phi.data()[eN_k_nSpace+I]; */
                /* norm = sqrt(norm); */
                /* for (int I=0;I<nSpace;I++) */
                /*        dir[I] = q_normal_phi.data()[eN_k_nSpace+I]/norm; */

                /* ck.calculateGScale(G,dir,h_phi); */
                epsHeaviside=epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                q_H.data()[eN_k] = gf.H(epsHeaviside,q_phi.data()[eN_k]);
              }//k
            // distribute rhs for mass correction
            for (int i=0;i<nDOF_trial_element;i++)
              {
                gf_nodes.set_quad(i);
                int eN_i = eN*nDOF_trial_element + i;
                int gi = phi_l2g.data()[eN_i];
                epsHeaviside = epsFactHeaviside*nodeDiametersArray.data()[mesh_l2g.data()[eN_i]];//cek hack, only works if isoparametric, but we can fix by including interpolation points
                H_dof.data() [gi] = gf_nodes.H(epsHeaviside,phi_dof.data()[gi]);
              }
          }//elements
      }

      void FCTStep(arguments_dict& args)
      {
        int NNZ = args.scalar<int>("NNZ");
        int numDOFs = args.scalar<int>("numDOFs");
        xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
        xt::pyarray<double>& solH = args.array<double>("solH");
        xt::pyarray<double>& solL = args.array<double>("solL");
        xt::pyarray<double>& limited_solution = args.array<double>("limited_solution");
        xt::pyarray<int>& csrRowIndeces_DofLoops = args.array<int>("csrRowIndeces_DofLoops");
        xt::pyarray<int>& csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops");
        xt::pyarray<double>& MassMatrix = args.array<double>("matrix");
        Rpos.resize(numDOFs,0.0), Rneg.resize(numDOFs,0.0);
        FluxCorrectionMatrix.resize(NNZ,0.0);
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

            double mini=0., maxi=1.0;
            double Pposi=0, Pnegi=0;
            // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
            for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
              {
                int j = csrColumnOffsets_DofLoops.data()[offset];
                // i-th row of flux correction matrix
                FluxCorrectionMatrix[ij] = ((i==j ? 1. : 0.)*mi - MassMatrix.data()[ij]) * (solH.data()[j]-solHi);

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
            limited_solution.data()[i] = fmax(0.0,solL.data()[i] + 1./lumped_mass_matrix.data()[i]*ith_Limiter_times_FluxCorrectionMatrix);
          }
      }

      // mql. copied from calculateElementJacobian. NOTE: there are some not necessary computations!!!
      inline void calculateElementMassMatrix(//element
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
                                             double* q_porosity,
                                             double* elementMassMatrix,
                                             double* elementLumpedMassMatrix,
                                             double* element_u,
                                             int eN)
      {
        for (int i=0;i<nDOF_test_element;i++)
          {
            elementLumpedMassMatrix[i] = 0.0;
            for (int j=0;j<nDOF_trial_element;j++)
              {
                elementMassMatrix[i*nDOF_trial_element+j]=0.0;
              }
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
                                 q_porosity[eN_k],
                                 r,
                                 dr);
            for(int i=0;i<nDOF_test_element;i++)
              {
                //int eN_k_i=eN_k*nDOF_test_element+i;
                //int eN_k_i_nSpace=eN_k_i*nSpace;
                elementLumpedMassMatrix[i] += u_test_dV[i];
                for(int j=0;j<nDOF_trial_element;j++)
                  {
                    elementMassMatrix[i*nDOF_trial_element+j] += u_trial_ref[k*nDOF_trial_element+j]*u_test_dV[i];
                  }//j
              }//i
          }//k
      }

      void calculateMassMatrix(arguments_dict& args)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& globalMassMatrix = args.array<double>("globalMassMatrix");
        xt::pyarray<double>& globalLumpedMassMatrix = args.array<double>("globalLumpedMassMatrix");
        //
        //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
        //
        for(int eN=0;eN<nElements_global;eN++)
          {
            register double  elementMassMatrix[nDOF_test_element*nDOF_trial_element],element_u[nDOF_trial_element], elementLumpedMassMatrix[nDOF_trial_element];
            for (int j=0;j<nDOF_trial_element;j++)
              {
                register int eN_j = eN*nDOF_trial_element+j;
                element_u[j] = u_dof.data()[u_l2g.data()[eN_j]];
              }
            calculateElementMassMatrix(mesh_trial_ref.data(),
                                       mesh_grad_trial_ref.data(),
                                       mesh_dof.data(),
                                       mesh_l2g.data(),
                                       dV_ref.data(),
                                       u_trial_ref.data(),
                                       u_grad_trial_ref.data(),
                                       u_test_ref.data(),
                                       u_grad_test_ref.data(),
                                       mesh_trial_trace_ref.data(),
                                       mesh_grad_trial_trace_ref.data(),
                                       dS_ref.data(),
                                       u_trial_trace_ref.data(),
                                       u_grad_trial_trace_ref.data(),
                                       u_test_trace_ref.data(),
                                       u_grad_test_trace_ref.data(),
                                       normal_ref.data(),
                                       boundaryJac_ref.data(),
                                       nElements_global,
                                       useMetrics,
                                       epsFactHeaviside,
                                       epsFactDirac,
                                       epsFactDiffusion,
                                       u_l2g.data(),
                                       elementDiameter.data(),
                                       nodeDiametersArray.data(),
                                       u_dof.data(),
                                       q_phi.data(),
                                       q_normal_phi.data(),
                                       q_H.data(),
                                       q_porosity.data(),
                                       elementMassMatrix,
                                       elementLumpedMassMatrix,
                                       element_u,
                                       eN);
            //
            //load into element Jacobian into global Jacobian
            //
            for (int i=0;i<nDOF_test_element;i++)
              {
                int eN_i = eN*nDOF_test_element+i;
                int gi = u_l2g.data()[eN_i];
                globalLumpedMassMatrix.data()[gi] += elementLumpedMassMatrix[i];
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    int eN_i_j = eN_i*nDOF_trial_element+j;
                    globalMassMatrix.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] +=
                      elementMassMatrix[i*nDOF_trial_element+j];
                  }//j
              }//i
          }//elements
      }//calculate mass matrix

      void setMassQuadratureEdgeBasedStabilizationMethods(arguments_dict& args,
                                                          bool useExact)
      {
        xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
        xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
        xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
        xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
        xt::pyarray<double>& x_ref = args.array<double>("x_ref");
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
        double epsFactHeaviside = args.scalar<double>("epsFactHeaviside");
        double epsFactDirac = args.scalar<double>("epsFactDirac");
        double epsFactDiffusion = args.scalar<double>("epsFactDiffusion");
        xt::pyarray<int>& phi_l2g = args.array<int>("phi_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& nodeDiametersArray = args.array<double>("nodeDiametersArray");
        xt::pyarray<double>& phi_dof = args.array<double>("phi_dof");
        xt::pyarray<double>& q_phi = args.array<double>("q_phi");
        xt::pyarray<double>& q_normal_phi = args.array<double>("q_normal_phi");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        xt::pyarray<double>& ebqe_normal_phi = args.array<double>("ebqe_normal_phi");
        xt::pyarray<double>& q_H = args.array<double>("q_H");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_n = args.array<double>("q_n");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_n = args.array<double>("ebqe_n");
        xt::pyarray<double>& q_r = args.array<double>("q_r");
        xt::pyarray<double>& q_porosity = args.array<double>("q_porosity");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& rhs_mass_correction = args.array<double>("rhs_mass_correction");
        xt::pyarray<double>& lumped_L2p_vof_mass_correction = args.array<double>("lumped_L2p_vof_mass_correction");
        xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
        int numDOFs = args.scalar<int>("numDOFs");
        gf.useExact=useExact;
        for(int eN=0;eN<nElements_global;eN++)
          {
            register double element_rhs_mass_correction[nDOF_test_element];
            for (int i=0;i<nDOF_test_element;i++)
              element_rhs_mass_correction[i] = 0.;
            double epsHeaviside;
            //loop over quadrature points and compute integrands
            double element_phi[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i=eN*nDOF_test_element+i;
                element_phi[i] = phi_dof.data()[phi_l2g.data()[eN_i]];
              }//i
            double element_nodes[nDOF_mesh_trial_element*3];
            for (int i=0;i<nDOF_mesh_trial_element;i++)
              {
                register int eN_i=eN*nDOF_mesh_trial_element+i;
                for(int I=0;I<3;I++)
                  element_nodes[i*3 + I] = mesh_dof.data()[mesh_l2g.data()[eN_i]*3 + I];
	      }//i
            gf.calculate(element_phi, element_nodes, x_ref.data(),false);
	    for  (int k=0;k<nQuadraturePoints_element;k++)
	      {
		//compute indeces and declare local storage
		register int eN_k = eN*nQuadraturePoints_element+k,
		  eN_k_nSpace = eN_k*nSpace;
		//eN_nDOF_trial_element = eN*nDOF_trial_element;
		//register double u=0.0,grad_u[nSpace],r=0.0,dr=0.0;
		register double jac[nSpace*nSpace],
		  jacDet,
		  jacInv[nSpace*nSpace],
		  //u_grad_trial[nDOF_trial_element*nSpace],
		  //u_test_dV[nDOF_trial_element],
		  //u_grad_test_dV[nDOF_test_element*nSpace],
		  dV,x,y,z,
		  u_test_dV[nDOF_test_element],
		  G[nSpace*nSpace],G_dd_G,tr_G,h_phi;
		gf.set_quad(k);
                //
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
                //get the physical integration weight
                dV = fabs(jacDet)*dV_ref.data()[k];
                ck.calculateG(jacInv,G,G_dd_G,tr_G);

                // precalculate test function times integration weight
                for (int j=0;j<nDOF_trial_element;j++)
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;

                /* double dir[nSpace]; */
                /* double norm = 1.0e-8; */
                /* for (int I=0;I<nSpace;I++) */
                /*        norm += q_normal_phi.data()[eN_k_nSpace+I]*q_normal_phi.data()[eN_k_nSpace+I]; */
                /* norm = sqrt(norm); */
                /* for (int I=0;I<nSpace;I++) */
                /*        dir[I] = q_normal_phi.data()[eN_k_nSpace+I]/norm; */

                /* ck.calculateGScale(G,dir,h_phi); */
                epsHeaviside=epsFactHeaviside*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter.data()[eN]);
                q_H.data()[eN_k] = gf.H(epsHeaviside,q_phi.data()[eN_k]);

                for (int i=0;i<nDOF_trial_element;i++)
                  element_rhs_mass_correction [i] += q_porosity.data()[eN_k]*q_H.data()[eN_k]*u_test_dV[i];
              }//k
            // distribute rhs for mass correction
            for (int i=0;i<nDOF_trial_element;i++)
              {
                int eN_i = eN*nDOF_trial_element + i;
                int gi = phi_l2g.data()[eN_i];
                rhs_mass_correction.data()[gi] += element_rhs_mass_correction[i];
              }
          }//elements
        // COMPUTE LUMPED L2 PROYJECTION
        for (int i=0; i<numDOFs; i++)
          {
            double mi = lumped_mass_matrix.data()[i];
            lumped_L2p_vof_mass_correction.data()[i] = 1./mi*rhs_mass_correction.data()[i];
          }
      }
    };//MCorr

  inline MCorr_base* newMCorr(int nSpaceIn,
                              int nQuadraturePoints_elementIn,
                              int nDOF_mesh_trial_elementIn,
                              int nDOF_trial_elementIn,
                              int nDOF_test_elementIn,
                              int nQuadraturePoints_elementBoundaryIn,
                              int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<MCorr_base,MCorr,CompKernel>(nSpaceIn,
                                                                                     nQuadraturePoints_elementIn,
                                                                                     nDOF_mesh_trial_elementIn,
                                                                                     nDOF_trial_elementIn,
                                                                                     nDOF_test_elementIn,
                                                                                     nQuadraturePoints_elementBoundaryIn,
                                                                                     CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<MCorr_base,MCorr,CompKernel>(nSpaceIn,
                                                                                   nQuadraturePoints_elementIn,
                                                                                   nDOF_mesh_trial_elementIn,
                                                                                   nDOF_trial_elementIn,
                                                                                   nDOF_test_elementIn,
                                                                                   nQuadraturePoints_elementBoundaryIn,
                                                                                   CompKernelFlag);
  }
}//proteus
#endif
