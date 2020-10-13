#ifndef Richards_H
#define Richards_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "../mprans/ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"
#define nnz nSpace

namespace py = pybind11;

namespace proteus
{
    class Richards_base
    {
        //The base class defining the interface
        public:
        virtual ~Richards_base() {}
        virtual void calculateResidual(arguments_dict& args)=0;
        virtual void calculateJacobian(arguments_dict& args)=0;
    };

    template<class CompKernelType,
             int nSpace,
             int nQuadraturePoints_element,
             int nDOF_mesh_trial_element,
             int nDOF_trial_element,
             int nDOF_test_element,
             int nQuadraturePoints_elementBoundary>
    class Richards : public Richards_base
    {
        public:
        const int nDOF_test_X_trial_element;
        CompKernelType ck;
        Richards():
            nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
            ck()
        {}
        inline
        void evaluateCoefficients(const int rowptr[nSpace],
                                  const int colind[nnz],
                                  const double rho,
                                  const double beta,
                                  const double gravity[nSpace],
                                  const double alpha,
                                  const double n_vg,
                                  const double thetaR,
                                  const double thetaSR,
                                  const double KWs[nnz],
                                  const double& u,
                                  double& m,
                                  double& dm,
                                  double f[nSpace],
                                  double df[nSpace],
                                  double a[nnz],
                                  double da[nnz])
        {
            const int nSpace2 = nSpace * nSpace;
            double psiC;
            double pcBar;
            double pcBar_n;
            double pcBar_nM1;
            double pcBar_nM2;
            double onePlus_pcBar_n;
            double sBar;
            double sqrt_sBar;
            double DsBar_DpsiC;
            double thetaW;
            double DthetaW_DpsiC;
            double vBar;
            double vBar2;
            double DvBar_DpsiC;
            double KWr;
            double DKWr_DpsiC;
            double rho2 = rho * rho;
            double thetaS;
            double rhom;
            double drhom;
            double m_vg;
            double pcBarStar;
            double sqrt_sBarStar;

            psiC = -u;
            m_vg = 1.0 - 1.0 / n_vg;
            thetaS = thetaR + thetaSR;
            if (psiC > 0.0)
            {
                pcBar = alpha * psiC;
                pcBarStar = pcBar;
                if (pcBar < 1.0e-8)
                    pcBarStar = 1.0e-8;
                pcBar_nM2 = pow(pcBarStar, n_vg - 2);
                pcBar_nM1 = pcBar_nM2 * pcBar;
                pcBar_n   = pcBar_nM1 * pcBar;
                onePlus_pcBar_n = 1.0 + pcBar_n;

                sBar = pow(onePlus_pcBar_n, -m_vg);
                /* using -mn = 1-n */
                DsBar_DpsiC = alpha * (1.0 - n_vg) * (sBar/onePlus_pcBar_n)*pcBar_nM1;

                vBar = 1.0-pcBar_nM1*sBar;
                vBar2 = vBar*vBar;
                DvBar_DpsiC = -alpha*(n_vg-1.0)*pcBar_nM2*sBar - pcBar_nM1*DsBar_DpsiC;

                thetaW = thetaSR*sBar + thetaR;
                DthetaW_DpsiC = thetaSR * DsBar_DpsiC;

                sqrt_sBar = sqrt(sBar);
                sqrt_sBarStar = sqrt_sBar;
                if (sqrt_sBar < 1.0e-8)
                    sqrt_sBarStar = 1.0e-8;
                KWr= sqrt_sBar*vBar2;
                DKWr_DpsiC= ((0.5/sqrt_sBarStar)*DsBar_DpsiC*vBar2
                        +
                        2.0*sqrt_sBar*vBar*DvBar_DpsiC);
            }
            else
            {
                thetaW        = thetaS;
                DthetaW_DpsiC = 0.0;
                KWr           = 1.0;
                DKWr_DpsiC    = 0.0;
            }
            //slight compressibility
            rhom = rho*exp(beta*u);
            drhom = beta*rhom;
            m = rhom*thetaW;
            dm = -rhom*DthetaW_DpsiC+drhom*thetaW;
            for (int I=0;I<nSpace;I++)
            {
                f[I] = 0.0;
                df[I] = 0.0;
                for (int ii=rowptr[I]; ii < rowptr[I+1]; ii++)
                {
                    f[I]  += rho2*KWr*KWs[ii]*gravity[colind[ii]];
                    df[I] += -rho2*DKWr_DpsiC*KWs[ii]*gravity[colind[ii]];
                    a[ii]  = rho*KWr*KWs[ii];
                    da[ii] = -rho*DKWr_DpsiC*KWs[ii];
                }
            }
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
        void calculateSubgridError_tau(const double&  Ct_sge,
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
        void exteriorNumericalFlux(const double& bc_flux,
                                   int rowptr[nSpace],
                                   int colind[nnz],
                                   int isSeepageFace,
                                   int& isDOFBoundary,
                                   double n[nSpace],
                                   double bc_u,
                                   double K[nnz],
                                   double grad_psi[nSpace],
                                   double u,
                                   double K_rho_g[nSpace],
                                   double penalty,
                                   double& flux)
        {
            double v_I,bc_u_seepage=0.0;
            if (isSeepageFace || isDOFBoundary)
            {
                flux = 0.0;
                for(int I=0;I<nSpace;I++)
                {
                    //gravity
                    v_I = K_rho_g[I];
                    //pressure head
                    for(int m=rowptr[I];m<rowptr[I+1];m++)
                    {
                        v_I -= K[m]*grad_psi[colind[m]];
                    }
                    flux += v_I*n[I];
                }
                if (isSeepageFace)
                    bc_u = bc_u_seepage;
                flux += penalty*(u-bc_u);
                if (isSeepageFace)
                {
                    if (flux > 0.0)
                    {
                        isDOFBoundary = 1;
                        bc_u = bc_u_seepage;
                    }
                    else
                    {
                        isDOFBoundary = 0;
                        flux = 0.0;
                    }
                }
                /* //set DOF flag and flux correctly if seepage face */
                /* if (isSeepageFace) */
                /*   { */
                /*     if (flux < 0.0 || u < bc_u_seepage) */
                /*     { */
                /*       isDOFBoundary = 0; */
                /*       flux = 0.0; */
                /*     } */
                /*     else */
                /*     { */
                /*       isDOFBoundary = 1; */
                /*       bc_u = bc_u_seepage; */
                /*     } */
                /*   } */
                /* //Dirichlet penalty */
                /* if (isDOFBoundary) */
                /*   flux += penalty*(u-bc_u); */
            }
            else
                flux = bc_flux;
        }

        void exteriorNumericalFluxJacobian(const int rowptr[nSpace],
                                           const int colind[nnz],
                                           const int isDOFBoundary,
                                           const double n[nSpace],
                                           const double K[nnz],
                                           const double dK[nnz],
                                           const double grad_psi[nSpace],
                                           const double grad_v[nSpace],
                                           const double dK_rho_g[nSpace],
                                           const double v,
                                           const double penalty,
                                           double& fluxJacobian)
        {
            if (isDOFBoundary)
            {
                fluxJacobian = 0.0;
                for(int I=0;I<nSpace;I++)
                {
                    //gravity
                    fluxJacobian += dK_rho_g[I]*v*n[I];
                    //pressure head
                    for(int m=rowptr[I]; m<rowptr[I+1]; m++)
                    {
                        fluxJacobian -= (K[m]*grad_v[colind[m]] + dK[m]*v*grad_psi[colind[m]])*n[I];
                    }
                }
                //Dirichlet penalty
                fluxJacobian += penalty*v;
            }
            else
                fluxJacobian = 0.0;
        }

        void calculateResidual(arguments_dict& args)
        {
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
        xt::pyarray<double>& ebqe_penalty_ext = args.array<double>("ebqe_penalty_ext");
        xt::pyarray<int>& elementMaterialTypes = args.array<int>("elementMaterialTypes");
        xt::pyarray<int>& isSeepageFace = args.array<int>("isSeepageFace");
        xt::pyarray<int>& a_rowptr = args.array<int>("a_rowptr");
        xt::pyarray<int>& a_colind = args.array<int>("a_colind");
        double rho = args.scalar<double>("rho");
        double beta = args.scalar<double>("beta");
        xt::pyarray<double>& gravity = args.array<double>("gravity");
        xt::pyarray<double>& alpha = args.array<double>("alpha");
        xt::pyarray<double>& n = args.array<double>("n");
        xt::pyarray<double>& thetaR = args.array<double>("thetaR");
        xt::pyarray<double>& thetaSR = args.array<double>("thetaSR");
        xt::pyarray<double>& KWs = args.array<double>("KWs");
        double useMetrics = args.scalar<double>("useMetrics");
        double alphaBDF = args.scalar<double>("alphaBDF");
        int lag_shockCapturing = args.scalar<int>("lag_shockCapturing");
        double shockCapturingDiffusion = args.scalar<double>("shockCapturingDiffusion");
        double sc_uref = args.scalar<double>("sc_uref");
        double sc_alpha = args.scalar<double>("sc_alpha");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& u_dof_old = args.array<double>("u_dof_old");
        xt::pyarray<double>& velocity = args.array<double>("velocity");
        xt::pyarray<double>& q_m = args.array<double>("q_m");
        xt::pyarray<double>& q_u = args.array<double>("q_u");
        xt::pyarray<double>& q_dV = args.array<double>("q_dV");
        xt::pyarray<double>& q_m_betaBDF = args.array<double>("q_m_betaBDF");
        xt::pyarray<double>& cfl = args.array<double>("cfl");
        xt::pyarray<double>& q_numDiff_u = args.array<double>("q_numDiff_u");
        xt::pyarray<double>& q_numDiff_u_last = args.array<double>("q_numDiff_u_last");
        int offset_u = args.scalar<int>("offset_u");
        int stride_u = args.scalar<int>("stride_u");
        xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
        xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
        xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
        xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
        xt::pyarray<double>& ebqe_bc_flux_ext = args.array<double>("ebqe_bc_flux_ext");
        xt::pyarray<double>& ebqe_phi = args.array<double>("ebqe_phi");
        double epsFact = args.scalar<double>("epsFact");
        xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
        xt::pyarray<double>& ebqe_flux = args.array<double>("ebqe_flux");
            assert(a_rowptr.data()[nSpace] == nnz);
            assert(a_rowptr.data()[nSpace] == nSpace);
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
                    register double u=0.0,grad_u[nSpace],grad_u_old[nSpace],
                    m=0.0,dm=0.0,
                    f[nSpace],df[nSpace],
                    a[nnz],da[nnz],
                    m_t=0.0,dm_t=0.0,
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
                    ck.calculateMappingVelocity_element(eN,
                            k,
                            mesh_velocity_dof.data(),
                            mesh_l2g.data(),
                            mesh_trial_ref.data(),
                            xt,yt,zt);
                    //get the physical integration weight
                    dV = fabs(jacDet)*dV_ref.data()[k];
                    q_dV.data()[eN_k] = dV;
                    ck.calculateG(jacInv,G,G_dd_G,tr_G);
                    //get the trial function gradients
                    ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                    //get the solution
                    ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],u);
                    //get the solution gradients
                    ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_u);
                    //precalculate test function products with integration weights
                    for (int j=0;j<nDOF_trial_element;j++)
                    {
                        u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                        for (int I=0;I<nSpace;I++)
                        {
                            u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                        }
                    }
                    //
                    //calculate pde coefficients at quadrature points
                    //
                    evaluateCoefficients(a_rowptr.data(),
                            a_colind.data(),
                            rho,
                            beta,
                            gravity.data(),
                            alpha.data()[elementMaterialTypes.data()[eN]],
                            n.data()[elementMaterialTypes.data()[eN]],
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]],
                            &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                            u,
                            m,
                            dm,
                            f,
                            df,
                            a,
                            da);
                    //
                    //calculate time derivative at quadrature points
                    //
                    ck.bdf(alphaBDF,
                            q_m_betaBDF.data()[eN_k],
                            m,
                            dm,
                            m_t,
                            dm_t);
                    /* // */
                    /* //calculate subgrid error (strong residual and adjoint) */
                    /* // */
                    /* //calculate strong residual */
                    /* pdeResidual_u = ck.Mass_strong(m_t) + ck.Advection_strong(df,grad_u); */
                    /* //calculate adjoint */
                    /* for (int i=0;i<nDOF_test_element;i++) */
                    /*     { */
                    /*       // register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace; */
                    /*       // Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]); */
                    /*       register int i_nSpace = i*nSpace; */
                    /*       Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]); */
                    /*     } */
                    /* //calculate tau and tau*Res */
                    /* calculateSubgridError_tau(elementDiameter[eN],dm_t,df,cfl[eN_k],tau0); */
                    /* calculateSubgridError_tau(Ct_sge, */
                    /*                           G, */
                    /*                           dm_t, */
                    /*                           df, */
                    /*                           tau1, */
                    /*                           cfl[eN_k]); */

                    /* tau = useMetrics*tau1+(1.0-useMetrics)*tau0; */

                    /* subgridError_u = -tau*pdeResidual_u; */
                    /* // */
                    /* //calculate shock capturing diffusion */
                    /* // */


                    /* ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u,grad_u,numDiff0);       */
                    /* //ck.calculateNumericalDiffusion(shockCapturingDiffusion,G,pdeResidual_u,grad_u_old,numDiff1); */
                    /* ck.calculateNumericalDiffusion(shockCapturingDiffusion,sc_uref, sc_alpha,G,G_dd_G,pdeResidual_u,grad_u,numDiff1); */
                    /* q_numDiff_u[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0; */
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
                            ck.Diffusion_weak(a_rowptr.data(),a_colind.data(),a,grad_u,&u_grad_test_dV[i_nSpace]);
                        /* +  */
                        /*   ck.SubgridError(subgridError_u,Lstar_u[i]) +  */
                        /*   ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[i_nSpace]);  */
                    }//i
                    //
                    q_m.data()[eN_k] = m;
                    q_u.data()[eN_k] = u;
                }
                //
                //load element into global residual and save element residual
                //
                for(int i=0;i<nDOF_test_element;i++)
                {
                    register int eN_i=eN*nDOF_test_element+i;

                    globalResidual.data()[offset_u+stride_u*u_l2g.data()[eN_i]] += elementResidual_u[i];
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
                             register double u_ext=0.0,
                                      grad_u_ext[nSpace],
                                      m_ext=0.0,
                                      dm_ext=0.0,
                                      f_ext[nSpace],
                                      df_ext[nSpace],
                                      a_ext[nnz],
                                      da_ext[nnz],
                                      flux_ext=0.0,
                                      bc_u_ext=0.0,
                                      bc_grad_u_ext[nSpace],
                                      bc_m_ext=0.0,
                                      bc_dm_ext=0.0,
                                      bc_f_ext[nSpace],
                                      bc_df_ext[nSpace],
                                      bc_a_ext[nnz],
                                      bc_da_ext[nnz],
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
                                      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
                                      //get the metric tensor
                                      //cek todo use symmetry
                                      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                                      //compute shape and solution information
                                      //shape
                                      ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                                      //solution and gradient
                                      ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],u_ext);
                                      ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
                                      //precalculate test function products with integration weights
                                      for (int j=0;j<nDOF_trial_element;j++)
                                      {
                                          u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                                      }
                                      //
                                      //load the boundary values
                                      //
                                      bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*ebqe_bc_u_ext.data()[ebNE_kb]+(1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;
                                      //
                                      //calculate the pde coefficients using the solution and the boundary values for the solution
                                      //
                                      evaluateCoefficients(a_rowptr.data(),
                                              a_colind.data(),
                                              rho,
                                              beta,
                                              gravity.data(),
                                              alpha.data()[elementMaterialTypes.data()[eN]],
                                              n.data()[elementMaterialTypes.data()[eN]],
                                              thetaR.data()[elementMaterialTypes.data()[eN]],
                                              thetaSR.data()[elementMaterialTypes.data()[eN]],
                                              &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                                              u_ext,
                                              m_ext,
                                              dm_ext,
                                              f_ext,
                                              df_ext,
                                              a_ext,
                                              da_ext);
                                      evaluateCoefficients(a_rowptr.data(),
                                              a_colind.data(),
                                              rho,
                                              beta,
                                              gravity.data(),
                                              alpha.data()[elementMaterialTypes.data()[eN]],
                                              n.data()[elementMaterialTypes.data()[eN]],
                                              thetaR.data()[elementMaterialTypes.data()[eN]],
                                              thetaSR.data()[elementMaterialTypes.data()[eN]],
                                              &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                                              bc_u_ext,
                                              bc_m_ext,
                                              bc_dm_ext,
                                              bc_f_ext,
                                              bc_df_ext,
                                              bc_a_ext,
                                              bc_da_ext);
                                      //
                                      //calculate the numerical fluxes
                                      //
                                      exteriorNumericalFlux(ebqe_bc_flux_ext[ebNE_kb],
                                              a_rowptr.data(),
                                              a_colind.data(),
                                              isSeepageFace.data()[ebNE],//tricky, this is a face flag not face quad
                                              isDOFBoundary_u.data()[ebNE_kb],
                                              normal,
                                              bc_u_ext,
                                              a_ext,
                                              grad_u_ext,
                                              u_ext,
                                              f_ext,
                                              ebqe_penalty_ext.data()[ebNE_kb],// penalty,
                                              flux_ext);
                                      ebqe_flux.data()[ebNE_kb] = flux_ext;
                                      ebqe_u.data()[ebNE_kb] = u_ext;
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

                             globalResidual.data()[offset_u+stride_u*u_l2g.data()[eN_i]] += elementResidual_u[i];
                         }//i
            }//ebNE
        }

        void calculateJacobian(arguments_dict& args)
        {
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
        xt::pyarray<double>& ebqe_penalty_ext = args.array<double>("ebqe_penalty_ext");
        xt::pyarray<int>& elementMaterialTypes = args.array<int>("elementMaterialTypes");
        xt::pyarray<int>& isSeepageFace = args.array<int>("isSeepageFace");
        xt::pyarray<int>& a_rowptr = args.array<int>("a_rowptr");
        xt::pyarray<int>& a_colind = args.array<int>("a_colind");
        double rho = args.scalar<double>("rho");
        double beta = args.scalar<double>("beta");
        xt::pyarray<double>& gravity = args.array<double>("gravity");
        xt::pyarray<double>& alpha = args.array<double>("alpha");
        xt::pyarray<double>& n = args.array<double>("n");
        xt::pyarray<double>& thetaR = args.array<double>("thetaR");
        xt::pyarray<double>& thetaSR = args.array<double>("thetaSR");
        xt::pyarray<double>& KWs = args.array<double>("KWs");
        double useMetrics = args.scalar<double>("useMetrics");
        double alphaBDF = args.scalar<double>("alphaBDF");
        int lag_shockCapturing = args.scalar<int>("lag_shockCapturing");
        double shockCapturingDiffusion = args.scalar<double>("shockCapturingDiffusion");
        xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
        xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
        xt::pyarray<double>& u_dof = args.array<double>("u_dof");
        xt::pyarray<double>& velocity = args.array<double>("velocity");
        xt::pyarray<double>& q_m_betaBDF = args.array<double>("q_m_betaBDF");
        xt::pyarray<double>& cfl = args.array<double>("cfl");
        xt::pyarray<double>& q_numDiff_u_last = args.array<double>("q_numDiff_u_last");
        xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
        xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
        xt::pyarray<double>& globalJacobian = args.array<double>("globalJacobian");
        int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
        xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
        xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
        xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
        xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
        xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
        xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
        xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
        xt::pyarray<double>& ebqe_bc_flux_ext = args.array<double>("ebqe_bc_flux_ext");
        xt::pyarray<int>& csrColumnOffsets_eb_u_u = args.array<int>("csrColumnOffsets_eb_u_u");
            assert(a_rowptr.data()[nSpace] == nnz);
            assert(a_rowptr.data()[nSpace] == nSpace);
            double Ct_sge = 4.0;

            //
            //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
            //
            for(int eN=0;eN<nElements_global;eN++)
            {
                register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
                for (int i=0;i<nDOF_test_element;i++)
                {
                    for (int j=0;j<nDOF_trial_element;j++)
                    {
                        elementJacobian_u_u[i][j]=0.0;
                    }
                }
                for (int k=0;k<nQuadraturePoints_element;k++)
                {
                    int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                        eN_k_nSpace = eN_k*nSpace,
                        eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

                    //declare local storage
                    register double u=0.0,
                             grad_u[nSpace],
                             m=0.0,dm=0.0,
                             f[nSpace],df[nSpace],
                             a[nnz],da[nnz],
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
                            mesh_dof.data(),
                            mesh_l2g.data(),
                            mesh_trial_ref.data(),
                            mesh_grad_trial_ref.data(),
                            jac,
                            jacDet,
                            jacInv,
                            x,y,z);
                    ck.calculateMappingVelocity_element(eN,
                            k,
                            mesh_velocity_dof.data(),
                            mesh_l2g.data(),
                            mesh_trial_ref.data(),
                            xt,yt,zt);
                    //get the physical integration weight
                    dV = fabs(jacDet)*dV_ref.data()[k];
                    ck.calculateG(jacInv,G,G_dd_G,tr_G);
                    //get the trial function gradients
                    ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
                    //get the solution
                    ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],u);
                    //get the solution gradients
                    ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_u);
                    //precalculate test function products with integration weights
                    for (int j=0;j<nDOF_trial_element;j++)
                    {
                        u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                        for (int I=0;I<nSpace;I++)
                        {
                            u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                        }
                    }
                    //
                    //calculate pde coefficients and derivatives at quadrature points
                    //
                    evaluateCoefficients(a_rowptr.data(),
                            a_colind.data(),
                            rho,
                            beta,
                            gravity.data(),
                            alpha.data()[elementMaterialTypes.data()[eN]],
                            n.data()[elementMaterialTypes.data()[eN]],
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]],
                            &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                            u,
                            m,
                            dm,
                            f,
                            df,
                            a,
                            da);
                    //
                    //calculate time derivatives
                    //
                    ck.bdf(alphaBDF,
                            q_m_betaBDF.data()[eN_k],
                            m,
                            dm,
                            m_t,
                            dm_t);
                    // //
                    // //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
                    // //
                    // //calculate the adjoint times the test functions
                    // for (int i=0;i<nDOF_test_element;i++)
                    //     {
                    //       // int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
                    //       // Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
                    //       register int i_nSpace = i*nSpace;
                    //       Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
                    //     }
                    // //calculate the Jacobian of strong residual
                    // for (int j=0;j<nDOF_trial_element;j++)
                    //     {
                    //       //int eN_k_j=eN_k*nDOF_trial_element+j;
                    //       //int eN_k_j_nSpace = eN_k_j*nSpace;
                    //       int j_nSpace = j*nSpace;
                    //       dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
                    //         ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);
                    //     }
                    // //tau and tau*Res
                    // calculateSubgridError_tau(elementDiameter[eN],
                    //                 dm_t,
                    //                 df,
                    //                 cfl[eN_k],
                    //                 tau0);
                    // calculateSubgridError_tau(Ct_sge,
                    //                           G,
                    //                 dm_t,
                    //                 df,
                    //                 tau1,
                    //                     cfl[eN_k]);
                    // tau = useMetrics*tau1+(1.0-useMetrics)*tau0;
                    // for(int j=0;j<nDOF_trial_element;j++)
                    //     dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
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
                            elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dm_t,u_trial_ref.data()[k*nDOF_trial_element+j],u_test_dV[i]) +
                                ck.AdvectionJacobian_weak(df,u_trial_ref.data()[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
                                ck.DiffusionJacobian_weak(a_rowptr.data(),a_colind.data(),a,da,
                                        grad_u,&u_grad_test_dV[i_nSpace],1.0,
                                        u_trial_ref.data()[k*nDOF_trial_element+j],&u_grad_trial[j_nSpace]);
                            // +
                            //     ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) +
                            //     ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]);
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
                        globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
                    }//j
                }//i
            }//elements
            //
            //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
            //
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
                             a_ext[nnz],
                             da_ext[nnz],
                             dflux_u_u_ext=0.0,
                             bc_u_ext=0.0,
                             //bc_grad_u_ext[nSpace],
                             bc_m_ext=0.0,
                             bc_dm_ext=0.0,
                             bc_f_ext[nSpace],
                             bc_df_ext[nSpace],
                             bc_a_ext[nnz],
                             bc_da_ext[nnz],
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
                             normal[3],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                             G[nSpace*nSpace],G_dd_G,tr_G;
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
                    dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
                    //dS = metricTensorDetSqrt*dS_ref.data()[kb];
                    ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                    //compute shape and solution information
                    //shape
                    ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                    //solution and gradients
                    ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],u_ext);
                    ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
                    //precalculate test function products with integration weights
                    for (int j=0;j<nDOF_trial_element;j++)
                    {
                        u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                    }
                    //
                    //load the boundary values
                    //
                    bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*ebqe_bc_u_ext.data()[ebNE_kb]+(1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;
                    //
                    //calculate the internal and external trace of the pde coefficients
                    //
                    evaluateCoefficients(a_rowptr.data(),
                            a_colind.data(),
                            rho,
                            beta,
                            gravity.data(),
                            alpha.data()[elementMaterialTypes.data()[eN]],
                            n.data()[elementMaterialTypes.data()[eN]],
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]],
                            &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                            u_ext,
                            m_ext,
                            dm_ext,
                            f_ext,
                            df_ext,
                            a_ext,
                            da_ext);
                    evaluateCoefficients(a_rowptr.data(),
                            a_colind.data(),
                            rho,
                            beta,
                            gravity.data(),
                            alpha.data()[elementMaterialTypes.data()[eN]],
                            n.data()[elementMaterialTypes.data()[eN]],
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]],
                            &KWs.data()[elementMaterialTypes.data()[eN]*nnz],
                            bc_u_ext,
                            bc_m_ext,
                            bc_dm_ext,
                            bc_f_ext,
                            bc_df_ext,
                            bc_a_ext,
                            bc_da_ext);
                    //
                    //calculate the flux jacobian
                    //
                    for (int j=0;j<nDOF_trial_element;j++)
                    {
                        //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
                        register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
                        exteriorNumericalFluxJacobian(a_rowptr.data(),
                                a_colind.data(),
                                isDOFBoundary_u.data()[ebNE_kb],
                                normal,
                                a_ext,
                                da_ext,
                                grad_u_ext,
                                &u_grad_trial_trace[j*nSpace],
                                df_ext,
                                u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element+j],
                                ebqe_penalty_ext.data()[ebNE_kb],//penalty,
                                fluxJacobian_u_u[j]);
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
                            globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS[i];
                        }//j
                    }//i
                }//kb
            }//ebNE
        }//computeJacobian
    };//Richards

    inline Richards_base* newRichards(int nSpaceIn,
                                      int nQuadraturePoints_elementIn,
                                      int nDOF_mesh_trial_elementIn,
                                      int nDOF_trial_elementIn,
                                      int nDOF_test_elementIn,
                                      int nQuadraturePoints_elementBoundaryIn,
                                      int CompKernelFlag)
    {
        return proteus::chooseAndAllocateDiscretization<Richards_base,Richards,CompKernel>(nSpaceIn,
                                                                                           nQuadraturePoints_elementIn,
                                                                                           nDOF_mesh_trial_elementIn,
                                                                                           nDOF_trial_elementIn,
                                                                                           nDOF_test_elementIn,
                                                                                           nQuadraturePoints_elementBoundaryIn,
                                                                                           CompKernelFlag);
    }
}//proteus
#endif
