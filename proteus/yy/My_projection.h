#ifndef PRESINC_H
#define PRESINC_H
#include <cmath>
#include <iostream>
#include <cstring>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
class cppMy_projection_base
{
public:
    virtual ~cppMy_projection_base(){}
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
            double* momentum_trial_ref,
            double* momentum_grad_trial_ref,
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
            int* velocity_l2g,
            int nDofs_velocity_per_ele,
            double* velocity_x_dof,
            double* velocity_y_dof,
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
            double* momentum_trial_ref,
            double* momentum_grad_trial_ref,
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
            int* velocity_l2g,
            int nDofs_velocity_per_ele,
            double* velocity_x_dof,
            double* velocity_y_dof,
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
class cppMy_projection : public cppMy_projection_base
{
public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    cppMy_projection():
        nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
        ck()
    {}

    inline void get_gradTrial_fromRef(const int ndofs, const int nd, const double* grad_trial_ref, const double* jacInv, double* grad_trial)
    {
        std::memset(grad_trial,0,ndofs*nd*sizeof(double));

        for (int j=0;j<ndofs;j++)
            for(int I=0;I<nd;I++)
                for(int J=0;J<nd;J++)
                    grad_trial[j*nd+I] += jacInv[J*nd+I]*grad_trial_ref[j*nd+J];
    }

    inline void get_gradFromDOF(int ndofs, int nd, int* l2g_element,double* dof,double* grad_trial, double* grad)
    {
        std::memset(grad,0,nd*sizeof(double));
        for (int j=0;j<ndofs;j++)
            for(int I=0;I<nd;I++)
                grad[I] += dof[l2g_element[j]]*grad_trial[j*2+I];
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
            double* momentum_trial_ref,
            double* momentum_grad_trial_ref,
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
            int* velocity_l2g,
            int nDofs_velocity_per_ele,
            double* velocity_x_dof,
            double* velocity_y_dof,
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
            int* elementBoundaryLocalElementBoundariesArray)
    {
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
            for(int i=0;i<nDOF_test_element;i++)
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
                register double u=0.0,
                        grad_u[nSpace],grad_velocity_x[nSpace],grad_velocity_y[nSpace],
                        a=0.0,
                        f[nSpace],
                        jac[nSpace*nSpace],
                        jacDet,
                        jacInv[nSpace*nSpace],
                        u_grad_trial[nDOF_trial_element*nSpace],
                        velocity_grad_trial[nDofs_velocity_per_ele*nSpace],
                        phi_dV[nDOF_trial_element],
                        grad_phi_dV[nDOF_test_element*nSpace],
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

                //compute divergence of velocity
                get_gradTrial_fromRef(nDofs_velocity_per_ele,
                                      nSpace,
                                      &momentum_grad_trial_ref[k*nDofs_velocity_per_ele*nSpace],
                                      jacInv,
                                      velocity_grad_trial);
                get_gradFromDOF(nDofs_velocity_per_ele,
                                nSpace,
                                &velocity_l2g[eN*nDofs_velocity_per_ele],
                                velocity_x_dof,
                                velocity_grad_trial,
                                grad_velocity_x);

                get_gradFromDOF(nDofs_velocity_per_ele,
                                nSpace,
                                &velocity_l2g[eN*nDofs_velocity_per_ele],
                                velocity_y_dof,
                                velocity_grad_trial,
                                grad_velocity_y);

                for (int j=0;j<nDOF_trial_element;j++)
                {
                    phi_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                    {
                        grad_phi_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;
                    }
                }

                for(int i=0;i<nDOF_test_element;i++)
                {
                    register int i_nSpace=i*nSpace;

                    elementResidual_u[i] += ck.Hamiltonian_strong(grad_u, &grad_phi_dV[i_nSpace])
                            +(grad_velocity_x[0]+grad_velocity_y[1])*phi_dV[i];
                }
            }//k

            for(int i=0;i<nDOF_test_element;i++)
            {
                register int eN_i=eN*nDOF_test_element+i;
                globalResidual[offset_u+stride_u*u_l2g[eN_i]]+=elementResidual_u[i];
            }//i
        }//elements
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
            double* momentum_trial_ref,
            double* momentum_grad_trial_ref,
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
            int* isDOFBoundary,
            int* isFluxBoundary,
            int* u_l2g,
            double* u_dof,
            int* velocity_l2g,
            int nDofs_velocity_per_ele,
            double* velocity_x_dof,
            double* velocity_y_dof,
            int* csrRowIndeces_u_u,
            int* csrColumnOffsets_u_u,
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
            register double  elementJacobian_u_u[nDOF_test_element*nDOF_trial_element];

            std::memset(elementJacobian_u_u,0,nDOF_trial_element*nDOF_test_element*sizeof(double));

            for  (int k=0;k<nQuadraturePoints_element;k++)
            {
                int eN_k = eN*nQuadraturePoints_element+k,
                    eN_k_nSpace = eN_k*nSpace;

                //declare local storage
                register double u=0.0,
                        grad_u[nSpace],
                        jac[nSpace*nSpace],
                        jacDet,
                        jacInv[nSpace*nSpace],
                        u_grad_trial[nDOF_trial_element*nSpace],
                        dV,
                        phi_dV[nDOF_test_element],
                        grad_phi_dV[nDOF_test_element*nSpace],
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

                for (int j=0;j<nDOF_trial_element;j++)
                {
                    phi_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                    {
                        grad_phi_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
                    }
                }
                for(int i=0;i<nDOF_test_element;i++)
                {
                    int i_nSpace=i*nSpace;
                    for(int j=0;j<nDOF_trial_element;j++)
                    {
                        int j_nSpace = j*nSpace;
                        elementJacobian_u_u[i*nDOF_trial_element+j] +=
                                ck.HamiltonianJacobian_strong(&u_grad_trial[j_nSpace],&grad_phi_dV[i_nSpace]);
                    }//j
                }//i
            }//k
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
};

inline cppMy_projection_base* newMy_projection(int nSpaceIn,
        int nQuadraturePoints_elementIn,
        int nDOF_mesh_trial_elementIn,
        int nDOF_trial_elementIn,
        int nDOF_test_elementIn,
        int nQuadraturePoints_elementBoundaryIn,
        int CompKernelFlag)
{
    if (nSpaceIn == 2)
        return proteus::chooseAndAllocateDiscretization2D<cppMy_projection_base,cppMy_projection,CompKernel>(nSpaceIn,
                nQuadraturePoints_elementIn,
                nDOF_mesh_trial_elementIn,
                nDOF_trial_elementIn,
                nDOF_test_elementIn,
                nQuadraturePoints_elementBoundaryIn,
                CompKernelFlag);
    else
        return proteus::chooseAndAllocateDiscretization<cppMy_projection_base,cppMy_projection,CompKernel>(nSpaceIn,
                nQuadraturePoints_elementIn,
                nDOF_mesh_trial_elementIn,
                nDOF_trial_elementIn,
                nDOF_test_elementIn,
                nQuadraturePoints_elementBoundaryIn,
                CompKernelFlag);
}
}//proteus
#endif
