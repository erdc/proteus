#ifndef PRESINC_H
#define PRESINC_H
#include <cmath>
#include <iostream>
#include <cstring>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
class cppMy_correction_base
{
public:
    virtual ~cppMy_correction_base(){}
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
            double* p_trial_ref,
            double* p_grad_trial_ref,
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
            double* u_dof,
            double* v_dof,
            double* u_old_dof,
            double* v_old_dof,
            int* p_l2g,
            int nDofs_p_per_ele,
            double* p_dof,
            int offset_u,
            int stride_u,
            int offset_v,
            int stride_v,
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
            double* p_trial_ref,
            double* p_grad_trial_ref,
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
            double* u_dof,
            double* v_dof,
            double* u_old_dof,
            double* v_old_dof,
            int* p_l2g,
            int nDofs_p_per_ele,
            double* p_dof,
            int* csrRowIndeces_u_u,
            int* csrColumnOffsets_u_u,
            int* csrRowIndeces_u_v,
            int* csrColumnOffsets_u_v,
            int* csrRowIndeces_v_u,
            int* csrColumnOffsets_v_u,
            int* csrRowIndeces_v_v,
            int* csrColumnOffsets_v_v,
            double* globalJacobian,
            int nExteriorElementBoundaries_global,
            int* exteriorElementBoundariesArray,
            int* elementBoundaryElementsArray,
            int* elementBoundaryLocalElementBoundariesArray)=0;
};

template<class CompKernelType,
int nSpace,
int nQuadraturePoints_element,
int nDOF_mesh_trial_element,
int nDOF_trial_element,
int nDOF_test_element,
int nQuadraturePoints_elementBoundary>
class cppMy_correction : public cppMy_correction_base
{
public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    cppMy_correction():
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
            double* p_trial_ref,
            double* p_grad_trial_ref,
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
            double* u_dof,
            double* v_dof,
            double* u_old_dof,
            double* v_old_dof,
            int* p_l2g,
            int nDofs_p_per_ele,
            double* p_dof,
            int offset_u,
            int stride_u,
            int offset_v,
            int stride_v,
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
            register double elementResidual_u[nDOF_test_element],elementResidual_v[nDOF_test_element],
                            element_u[nDOF_trial_element],element_v[nDOF_trial_element],
                            element_u_old[nDOF_trial_element],element_v_old[nDOF_trial_element];
            for (int i=0;i<nDOF_test_element;i++)
            {
                register int eN_i=eN*nDOF_test_element+i;
                element_u[i] = u_dof[u_l2g[eN_i]];
                element_v[i] = v_dof[u_l2g[eN_i]];
                element_u_old[i] = u_old_dof[u_l2g[eN_i]];
                element_v_old[i] = v_old_dof[u_l2g[eN_i]];
            }//i
            for(int i=0;i<nDOF_test_element;i++)
            {
                elementResidual_u[i]=0.0;
                elementResidual_v[i]=0.0;
            }//i
            //loop over quadrature points and compute integrands
            for  (int k=0;k<nQuadraturePoints_element;k++)
            {
                //compute indeces and declare local storage
                register int eN_k = eN*nQuadraturePoints_element+k,
                             eN_k_nSpace = eN_k*nSpace;
                //eN_nDOF_trial_element = eN*nDOF_trial_element;
                register double u,v,u_old,v_old,
                        grad_p[nSpace],
                        jac[nSpace*nSpace],
                        jacDet,
                        jacInv[nSpace*nSpace],
                        u_grad_trial[nDOF_trial_element*nSpace],
                        p_grad_trial[nDofs_p_per_ele*nSpace],
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
                ck.valFromElementDOF(element_v,&u_trial_ref[k*nDOF_trial_element],v);
                ck.valFromElementDOF(element_u_old,&u_trial_ref[k*nDOF_trial_element],u_old);
                ck.valFromElementDOF(element_v_old,&u_trial_ref[k*nDOF_trial_element],v_old);

                //compute divergence of velocity
                get_gradTrial_fromRef(nDofs_p_per_ele,
                                      nSpace,
                                      &p_grad_trial_ref[k*nDofs_p_per_ele*nSpace],
                                      jacInv,
                                      p_grad_trial);
                get_gradFromDOF(nDofs_p_per_ele,
                                nSpace,
                                &p_l2g[eN*nDofs_p_per_ele],
                                p_dof,
                                p_grad_trial,
                                grad_p);

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
                    elementResidual_u[i] += (u-u_old+grad_p[0])*phi_dV[i];
                    elementResidual_v[i] += (v-v_old+grad_p[1])*phi_dV[i];
                }
            }//k

            for(int i=0;i<nDOF_test_element;i++)
            {
                register int eN_i=eN*nDOF_test_element+i;
                globalResidual[offset_u+stride_u*u_l2g[eN_i]]+=elementResidual_u[i];
                globalResidual[offset_v+stride_v*u_l2g[eN_i]]+=elementResidual_v[i];

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
            double* p_trial_ref,
            double* p_grad_trial_ref,
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
            double* u_dof,
            double* v_dof,
            double* u_old_dof,
            double* v_old_dof,
            int* p_l2g,
            int nDofs_p_per_ele,
            double* p_dof,
            int* csrRowIndeces_u_u,
            int* csrColumnOffsets_u_u,
            int* csrRowIndeces_u_v,
            int* csrColumnOffsets_u_v,
            int* csrRowIndeces_v_u,
            int* csrColumnOffsets_v_u,
            int* csrRowIndeces_v_v,
            int* csrColumnOffsets_v_v,
            double* globalJacobian,
            int nExteriorElementBoundaries_global,
            int* exteriorElementBoundariesArray,
            int* elementBoundaryElementsArray,
            int* elementBoundaryLocalElementBoundariesArray)
    {
        //
        //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
        //
        for(int eN=0;eN<nElements_global;eN++)
        {
            register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
//            elementJacobian_u_v[nDOF_test_element*nDOF_trial_element],
//            elementJacobian_v_u[nDOF_test_element*nDOF_trial_element],
            elementJacobian_v_v[nDOF_test_element][nDOF_trial_element];

            std::memset(elementJacobian_u_u,0,nDOF_trial_element*nDOF_test_element*sizeof(double));
//            std::memset(elementJacobian_u_v,0,nDOF_trial_element*nDOF_test_element*sizeof(double));
//            std::memset(elementJacobian_v_u,0,nDOF_trial_element*nDOF_test_element*sizeof(double));
            std::memset(elementJacobian_v_v,0,nDOF_trial_element*nDOF_test_element*sizeof(double));

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
                        phi[nDOF_test_element],
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
                    phi[j]    = u_test_ref[k*nDOF_trial_element+j];
                    phi_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
                    for (int I=0;I<nSpace;I++)
                    {
                        grad_phi_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;
                    }
                }
                for(int i=0;i<nDOF_test_element;i++)
                {
                    for(int j=0;j<nDOF_trial_element;j++)
                    {
                        elementJacobian_u_u[i][j] += phi[j]*phi_dV[i];
                        elementJacobian_v_v[i][j] += phi[j]*phi_dV[i];
                    }//j
                }//i
            }//k
            for (int i=0;i<nDOF_test_element;i++)
              {
                register int eN_i = eN*nDOF_test_element+i;
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    register int eN_i_j = eN_i*nDOF_trial_element+j;
                    globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
                    //globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += elementJacobian_u_v[i][j];
                    //globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += elementJacobian_v_u[i][j];
                    globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += elementJacobian_v_v[i][j];
                  }//j
              }//i
        }//elements

    }//computeJacobian
};

inline cppMy_correction_base* newMy_correction(int nSpaceIn,
        int nQuadraturePoints_elementIn,
        int nDOF_mesh_trial_elementIn,
        int nDOF_trial_elementIn,
        int nDOF_test_elementIn,
        int nQuadraturePoints_elementBoundaryIn,
        int CompKernelFlag)
{
    if (nSpaceIn == 2)
        return proteus::chooseAndAllocateDiscretization2D<cppMy_correction_base,cppMy_correction,CompKernel>(nSpaceIn,
                nQuadraturePoints_elementIn,
                nDOF_mesh_trial_elementIn,
                nDOF_trial_elementIn,
                nDOF_test_elementIn,
                nQuadraturePoints_elementBoundaryIn,
                CompKernelFlag);
    else
        std::cout<<"not work for nSpaceIn=3"<<std::endl;
}
}//proteus
#endif
