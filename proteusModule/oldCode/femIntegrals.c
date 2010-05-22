#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include "femIntegrals.h"

/**
   \defgroup femIntegrals femIntegrals

   A library of functions for computing the discrete finite element formulations.
   @{
*/
   
/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the mass integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param m (nElements_global x nQuadraturePoints_element). The mass associated with each quadrature point (includes integration weight).
   @param w (nElements_global x  nQuadraturePoints_element  x nDOF_element). The finite element test function values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} + \int_{\Omega_e} \rho w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \rho w_i dV =  \int_{\Omega_r} \rho w_i |J_e| d\hat{V} \approx  \sum_k m_k w_i 
   \f]
   \f[
   m_k = \rho(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateMass(int nElements_global,
                int nQuadraturePoints_element,
                int nDOF_element,
                double* m,
                double* w,
                double* residual)
{
  int eN,i,k;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        residual[eN*nDOF_element + 
                 i] 
          += 
          m[eN*nQuadraturePoints_element + 
            k]
          *
          w[eN*nQuadraturePoints_element*nDOF_element + 
            k*nDOF_element + 
            i];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the mass integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param dm (nElements_global x nQuadraturePoints_element). The derivate of the mass associated with each quadrature point (includes integration weight).
   @param v_x_w (nElements_global x  nQuadraturePoints_element  x nDOF_element). The tensor product of the  finite element trial and test function values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} +\int_{\Omega_e} \frac{\partial \rho}{\partial u} v_j w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \frac{\partial \rho}{\partial u} v_j w_i dV =  \int_{\Omega_r} \frac{\partial \rho}{\partial u} v_j w_i |J_e| d\hat{V} \approx  \sum_k dm_k (v \otimes w)_{j,i} 
   \f]
   \f[
   dm_k = \frac{\partial \rho}{\partial u}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateMassJacobian(int nElements_global,
                        int nQuadraturePoints_element,
                        int nDOF_element,
                        double* dm,
                        double* v_x_w,
                        double* jacobian)
{
  int eN,i,j,k,nDOF_element2=nDOF_element*nDOF_element;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (j=0;j<nDOF_element;j++)
          jacobian[eN*nDOF_element2 + 
                   i*nDOF_element+
                   j] 
            += 
            dm[eN*nQuadraturePoints_element + 
               k]
            *
            v_x_w[eN*nQuadraturePoints_element*nDOF_element2 + 
                  k*nDOF_element2 + 
                  j*nDOF_element+
                  i];
}

/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the advection integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param f (nElements_global x nQuadraturePoints_element x nSpace). The advection associated with each quadrature point (includes integration weight).
   @param grad_w (nElements_global x  nQuadraturePoints_element  x nDOF_element x nSpace). The finite element test function gradient values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} - \int_{\Omega_e} \mathbf{f} \cdot \nabla w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \mathbf{f} \cdot \nabla w_i dV =  \int_{\Omega_r} \mathbf{f} \cdot \nabla w_i |J_e| d\hat{V} \int_{\Omega_r} \sum_k \mathbf{f}_k \cdot \nabla w_i 
   \f]
   \f[
   \mathbf{f}_k = \mathbf{f}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateAdvection(int nElements_global,
                     int nQuadraturePoints_element,
                     int nDOF_element,
                     int nSpace,
                     double* f,
                     double* grad_w,
                     double* residual)
{
  int eN,i,k,I;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (I=0;I<nSpace;I++)
          residual[eN*nDOF_element + 
                   i] 
            -=
            f[eN*nQuadraturePoints_element*nSpace + 
              k*nSpace + 
              I]
            *
            grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                   k*nDOF_element*nSpace + 
                   i*nSpace + 
                   I];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the advection integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param df (nElements_global x nQuadraturePoints_element x nSpace). The derivative of the advection associated with each quadrature point (includes integration weight).
   @param v_x_grad_w (nElements_global x  nQuadraturePoints_element  x nDOF_element x nSpace). The tensor product of the  finite element trial and test function gradient values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} - \int_{\Omega_e} \frac{\partial \mathbf{f}}{\partial u} v_j \cdot \nabla w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \frac{\partial \mathbf{f}}{\partial u} v_j \cdot \nabla w_i dV =  \int_{\Omega_r} \frac{\partial \mathbf{f}}{\partial u} v_j \cdot \nabla w_i |J_e| d\hat{V} \approx  \sum_k df_k \cdot (v \otimes \nabla w)_{j,i} 
   \f]
   \f[
   df_k = \frac{\partial \mathbf{f}}{\partial u}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateAdvectionJacobian(int nElements_global,
                             int nQuadraturePoints_element,
                             int nDOF_element,
                             int nSpace,
                             double* df,
                             double* v_x_grad_w,
                             double* jacobian)
{
  int eN,i,j,k,I,nDOF_element2=nDOF_element*nDOF_element;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            jacobian[eN*nDOF_element2 + 
                     i*nDOF_element+
                     j] 
              -=
              df[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace + 
                         k*nDOF_element2*nSpace + 
                         j*nDOF_element*nSpace + 
                         i*nSpace + 
                         I];
}

/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the diffusion integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param a (nElements_global x nQuadraturePoints_element x nSpace x nSpace). The diffusion coefficient associated with each quadrature point (includes integration weight).
   @param grad_phi_x_grad_w (nElements_global x  nQuadraturePoints_element  x nDOF_element x nSpace x  nSpace). The tensor product of potential gradient values and finite element test function gradient values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} + \int_{\Omega_e} \bar{\mathbf{a}} \nabla  \phi \cdot \nabla w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \bar{\mathbf{a}} \nabla \phi \cdot \nabla w_i dV =  \int_{\Omega_r} \bar{\mathbf{a}} \nabla \phi  \cdot \nabla w_i |J_e| d\hat{V} = \sum_k \bar{\mathbf{a}}_k \cdot (\nabla  \phi \otimes \nabla w)_i 
   \f]
   \f[
   \bar{\mathbf{a}}_k = \bar{\mathbf{a}}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateDiffusion(int nElements_global,
                     int nQuadraturePoints_element,
                     int nDOF_element,
                     int nSpace,
                     double* a,
                     double* grad_phi_x_grad_w,
                     double* residual)
{
  int eN,i,k,I,J,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (I=0;I<nSpace;I++)
          for (J=0;J<nSpace;J++)
            residual[eN*nDOF_element + i] 
              +=
              a[eN*nQuadraturePoints_element*nSpace2 + 
                k*nSpace2 + 
                I*nSpace + 
                J]
              *
              grad_phi_x_grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace2 + 
                                k*nDOF_element*nSpace2 + 
                                i*nSpace2+
                                J*nSpace + 
                                I];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the diffusion integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param l2g (nElements_global x nDOF_element) The mapping between element degrees of  freedom and global degrees of freedom.
   @param nSpace The number of spatial dimensions.
   @param a (nElements_global x nQuadraturePoints_element x nSpace x nSpace). The diffusion coefficient associated with each quadrature point (includes integration weight).
   @param da (nElements_global x nQuadraturePoints_element x nSpace x nSpace). The derivative of the diffusion coefficient associated with each quadrature point (includes integration weight).
   @param grad_phi_x_grad_w (nElements_global x  nQuadraturePoints_element  x nDOF_element x nSpace x nSpace). The tensor product of the  finite element trial and test function gradient values.
   @param dphi (nDOF_global). The global degrees of  freedom of phi.
   @param v (nElements_global x nQuadraturePoints_element x nDOF_element). The trial function values at each quadrature point.
   @param grad_v_x_grad_w (nElements_global x nQuadraturePoints_element x nDOF_element x nDOF_element x nSpace x  nSpace). The tensor product  of trial function gradient values and  test function gradient values at each quadrature point.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} + \int_{\Omega_e} (\frac{\partial \bar{\mathbf{a}}}{\partial u} v_j \nabla \phi + \bar{\mathbf{a}} \frac{\partial \phi}{\partial u} \nabla v_j ) \cdot \nabla w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \frac{\partial \bar{\mathbf{a}}}{\partial u} v_j \nabla \phi \cdot \nabla w_i dV =  \int_{\Omega_r} \frac{\partial \bar{\mathbf{a}}}{\partial u} v_j \nabla \phi \cdot \nabla w_i |J_e| d\hat{V} \approx  \sum_k da_k \cdot (\nabla \phi \otimes \nabla w)_{j,i} 
   \f]
   \f[
   \int_{\Omega_e} \bar{\mathbf{a}}\frac{\partial \phi}{\partial u} \nabla v_j \cdot \nabla w_i dV =  \int_{\Omega_r} \bar{\mathbf{a}} \frac{\partial \phi}{\partial u} \nabla v_j \cdot \nabla w_i |J_e| d\hat{V} \approx  \sum_k a_k dphi_{j} \cdot (\nabla v \otimes \nabla w)_{j,i} 
   \f]
   \f[
   a_k = \bar{\mathbf{a}}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
   \f[
   da_k = \frac{\partial \bar{\mathbf{a}}}{\partial u}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
   \f[
   dphi_j = dphi[l2g[e,j]]
   \f]
*/
void updateDiffusionJacobian(int nElements_global,
                             int nQuadraturePoints_element,
                             int nDOF_element,
                             int nSpace,
                             int* l2g,
                             double* a,
                             double* da,
                             double* grad_phi_x_grad_w,
                             double* dphi,
                             double* v,
                             double* grad_v_x_grad_w,
                             double* jacobian)
{
  int eN,i,j,k,I,J,nSpace2=nSpace*nSpace,nDOF_element2=nDOF_element*nDOF_element;
  double daProduct,dphiProduct;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          daProduct=0.0;
          for (I=0;I<nSpace;I++)
            for (J=0;J<nSpace;J++)
              daProduct 
                += 
                da[eN*nQuadraturePoints_element*nSpace2 +
                   k*nSpace2 +
                   I*nSpace +
                   J]
                *
                grad_phi_x_grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace2 +
                                  k*nDOF_element*nSpace2 +
                                  i*nSpace2+
                                  J*nSpace +
                                  I];
          for (j=0;j<nDOF_element;j++)
            {
              dphiProduct=0.0;
              for (I=0;I<nSpace;I++)
                for (J=0;J<nSpace;J++)
                  dphiProduct 
                    +=
                    a[eN*nQuadraturePoints_element*nSpace2 + 
                      k*nSpace2+
                      I*nSpace + 
                      J]
                    *
                    grad_v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace2 +
                                    k*nDOF_element2*nSpace2 +
                                    j*nDOF_element*nSpace2 + 
                                    i*nSpace2 + 
                                    J*nSpace + 
                                    I];
              jacobian[eN*nDOF_element2 + 
                       i*nDOF_element + 
                       j] 
                += 
                daProduct
                *
                v[eN*nQuadraturePoints_element*nDOF_element+
                  k*nDOF_element+
                  j] 
                +
                dphiProduct
                *
                dphi[l2g[eN*nDOF_element + 
                         j]];
            }
        }
}

/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the reaction integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param r (nElements_global x nQuadraturePoints_element). The reaction associated with each quadrature point (includes integration weight).
   @param w (nElements_global x  nQuadraturePoints_element  x nDOF_element). The finite element test function values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} + \int_{\Omega_e} r w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} r w_i dV =  \int_{\Omega_r} r w_i |J_e| d\hat{V} \approx  \sum_k r_k w_i 
   \f]
   \f[
   r_k = \rho(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateReaction(int nElements_global,
                    int nQuadraturePoints_element,
                    int nDOF_element,
                    double* r,
                    double* w,
                    double* residual)
{
  updateMass(nElements_global,nQuadraturePoints_element,nDOF_element,r,w,residual);
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the reaction integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param dr (nElements_global x nQuadraturePoints_element). The derivate of the reaction associated with each quadrature point (includes integration weight).
   @param v_x_w (nElements_global x  nQuadraturePoints_element  x nDOF_element). The tensor product of the  finite element trial and test function values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} +\int_{\Omega_e} \frac{\partial r}{\partial u} v_j w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \frac{\partial r}{\partial u} v_j w_i dV =  \int_{\Omega_r} \frac{\partial r}{\partial u} v_j w_i |J_e| d\hat{V} \approx  \sum_k dr_k (v \otimes w)_{j,i} 
   \f]
   \f[
   dr_k = \frac{\partial r}{\partial u}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateReactionJacobian(int nElements_global,
                            int nQuadraturePoints_element,
                            int nDOF_element,
                            double* dr,
                            double* v_x_w,
                            double* jacobian)
{
  updateMassJacobian(nElements_global,nQuadraturePoints_element,nDOF_element,dr,v_x_w,jacobian);
}
/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the stabilization integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param tau (nElements_global x nQuadraturePoints_element). The stabilization parameter associated with each quadrature point (includes integration weight).
   @param pdeResidual (nElements_global x nQuadraturePoints_element). The pde residual associated with each quadrature point.
   @param LstarW (nElements_global x  nQuadraturePoints_element  x nDOF_element). The linearized adjoint applied to the finite element test function values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} - \int_{\Omega_e} \tau \mathcal{R} \mathcal{L^*}(w_i) dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \tau \mathcal{R} \mathcal{L^*}(w_i) dV =  \int_{\Omega_r} \tau \mathcal{R} \mathcal{L^*}(w_i) |J_e| d\hat{V} \approx  \sum_k \tau_k \mathcal{R}_k \mathcal{L^*}(w_i)  
   \f]
   \f[
   \tau_k = \tau(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateStabilization(int nElements_global,
                         int nQuadraturePoints_element,
                         int nDOF_element,
                         double* tau,
                         double* pdeResidual,
                         double* LstarW,
                         double* residual)
{
  int eN,i,k;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        residual[eN*nDOF_element + i] 
          -= 
          tau[eN*nQuadraturePoints_element + 
              k]
          *
          pdeResidual[eN*nQuadraturePoints_element+
                      k]
          *
          LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                 k*nDOF_element + 
                 i];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the stabilization integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param tau (nElements_global x nQuadraturePoints_element). The stabilization parameter associated with each quadrature point (includes integration weight).
   @param dpdeResidual (nElements_global x nQuadraturePoints_element x nDOF_element). The Jacobian of the pde residual associated with each quadrature point.
   @param LstarW (nElements_global x  nQuadraturePoints_element  x nDOF_element). The linearized adjoint applied to the finite element test function values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} - \int_{\Omega_e} \tau \frac{\partial \mathcal{R}}{\partial u_j} \mathcal{L^*}(w_i) dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \tau \frac{\partial \mathcal{R}}{\partial u_j} \mathcal{L^*}(w_i) dV =  \int_{\Omega_r} \tau \frac{\partial \mathcal{R}}{\partial u_j} \mathcal{L^*}(w_i) |J_e| d\hat{V} \approx  \sum_k \tau_k (\frac{\partial \mathcal{R}}{\partial u_j})_k \mathcal{L^*}(w_i)  
   \f]
   \f[
   \tau_k = \tau(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateStabilizationJacobian(int nElements_global,
                                 int nQuadraturePoints_element,
                                 int nDOF_element,
                                 double* tau,
                                 double* dpdeResidual,
                                 double* LstarW,
                                 double* jacobian)
{
  int eN,i,j,k,nDOF_element2=nDOF_element*nDOF_element;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (j=0;j<nDOF_element;j++)
          jacobian[eN*nDOF_element2 + 
                   i*nDOF_element + 
                   j] 
            -= 
            tau[eN*nQuadraturePoints_element + 
                k]
            *
            dpdeResidual[eN*nQuadraturePoints_element*nDOF_element+
                         k*nDOF_element+
                         j]
            *
            LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                   k*nDOF_element + 
                   i];
}

/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the shock capturing integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param numDiff (nElements_global x nQuadraturePoints_element). The shock capturing diffusion associated with each quadrature point (includes integration weight).
   @param grad_u_x_grad_w (nElements_global x nQuadraturePoints_element x nDOF_element x nSpace  x nSpace). The tensor product of the solution gradient values and the test function gradient values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} + \int_{\Omega_e} \epsilon \nabla u \cdot  \nabla w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} \epsilon \nabla u \cdot  \nabla w_i dV =  \int_{\Omega_r} \epsilon \nabla u \cdot  \nabla w_i |J_e| d\hat{V} \approx  \sum_k numDiff_k \bar{\mathbf{I}} \cdot (\nabla u \otimes \nabla  w)_i  
   \f]
   \f[
   numDiff_k = \epsilon(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateShockCapturing(int nElements_global,
                          int nQuadraturePoints_element,
                          int nDOF_element,
                          int nSpace,
                          double* numDiff,
                          double* grad_u_x_grad_w,
                          double* residual)
{
  int eN,i,k,I,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (I=0;I<nSpace;I++)
          residual[eN*nDOF_element + i] 
            += 
            numDiff[eN*nQuadraturePoints_element + 
                    k]
            *
            grad_u_x_grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace2 + 
                            k*nDOF_element*nSpace2 + 
                            i*nSpace2 + 
                            I*nSpace + 
                            I];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the shock capturing integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param numDiff (nElements_global x nQuadraturePoints_element). The shock capturing diffusion associated with each quadrature point (includes integration weight).
   @param grad_v_x_grad_w (nElements_global x nQuadraturePoints_element x nDOF_element x nDOF_element x nSpace  x nSpace). The tensor product of the trial function gradient values and the test function gradient values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} + \int_{\Omega_e} \epsilon \nabla v_j \cdot  \nabla w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \epsilon \nabla v_j \cdot  \nabla w_i dV =  \int_{\Omega_r} \epsilon \nabla v_j \cdot  \nabla w_i |J_e| d\hat{V} \approx  \sum_k numDiff_k \bar{\mathbf{I}} \cdot (\nabla v \otimes \nabla  w)_{j,i}  
   \f]
   \f[
   numDiff_k = \epsilon(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateShockCapturingJacobian(int nElements_global,
                                  int nQuadraturePoints_element,
                                  int nDOF_element,
                                  int nSpace,
                                  double* numDiff,
                                  double* grad_v_x_grad_w,
                                  double* jacobian)
{
  int eN,i,j,k,I,nSpace2=nSpace*nSpace,nDOF_element2=nDOF_element*nDOF_element;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            jacobian[eN*nDOF_element2 + 
                     i*nDOF_element + 
                     j] 
              += 
              numDiff[eN*nQuadraturePoints_element + 
                      k]
              *
              grad_v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace2 + 
                              k*nDOF_element2*nSpace2 + 
                              j*nDOF_element*nSpace2 + 
                              i*nSpace2 + 
                              I*nSpace + 
                              I];
}

/**
   \brief Calculate the linearized adjoint of the scalar advection-diffusion-reaction equation
*/
void calculateAdjointADR(int nElements_global,
                         int nQuadraturePoints_element,
                         int nDOF_element,
                         int nSpace,
                         double* w,
                         double* grad_w,
                         double* df,
                         double* da,
                         double* grad_phi,
                         double* dr,
                         double* LstarW)
{
  int eN,i,k,I,J,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                 k*nDOF_element + 
                 i] 
            = 
            dr[eN*nQuadraturePoints_element + 
               k]
            *
            w[eN*nQuadraturePoints_element*nDOF_element + 
              k*nDOF_element + 
              i];
          for (I=0;I<nSpace;I++)
            LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                   k*nDOF_element + 
                   i] 
              -=  
              df[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                     k*nDOF_element*nSpace + 
                     i*nSpace + 
                     I];
          for (I=0;I<nSpace;I++)
            for(J=0;J<nSpace;J++)
              LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                     k*nDOF_element + 
                     i] 
                +=  
                da[eN*nQuadraturePoints_element*nSpace2 + 
                   k*nSpace2 + 
                   I*nSpace+
                   J]
                *
                grad_phi[eN*nQuadraturePoints_element*nSpace+
                         k*nSpace+
                         J]
                *
                grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                       k*nDOF_element*nSpace + 
                       i*nSpace + 
                       I];
        }
}

/** 
    \brief  Calculate the product of two scalars at the quadrature points
*/
void calculateScalarScalarProduct(int nElements_global,
                                  int nQuadraturePoints_element,
                                  double* s1, 
                                  double* s2, 
                                  double* sResult)
{
  int eN,k;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      sResult[eN*nQuadraturePoints_element + 
              k] 
        = 
        s1[eN*nQuadraturePoints_element + 
           k]
        *
        s2[eN*nQuadraturePoints_element + 
           k];
}

/**
   \brief Calculate the product of a vector and a scalar at the quadrature points.
*/
void calculateVectorScalarProduct(int nElements_global,
                                  int nQuadraturePoints_element,
                                  int nSpace,
                                  double* v, 
                                  double* s, 
                                  double* vResult)
{
  int eN,k,I;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (I=0;I<nSpace;I++)
        vResult[eN*nQuadraturePoints_element*nSpace + 
                k*nSpace + 
                I] 
          = 
          v[eN*nQuadraturePoints_element*nSpace + 
            k*nSpace + 
            I]
          *
          s[eN*nQuadraturePoints_element + 
            k];
}

/**
   \brief Calculate the product of a tensor and scalar at the quadrature points
*/
void calculateTensorScalarProduct(int nElements_global,
                                  int nQuadraturePoints_element,
                                  int nSpace,
                                  double* t, 
                                  double* s, 
                                  double* tResult)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (I=0;I<nSpace;I++)
        for (J=0;J<nSpace;J++)
          tResult[eN*nQuadraturePoints_element*nSpace2 + 
                  k*nSpace2 + 
                  I*nSpace + 
                  J] 
            = 
            t[eN*nQuadraturePoints_element*nSpace2 + 
              k*nSpace2 + 
              I*nSpace + 
              J]
            *
            s[eN*nQuadraturePoints_element + 
              k];
}

/**
   \brief  Calculate the strong form of the residual of the scalar advection-diffusion-reaction equation
*/
void calculatePDEResidualADR(int nElements_global,
                             int nQuadraturePoints_element,
                             double* div_f,
                             double* div_a,
                             double* r,
                             double* mt,
                             double* pdeResidual)
{
  int eN,k;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      {
        pdeResidual[eN*nQuadraturePoints_element+
                    k] 
          = 
          div_f[eN*nQuadraturePoints_element+
                k] 
          +
          div_a[eN*nQuadraturePoints_element+
                k] 
          +
          r[eN*nQuadraturePoints_element+
            k] 
          +
          mt[eN*nQuadraturePoints_element+
             k];
      }
}

/**
   \brief Calculate the Jacobian of the strong form of the residual of the scalar advection-diffusion-reaction equation.
*/
void calculatePDEResidualJacobianADR(int nElements_global,
                                     int nQuadraturePoints_element,
                                     int nDOF_element,
                                     double* ddiv_f,
                                     double* ddiv_a,
                                     double* dr,
                                     double* dmt,
                                     double* v,
                                     double* dpdeResidual)
{
  int eN,k,j;
  for(eN=0;eN<nElements_global;eN++)
    for(j=0;j<nDOF_element;j++)
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          dpdeResidual[eN*nQuadraturePoints_element*nDOF_element+
                       k*nDOF_element + 
                       j] 
            = 
            ddiv_f[eN*nQuadraturePoints_element*nDOF_element+
                   k*nDOF_element + 
                   j]
            +
            ddiv_a[eN*nQuadraturePoints_element*nDOF_element+
                   k*nDOF_element + 
                   j] 
            +
            (dr[eN*nQuadraturePoints_element+
                k] 
             +
             dmt[eN*nQuadraturePoints_element+
                 k])
            *
            v[eN*nQuadraturePoints_element*nDOF_element+
              k*nDOF_element + 
              j];
        }
}

/**
   \brief  Calculate the divergence of a nonlinear vector field via the chain rule
*/
void calculateDiv_f(int nElements_global,
                    int nQuadraturePoints_element,
                    int nDOF_element,
                    int nSpace,
                    double* df,
                    double* grad_u,
                    double* grad_v,
                    double* div_f,
                    double* ddiv_f)
{
  int eN,k,j,I;
  memset(div_f, 0, sizeof(double) * nElements_global * nQuadraturePoints_element);
  memset(ddiv_f, 0, sizeof(double) * nElements_global * nQuadraturePoints_element * nDOF_element);
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      {
        for(I=0;I<nSpace;I++)
          div_f[eN*nQuadraturePoints_element + 
                k]
            +=
            df[eN*nQuadraturePoints_element*nSpace + 
               k*nSpace + 
               I]
            *
            grad_u[eN*nQuadraturePoints_element*nSpace + 
                   k*nSpace + 
                   I];
        for(j=0;j<nDOF_element;j++)
          {
            for(I=0;I<nSpace;I++)
              {
                ddiv_f[eN*nQuadraturePoints_element*nDOF_element + 
                       k*nDOF_element + 
                       j]
                  +=
                  df[eN*nQuadraturePoints_element*nSpace + 
                     k*nSpace + 
                     I]
                  *
                  grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_element + 
                         k*nSpace*nDOF_element + 
                         j*nSpace + 
                         I];
              }
          }
      }
}

/**
   \brief  Calculate the divergence of a nonlinear diffusion via the chain rule
*/
void calculateDiv_a(int nElements_global,
                    int nQuadraturePoints_element,
                    int nDOF_element,
                    int nSpace,
                    int* l2g,
                    double* da,
                    double* dphi,
                    double* grad_phi,
                    double* grad_u,
                    double* grad_v,
                    double* div_a,
                    double* ddiv_a)
{
  int eN,k,j,I,J,nSpace2=nSpace;
  memset(div_a, 0, sizeof(double) * nElements_global * nQuadraturePoints_element);
  memset(ddiv_a, 0, sizeof(double) * nElements_global * nQuadraturePoints_element * nDOF_element);
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      {
        for(I=0;I<nSpace;I++)
          for (J=0;J<nSpace;J++)
            div_a[eN*nQuadraturePoints_element + 
                  k]
              -=
              da[eN*nQuadraturePoints_element*nSpace2 + 
                 k*nSpace2 + 
                 I*nSpace+
                 J]
              *
              grad_phi[eN*nQuadraturePoints_element*nSpace + 
                       k*nSpace + 
                       J]
              *grad_u[eN*nQuadraturePoints_element*nSpace + 
                      k*nSpace + 
                      I];
        for(j=0;j<nDOF_element;j++)
          {
            for(I=0;I<nSpace;I++)
              for(J=0;J<nSpace;J++)
                ddiv_a[eN*nQuadraturePoints_element*nDOF_element + 
                       k*nDOF_element + 
                       j]
                  -=
                  da[eN*nQuadraturePoints_element*nSpace2 + 
                     k*nSpace2 + 
                     I*nSpace+
                     J]
                  *
                  (grad_phi[eN*nQuadraturePoints_element*nSpace + 
                            k*nSpace + 
                            J]
                   *
                   grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_element + 
                          k*nSpace*nDOF_element + 
                          j*nSpace + 
                          I]
                   +
                   dphi[l2g[eN*nDOF_element + 
                            j]]*
                   grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_element + 
                          k*nSpace*nDOF_element + 
                          j*nSpace + 
                          J]
                   *
                   grad_u[eN*nQuadraturePoints_element*nSpace+ 
                          k*nSpace+ 
                          I]);
          }
      }
}


/**
   \brief Calculate the stabilization parameter for the scalar advection-diffusion-reaction equation
*/
void calculateStabilizationADR(int nElements_global,
                               int nQuadraturePoints_element,
                               int nSpace,
                               char stabilization,
                               double* elementDiameter,
                               double* df,
                               double* a,
                               double* da,
                               double* grad_phi,
                               double* dphi,
                               double* dr,
                               double* dmt,
                               double* pe,
                               double* cfl,
                               double* tau)
{
  if(stabilization == '2')
    calculateStabilizationADR_2(nElements_global, nQuadraturePoints_element, nSpace, elementDiameter, df, a, da, grad_phi, dphi, dr, dmt, pe, cfl, tau);
  else if(stabilization == '1')
    calculateStabilizationADR_1(nElements_global, nQuadraturePoints_element, nSpace, elementDiameter, df, a, da, grad_phi, dphi, dr, dmt, pe, cfl, tau);
  else if(stabilization == 'p')
    calculateStabilizationADR_p(nElements_global, nQuadraturePoints_element, nSpace, elementDiameter, df, a, da, grad_phi, dphi, dr, dmt, pe, cfl, tau);
	
}

/**
   \brief Calculate the stabilization parameter for the scalar advection-diffusion-reaction equation using the peclet  number formula
*/
void calculateStabilizationADR_p(int nElements_global,
                                 int nQuadraturePoints_element,
                                 int nSpace,
                                 double* elementDiameter,
                                 double* df,
                                 double* a,
                                 double* da,
                                 double* grad_phi,
                                 double* dphi,
                                 double* dr,
                                 double* dmt,
                                 double* pe,
                                 double* cfl,
                                 double* tau)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  double h,Vlin,Alin,A_II,num,den,cfl1,cfl2,vlin[3];
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          Vlin = 0.0;
          Alin = 0.0;
          for(I=0;I<nSpace;I++)
            {
              vlin[I] = df[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           I];
              for(J=0;J<nSpace;J++)
                vlin[I] 
                  += 
                  da[eN*nQuadraturePoints_element*nSpace2 + 
                     k*nSpace2 + 
                     I*nSpace + 
                     J]
                  *
                  grad_phi[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           J];
            }
          for(I=0;I<nSpace;I++)
            Vlin += vlin[I]*vlin[I];
          Vlin = sqrt(Vlin);
          /*max diagonal entry for A. This choice needs to be looked into*/
          for(I=0;I<nSpace;I++)
            {
              A_II = a[eN*nQuadraturePoints_element*nSpace2 + 
                       k*nSpace2 + 
                       I*nSpace + 
                       I];
              Alin = (A_II > Alin) ? A_II : Alin;
            }
          Alin*=dphi[eN*nQuadraturePoints_element + 
                     k];
          cfl1 = Vlin/h;
          cfl2 = Alin/(h*h);
          if (cfl1 > cfl2)
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl1;
          else
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl2;
          num = 0.5*Vlin*h + 1.0e-8;
          den = Alin + num*1.0e-8;
          pe[eN*nQuadraturePoints_element + 
             k] = num/den;
          tau[eN*nQuadraturePoints_element + k]=0.5*h*(
                                                       1.0/tanh(pe[eN*nQuadraturePoints_element+
                                                                   k]) -
                                                       1.0/pe[eN*nQuadraturePoints_element+
                                                              k])/(Vlin+1.0e-8);
        }
    }
}

/**
   \brief Calculate the stabilization parameter for the scalar advection-diffusion-reaction equation using the "l_1 norm formula"
*/
void calculateStabilizationADR_1(int nElements_global,
                                 int nQuadraturePoints_element,
                                 int nSpace,
                                 double* elementDiameter,
                                 double* df,
                                 double* a,
                                 double* da,
                                 double* grad_phi,
                                 double* dphi,
                                 double* dr,
                                 double* dmt,
                                 double* pe,
                                 double* cfl,
                                 double* tau)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  double h,Vlin,Alin,A_II,num,den,cfl1,cfl2,vlin[3];
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          Vlin = 0.0;
          Alin = 0.0;
          for(I=0;I<nSpace;I++)
            {
              vlin[I] = df[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           I];
              for(J=0;J<nSpace;J++)
                vlin[I] 
                  += 
                  da[eN*nQuadraturePoints_element*nSpace2 + 
                     k*nSpace2 + 
                     I*nSpace + 
                     J]
                  *
                  grad_phi[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           J];
            }
          for(I=0;I<nSpace;I++)
            Vlin += vlin[I]*vlin[I];
          Vlin = sqrt(Vlin);
          /*max diagonal entry for A. This choice needs to be looked into*/
          for(I=0;I<nSpace;I++)
            {
              A_II = a[eN*nQuadraturePoints_element*nSpace2 + 
                       k*nSpace2 + 
                       I*nSpace + 
                       I];
              Alin = (A_II > Alin) ? A_II : Alin;
            }
          Alin*=dphi[eN*nQuadraturePoints_element + 
                     k];
          cfl1 = Vlin/h;
          cfl2 = Alin/(h*h);
          if (cfl1 > cfl2)
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl1;
          else
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl2;
          num = 0.5*Vlin*h + 1.0e-8;
          den = Alin + num*1.0e-8;
          pe[eN*nQuadraturePoints_element + 
             k] = num/den;
          tau[eN*nQuadraturePoints_element + k]=1.0/((2.0*Vlin/h)+ 
                                                     (4.0*Alin/(h*h)) +
                                                     fabs(dr[eN*nQuadraturePoints_element + 
                                                             k])+
                                                     fabs(dmt[eN*nQuadraturePoints_element + 
                                                              k]));
        }
    }
}

/**
   \brief Calculate the stabilization parameter for the scalar advection-diffusion-reaction equation using the "l_2 norm" formula
*/
void calculateStabilizationADR_2(int nElements_global,
                                 int nQuadraturePoints_element,
                                 int nSpace,
                                 double* elementDiameter,
                                 double* df,
                                 double* a,
                                 double* da,
                                 double* grad_phi,
                                 double* dphi,
                                 double* dr,
                                 double* dmt,
                                 double* pe,
                                 double* cfl,
                                 double* tau)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  double h,Vlin,Alin,A_II,num,den,cfl1,cfl2,vlin[3];
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          Vlin = 0.0;
          Alin = 0.0;
          for(I=0;I<nSpace;I++)
            {
              vlin[I] = df[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           I];
              for(J=0;J<nSpace;J++)
                vlin[I] 
                  += 
                  da[eN*nQuadraturePoints_element*nSpace2 + 
                     k*nSpace2 + 
                     I*nSpace + 
                     J]
                  *
                  grad_phi[eN*nQuadraturePoints_element*nSpace + 
                           k*nSpace + 
                           J];
            }
          for(I=0;I<nSpace;I++)
            Vlin += vlin[I]*vlin[I];
          Vlin = sqrt(Vlin);
          /*max diagonal entry for A. This choice needs to be looked into*/
          for(I=0;I<nSpace;I++)
            {
              A_II = a[eN*nQuadraturePoints_element*nSpace2 + 
                       k*nSpace2 + 
                       I*nSpace + 
                       I];
              Alin = (A_II > Alin) ? A_II : Alin;
            }
          Alin*=dphi[eN*nQuadraturePoints_element + 
                     k];
          cfl1 = Vlin/h;
          cfl2 = Alin/(h*h);
          if (cfl1 > cfl2)
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl1;
          else
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl2;
          num = 0.5*Vlin*h + 1.0e-8;
          den = Alin + num*1.0e-8;
          tau[eN*nQuadraturePoints_element + 
              k]
            =1.0/sqrt((2.0*Vlin/h)*(2.0*Vlin/h)+ 
                      9.0*(4.0*Alin/(h*h))*(4.0*Alin/(h*h)) +
                      dr[eN*nQuadraturePoints_element + 
                         k]
                      *
                      dr[eN*nQuadraturePoints_element + 
                         k]
                      +
                      dmt[eN*nQuadraturePoints_element + 
                          k]
                      *
                      dmt[eN*nQuadraturePoints_element + 
                          k]);
          pe[eN*nQuadraturePoints_element + 
             k] = num/den;
        }
    }
}

/**
   \brief Calculate the Peclet and Courant-Friedrichs-Lewy numbers for the scalar advection-diffusion-reaction equation
*/
void calculateDimensionlessNumbersADR(int nElements_global,
                                      int nQuadraturePoints_element,
                                      int nSpace,
                                      double* elementDiameter,
                                      double* df,
                                      double* a,
                                      double* dphi,
                                      double* dr,
                                      double* dmt,
                                      double* pe,
                                      double* cfl)
{
  int eN,k,I,nSpace2=nSpace*nSpace;
  double h,Vlin,Alin,A_II,num,den,cfl1,cfl2;
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          Vlin = 0.0;
          Alin = 0.0;
          for(I=0;I<nSpace;I++)
            Vlin 
              += 
              df[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              df[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I];
          Vlin = sqrt(Vlin);
          /*max diagonal entry for A. This choice needs to be looked into*/
          for(I=0;I<nSpace;I++)
            {
              A_II = a[eN*nQuadraturePoints_element*nSpace2 + 
                       k*nSpace2 + 
                       I*nSpace + 
                       I];
              Alin = (A_II > Alin) ? A_II : Alin;
            }
          Alin*=dphi[eN*nQuadraturePoints_element + 
                     k];
          cfl1 = Vlin/h;
          cfl2 = Alin/(h*h);
          if (cfl1 > cfl2)
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl1;
          else
            cfl[eN*nQuadraturePoints_element + 
                k] = cfl2;
          num = 0.5*Vlin*h + 1.0e-8;
          den = Alin + num*1.0e-8;
          pe[eN*nQuadraturePoints_element + k] = num/den;
        }
    }
}

/**
   \brief Calculate the shock capturing diffusion for  the scalar advection-diffusion-reaction equation
*/
void calculateShockCapturingADR(int nElements_global,
                                int nQuadraturePoints_element,
                                int nSpace,
                                char shockCapturing,
                                double shockCapturingDiffusion,
                                double* elementDiameter,
                                double* pdeResidual,
                                double* mt,
                                double* grad_u,
                                double* numDiff)
{
  int eN,k,I;
  double h,num,den=0.0,enorm_grad_u;
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          num = shockCapturingDiffusion*h*fabs(pdeResidual[eN*nQuadraturePoints_element+
                                                           k]); 
          enorm_grad_u = 0.0;
          for(I=0;I<nSpace;I++)
            {
              enorm_grad_u
                +=
                grad_u[eN*nQuadraturePoints_element*nSpace+
                       k*nSpace+
                       I]
                *
                grad_u[eN*nQuadraturePoints_element*nSpace+
                       k*nSpace+
                       I];
            }
          enorm_grad_u = sqrt(enorm_grad_u);
          if(shockCapturing == '1')
            {
              den = (fabs(mt[eN*nQuadraturePoints_element+
                             k]) +
                     enorm_grad_u +
                     1.0e-8);
            }
          else if(shockCapturing == '2')
            {
              den = sqrt(mt[eN*nQuadraturePoints_element+
                            k]
                         *
                         mt[eN*nQuadraturePoints_element+
                            k] 
                         +
                         enorm_grad_u*enorm_grad_u +
                         1.0e-16);
            }
          numDiff[eN*nQuadraturePoints_element+k] = num/den;
        }
    }
}

/**
   \brief Loop over all the elements and update the element residual with the numerical quadrature approximation of the Hamiltonian integral.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param h (nElements_global x nQuadraturePoints_element). The Hamiltonian associated with each quadrature point (includes integration weight).
   @param w (nElements_global x  nQuadraturePoints_element  x nDOF_element). The finite element test function values.
   @param residual (nElements_global x nDOF_element). The element residual, which is  updated  upon return.

   The result of calling this function is
   \f[
   residual_{e,i} \mapsto residual_{e,i} - \int_{\Omega_e} h w_i dV \quad \forall e,i 
   \f]

   where

   \f[
   \int_{\Omega_e} h w_i dV =  \int_{\Omega_r} h w_i |J_e| d\hat{V} \int_{\Omega_r} \sum_k h_k w_i 
   \f]
   \f[
   h_k = h(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateHamiltonian(int nElements_global,
                       int nQuadraturePoints_element,
                       int nDOF_element,
                       double* h,
                       double* w,
                       double* residual)
{
  int eN,i,k;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        residual[eN*nDOF_element + i] 
          +=
          h[eN*nQuadraturePoints_element + 
            k]
          *
          w[eN*nQuadraturePoints_element*nDOF_element + 
            k*nDOF_element + 
            i];
}

/**
   \brief Loop over all the elements and update the element Jacobian with the numerical quadrature approximation of the Hamiltonian integral Jacobian.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param dh (nElements_global x nQuadraturePoints_element x nSpace). The derivative of the Hamiltonian (with respect to  the solution gradient) associated with each quadrature point (includes integration weight).
   @param grad_v_x_w (nElements_global x  nQuadraturePoints_element  x nDOF_element x nSpace). The tensor product of the  finite element trial and test function gradient values.
   @param jacobian (nElements_global x nDOF_element x nDOF_element). The element jacobian, which is  updated  upon return.

   The result of calling this function is
   \f[
   jacobian_{e,i,j} \mapsto jacobian_{e,i,j} - \int_{\Omega_e} \frac{\partial \mathbf{h}}{\partial \nabla u} \cdot \nabla v_j w_i dV \quad \forall e,i,j 
   \f]

   where

   \f[
   \int_{\Omega_e} \frac{\partial \mathbf{h}}{\partial \nabla  u} \cdot \nabla v_j \nabla w_i dV =  \int_{\Omega_r} \frac{\partial \mathbf{h}}{\partial \nabla u} \nabla v_j \cdot w_i |J_e| d\hat{V} \approx  \sum_k dh_k \cdot (\nabla v \otimes w)_{j,i} 
   \f]
   \f[
   dh_k = \frac{\partial \mathbf{h}}{\partial \nabla u}(\mathbf{\hat{x}}_k) |J_e(\mathbf{\hat{x}}_k)| \omega_k
   \f]
*/
void updateHamiltonianJacobian(int nElements_global,
                               int nQuadraturePoints_element,
                               int nDOF_element,
                               int nSpace,
                               double* dh,
                               double* grad_v_x_w,
                               double* jacobian)
{
  int eN,i,j,k,I,nDOF_element2=nDOF_element*nDOF_element;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            jacobian[eN*nDOF_element2 + 
                     i*nDOF_element +
                     j]
              +=
              dh[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace+
                 I]
              *
              grad_v_x_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace + 
                         k*nDOF_element2*nSpace +
                         j*nDOF_element*nSpace +
                         i*nSpace+
                         I];
}

/**
   \brief Calculate the linearized adjoint of a Hamilton-Jacobi equation at the quadrature points.
*/
void calculateAdjointHJ(int nElements_global,
                        int nQuadraturePoints_element,
                        int nDOF_element,
                        int nSpace,
                        double* w,
                        double* grad_w,
                        double* dh,
                        double* rh,
                        double* LstarW)
{
  int eN,i,k,I;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                 k*nDOF_element + 
                 i] 
            = 0.0;
          for (I=0;I<nSpace;I++)
            LstarW[eN*nQuadraturePoints_element*nDOF_element + 
                   k*nDOF_element + 
                   i] 
              -=  
              dh[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                     k*nDOF_element*nSpace + 
                     i*nSpace + 
                     I];
        }
}

/**
   \brief Calculate the strong form of the residual of a Hamilton-Jacobi equation at the quadrature points
*/
void calculatePDEResidualHJ(int nElements_global,
                            int nQuadraturePoints_element,
                            int nSpace,
                            double* dh,
                            double* grad_u,
                            double* rh,
                            double* mt,
                            double* pdeResidual)
{
  int eN,k,I;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      {
        pdeResidual[eN*nQuadraturePoints_element+k] =  mt[eN*nQuadraturePoints_element+k] +
          rh[eN*nQuadraturePoints_element+k];
        for (I=0;I<nSpace;I++)
          pdeResidual[eN*nQuadraturePoints_element+k] += dh[eN*nQuadraturePoints_element*nSpace+k*nSpace+I]*
            grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace+I];
      }
}

/**
   \brief Calculate the Jacobian of the strong form of the residual of a Hamilton-Jacobi equation at the quadrature points.
*/
void calculatePDEResidualJacobianHJ(int nElements_global,
                                    int nQuadraturePoints_element,
                                    int nDOF_element,
                                    int nSpace,
                                    double* dh,
                                    double* grad_v,
                                    double* dmt,
                                    double* v,
                                    double* dpdeResidual)
{
  int eN,k,j,I;
  for(eN=0;eN<nElements_global;eN++)
    for(j=0;j<nDOF_element;j++)
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          dpdeResidual[eN*nQuadraturePoints_element*nDOF_element+
                       k*nDOF_element + 
                       j] 
            = 
            dmt[eN*nQuadraturePoints_element+
                k]
            *
            v[eN*nQuadraturePoints_element*nDOF_element+
              k*nDOF_element + 
              j];
          for (I=0;I<nSpace;I++)
            dpdeResidual[eN*nQuadraturePoints_element*nDOF_element+
                         k*nDOF_element + 
                         j] 
              += 
              dh[eN*nQuadraturePoints_element*nSpace+
                 k*nSpace+
                 I]
              *
              grad_v[eN*nQuadraturePoints_element*nDOF_element*nSpace+
                     k*nDOF_element*nSpace + 
                     j*nSpace + 
                     I];
        }
}

/**
   \brief Calculate the stabilization parameter for a Hamilton-Jacobi equation at the quadrature points
*/
void calculateStabilizationHJ(int nElements_global,
                              int nQuadraturePoints_element,
                              int nSpace,
                              char stabilization,
                              double* elementDiameter,
                              double* dh,
                              double* dmt,
                              double* cfl,
                              double* tau)
{
  int eN,k,I;
  double h,Vlin;
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          Vlin = 0.0;
          for(I=0;I<nSpace;I++)
            Vlin 
              += 
              dh[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              dh[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I];
          Vlin = sqrt(Vlin);
          cfl[eN*nQuadraturePoints_element + 
              k] 
            = Vlin/h;
          if(stabilization == '2')
            tau[eN*nQuadraturePoints_element + 
                k]=1.0/sqrt((2.0*Vlin/h)*(2.0*Vlin/h)
                            + 
                            dmt[eN*nQuadraturePoints_element + 
                                k]
                            *
                            dmt[eN*nQuadraturePoints_element + 
                                k]);
          else if (stabilization == '1')
            tau[eN*nQuadraturePoints_element + k]=1.0/((2.0*Vlin/h)+ 
                                                       fabs(dmt[eN*nQuadraturePoints_element + 
                                                                k]));
        }
    }
}

/**
   \brief Calculate the shock capturing diffusion for a Hamilton-Jacobi equation at the quadrature  points
*/
void calculateShockCapturingHJ(int nElements_global,
                               int nQuadraturePoints_element,
                               int nSpace,
                               char shockCapturing,
                               double shockCapturingDiffusion,
                               double* elementDiameter,
                               double* pdeResidual,
                               double* mt,
                               double* dh,
                               double* numDiff)
{
  int eN,k,I;
  double h,num,den=0.0,enorm_dh;
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
          num = shockCapturingDiffusion*h*fabs(pdeResidual[eN*nQuadraturePoints_element+
                                                           k]); 
          enorm_dh = 0.0;
          for(I=0;I<nSpace;I++)
            {
              enorm_dh
                +=
                dh[eN*nQuadraturePoints_element*nSpace+
                   k*nSpace+
                   I]
                *
                dh[eN*nQuadraturePoints_element*nSpace+
                   k*nSpace+
                   I];
            }
          enorm_dh = sqrt(enorm_dh);
          if(shockCapturing == '1')
            {
              den = (fabs(mt[eN*nQuadraturePoints_element+
                             k]) +
                     enorm_dh)+num*1.0e-8;
            }
          else if(shockCapturing == '2')
            {
              den = sqrt(mt[eN*nQuadraturePoints_element+
                            k]
                         *
                         mt[eN*nQuadraturePoints_element+
                            k] 
                         +
                         enorm_dh*enorm_dh) + num*1.0e-8;
            }
          numDiff[eN*nQuadraturePoints_element+k] = num/den;
        }
    }
}

/**
   \brief Calcualte the tensor product of trial and test functions at the quadrature  points
*/
void calculateShape_x_Shape(int nElements_global,
                            int nQuadraturePoints_element,
                            int nDOF_element,
                            double* v,
                            double* w,
                            double* v_x_w)
{
  int eN,k,i,j,nDOF_element2=nDOF_element*nDOF_element;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (i=0;i<nDOF_element;i++)
        for (j=0;j<nDOF_element;j++)
          v_x_w[eN*nQuadraturePoints_element*nDOF_element2 + 
                k*nDOF_element2+
                j*nDOF_element+
                i] 
            = 
            v[eN*nQuadraturePoints_element*nDOF_element + 
              k*nDOF_element+
              j]
            *
            w[eN*nQuadraturePoints_element*nDOF_element + 
              k*nDOF_element+
              i];
}

/**
   \brief Calculate  the tensor  product of trial functions and test function gradients at the quadrature  points.
*/
void calculateShape_x_GradShape(int nElements_global,
                                int nQuadraturePoints_element,
                                int nDOF_element,
                                int nSpace,
                                double* v,
                                double* grad_w,
                                double* v_x_grad_w)
{
  int eN,k,i,j,I,nDOF_element2=nDOF_element*nDOF_element;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (i=0;i<nDOF_element;i++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace + 
                       k*nDOF_element2*nSpace+
                       j*nDOF_element*nSpace+
                       i*nSpace+I] 
              = 
              v[eN*nQuadraturePoints_element*nDOF_element + 
                k*nDOF_element+
                j]
              *
              grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                     k*nDOF_element*nSpace+
                     i*nSpace+
                     I];
}

/**
   \brief Calculate  the tensor  product of trial function gradients and test functions at the quadrature points.
*/
void calculateGradShape_x_Shape(int nElements_global,
                                int nQuadraturePoints_element,
                                int nDOF_element,
                                int nSpace,
                                double* grad_v,
                                double* w,
                                double* grad_v_x_w)
{
  int eN,k,i,j,I,nDOF_element2=nDOF_element*nDOF_element;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (i=0;i<nDOF_element;i++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            grad_v_x_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace + 
                       k*nDOF_element2*nSpace+
                       j*nDOF_element*nSpace+
                       i*nSpace+
                       I] 
              = 
              grad_v[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                     k*nDOF_element*nSpace+
                     j*nSpace+
                     I]
              *
              w[eN*nQuadraturePoints_element*nDOF_element + 
                k*nDOF_element+
                i];
}

/**
   \brief Calculate  the tensor  product of trial function gradients and test function gradients at the quadrature points.
*/
void calculateGradShape_x_GradShape(int nElements_global,
                                    int nQuadraturePoints_element,
                                    int nDOF_element,
                                    int nSpace,
                                    double* grad_v,
                                    double* grad_w,
                                    double* grad_v_x_grad_w)
{
  int eN,k,i,j,I,J,nDOF_element2=nDOF_element*nDOF_element,nSpace2=nSpace*nSpace;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (i=0;i<nDOF_element;i++)
        for (j=0;j<nDOF_element;j++)
          for (I=0;I<nSpace;I++)
            for (J=0;J<nSpace;J++)
              grad_v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace2 + 
                              k*nDOF_element2*nSpace2+
                              j*nDOF_element*nSpace2+
                              i*nSpace2+
                              I*nSpace+
                              J] 
                = 
                grad_v[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                       k*nDOF_element*nSpace+
                       j*nSpace+
                       I]
                *
                grad_w[eN*nQuadraturePoints_element*nDOF_element*nSpace + 
                       k*nDOF_element*nSpace+
                       i*nSpace+
                       J];
}

/**
   \brief Calculate the physical space integration weights from the reference element weights and Jacobian determinants.
*/
void calculateIntegrationWeights(int nElements_global,
                                 int nQuadraturePoints_element,
                                 double* abs_det_J,
                                 double* referenceWeights,
                                 double* weights)
{
  int eN,k;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      weights[eN*nQuadraturePoints_element+
              k]
        =
        abs_det_J[eN*nQuadraturePoints_element+
                  k]
        *
        referenceWeights[k];
}

/**
   \brief Calculate the values of a multicomponent finite element function at the quadrature  points from the degrees of freedom and the test function values at the quadrature points
*/ 
void calculateFiniteElementFunctionValues(int nElements_global,
                                          int nQuadraturePoints_element,
                                          int nDOF_element,
                                          int  nComponents,
                                          int* l2g,
                                          double* dof,
                                          double* v,
                                          double* u)
{
  int eN,k,j,t;
  for (eN=0;eN<nElements_global;eN++)
    for  (k=0;k<nQuadraturePoints_element;k++)
      for (t=0;t<nComponents;t++)
        u[eN*nQuadraturePoints_element*nComponents+
          k*nComponents+
          t] = 0.0;
  for (eN=0;eN<nElements_global;eN++)
    for  (k=0;k<nQuadraturePoints_element;k++)
      for (j=0;j<nDOF_element;j++)
        for (t=0;t<nComponents;t++)
          u[eN*nQuadraturePoints_element*nComponents+
            k*nComponents+
            t] 
            += 
            dof[l2g[eN*nDOF_element+
                    j]*nComponents+
                t]
            *
            v[eN*nQuadraturePoints_element*nDOF_element+
              k*nDOF_element+
              j];
}

/**
   \brief Calculate the gradient values of a multicomponent finite element function at the quadrature  points from the degrees of freedom and the test function gradient values at the quadrature points
*/ 
void calculateFiniteElementFunctionGradientValues(int nElements_global,
                                                  int nQuadraturePoints_element,
                                                  int nDOF_element,
                                                  int nComponents,
                                                  int nSpace,
                                                  int* l2g,
                                                  double* dof,
                                                  double* grad_v,
                                                  double* grad_u)
{
  int eN,k,j,t,I;
  for (eN=0;eN<nElements_global;eN++)
    for  (k=0;k<nQuadraturePoints_element;k++)
      for (t=0;t<nComponents;t++)
        for  (I=0;I<nSpace;I++)
          grad_u[eN*nQuadraturePoints_element*nComponents*nSpace+
                 k*nComponents*nSpace+
                 t*nSpace+
                 I] = 0.0;
  for (eN=0;eN<nElements_global;eN++)
    for  (k=0;k<nQuadraturePoints_element;k++)
      for (j=0;j<nDOF_element;j++)
        for (t=0;t<nComponents;t++)
          for  (I=0;I<nSpace;I++)
            grad_u[eN*nQuadraturePoints_element*nComponents*nSpace+
                   k*nComponents*nSpace+
                   t*nSpace+
                   I]
              +=
              dof[l2g[eN*nDOF_element+
                      j]*nComponents+
                  t]
              *
              grad_v[eN*nQuadraturePoints_element*nDOF_element*nSpace+
                     k*nDOF_element*nSpace+
                     j*nSpace+
                     I];
}

/**
   \brief  Loop over all the quadrature points and calculate the tensor product of the solution gradient with the test functions.

   @param nElements_global The number of elements in the mesh
   @param nQuadraturePoints_element The number of quadrature points per element.
   @param nDOF_element The number of degrees of freedom per element.
   @param nSpace The number of spatial dimensions.
   @param l2g (nElements_global x nDOF_element) The mapping between element degrees of  freedom and global degrees of freedom.
   @param dof (nDOF_global). The global degrees of  freedom of u.
   @param grad_v_x_grad_w (nElements_global x nQuadraturePoints_element x x nDOF_element x nSpace x nSpace) The tensor product of the solution gradient values and  trial function gradient values at the quadrature points.
   @param grad_u_x_grad_w (nElements_global x nQuadraturePoints_element x x nDOF_element x nSpace x nSpace) The tensor product of the solution gradient values and  trial function gradient values at the quadrature points.

   \f[
   u_{e,q} \mapsto \sum_j u_j v_{e,q,j}
   \f]
   \f[
   \nabla u_{e,q} \mapsto \sum_j u_j \nabla v_{e,q,j}
   \f]
   \f[
   (\nabla u \otimes \nabla w)_{e,q,j,i} \mapsto \sum_j u_j (\nabla v \otimes \nabla w)_{e,q,j,i}
   \f]

   where

   \f[
   u_j = u[l2g[e,j]]
   \f]
*/
void calculateFiniteElementFunctionGradientTensorValues(int nElements_global,
                                                        int nQuadraturePoints_element,
                                                        int nDOF_element,
                                                        int nComponents,
                                                        int nSpace,
                                                        int* l2g,
                                                        double* dof,
                                                        double* grad_v_x_grad_w,
                                                        double* grad_u_x_grad_w)
{
  int eN,k,j,i,I,J,t,nSpace2=nSpace*nSpace;
  int nDOF_element2 = nDOF_element*nDOF_element;
  /* initialize  to 0.0*/
  for(eN=0;eN<nElements_global;eN++)
    for(k=0;k<nQuadraturePoints_element;k++)
      for(i=0;i<nDOF_element;i++)
        for (t=0;t<nComponents;t++)
          for(I=0;I<nSpace;I++)
            for(J=0;J<nSpace;J++)
              grad_u_x_grad_w[eN*nQuadraturePoints_element*nDOF_element*nComponents*nSpace2+
                              k*nDOF_element*nComponents*nSpace2 + 
                              i*nComponents*nSpace2 +
                              t*nSpace2 +
                              I*nSpace + 
                              J] = 0.0;
  for(eN=0;eN<nElements_global;eN++)
    for(k=0;k<nQuadraturePoints_element;k++)
      for(j=0;j<nDOF_element;j++)
        for(i=0;i<nDOF_element;i++)
          for (t=0;t<nComponents;t++)
            for(I=0;I<nSpace;I++)
              for(J=0;J<nSpace;J++)
                grad_u_x_grad_w[eN*nQuadraturePoints_element*nDOF_element*nComponents*nSpace2+
                                k*nDOF_element*nComponents*nSpace2 + 
                                i*nComponents*nSpace2 +
                                t*nSpace2 +
                                I*nSpace + 
                                J] 
                  +=
                  dof[l2g[eN*nDOF_element+
                          j]*nComponents+
                      t]
                  *
                  grad_v_x_grad_w[eN*nQuadraturePoints_element*nDOF_element2*nSpace2 + 
                                  k*nDOF_element2*nSpace2 + 
                                  j*nDOF_element*nSpace2 + 
                                  i*nSpace2 + 
                                  I*nSpace + 
                                  J];
}

/**
   \brief Calculate the values of a multi-component finite elment function at element boundary quadrature points from the degrees of freedom and the trace of the trial functions at the element boundary quadrature points
*/
void calculateFiniteElementFunctionValuesTrace(int nElements_global,
                                               int nElementBoundaries_element,
                                               int nQuadraturePoints_elementBoundary,
                                               int nDOF_element,
                                               int nComponents,
                                               int* l2g,
                                               double* dof,
                                               double* v,
                                               double* u)
{
  int eN,ebN,k,j,t;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(t=0;t<nComponents;t++)
          u[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nComponents+
            ebN*nQuadraturePoints_elementBoundary*nComponents+
            k*nComponents+
            t] = 0.0;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(j=0;j<nDOF_element;j++)
          for(t=0;t<nComponents;t++)
            u[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nComponents+
              ebN*nQuadraturePoints_elementBoundary*nComponents+
              k*nComponents+
              t] 
              += 
              dof[l2g[eN*nDOF_element+j]*nComponents+
                  t]
              *
              v[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                ebN*nQuadraturePoints_elementBoundary*nDOF_element+
                k*nDOF_element+
                j];
}

/**
   \brief Calculate the gradients of a multi-component finite element function at the element boundary quadrature points from the degress of freedom and the trace of the trial functions at the element boundary quadrature  points.
*/
void calculateFiniteElementFunctionGradientValuesTrace(int nElements_global,
                                                       int nElementBoundaries_element,
                                                       int nQuadraturePoints_elementBoundary,
                                                       int nDOF_element,
                                                       int nComponents,
                                                       int nSpace,
                                                       int* l2g,
                                                       double* dof,
                                                       double* grad_v,
                                                       double* grad_u)
{
  int eN,ebN,k,j,t,I;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(t=0;t<nComponents;t++)
          for(I=0;I<nSpace;I++)
            grad_u[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nComponents*nSpace+
                   ebN*nQuadraturePoints_elementBoundary*nComponents*nSpace+
                   k*nComponents*nSpace+
                   t*nSpace+
                   I] = 0.0;
  for(eN=0;eN<nElements_global;eN++)
    for(ebN=0;ebN<nElementBoundaries_element;ebN++)
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(j=0;j<nDOF_element;j++)
          for(t=0;t<nComponents;t++)
            for(I=0;I<nSpace;I++)
              grad_u[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nComponents*nSpace+
                     ebN*nQuadraturePoints_elementBoundary*nComponents*nSpace+
                     k*nComponents*nSpace+
                     t*nSpace+
                     I] 
                += 
                dof[l2g[eN*nDOF_element+j]*nComponents+
                    t]
                *
                grad_v[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                       ebN*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                       k*nDOF_element*nSpace+
                       j*nSpace+
                       I];
}

/**
   \brief Update the global residuals from the element residuals.
*/
void updateGlobalResidualFromElementResidual(int nElements_global,
                                             int nDOF_element,
                                             int* nFreeDOF_element,
                                             int* freeLocal,
                                             int* freeGlobal,
                                             double* elementResidual,
                                             double* globalResidual)
{
  int eN,ii;
  for (eN=0;eN<nElements_global;eN++)
    for (ii=0;ii<nFreeDOF_element[eN];ii++)
      globalResidual[freeGlobal[eN*nDOF_element+
                                ii]] 
        +=
        elementResidual[eN*nDOF_element + 
                        freeLocal[eN*nDOF_element+
                                  ii]];
}

/**
   \brief  Update the global CSR jacobian  from the element Jacobians.
*/
void updateGlobalJacobianFromElementJacobian_CSR(int nElements_global,
                                                 int nDOF_element,
                                                 int* nFreeDOF_element,
                                                 int* freeLocal,
                                                 int* csrRowIndeces,
                                                 int* csrColumnOffsets,
                                                 double* elementJacobian,
                                                 double* globalJacobian)
{
  int eN,ii,jj,i,j,jacIndex,nDOF_element2=nDOF_element*nDOF_element;
  for (eN=0;eN<nElements_global;eN++)
    for (ii=0;ii<nFreeDOF_element[eN];ii++)
      {
	i = freeLocal[eN*nDOF_element+
		      ii];
	for (jj=0;jj<nFreeDOF_element[eN];jj++)
	  {
	    j = freeLocal[eN*nDOF_element+
			  jj];
	    jacIndex = csrRowIndeces[eN*nDOF_element+
				     ii]
	      +
	      csrColumnOffsets[eN*nDOF_element2+
			       ii*nDOF_element+
			       jj];
	    globalJacobian[jacIndex]
	      +=
	      elementJacobian[eN*nDOF_element2 + 
			      i*nDOF_element + 
			      j];
	  }
      }
}

/**
   \brief Update the global CSR Jacobian from the element boundary flux Jacobians at interior boundaries
*/    
void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(int nInteriorElementBoundaries_global,
                                                                     int nDOF_element,
                                                                     int nQuadraturePoints_elementBoundary,
                                                                     int nElementBoundaries_element,
                                                                     int* interiorElementBoundaries,
                                                                     int* elementBoundaryElements,
                                                                     int* elementBoundaryLocalElementBoundaries,
                                                                     int* nFreeDOF_element,
                                                                     int* freeLocal,
                                                                     int* csrRowIndeces,
                                                                     int* csrColumnOffsets_eb,
                                                                     double* elementBoundaryFluxJacobian,
                                                                     double* w,
                                                                     double* jac)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,ii,i,k,jj,j,nDOF_element2=nDOF_element*nDOF_element,jacIndex;
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global   = elementBoundaryElements[ebN*2+0];
      right_eN_global  = elementBoundaryElements[ebN*2+1];
      left_ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(ii=0;ii<nFreeDOF_element[left_eN_global];ii++)
        {
          i = freeLocal[left_eN_global*nDOF_element+ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[left_eN_global];jj++)
                {
                  j = freeLocal[left_eN_global*nDOF_element+
                                jj];
                  jacIndex = csrRowIndeces[left_eN_global*nDOF_element+
                                           ii] 
                    + 
                    csrColumnOffsets_eb[ebN*4*nDOF_element2 + 
                                        0*2*nDOF_element2 + 
                                        0*nDOF_element2 + 
                                        ii*nDOF_element + 
                                        jj];
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
              for(jj=0;jj<nFreeDOF_element[right_eN_global];jj++)
                {
                  j = freeLocal[right_eN_global*nDOF_element+
                                jj];
                  jacIndex = csrRowIndeces[left_eN_global*nDOF_element+
                                           ii] 
                    + 
                    csrColumnOffsets_eb[ebN*4*nDOF_element2+
                                        0*2*nDOF_element2+
                                        1*nDOF_element2+
                                        ii*nDOF_element+
                                        jj];
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                1*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
            }
        }
      for(ii=0;ii<nFreeDOF_element[right_eN_global];ii++)
        {
          i = freeLocal[right_eN_global*nDOF_element+
                        ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[left_eN_global];jj++)
                {
                  j = freeLocal[left_eN_global*nDOF_element+
                                jj];
                  jacIndex = csrRowIndeces[right_eN_global*nDOF_element+
                                           ii] 
                    + 
                    csrColumnOffsets_eb[ebN*4*nDOF_element2+
                                        1*2*nDOF_element2+
                                        0*nDOF_element2+
                                        ii*nDOF_element+
                                        jj];
                  jac[jacIndex] 
                    -= 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
              for(jj=0;jj<nFreeDOF_element[right_eN_global];jj++)
                {
                  j = freeLocal[right_eN_global*nDOF_element+
                                jj];
                  jacIndex = csrRowIndeces[right_eN_global*nDOF_element+
                                           ii] 
                    + 
                    csrColumnOffsets_eb[ebN*4*nDOF_element2+
                                        1*2*nDOF_element2+
                                        1*nDOF_element2+
                                        ii*nDOF_element+
                                        jj];
                  jac[jacIndex] 
                    -= 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                                                1*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]*
                    w[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
            }
        }
    }
}

/**
   \brief Update the global CSR Jacobian from the element boundary flux Jacobians at exterior boundaries
*/    
void updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(int nExteriorElementBoundaries_global,
                                                                     int nDOF_element,
                                                                     int nQuadraturePoints_elementBoundary,
                                                                     int nElementBoundaries_element,
                                                                     int* exteriorElementBoundaries,
                                                                     int* elementBoundaryElements,
                                                                     int* elementBoundaryLocalElementBoundaries,
                                                                     int* nFreeDOF_element,
                                                                     int* freeLocal,
                                                                     int* csrRowIndeces,
                                                                     int* csrColumnOffsets_eb,
                                                                     double* elementBoundaryFluxJacobian,
                                                                     double* w,
                                                                     double* jac)
{
  int ebNE,ebN,eN_global,ebN_element,ii,i,k,jj,j,nDOF_element2=nDOF_element*nDOF_element,jacIndex;
  for (ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global   = elementBoundaryElements[ebN*2+0];
      ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(ii=0;ii<nFreeDOF_element[eN_global];ii++)
        {
          i = freeLocal[eN_global*nDOF_element+ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[eN_global];jj++)
                {
                  j = freeLocal[eN_global*nDOF_element+
                                jj];
                  jacIndex = csrRowIndeces[eN_global*nDOF_element+
                                           ii] 
                    + 
                    csrColumnOffsets_eb[ebN*4*nDOF_element2 + 
                                        0*2*nDOF_element2 + 
                                        0*nDOF_element2 + 
                                        ii*nDOF_element + 
                                        jj];
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
            }
        }
    }
}

/**
   \brief  Update the global dense jacobian  from the element Jacobians.
*/
void updateGlobalJacobianFromElementJacobian_dense(int nElements_global,
						   int nDOF_element,
						   int nFreeDOF_global,
						   int* nFreeDOF_element,
						   int* freeLocal,
						   int* freeGlobal,
						   double* elementJacobian,
						   double* globalJacobian)
{
  int eN,ii,jj,nDOF_element2=nDOF_element*nDOF_element,i,j,jacIndex,I,J;
  for (eN=0;eN<nElements_global;eN++)
    for (ii=0;ii<nFreeDOF_element[eN];ii++)
      {
	i = freeLocal[eN*nDOF_element+
		      ii];
	I = freeGlobal[eN*nDOF_element+
		       ii];
	for (jj=0;jj<nFreeDOF_element[eN];jj++)
	  {
	    j = freeLocal[eN*nDOF_element+
			  jj];
	    J = freeGlobal[eN*nDOF_element+
			   jj];
	    /*jacIndex = I*nFreeDOF_global + J;*/
            jacIndex = I + J*nFreeDOF_global;
	    globalJacobian[jacIndex]
	      +=
	      elementJacobian[eN*nDOF_element2 +
			      i*nDOF_element+
			      j];
	  }
      }
}

/**
   \brief Update the global dense Jacobian from the element boundary flux Jacobians at interior boundaries
*/    
void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(int nInteriorElementBoundaries_global,
								       int nDOF_element,
								       int nQuadraturePoints_elementBoundary,
								       int nElementBoundaries_element,
								       int nFreeDOF_global,
								       int* interiorElementBoundaries,
								       int* elementBoundaryElements,
								       int* elementBoundaryLocalElementBoundaries,
								       int* nFreeDOF_element,
								       int* freeLocal,
								       int* freeGlobal,
								       double* elementBoundaryFluxJacobian,
								       double* w,
								       double* jac)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,ii,i,k,jj,j,jacIndex,I,J;
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global   = elementBoundaryElements[ebN*2+0];
      right_eN_global  = elementBoundaryElements[ebN*2+1];
      left_ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(ii=0;ii<nFreeDOF_element[left_eN_global];ii++)
        {
          i = freeLocal[left_eN_global*nDOF_element+
			ii];
	  I = freeGlobal[left_eN_global*nDOF_element+
			 ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[left_eN_global];jj++)
                {
                  j = freeLocal[left_eN_global*nDOF_element+
                                jj];
		  J = freeGlobal[left_eN_global*nDOF_element+
				 jj];
                  jacIndex = I*nFreeDOF_global+J;
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
              for(jj=0;jj<nFreeDOF_element[right_eN_global];jj++)
                {
                  j = freeLocal[right_eN_global*nDOF_element+
                                jj];
		  J = freeGlobal[right_eN_global*nDOF_element+
				 jj];
                  jacIndex = I*nFreeDOF_global+J;
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                1*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]*
                    w[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+i];
                }
            }
        }
      for(ii=0;ii<nFreeDOF_element[right_eN_global];ii++)
        {
          i = freeLocal[right_eN_global*nDOF_element+
                        ii];
	  I = freeGlobal[right_eN_global*nDOF_element+
			 ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[left_eN_global];jj++)
                {
                  j = freeLocal[left_eN_global*nDOF_element+
                                jj];
		  J = freeGlobal[left_eN_global*nDOF_element+
				 jj];
                  jacIndex = I*nFreeDOF_global+J;
                  jac[jacIndex] 
                    -= 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
              for(jj=0;jj<nFreeDOF_element[right_eN_global];jj++)
                {
                  j = freeLocal[right_eN_global*nDOF_element+
                                jj];
		  J = freeGlobal[right_eN_global*nDOF_element+
				 jj];
                  jacIndex = I*nFreeDOF_global+J;
                  jac[jacIndex] 
                    -= 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                                                1*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]*
                    w[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
            }
        }
    }
}

/**
   \brief Update the global dense Jacobian from the element boundary flux Jacobians at exterior boundaries
*/    
void updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(int nExteriorElementBoundaries_global,
								       int nDOF_element,
								       int nQuadraturePoints_elementBoundary,
								       int nElementBoundaries_element,
								       int nFreeDOF_global,
								       int* exteriorElementBoundaries,
								       int* elementBoundaryElements,
								       int* elementBoundaryLocalElementBoundaries,
								       int* nFreeDOF_element,
								       int* freeLocal,
								       int* freeGlobal,
								       double* elementBoundaryFluxJacobian,
								       double* w,
								       double* jac)
{
  int ebNE,ebN,eN_global,ebN_element,ii,i,k,jj,j,I,J,jacIndex;
  for (ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global   = elementBoundaryElements[ebN*2+0];
      ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(ii=0;ii<nFreeDOF_element[eN_global];ii++)
        {
          i = freeLocal[eN_global*nDOF_element+
			ii];
	  I = freeGlobal[eN_global*nDOF_element+
			 ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element[eN_global];jj++)
                {
                  j = freeLocal[eN_global*nDOF_element+
                                jj];
		  J = freeGlobal[eN_global*nDOF_element+
				 jj];
		  jacIndex = I*nFreeDOF_global + J;
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element + 
                                                0*nQuadraturePoints_elementBoundary*nDOF_element+
                                                k*nDOF_element+
                                                j]
                    *
                    w[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element + 
                      ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                      k*nDOF_element+
                      i];
                }
            }
        }
    }
}

/**
   \brief Calculate the total (advective + diffusive) flow  velocity.
*/
void calculateFlowVelocity(int nElements_global,
                           int nQuadraturePoints_element,
                           int nSpace,
                           double* f,
                           double* a,
                           double* grad_phi,
                           double* v)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  for  (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (I=0;I<nSpace;I++)
        {
          v[eN*nQuadraturePoints_element*nSpace+
            k*nSpace+
            I]
            =
            f[eN*nQuadraturePoints_element*nSpace+
              k*nSpace+
              I];
          for (J=0;J<nSpace;J++)
            v[eN*nQuadraturePoints_element*nSpace+
              k*nSpace+
              I]
              -=
              a[eN*nQuadraturePoints_element*nSpace2+
                k*nSpace2+
                I*nSpace+
                J]
              *
              grad_phi[eN*nQuadraturePoints_element*nSpace+
                       k*nSpace+
                       J];
        }
}

/**
   \brief Update a single  element of the Jacobian
*/
void updateAddJacobian_CSR(int jacIndex, 
                           double val, 
                           double* jac)
{
  jac[jacIndex] += val;
}

/**
   \brief  Set all the Jacobian entries  to 0.0
*/
void zeroJacobian_CSR(int nNonzeros, 
                      double* jac)
{
  int  i;
  for (i=0;i<nNonzeros;i++)
    jac[i] = 0.0;
}

/**
   \brief Update the element boundary flux on interior element boundaries
*/
void updateInteriorElementBoundaryFlux(int nInteriorElementBoundaries_global,
                                       int nQuadraturePoints_elementBoundary,
                                       int nElementBoundaries_element,
                                       int nDOF_element,
                                       int* interiorElementBoundaries,
                                       int* elementBoundaryElements,
                                       int* elementBoundaryLocalElementBoundaries,
                                       double* flux,
                                       double* w,
                                       double* residual)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,i,k;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(i=0;i<nDOF_element;i++)
        for(k=0;k<nQuadraturePoints_elementBoundary;k++)
          {
            residual[left_eN_global*nDOF_element+
		     i] 
              += 
              flux[ebN*nQuadraturePoints_elementBoundary+
                   k]*
              w[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                k*nDOF_element+
                i];
            residual[right_eN_global*nDOF_element+
		     i] 
              -= 
              flux[ebN*nQuadraturePoints_elementBoundary+
                   k]*
              w[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                k*nDOF_element+
                i];
          }
    }
}

/**
   \brief Update the element boundary flux on exterior element boundaries
*/
void updateExteriorElementBoundaryFlux(int nExteriorElementBoundaries_global,
                                       int nQuadraturePoints_elementBoundary,
                                       int nElementBoundaries_element,
                                       int nDOF_element,
                                       int* exteriorElementBoundaries,
                                       int* elementBoundaryElements,
                                       int* elementBoundaryLocalElementBoundaries,
                                       double* flux,
                                       double* w,
                                       double* residual)
{
  int ebNE,ebN,eN_global,ebN_element,i,k;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(i=0;i<nDOF_element;i++)
        for(k=0;k<nQuadraturePoints_elementBoundary;k++)
          {
            residual[eN_global*nDOF_element+
                     i] 
              += 
              flux[ebN*nQuadraturePoints_elementBoundary+
                   k]
              *
              w[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                k*nDOF_element+
                i];
          }
    }
}

/**
   \brief Calculate the advective flux at at interior element boundaries
*/
void calculateInteriorNumericalAdvectiveFlux(int nInteriorElementBoundaries_global,
                                             int nQuadraturePoints_elementBoundary,
                                             int nElementBoundaries_element,
                                             int nSpace,
                                             int* interiorElementBoundaries,
                                             int* elementBoundaryElements,
                                             int* elementBoundaryLocalElementBoundaries,
                                             double* n,
                                             double* f,
                                             double* df,
                                             double* dx_f,
                                             double* flux)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,k,J;
  double left_speed,right_speed,max_speed;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] = 0.0;
          left_speed=0.0;
          right_speed=0.0;
          for(J=0;J<nSpace;J++)
            {
              left_speed 
                += 
                n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
		  J]
                *
                df[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                   left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                   k*nSpace+
		   J];
              right_speed 
                += 
                n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
		  J]
                *
                df[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                   right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                   k*nSpace+
		   J];
            }
          max_speed = right_speed;
          if (fmax(fabs(left_speed),fabs(right_speed)) == fabs(left_speed))
            max_speed = left_speed;
          if (max_speed >= 0.0)
            {
              for(J=0;J<nSpace;J++)
                flux[ebN*nQuadraturePoints_elementBoundary+
                     k] 
                  += 
                  n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                    left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                    k*nSpace+
                    J]
                  *
                  f[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                    left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                    k*nSpace+
                    J];
            }
          else
            {
              for(J=0;J<nSpace;J++)
                flux[ebN*nQuadraturePoints_elementBoundary+
                     k] 
                  += 
                  n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                    left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                    k*nSpace+J]
                  *
                  f[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                    right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                    k*nSpace+J];
            }
	  flux[ebN*nQuadraturePoints_elementBoundary+
	       k]
	    *=
	    dx_f[ebN*nQuadraturePoints_elementBoundary+
		 k];
        }
    }
}

/**
   \brief Update the advective flux Jacobian at at interior element boundaries.
*/
void updateInteriorNumericalAdvectiveFluxJacobian(int nInteriorElementBoundaries_global,
                                                  int nQuadraturePoints_elementBoundary,
                                                  int nElementBoundaries_element,
                                                  int nDOF_element,
                                                  int nSpace,
                                                  int* interiorElementBoundaries,
                                                  int* elementBoundaryElements,
                                                  int* elementBoundaryLocalElementBoundaries,
                                                  double* n,
                                                  double* df,
						  double* v,
                                                  double* dx_f,
                                                  double* fluxJacobian)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,k,J,j;
  double left_speed,right_speed,max_speed,leftJacobian,rightJacobian;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          left_speed=0.0;
          right_speed=0.0;
          for(J=0;J<nSpace;J++)
            {
              left_speed 
                += 
                n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+J]
                *
                df[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                   left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                   k*nSpace+
		   J];
              right_speed 
                += 
                n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
		  J]
                *
                df[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                   right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                   k*nSpace+
		   J];
            }
          if (fmax(fabs(left_speed),fabs(right_speed)) == fabs(left_speed))
            max_speed = left_speed;
          else
            max_speed = right_speed;
          for(j=0;j<nDOF_element;j++)
            {
              leftJacobian=0.0;
              rightJacobian=0.0;
              {
                if (max_speed >= 0.0)
                  {
                    for(J=0;J<nSpace;J++)
                      {
                        leftJacobian 
                          += 
                          n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                            left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
			    J]
                          *
                          df[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                             left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                             k*nSpace+
			     J];
                      }
                    leftJacobian *= 
		      dx_f[ebN*nQuadraturePoints_elementBoundary+
			   k]
		      *
		      v[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
			left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
			k*nDOF_element+
			j];
		    fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
				 0*nQuadraturePoints_elementBoundary*nDOF_element+
				 k*nDOF_element+
				 j] 
		      += leftJacobian;
                  }
                else
                  {
                    for(J=0;J<nSpace;J++)
                      {
                        rightJacobian 
                          += 
                          n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                            left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
			    J]
                          *
                          df[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                             right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                             k*nSpace+
			     J];
                      }
                    rightJacobian *= 
		      dx_f[ebN*nQuadraturePoints_elementBoundary+
			   k]
		      *
		      v[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
			right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
			k*nDOF_element+
			j];
		    fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
				 1*nQuadraturePoints_elementBoundary*nDOF_element+
				 k*nDOF_element+
				 j]
		      += rightJacobian;
                  }
              }
            }
        }
    }
}

/**
   \brief Update the advective flux at exterior element boundaries.
*/
void calculateExteriorNumericalAdvectiveFlux(int nExteriorElementBoundaries_global,
                                             int nQuadraturePoints_elementBoundary,
                                             int nElementBoundaries_element,
                                             int nSpace,
                                             int* exteriorElementBoundaries,
                                             int* elementBoundaryElements,
                                             int* elementBoundaryLocalElementBoundaries,
                                             int* inflowBoundary,
                                             double* inflowFlux,
                                             double* n,
                                             double* f,
                                             double* df,
                                             double* dx_f,
                                             double* flux)
{
  int ebNE,ebN,eN_global,ebN_element,k,J;
  register double speed;
  int inflowFlag;
  memset(inflowBoundary,0,sizeof(int)*nExteriorElementBoundaries_global);
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      inflowFlag=0;
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          flux[ebN*nQuadraturePoints_elementBoundary+k] = 0.0;
          speed = 0.0;
          for(J=0;J<nSpace;J++)
            {
              speed 
                += 
                n[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
		  J]
                *
                df[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                   ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                   k*nSpace+
		   J];
              flux[ebN*nQuadraturePoints_elementBoundary+k] 
                += 
                n[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
                  J]
                *
                f[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
                  J];
            }
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] 
            *= 
            dx_f[ebN*nQuadraturePoints_elementBoundary+
                 k];
          if(speed < 0)
            inflowFlag=1;
        }
      inflowBoundary[ebNE]=inflowFlag;
      if (inflowFlag == 1)
        for(k=0;k<nQuadraturePoints_elementBoundary;k++)
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] = 
            inflowFlux[ebNE*nQuadraturePoints_elementBoundary+
                       k];
    }
}

/**
   \brief Set the advective flux boundary condition at exterior element boundaries from the current exterior flux.
*/
void setInflowFlux(int nExteriorElementBoundaries_global,
                   int nQuadraturePoints_elementBoundary,
                   int* exteriorElementBoundaries,
                   double* inflowFlux,
                   double* flux)
{
  int ebNE,ebN,k;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        inflowFlux[ebNE*nQuadraturePoints_elementBoundary+
                   k]
          =
          flux[ebN*nQuadraturePoints_elementBoundary+
               k]; 
    }
}


/**
   \brief Update the advective flux Jacobian at at exterior element boundaries.
*/
void updateExteriorNumericalAdvectiveFluxJacobian(int nExteriorElementBoundaries_global,
                                                  int nQuadraturePoints_elementBoundary,
                                                  int nElementBoundaries_element,
                                                  int nDOF_element,
                                                  int nSpace,
                                                  int* exteriorElementBoundaries,
                                                  int* elementBoundaryElements,
                                                  int* elementBoundaryLocalElementBoundaries,
                                                  int* inflowBoundary,
                                                  double* n,
                                                  double* df,
						  double* v,
                                                  double* dx_f,
                                                  double* fluxJacobian)
{
  int ebNE,ebN,eN_global,ebN_element,k,j,J;
  register double jacobian;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(j=0;j<nDOF_element;j++)
        {
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
	      jacobian=0.0;
              for(J=0;J<nSpace;J++)
                {
                  jacobian 
                    += 
                    n[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                      ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                      k*nSpace+
		      J]
                    *
                    df[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                       ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                       k*nSpace+
                       J];
                }
              jacobian *= 
		dx_f[ebN*nQuadraturePoints_elementBoundary+
		     k]
		*
		v[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
		  ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
		  k*nDOF_element+
		  j];
              fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                           0*nQuadraturePoints_elementBoundary*nDOF_element+
                           k*nDOF_element+
                           j] 
                += jacobian;
            }      
          if (inflowBoundary[ebNE] == 1)
            for(k=0;k<nQuadraturePoints_elementBoundary;k++)
              fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                           0*nQuadraturePoints_elementBoundary*nDOF_element+
                           k*nDOF_element+
                           j] = 0.0;
        }
    }
}
 
/**
   \brief Calculate the diffusive flux at interior element boundary quadrature points
*/ 
void calculateInteriorNumericalDiffusiveFlux(int nInteriorElementBoundaries_global,
                                             int nQuadraturePoints_elementBoundary,
                                             int nElementBoundaries_element,
                                             int nSpace,
                                             int* interiorElementBoundaries,
                                             int* elementBoundaryElements,
                                             int* elementBoundaryLocalElementBoundaries,
                                             double* n,
                                             double* a,
                                             double* grad_phi,
                                             double* u,
                                             double* penalty,
                                             double* dx_a,
                                             double* flux)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,k,J,I,nSpace2=nSpace*nSpace;
  double diffusiveVelocityComponent_I;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          flux[ebN*nQuadraturePoints_elementBoundary+k] = 0.0;
          for(I=0;I<nSpace;I++)
            {
              diffusiveVelocityComponent_I=0.0;
              for(J=0;J<nSpace;J++)
                {
                  diffusiveVelocityComponent_I 
                    -= 
                    (a[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                       left_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                       k*nSpace2+
                       I*nSpace+
                       J]
                     *
                     grad_phi[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                              left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                              k*nSpace+J]
                     +
                     a[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                       right_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                       k*nSpace2+
                       I*nSpace+
                       J]
                     *
                     grad_phi[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                              right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                              k*nSpace+J]);
                }
              flux[ebN*nQuadraturePoints_elementBoundary+
                   k] 
                += 
                diffusiveVelocityComponent_I
                *
                n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+
                  I];
            }
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] *= 0.5;
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] 
            += 
            penalty[ebN*nQuadraturePoints_elementBoundary+
                    k]
            *
            (u[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
               left_ebN_element*nQuadraturePoints_elementBoundary+
               k]-
             u[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
               right_ebN_element*nQuadraturePoints_elementBoundary+
               k]);
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] 
            *= 
            dx_a[ebN*nQuadraturePoints_elementBoundary+
                 k];
        }
    }
}

/**
   \brief Calculate the diffusive flux Jacobian at interior element boundary quadrature points
*/ 
void updateInteriorNumericalDiffusiveFluxJacobian(int nInteriorElementBoundaries_global,
                                                  int nQuadraturePoints_elementBoundary,
                                                  int nElementBoundaries_element,
                                                  int nDOF_element,
                                                  int nSpace,
                                                  int* l2g,
                                                  int* interiorElementBoundaries,
                                                  int* elementBoundaryElements,
                                                  int* elementBoundaryLocalElementBoundaries,
                                                  double* n,
                                                  double* a,
                                                  double* da,
                                                  double* grad_phi,
                                                  double* dphi,
                                                  double* v,
                                                  double* grad_v,
                                                  double* penalty,
                                                  double* dx_a,
                                                  double* fluxJacobian)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,k,j,left_j_global,right_j_global,I,J,nSpace2=nSpace*nSpace;
  double leftJacobian,rightJacobian,diffusiveVelocityComponent_I_leftJacobian,diffusiveVelocityComponent_I_rightJacobian,diffusiveVelocityComponent_I_leftJacobian2,diffusiveVelocityComponent_I_rightJacobian2;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          for(j=0;j<nDOF_element;j++)
            {
              leftJacobian=0.0;
              rightJacobian=0.0;
              left_j_global = l2g[left_eN_global*nDOF_element+j];
              right_j_global= l2g[right_eN_global*nDOF_element+j];
              for(I=0;I<nSpace;I++)
                {
                  diffusiveVelocityComponent_I_leftJacobian=0.0;
                  diffusiveVelocityComponent_I_leftJacobian2=0.0;
                  diffusiveVelocityComponent_I_rightJacobian=0.0;
                  diffusiveVelocityComponent_I_rightJacobian2=0.0;
                  for(J=0;J<nSpace;J++)
                    {
                      diffusiveVelocityComponent_I_leftJacobian 
                        -= 
                        da[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                           left_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                           k*nSpace2+
                           I*nSpace+
                           J]
                        *
                        grad_phi[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                 left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                 k*nSpace+
                                 J];
                      diffusiveVelocityComponent_I_rightJacobian 
                        -= 
                        da[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                           right_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                           k*nSpace2+
                           I*nSpace+
                           J]
                        *
                        grad_phi[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                 right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                 k*nSpace+
                                 J];
                      diffusiveVelocityComponent_I_leftJacobian2 
                        -= 
                        a[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                          left_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                          k*nSpace2+
                          I*nSpace+
                          J]
                        *
                        grad_v[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               k*nDOF_element*nSpace+
                               j*nSpace+
                               J];
                      diffusiveVelocityComponent_I_rightJacobian2 
                        -= 
                        a[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                          right_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                          k*nSpace2+
                          I*nSpace+
                          J]
                        *
                        grad_v[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               k*nDOF_element*nSpace+
                               j*nSpace+
                               J];
                      
                    }
                  diffusiveVelocityComponent_I_leftJacobian=0.0;
                  diffusiveVelocityComponent_I_rightJacobian=0.0;
                  leftJacobian 
                    += 
                    (diffusiveVelocityComponent_I_leftJacobian
                     *
                     v[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       k*nDOF_element+
                       j] 
                     +
                     diffusiveVelocityComponent_I_leftJacobian2*
                     dphi[left_j_global])
                    *
                    n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                      left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                      k*nSpace+
                      I];
                  rightJacobian 
                    += 
                    (diffusiveVelocityComponent_I_rightJacobian
                     *
                     v[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       k*nDOF_element+
                       j] 
                     +
                     diffusiveVelocityComponent_I_rightJacobian2*dphi[right_j_global])
                    *
                    n[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                      left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                      k*nSpace+
                      I];
                }
              leftJacobian *= 0.5;
              rightJacobian *= 0.5;
              leftJacobian 
                += 
                penalty[ebN*nQuadraturePoints_elementBoundary+
                        k]
                *
                v[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                  left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                  k*nDOF_element+
                  j];
              rightJacobian 
                -= 
                penalty[ebN*nQuadraturePoints_elementBoundary+
                        k]
                *
                v[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                  right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                  k*nDOF_element+
                  j];
              leftJacobian *= dx_a[ebN*nQuadraturePoints_elementBoundary+
                                   k];
              rightJacobian *= dx_a[ebN*nQuadraturePoints_elementBoundary+
                                    k];
              fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                           0*nQuadraturePoints_elementBoundary*nDOF_element+
                           k*nDOF_element+
                           j] += leftJacobian;
              fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                           1*nQuadraturePoints_elementBoundary*nDOF_element+
                           k*nDOF_element+
                           j] += rightJacobian;
            }
        }
    }
}
 
/**
   \brief Calculate the diffusive flux at exterior element boundary quadrature points
*/ 
void calculateExteriorNumericalDiffusiveFlux(int nExteriorElementBoundaries_global,
                                             int nQuadraturePoints_elementBoundary,
                                             int nElementBoundaries_element,
                                             int nSpace,
                                             int* exteriorElementBoundaries,
                                             int* elementBoundaryElements,
                                             int* elementBoundaryLocalElementBoundaries,
                                             double* n,
                                             double* a,
                                             double* grad_phi,
                                             double* u,
                                             double* penalty,
                                             double* dx_a,
                                             double* flux)
{
  int ebNE,ebN,eN_global,ebN_element,k,J,I,nSpace2=nSpace*nSpace;
  double diffusiveVelocityComponent_I;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          flux[ebN*nQuadraturePoints_elementBoundary+k] = 0.0;
          for(I=0;I<nSpace;I++)
            {
              diffusiveVelocityComponent_I=0.0;
              for(J=0;J<nSpace;J++)
                {
                  diffusiveVelocityComponent_I 
                    -= 
                    a[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                      ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                      k*nSpace2+
                      I*nSpace+
                      J]
                    *
                    grad_phi[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                             ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                             k*nSpace+J];
                }
              flux[ebN*nQuadraturePoints_elementBoundary+k] 
                += 
                diffusiveVelocityComponent_I
                *
                n[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                  ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                  k*nSpace+I];
            }
          /*           flux[ebN*nQuadraturePoints_elementBoundary+k] += penalty[ebN*nQuadraturePoints_elementBoundary+k]* */
          /*             (u[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+ */
          /*                ebN_element*nQuadraturePoints_elementBoundary+ */
          /*                k]- */
          /*              u[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+ */
          /*                right_ebN_element*nQuadraturePoints_elementBoundary+ */
          /*                k]); */
          flux[ebN*nQuadraturePoints_elementBoundary+
               k] 
            *= 
            dx_a[ebN*nQuadraturePoints_elementBoundary+
                 k];
        }
    }
}

/**
   \brief Update the diffusive flux Jacobian at exterior element boundary quadrature points
*/ 
void updateExteriorNumericalDiffusiveFluxJacobian(int nExteriorElementBoundaries_global,
                                                  int nQuadraturePoints_elementBoundary,
                                                  int nElementBoundaries_element,
                                                  int nDOF_element,
                                                  int nSpace,
                                                  int* l2g,
                                                  int* exteriorElementBoundaries,
                                                  int* elementBoundaryElements,
                                                  int* elementBoundaryLocalElementBoundaries,
                                                  double* n,
                                                  double* a,
                                                  double* da,
                                                  double* grad_phi,
                                                  double* dphi,
                                                  double* v,
                                                  double* grad_v,
                                                  double* penalty,
                                                  double* dx_a,
                                                  double* fluxJacobian)
{
  int ebNE,ebN,eN_global,ebN_element,k,j,j_global,I,J,nSpace2=nSpace*nSpace;
  double Jacobian,diffusiveVelocityComponent_I_Jacobian,diffusiveVelocityComponent_I_Jacobian2;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          for(j=0;j<nDOF_element;j++)
            {
              Jacobian=0.0;
              j_global = l2g[eN_global*nDOF_element+j];
              for(I=0;I<nSpace;I++)
                {
                  diffusiveVelocityComponent_I_Jacobian=0.0;
                  diffusiveVelocityComponent_I_Jacobian2=0.0;
                  for(J=0;J<nSpace;J++)
                    {
                      diffusiveVelocityComponent_I_Jacobian 
                        -= 
                        da[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                           ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                           k*nSpace2+
                           I*nSpace+
                           J]
                        *
                        grad_phi[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                 ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                 k*nSpace+
                                 J];
                      diffusiveVelocityComponent_I_Jacobian2 
                        -= 
                        a[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
                          ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
                          k*nSpace2+
                          I*nSpace+
                          J]
                        *
                        grad_v[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               ebN_element*nQuadraturePoints_elementBoundary*nDOF_element*nSpace+
                               k*nDOF_element*nSpace+
                               j*nSpace+
                               J];
                    }
                  Jacobian 
                    += 
                    (diffusiveVelocityComponent_I_Jacobian*
                     v[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+
                       k*nDOF_element+
                       j] 
                     +
                     diffusiveVelocityComponent_I_Jacobian2*
                     dphi[j_global])
                    *
                    n[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                      ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                      k*nSpace+I];
                }
              Jacobian*= 0.5;
              /*               Jacobian += penalty[ebN*nQuadraturePoints_elementBoundary+k]*v[eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_element+ */
              /*                                         ebN_element*nQuadraturePoints_elementBoundary*nDOF_element+ */
              /*                                         k*nDOF_element+ */
              /*                                         j]; */
              Jacobian *=dx_a[ebN*nQuadraturePoints_elementBoundary+
                              k];
              fluxJacobian[ebN*2*nQuadraturePoints_elementBoundary*nDOF_element+
                           0*nQuadraturePoints_elementBoundary*nDOF_element+
                           k*nDOF_element+
                           j] 
                += Jacobian;
            }
        }
    }
}
 
void calculateInteriorElementBoundaryVelocities(int nInteriorElementBoundaries_global,
                                                int nQuadraturePoints_elementBoundary,
                                                int nElementBoundaries_element,
						int nSpace,
                                                int* interiorElementBoundaries,
                                                int* elementBoundaryElements,
                                                int* elementBoundaryLocalElementBoundaries,
						double* m,
                                                double* a,
						double* grad_phi,
						double* f,
                                                double* vAverage,
                                                double* vJump,
                                                double* mAverage,
                                                double* mJump)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,k,I,J,nSpace2=nSpace*nSpace;
  register double vLeft,vRight,mLeft,mRight;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
	  mLeft 
	    = 
	    m[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
	      left_ebN_element*nQuadraturePoints_elementBoundary+
	      k];
	  mRight
	    =
	    f[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
	      right_ebN_element*nQuadraturePoints_elementBoundary+
	      k];
	  for (I=0;I<nSpace;I++)
	    {
	      vLeft 
		= 
		f[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
		  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
		  k*nSpace+
		  I];
	      vRight
		=
		f[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
		  right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
		  k*nSpace+
		  I];
	      for (J=0;J<nSpace;J++)
		{
		  vLeft 
		    -= 
		    a[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
		      left_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
		      k*nSpace2+
		      I*nSpace+
		      J]
		    *
		    grad_phi[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			     left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
			     k*nSpace+
			     J];
		  vRight -=
		    a[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
		      right_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
		      k*nSpace2+
		      I*nSpace+
		      J]
		    *
		    grad_phi[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			     right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
			     k*nSpace+
			     J];
		}
	      vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
		       k*nSpace+
		       I] = 0.5*(vLeft + vRight);
	      vJump[ebN*nQuadraturePoints_elementBoundary*nSpace+
		    k*nSpace+
		    I] = (vLeft - vRight);
	    }
	  mAverage[ebN*nQuadraturePoints_elementBoundary+
		   k] = 0.5*(mLeft+mRight);
	  mJump[ebN*nQuadraturePoints_elementBoundary+
		k] = mLeft - mRight;
	}
    }
}

void calculateExteriorElementBoundaryVelocities(int nExteriorElementBoundaries_global,
                                                int nQuadraturePoints_elementBoundary,
                                                int nElementBoundaries_element,
						int nSpace,
                                                int* exteriorElementBoundaries,
                                                int* elementBoundaryElements,
                                                int* elementBoundaryLocalElementBoundaries,
						double* m,
                                                double* a,
						double* grad_phi,
						double* f,
                                                double* vAverage,
                                                double* vJump,
                                                double* mAverage,
                                                double* mJump)
{
  int ebNE,ebN,left_eN_global,left_ebN_element,k,I,J,nSpace2=nSpace*nSpace;
  register double vLeft,mLeft;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
	  mLeft 
	    = 
	    m[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
	      left_ebN_element*nQuadraturePoints_elementBoundary+
	      k];
	  for (I=0;I<nSpace;I++)
	    {
	      vLeft 
		= 
		f[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
		  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
		  k*nSpace+
		  I];
	      for (J=0;J<nSpace;J++)
		{
		  vLeft 
		    -= 
		    a[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace2+
		      left_ebN_element*nQuadraturePoints_elementBoundary*nSpace2+
		      k*nSpace2+
		      I*nSpace+
		      J]
		    *
		    grad_phi[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			     left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
			     k*nSpace+
			     J];
		}
	      vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
		       k*nSpace+
		       I] = vLeft;
	      vJump[ebN*nQuadraturePoints_elementBoundary*nSpace+
		    k*nSpace+
		    I] = vLeft;
	    }
	  mAverage[ebN*nQuadraturePoints_elementBoundary+
		   k] = mLeft;
	  mJump[ebN*nQuadraturePoints_elementBoundary+
		k] = mLeft;
	}
    }
}

void calculateConservationResidualPWL(int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nQuadraturePoints_elementBoundary,
				      int nElementBoundaries_element,
				      int nNodes_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nodeStarOffsets,
				      int* nElements_node,
				      double* elementResidual,
				      double* vAverage,
				      double* starU,
				      double* w,
				      double* normal,
				      double* dx,
				      double* conservationResidual,
				      double* starR,
				      double* vConservative,
                                      double* vConservative_element)
{
  int ebNI,ebNE,ebN,eN,eN_star,left_eN,right_eN,ebN_element,left_ebN_element,right_ebN_element,left_eN_star,right_eN_star,nN,nN_global,k,I;
  register double flux,fluxAverage,fluxCorrection;
  memset(conservationResidual,0,sizeof(double)*nElements_global);
  /*initial residual with element residual*/
  for (eN=0;eN<nElements_global;eN++)
    for (nN=0;nN<nNodes_element;nN++)
      {
	nN_global = elementNodes[eN*nNodes_element+
				nN];
	eN_star  = nodeStarElements[eN*nNodes_element+
				    nN];
	starR[nodeStarOffsets[nN_global] +
              eN_star] 
	  = 
	  elementResidual[eN*nNodes_element+
			  nN];
	conservationResidual[eN]
	  += 
	  elementResidual[eN*nNodes_element+
			  nN];
      }
  /*calculate interior element boundary fluxes and update residual*/
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN = elementBoundaryElements[ebN*2+0];
      right_eN = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
          fluxAverage = 0.0;
	  for (I=0;I<nSpace;I++)
	    {
              vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
                            I] 
                =
                vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
                         k*nSpace+
                         I];
              fluxAverage
                +=
                vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
                         k*nSpace+
                         I]
                *
                normal[ebN*nQuadraturePoints_elementBoundary*nSpace+
                       k*nSpace+
                       I];
            }
          for (nN=0;nN<nNodes_element;nN++)
	    {
	      nN_global = elementNodes[left_eN*nNodes_element+
				      nN];
	      left_eN_star = nodeStarElements[left_eN*nNodes_element+
					      nN];
	      /* check if node is opposite element boundary we're computing and ignore 0 contribution */
	      /* this shouldn't be necessary */
	      if (nN != left_ebN_element)
		{
		  right_eN_star = nodeStarElementNeighbors[left_eN*nNodes_element*nElementBoundaries_element+
							   nN*nElementBoundaries_element+
							   left_ebN_element];
                  fluxCorrection = (starU[nodeStarOffsets[nN_global]+
                                          left_eN_star]
                                    -
                                    starU[nodeStarOffsets[nN_global]+
                                          right_eN_star])
                    *
                    w[left_eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      left_ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      k*nNodes_element+
                      nN];
		  for (I=0;I<nSpace;I++)
		    {
		      vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
				    k*nSpace+
				    I] 
                        +=
			fluxCorrection
			*
			normal[ebN*nQuadraturePoints_elementBoundary*nSpace+
                               k*nSpace+
                               I];
		    }
		  flux = (fluxAverage*w[left_eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                                        left_ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                                        k*nNodes_element+
                                        nN]
                          + fluxCorrection)*dx[ebN*nQuadraturePoints_elementBoundary+
                                               k];
		  starR[nodeStarOffsets[nN_global]+
                        left_eN_star]
		    += flux;
		  starR[nodeStarOffsets[nN_global]+
                        right_eN_star]
		    -= flux;
		  conservationResidual[left_eN] += flux;
		  conservationResidual[right_eN] -= flux;
		}
	    }
	}
    }
  /*calculate fluxes on exterior element boundaries and update residual*/
  for (ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
          fluxAverage=0.0;
	  for (I=0;I<nSpace;I++)
            {
              vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
                            I] 
                =
                vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
                         k*nSpace+
                         I];
              fluxAverage += 
                vAverage[ebN*nQuadraturePoints_elementBoundary*nSpace+
                         k*nSpace+
                         I]
                *
                normal[ebN*nQuadraturePoints_elementBoundary*nSpace+
                       k*nSpace+
                       I];
            }
	  for (nN=0;nN<nNodes_element;nN++)
	    {
	      nN_global = elementNodes[eN*nNodes_element+
				       nN];
	      eN_star = nodeStarElements[eN*nNodes_element+
					      nN];
	      /* check if this node lies opposite the element boundary whose contribution we're computing */
	      /* in that case there is no flux contribution because the test function is zero*/
	      /* however, we may be able to remove this conditional for speed */
	      if (nN != ebN_element)
		{
                  fluxCorrection = starU[nodeStarOffsets[nN_global]+
                                         eN_star]
                    *
                    w[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      k*nNodes_element+
                      nN];
		  for (I=0;I<nSpace;I++)
                    vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                                  k*nSpace+
                                  I] +=
                      fluxCorrection
                      *
                      normal[ebN*nQuadraturePoints_elementBoundary*nSpace+
                             k*nSpace+
                             I];
		  flux = (fluxAverage
                          *
                          w[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                            ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                            k*nNodes_element+
                            nN]
                          +
                          fluxCorrection)
                    *
                    dx[ebN*nQuadraturePoints_elementBoundary+
                       k];
		  starR[nodeStarOffsets[nN_global]+
                        eN_star]
		    += flux;
		  conservationResidual[eN] += flux;
		}
	    }
	}
    }
  /*copy the global velocity onto the elements*/
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN = elementBoundaryElements[ebN*2+0];
      right_eN = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(I=0;I<nSpace;I++)
          {
            vConservative_element[left_eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                  left_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                  k*nSpace+
                                  I]
              =
              vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
                            I];
            vConservative_element[right_eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                  right_ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                  k*nSpace+
                                  I]
              =
              vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                            k*nSpace+
                            I];
          }
    }
  /*copy the global velocity onto the elements*/
  for (ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        for(I=0;I<nSpace;I++)
          vConservative_element[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
                                ebN_element*nQuadraturePoints_elementBoundary*nSpace+
                                k*nSpace+
                                I]
            =
            vConservative[ebN*nQuadraturePoints_elementBoundary*nSpace+
                          k*nSpace+
                          I];
    }
}

#ifndef PYADH_CLAPACK_H
#define PYADH_CLAPACK_H "vecLib/clapack.h"
#endif

#include PYADH_CLAPACK_H

static __CLPK_integer* starJacobianPivots=NULL;

void calculateConservationJacobianPWL(int nNodes_global,
                                      int nNodes_internal,
                                      int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nQuadraturePoints_elementBoundary,
				      int nElementBoundaries_element,
				      int nNodes_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nodeStarOffsets,
				      int* nodeStarJacobianOffsets,
                                      int* nElements_node,
                                      int* internalNodes,
				      double* w,
				      double* normal,
				      double* dx,
				      double* starJacobian)
{
  int eN,ebNI,ebNE,ebN,eN_star,left_eN,right_eN,ebN_element,left_ebN_element,right_ebN_element,left_eN_star,right_eN_star,nN,nN_global,nNI,k;
  __CLPK_integer INFO=0;
  register double wflux;
  if (starJacobianPivots != NULL)
    free(starJacobianPivots);
  starJacobianPivots = calloc(sizeof(__CLPK_integer*),nodeStarOffsets[nNodes_global-1]+nElements_node[nNodes_global-1]);
  memset(starJacobian,0,sizeof(double)*(nodeStarJacobianOffsets[nNodes_global-1] + nElements_node[nNodes_global-1]*nElements_node[nNodes_global-1]));
  /*Load Jacobian entries arising from iterior element boundaries*/
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN = elementBoundaryElements[ebN*2+0];
      right_eN = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
	  for (nN=0;nN<nNodes_element;nN++)
	    {
	      nN_global = elementNodes[left_eN*nNodes_element+
				      nN];
	      left_eN_star = nodeStarElements[left_eN*nNodes_element+
					      nN];
	      /* check if node is opposite element boundary we're computing and ignore 0 contribution */
	      /* this shouldn't be necessary */
	      if (nN != left_ebN_element)
		{
		  right_eN_star = nodeStarElementNeighbors[left_eN*nNodes_element*nElementBoundaries_element+
							   nN*nElementBoundaries_element+
							   left_ebN_element];
		  wflux = w[left_eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                            left_ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                            k*nNodes_element+
                            nN]
                    *
                    dx[ebN*nQuadraturePoints_elementBoundary+
                       k];
                  
		  starJacobian[nodeStarJacobianOffsets[nN_global]+
			       left_eN_star+
			       left_eN_star*nElements_node[nN_global]]
		    += wflux;
		  starJacobian[nodeStarJacobianOffsets[nN_global]+
			       left_eN_star+
			       right_eN_star*nElements_node[nN_global]]
		    -= wflux;
		  starJacobian[nodeStarJacobianOffsets[nN_global]+
			       right_eN_star+
			       left_eN_star*nElements_node[nN_global]]
		    -= wflux;
		  starJacobian[nodeStarJacobianOffsets[nN_global]+
			       right_eN_star+
			       right_eN_star*nElements_node[nN_global]]
		    += wflux;
		}
	    }
	}
    }  
  /*Load Jacobian entries arising from exterior element boundaries*/
  for (ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	{
	  for (nN=0;nN<nNodes_element;nN++)
	    {
	      nN_global = elementNodes[eN*nNodes_element+
				       nN];
	      eN_star = nodeStarElements[eN*nNodes_element+
					 nN];
	      /* check if this node lies opposite the element boundary whose contribution we're computing */
	      /* in that case there is no flux contribution because the test function is zero*/
	      /* however, we may be able to remove this conditional for speed */
	      if (nN != ebN_element)
		{
		  starJacobian[nodeStarJacobianOffsets[nN_global]+
                               eN_star+
                               eN_star*nElements_node[nN_global]]
		    += 
                    w[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      ebN_element*nQuadraturePoints_elementBoundary*nNodes_element+
                      k*nNodes_element+
                      nN]
                    *
                    dx[ebN*nQuadraturePoints_elementBoundary+
                       k];
		}
	    }
	}
    }
  /*set Dirichlet boundary conditions at one element on interior node stars*/
  for (nNI=0;nNI<nNodes_internal;nNI++)
    {
      nN = internalNodes[nNI];
      starJacobian[nodeStarJacobianOffsets[nN]]=1.0;
      for (eN=1;eN<nElements_node[nN];eN++)
        starJacobian[nodeStarJacobianOffsets[nN]+
		     eN*nElements_node[nN]]
	  =0.0;
    }
  /*factor with lapack*/
  for (nN=0;nN<nNodes_global;nN++)
    dgetrf_(&nElements_node[nN],
            &nElements_node[nN],
            &starJacobian[nodeStarJacobianOffsets[nN]],
            &nElements_node[nN],
            &starJacobianPivots[nodeStarOffsets[nN]],
            &INFO);
}

void calculateConservationFluxPWL(int nNodes_global,
                                  int nNodes_internal,
                                  int* nElements_node,
                                  int* nodeStarOffsets,
                                  int* nodeStarJacobianOffsets,
                                  int* internalNodes,
                                  double* starR,
                                  double* starJ,
                                  double* starU)
{
  int nN,nNI,eN;
  __CLPK_integer NRHS=1,INFO=0;
  char TRANS='N';
  /*load -R into U*/
  for (nN=0;nN<nNodes_global;nN++)
    for(eN=0;eN<nElements_node[nN];eN++)
      starU[nodeStarOffsets[nN]+eN] = -starR[nodeStarOffsets[nN]+eN];
  /*set Dirichlet boundary conditions on interior node stars*/
  for (nNI=0;nNI<nNodes_internal;nNI++)
    {
      nN = internalNodes[nNI];
      starU[nodeStarOffsets[nN]]=0.0;
    }
  /*solve with lapack*/
  for  (nN=0;nN<nNodes_global;nN++)
    dgetrs_(&TRANS,
            &nElements_node[nN],
            &NRHS,
            &starJ[nodeStarJacobianOffsets[nN]],
            &nElements_node[nN],
            &starJacobianPivots[nodeStarOffsets[nN]],
            &starU[nodeStarOffsets[nN]],
            &nElements_node[nN],
            &INFO);
}

/** @} */
