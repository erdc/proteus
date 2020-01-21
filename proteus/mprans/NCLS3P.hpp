#ifndef NCLS3P_HPP
#define NCLS3P_HPP

#include "xtensor-python/pyarray.hpp"

namespace proteus
{
    class cppNCLS3P_base
    {
    public:

        virtual ~cppNCLS3P_base() = default;
        virtual void calculateResidual(//element
                                       xt::pyarray<double> mesh_trial_ref,
                                       xt::pyarray<double> mesh_grad_trial_ref,
                                       xt::pyarray<double> mesh_dof,
                                       xt::pyarray<double> meshVelocity_dof,
                                       double MOVING_DOMAIN,
                                       xt::pyarray<int> mesh_l2g,
                                       xt::pyarray<double> dV_ref,
                                       xt::pyarray<double> u_trial_ref,
                                       xt::pyarray<double> u_grad_trial_ref,
                                       xt::pyarray<double> u_test_ref,
                                       xt::pyarray<double> u_grad_test_ref,
                                       //element boundary
                                       xt::pyarray<double> mesh_trial_trace_ref,
                                       xt::pyarray<double> mesh_grad_trial_trace_ref,
                                       xt::pyarray<double> dS_ref,
                                       xt::pyarray<double> u_trial_trace_ref,
                                       xt::pyarray<double> u_grad_trial_trace_ref,
                                       xt::pyarray<double> u_test_trace_ref,
                                       xt::pyarray<double> u_grad_test_trace_ref,
                                       xt::pyarray<double> normal_ref,
                                       xt::pyarray<double> boundaryJac_ref,
                                       //physics
                                       int nElements_global,
                                       double useMetrics,
                                       double alphaBDF,
                                       int lag_shockCapturing, /*mwf not used yet*/
                                       double shockCapturingDiffusion,
                                       double sc_uref, double sc_alpha,
                                       xt::pyarray<int> u_l2g,
                                       xt::pyarray<double> elementDiameter,
                                       xt::pyarray<double> u_dof,xt::pyarray<double> u_dof_old,
                                       xt::pyarray<double> velocity,
                                       xt::pyarray<double> q_m,
                                       xt::pyarray<double> q_u,
                                       xt::pyarray<double> q_n,
                                       xt::pyarray<double> q_dH,
                                       xt::pyarray<double> q_m_betaBDF,
                                       xt::pyarray<double> q_dV,
                                       xt::pyarray<double> q_dV_last,
                                       xt::pyarray<double> cfl,
                                       xt::pyarray<double> q_numDiff_u,
                                       xt::pyarray<double> q_numDiff_u_last,
                                       int offset_u, int stride_u,
                                       xt::pyarray<double> globalResidual,
                                       int nExteriorElementBoundaries_global,
                                       xt::pyarray<int> exteriorElementBoundariesArray,
                                       xt::pyarray<int> elementBoundaryElementsArray,
                                       xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                       xt::pyarray<double> ebqe_velocity_ext,
                                       xt::pyarray<int> isDOFBoundary_u,
                                       xt::pyarray<double> ebqe_rd_u_ext,
                                       xt::pyarray<double> ebqe_bc_u_ext,
                                       xt::pyarray<double> ebqe_u,
                                       xt::pyarray<double> cell_interface_locator,
				       xt::pyarray<double> interface_locator,
                                       // FOR TAYLOR GALERKIN METHODS                                   
                                       int EXPLICIT_METHOD,
				       double degree_polynomial,
                                       int stage,
				       xt::pyarray<double> uTilde_dof,
                                       double dt,
                                       // TO KILL SUPG AND SHOCK CAPTURING
                                       int PURE_BDF)=0;
        virtual void calculateJacobian(//element
                                       xt::pyarray<double> mesh_trial_ref,
                                       xt::pyarray<double> mesh_grad_trial_ref,
                                       xt::pyarray<double> mesh_dof,
                                       xt::pyarray<double> mesh_velocity_dof,
                                       double MOVING_DOMAIN,
                                       xt::pyarray<int> mesh_l2g,
                                       xt::pyarray<double> dV_ref,
                                       xt::pyarray<double> u_trial_ref,
                                       xt::pyarray<double> u_grad_trial_ref,
                                       xt::pyarray<double> u_test_ref,
                                       xt::pyarray<double> u_grad_test_ref,
                                       //element boundary
                                       xt::pyarray<double> mesh_trial_trace_ref,
                                       xt::pyarray<double> mesh_grad_trial_trace_ref,
                                       xt::pyarray<double> dS_ref,
                                       xt::pyarray<double> u_trial_trace_ref,
                                       xt::pyarray<double> u_grad_trial_trace_ref,
                                       xt::pyarray<double> u_test_trace_ref,
                                       xt::pyarray<double> u_grad_test_trace_ref,
                                       xt::pyarray<double> normal_ref,
                                       xt::pyarray<double> boundaryJac_ref,
                                       //physics
                                       int nElements_global,
                                       double useMetrics,
                                       double alphaBDF,
                                       int lag_shockCapturing,/*mwf not used yet*/
                                       double shockCapturingDiffusion,
                                       xt::pyarray<int> u_l2g,
                                       xt::pyarray<double> elementDiameter,
                                       xt::pyarray<double> u_dof,
                                       xt::pyarray<double> velocity,
                                       xt::pyarray<double> q_m_betaBDF,
                                       xt::pyarray<double> cfl,
                                       xt::pyarray<double> q_numDiff_u_last,
                                       xt::pyarray<int> csrRowIndeces_u_u,xt::pyarray<int> csrColumnOffsets_u_u,
                                       xt::pyarray<double> globalJacobian,
                                       int nExteriorElementBoundaries_global,
                                       xt::pyarray<int> exteriorElementBoundariesArray,
                                       xt::pyarray<int> elementBoundaryElementsArray,
                                       xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                       xt::pyarray<double> ebqe_velocity_ext,
                                       xt::pyarray<int> isDOFBoundary_u,
                                       xt::pyarray<double> ebqe_rd_u_ext,
                                       xt::pyarray<double> ebqe_bc_u_ext,
                                       xt::pyarray<int> csrColumnOffsets_eb_u_u,
                                       // FOR TAYLOR GALERKIN METHODS
                                       int EXPLICIT_METHOD,
                                       // TO KILL SUPG AND SHOCK CAPTURING
                                       int PURE_BDF) = 0;
        virtual void calculateWaterline(//element
                                        xt::pyarray<int> wlc,
                                        xt::pyarray<double> waterline,
                                        xt::pyarray<double> mesh_trial_ref,
                                        xt::pyarray<double> mesh_grad_trial_ref,
                                        xt::pyarray<double> mesh_dof,
                                        xt::pyarray<double> meshVelocity_dof,
                                        double MOVING_DOMAIN,
                                        xt::pyarray<int> mesh_l2g,
                                        xt::pyarray<double> dV_ref,
                                        xt::pyarray<double> u_trial_ref,
                                        xt::pyarray<double> u_grad_trial_ref,
                                        xt::pyarray<double> u_test_ref,
                                        xt::pyarray<double> u_grad_test_ref,
                                        //element boundary
                                        xt::pyarray<double> mesh_trial_trace_ref,
                                        xt::pyarray<double> mesh_grad_trial_trace_ref,
                                        xt::pyarray<double> dS_ref,
                                        xt::pyarray<double> u_trial_trace_ref,
                                        xt::pyarray<double> u_grad_trial_trace_ref,
                                        xt::pyarray<double> u_test_trace_ref,
                                        xt::pyarray<double> u_grad_test_trace_ref,
                                        xt::pyarray<double> normal_ref,
                                        xt::pyarray<double> boundaryJac_ref,
                                        //physics
                                        int nElements_global,
                                        double useMetrics,
                                        double alphaBDF,
                                        int lag_shockCapturing, /*mwf not used yet*/
                                        double shockCapturingDiffusion,
                                        double sc_uref, double sc_alpha,
                                        xt::pyarray<int> u_l2g,
                                        xt::pyarray<double> elementDiameter,
                                        xt::pyarray<double> u_dof,xt::pyarray<double> u_dof_old,
                                        xt::pyarray<double> velocity,
                                        xt::pyarray<double> q_m,
                                        xt::pyarray<double> q_u,
                                        xt::pyarray<double> q_n,
                                        xt::pyarray<double> q_dH,
                                        xt::pyarray<double> q_m_betaBDF,
                                        xt::pyarray<double> cfl,
                                        xt::pyarray<double> q_numDiff_u,
                                        xt::pyarray<double> q_numDiff_u_last,
                                        int offset_u, int stride_u,
                                        int nExteriorElementBoundaries_global,
                                        xt::pyarray<int> exteriorElementBoundariesArray,
                                        xt::pyarray<int> elementBoundaryElementsArray,
                                        xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                        xt::pyarray<int> elementBoundaryMaterialTypes,
                                        xt::pyarray<double> ebqe_velocity_ext,
                                        xt::pyarray<int> isDOFBoundary_u,
                                        xt::pyarray<double> ebqe_bc_u_ext,
                                        xt::pyarray<double> ebqe_u) = 0;
    };

    cppNCLS3P_base* newNCLS3P(int nSpaceIn,
                              int nQuadraturePoints_elementIn,
                              int nDOF_mesh_trial_elementIn,
                              int nDOF_trial_elementIn,
                              int nDOF_test_elementIn,
                              int nQuadraturePoints_elementBoundaryIn,
                              int CompKernelFlag);
}

#endif

