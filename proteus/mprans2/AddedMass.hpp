#ifndef ADDEDMASS_HPP
#define ADDEDMASS_HPP

#include "xtensor-python/pyarray.hpp"

namespace proteus
{

    class cppAddedMass_base
    {
    public:
    
        virtual ~cppAddedMass_base() = default;

        virtual void calculateResidual(//element
                                       xt::pyarray<double> mesh_trial_ref,
                                       xt::pyarray<double> mesh_grad_trial_ref,
                                       xt::pyarray<double> mesh_dof,
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
                                       int nElementBoundaries_owned,
                                       xt::pyarray<int> u_l2g,
                                       xt::pyarray<double> u_dof,
                                       xt::pyarray<double> q_rho,
                                       int offset_u,
                                       int stride_u,
                                       xt::pyarray<double> globalResidual,
                                       int nExteriorElementBoundaries_global,
                                       xt::pyarray<int> exteriorElementBoundariesArray,
                                       xt::pyarray<int> elementBoundaryElementsArray,
                                       xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                       xt::pyarray<int> elementBoundaryMaterialTypesArray,
                                       xt::pyarray<double> Aij,
                                       int added_mass_i,
                                       xt::pyarray<double> barycenters,
                                       xt::pyarray<int> flags_rigidbody) = 0;
    
        virtual void calculateJacobian(//element
                                       xt::pyarray<double> mesh_trial_ref,
                                       xt::pyarray<double> mesh_grad_trial_ref,
                                       xt::pyarray<double> mesh_dof,
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
                                       xt::pyarray<int> u_l2g,
                                       xt::pyarray<double> u_dof,
                                       xt::pyarray<double> q_rho,
                                       xt::pyarray<int> csrRowIndeces_u_u,
                                       xt::pyarray<int> csrColumnOffsets_u_u,
                                       xt::pyarray<double> globalJacobian,
                                       int nExteriorElementBoundaries_global,
                                       xt::pyarray<int> exteriorElementBoundariesArray,
                                       xt::pyarray<int> elementBoundaryElementsArray,
                                       xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                       xt::pyarray<int> csrColumnOffsets_eb_u_u) = 0;
    };

    cppAddedMass_base* newAddedMass(int nSpaceIn,
				    int nQuadraturePoints_elementIn,
				    int nDOF_mesh_trial_elementIn,
				    int nDOF_trial_elementIn,
				    int nDOF_test_elementIn,
				    int nQuadraturePoints_elementBoundaryIn,
				    int CompKernelFlag);

}

#endif

