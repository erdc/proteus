#ifndef VOF3P_H
#define VOF3P_H

#include "xtensor-python/pyarray.hpp"

namespace proteus
{
    class cppVOF3P_base
    {
    public:

        virtual ~cppVOF3P_base() = default;

        virtual void calculateResidualElementBased(//element
                                                   double dt,
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
                                                   int lag_shockCapturing,
                                                   double shockCapturingDiffusion,
                                                   double sc_uref,
                                                   double sc_alpha,
                                                   //VRANS
                                                   const xt::pyarray<double>& q_vos,
                                                   const xt::pyarray<double>& vos_dof,
                                                   //
                                                   xt::pyarray<int> u_l2g,
                                                   xt::pyarray<int> r_l2g,
                                                   xt::pyarray<double> elementDiameter,
                                                   double degree_polynomial,
                                                   xt::pyarray<double> u_dof,
                                                   xt::pyarray<double> u_dof_old,
                                                   xt::pyarray<double> velocity,
                                                   xt::pyarray<double> q_m,
                                                   xt::pyarray<double> q_u,
                                                   xt::pyarray<double> q_m_betaBDF,
                                                   xt::pyarray<double> q_dV,
                                                   xt::pyarray<double> q_dV_last,
                                                   xt::pyarray<double> cfl,
                                                   xt::pyarray<double> edge_based_cfl,
                                                   xt::pyarray<double> q_numDiff_u,
                                                   xt::pyarray<double> q_numDiff_u_last,
                                                   int offset_u, int stride_u,
                                                   xt::pyarray<double> globalResidual,
                                                   int nExteriorElementBoundaries_global,
                                                   xt::pyarray<int> exteriorElementBoundariesArray,
                                                   xt::pyarray<int> elementBoundaryElementsArray,
                                                   xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                                   xt::pyarray<double> ebqe_velocity_ext,
                                                   //VRANS
                                                   const xt::pyarray<double>& ebqe_vos_ext,
                                                   //
                                                   xt::pyarray<int> isDOFBoundary_u,
                                                   xt::pyarray<double> ebqe_bc_u_ext,
                                                   xt::pyarray<int> isFluxBoundary_u,
                                                   xt::pyarray<double> ebqe_bc_flux_u_ext,
                                                   xt::pyarray<double> ebqe_phi,double epsFact,
                                                   xt::pyarray<double> ebqe_u,
                                                   xt::pyarray<double> ebqe_flux,
                                                   // TAYLOR GALERKIN
                                                   int stage,
                                                   double * uTilde_dof,
                                                   // PARAMETERS FOR ENTROPY VISCOSITY
                                                   double cE,
                                                   double cMax,
                                                   double cK,
                                                   // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
                                                   double uL,
                                                   double uR,
                                                   // PARAMETERS FOR EDGE VISCOSITY
                                                   int numDOFs,
                                                   int NNZ,
                                                   xt::pyarray<int> csrRowIndeces_DofLoops,
                                                   xt::pyarray<int> csrColumnOffsets_DofLoops,
                                                   xt::pyarray<int> csrRowIndeces_CellLoops,
                                                   xt::pyarray<int> csrColumnOffsets_CellLoops,
                                                   xt::pyarray<int> csrColumnOffsets_eb_CellLoops,
                                                   xt::pyarray<double> ML,
                                                   // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
                                                   int LUMPED_MASS_MATRIX,
                                                   int STABILIZATION_TYPE,
                                                   int ENTROPY_TYPE,
                                                   // FOR FCT
                                                   xt::pyarray<double> uLow,
                                                   xt::pyarray<double> dLow,
                                                   xt::pyarray<double> dt_times_dH_minus_dL,
                                                   xt::pyarray<double> min_u_bc,
                                                   xt::pyarray<double> max_u_bc,
                                                   // AUX QUANTITIES OF INTEREST
                                                   xt::pyarray<double> quantDOFs) = 0;

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
                                       //VRANS
                                       const xt::pyarray<double>& q_vos,
                                       //
                                       xt::pyarray<int> u_l2g,
                                       xt::pyarray<int> r_l2g,
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
                                       //VRANS
                                       const xt::pyarray<double>& ebqe_vos_ext,
                                       //
                                       xt::pyarray<int> isDOFBoundary_u,
                                       xt::pyarray<double> ebqe_bc_u_ext,
                                       xt::pyarray<int> isFluxBoundary_u,
                                       xt::pyarray<double> ebqe_bc_flux_u_ext,
                                       xt::pyarray<int> csrColumnOffsets_eb_u_u,
                                       int STABILIZATION_TYPE) = 0;

        virtual void FCTStep(double dt,
                             int NNZ, //number on non-zero entries on sparsity pattern
                             int numDOFs, //number of DOFs
                             xt::pyarray<double> lumped_mass_matrix, //lumped mass matrix (as vector)
                             xt::pyarray<double> soln, //DOFs of solution at time tn
                             xt::pyarray<double> solH, //DOFs of high order solution at tnp1
                             xt::pyarray<double> uLow,
                             xt::pyarray<double> dLow,
                             xt::pyarray<double> limited_solution,
                             xt::pyarray<int> csrRowIndeces_DofLoops, //csr row indeces
                             xt::pyarray<int> csrColumnOffsets_DofLoops, //csr column offsets
                             xt::pyarray<double> MassMatrix, //mass matrix
                             xt::pyarray<double> dt_times_dH_minus_dL, //low minus high order dissipative matrices
                             xt::pyarray<double> min_u_bc, //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
                             xt::pyarray<double> max_u_bc,
                             int LUMPED_MASS_MATRIX,
                             int STABILIZATION_TYPE) = 0;

        virtual void calculateResidualEdgeBased(//element
                                                double dt,
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
                                                int lag_shockCapturing,
                                                double shockCapturingDiffusion,
                                                double sc_uref,
                                                double sc_alpha,
                                                //VRANS
                                                const xt::pyarray<double>& q_vos,
                                                const xt::pyarray<double>& vos_dof,
                                                //
                                                xt::pyarray<int> u_l2g,
                                                xt::pyarray<int> r_l2g,
                                                xt::pyarray<double> elementDiameter,
                                                double degree_polynomial,
                                                xt::pyarray<double> u_dof,
                                                xt::pyarray<double> u_dof_old,
                                                xt::pyarray<double> velocity,
                                                xt::pyarray<double> q_m,
                                                xt::pyarray<double> q_u,
                                                xt::pyarray<double> q_m_betaBDF,
                                                xt::pyarray<double> q_dV,
                                                xt::pyarray<double> q_dV_last,
                                                xt::pyarray<double> cfl,
                                                xt::pyarray<double> edge_based_cfl,
                                                xt::pyarray<double> q_numDiff_u,
                                                xt::pyarray<double> q_numDiff_u_last,
                                                int offset_u, int stride_u,
                                                xt::pyarray<double> globalResidual,
                                                int nExteriorElementBoundaries_global,
                                                xt::pyarray<int> exteriorElementBoundariesArray,
                                                xt::pyarray<int> elementBoundaryElementsArray,
                                                xt::pyarray<int> elementBoundaryLocalElementBoundariesArray,
                                                xt::pyarray<double> ebqe_velocity_ext,
                                                //VRANS
                                                const xt::pyarray<double>& ebqe_vos_ext,
                                                //
                                                xt::pyarray<int> isDOFBoundary_u,
                                                xt::pyarray<double> ebqe_bc_u_ext,
                                                xt::pyarray<int> isFluxBoundary_u,
                                                xt::pyarray<double> ebqe_bc_flux_u_ext,
                                                xt::pyarray<double> ebqe_phi,double epsFact,
                                                xt::pyarray<double> ebqe_u,
                                                xt::pyarray<double> ebqe_flux,
                                                //EXPLICIT METHODS
                                                int stage,
                                                double * uTilde_dof,
                                                // PARAMETERS FOR EDGE BASED STABILIZATION
                                                double cE,
                                                double cMax,
                                                double cK,
                                                // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
                                                double uL,
                                                double uR,
                                                // PARAMETERS FOR EDGE VISCOSITY
                                                int numDOFs,
                                                int NNZ,
                                                xt::pyarray<int> csrRowIndeces_DofLoops,
                                                xt::pyarray<int> csrColumnOffsets_DofLoops,
                                                xt::pyarray<int> csrRowIndeces_CellLoops,
                                                xt::pyarray<int> csrColumnOffsets_CellLoops,
                                                xt::pyarray<int> csrColumnOffsets_eb_CellLoops,
                                                xt::pyarray<double> ML,
                                                // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
                                                int LUMPED_MASS_MATRIX,
                                                int STABILIZATION_TYPE,
                                                int ENTROPY_TYPE,
                                                // FOR FCT
                                                xt::pyarray<double> uLow,
                                                xt::pyarray<double> dLow,
                                                xt::pyarray<double> dt_times_dH_minus_dL,
                                                xt::pyarray<double> min_u_bc,
                                                xt::pyarray<double> max_u_bc,
                                                // AUX QUANTITIES OF INTEREST
                                                xt::pyarray<double> quantDOFs) = 0;
    };

    cppVOF3P_base* newVOF3P(int nSpaceIn,
                            int nQuadraturePoints_elementIn,
                            int nDOF_mesh_trial_elementIn,
                            int nDOF_trial_elementIn,
                            int nDOF_test_elementIn,
                            int nQuadraturePoints_elementBoundaryIn,
                            int CompKernelFlag);

}

#endif //VOF3P_H
