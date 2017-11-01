from proteus import *
from proteus.default_p import *
from math import *
from proteus.mprans import NCLS
import numpy as np

run_times = 0
dij = []
Cx_T = []
Cy_T = []
ML_new = []


def getResidual(
    # double, self.timeIntegration.dt
    dt,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.psi,
    mesh_trial_ref,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.grad_psi
    mesh_grad_trial_ref,
    # numpy.ndarray, self.mesh.nodeArray
    mesh_dof,
    # numpy.ndarray, self.mesh.nodeVelocityArray
    mesh_velocity_dof,
    # double, self.MOVING_DOMAIN
    MOVING_DOMAIN,
    # numpy.ndarray, self.mesh.elementNodesArray
    mesh_l2g,
    # numpy.ndarray, self.elementQuadratureWeights[('u',0)]
    dV_ref,
    # numpy.ndarray, self.u[0].femSpace.psi
    u_trial_ref,
    # numpy.ndarray, self.u[0].femSpace.grad_psi
    u_grad_trial_ref,
    # numpy.ndarray, self.u[0].femSpace.psi
    u_test_ref,
    # numpy.ndarray, self.u[0].femSpace.grad_psi
    u_grad_test_ref,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.psi_trace
    mesh_trial_trace_ref,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.grad_psi_trace
    mesh_grad_trial_trace_ref,
    # numpy.ndarray, self.elementBoundaryQuadratureWeights[('u',0)]
    dS_ref,
    # numpy.ndarray, self.u[0].femSpace.psi_trace
    u_trial_trace_ref,
    # numpy.ndarray, self.u[0].femSpace.grad_psi_trace
    u_grad_trial_trace_ref,
    # numpy.ndarray, self.u[0].femSpace.psi_trace
    u_test_trace_ref,
    # numpy.ndarray, self.u[0].femSpace.grad_psi_trace
    u_grad_test_trace_ref,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.boundaryNormals
    normal_ref,
    # numpy.ndarray, self.u[0].femSpace.elementMaps.boundaryJacobians
    boundaryJac_ref,
    # int, self.mesh.nElements_global
    nElements_global,
    # double,  self.coefficients.useMetrics
    useMetrics,
    # double, self.timeIntegration.alpha_bdf
    alphaBDF,
    # int, self.shockCapturing.lag
    lag_shockCapturing,
    # double, self.shockCapturing.shockCapturingFactor
    shockCapturingDiffusion,
    # double, self.coefficients.sc_uref
    sc_uref,
    # double, self.coefficients.sc_beta
    sc_alpha,
    # numpy.ndarray, self.u[0].femSpace.dofMap.l2g
    u_l2g,
    # numpy.ndarray, self.mesh.elementDiametersArray,
    elementDiameter,
    # numpy.ndarray, self.mesh.nodeDiametersArray,
    nodeDiametersArray,
    # int, degree_polynomial,
    degree_polynomial,
    # numpy.ndarray, self.u[0].dof,
    u_dof,
    # numpy.ndarray, self.u_dof_old, #DOFs at lstage. Used only when
    # STABILIZATION_TYPE>0; i.e., EV
    u_dof_old,
    # numpy.ndarray, self.uStar_dof,
    uStar_dof,
    # numpy.ndarray, self.coefficients.q_v,
    velocity,
    # numpy.ndarray, self.timeIntegration.m_tmp[0],
    q_m,
    # numpy.ndarray, self.q[('u',0)],
    q_u,
    # numpy.ndarray, self.q[('grad(u)',0)],
    q_n,
    # numpy.ndarray, self.q[('dH_sge',0,0)],
    q_dH,
    # numpy.ndarray, self.timeIntegration.beta_bdf[0],#betaBDF. Used only when
    # STABILIZATION_TYPE=0
    q_m_betaBDF,
    # numpy.ndarray, self.q['dV'],
    q_dV,
    # numpy.ndarray, self.q['dV_last'],
    q_dV_last,
    # numpy.ndarray, self.q[('cfl',0)],
    cfl,
    # numpy.ndarray, self.edge_based_cfl,
    edge_based_cfl,
    # numpy.ndarray, self.shockCapturing.numDiff[0],
    q_numDiff_u,
    # numpy.ndarray, self.shockCapturing.numDiff_last[0],
    q_numDiff_u_last,
    # int, self.offset[0]
    offset_u,
    # int, self.stride[0],
    stride_u,
    # numpy.ndarray, r,
    globalResidual,
    # int, self.mesh.nExteriorElementBoundaries_global,
    nExteriorElementBoundaries_global,
    # numpy.ndarray, self.mesh.exteriorElementBoundariesArray,
    exteriorElementBoundariesArray,
    # numpy.ndarray, self.mesh.elementBoundaryElementsArray,
    elementBoundaryElementsArray,
    # numpy.ndarray, self.mesh.elementBoundaryLocalElementBoundariesArray,
    elementBoundaryLocalElementBoundariesArray,
    # numpy.ndarray, self.coefficients.ebqe_v,
    ebqe_velocity_ext,
    # numpy.ndarray, self.numericalFlux.isDOFBoundary[0],
    isDOFBoundary_u,
    # numpy.ndarray, self.coefficients.rdModel.ebqe[('u',0)],
    ebqe_rd_u_ext,
    # numpy.ndarray, self.numericalFlux.ebqe[('u',0)],
    ebqe_bc_u_ext,
    # numpy.ndarray, self.ebqe[('u',0)],
    ebqe_u,
    # int, len(rowptr)-1,
    numDOFs,
    # int, self.nnz, number of non-zero entries in sparse matrix
    NNZ,
    # numpy.ndarray, rowptr, #Row indices for Sparsity Pattern (convenient for
    # DOF loops)
    csrRowIndeces_DofLoops,
    # numpy.ndarray, colind, #Column indices for Sparsity Pattern (convenient
    # for DOF loops)
    csrColumnOffsets_DofLoops,
    # numpy.ndarray, self.csrRowIndeces[(0,0)], #row indices (convenient for
    # element loops)
    csrRowIndeces_CellLoops,
    # numpy.ndarray, self.csrColumnOffsets[(0,0)], #column indices (convenient
    # for element loops)
    csrColumnOffsets_CellLoops,
    # numpy.ndarray, self.csrColumnOffsets_eb[(0, 0)], #indices for boundary
    # terms
    csrColumnOffsets_eb_CellLoops,
    # int, self.coefficients.LUMPED_MASS_MATRIX,
    LUMPED_MASS_MATRIX,
    # numpy.ndarray, self.quantDOFs,
    quantDOFs,
    # double, self.coefficients.lambda_coupez,
    lambda_coupez,
    # double, self.coefficients.epsCoupez,
    epsCoupez,
    # double, self.coefficients.epsFactRedistancing*self.mesh.h,
    epsFactRedistancing,
    # int, self.coefficients.COUPEZ,
    COUPEZ,
    # int, self.coefficients.SATURATED_LEVEL_SET,
    SATURATED_LEVEL_SET,
    # numpy.ndarray, Cx, 1d array
    Cx,
    # numpy.ndarray, Cy
    Cy,
    # numpy.ndarray, Cz
    Cz,
    # numpy.ndarray, self.ML,
    ML,
    # int, self.coefficients.STABILIZATION_TYPE,
    STABILIZATION_TYPE,
    # int, self.coefficients.ENTROPY_TYPE,
    ENTROPY_TYPE,
    # double, self.coefficients.cE
    cE,
    ad_function,
        moving_function):

    global run_times, dij, Cx_T, Cy_T, ML_new

    Cx[:] = 0.0
    Cy[:] = 0.0

    if run_times == 0:
        Cx_T = np.zeros_like(Cx, 'd')
        Cy_T = np.zeros_like(Cy, 'd')
        dij = np.zeros_like(Cx, 'd')
        ML_new = np.zeros_like(ML, 'd')
    else:
        Cx_T[:] = 0.0
        Cy_T[:] = 0.0
        dij[:] = 0.0
        ML_new[:] = 0.0

    # Since the mesh is new
    for e in xrange(nElements_global):
        qpt, J, invJ, invJT, detJ = P1_calculateMapping_element(
            mesh_dof[mesh_l2g[e]], mesh_trial_ref, mesh_grad_trial_ref)

        dV = np.abs(detJ) * dV_ref

        n_pts = len(dV)

#             u = np.dot(u_trial_ref, u_dof[u_l2g[e]])
#             ux = np.dot(u_grad_trial_ref[:, :, 0], u_dof[u_l2g[e]])
#             uy = np.dot(u_grad_trial_ref[:, :, 1], u_dof[u_l2g[e]])

        v_basis = u_trial_ref
        v_grad_basis = [np.dot(u_grad_trial_ref[k, :, :], invJ[k, :, :])
                        for k in xrange(n_pts)]

        w_basis = u_test_ref
        w_grad_basis = [np.dot(u_grad_test_ref[k, :, :], invJ[k, :, :])
                        for k in xrange(n_pts)]

        # Computer cell cfl number
        ele_max_speed = 1e-10
        for i in range(3):
            dof_i = u_l2g[e, i]
            vi = velocity[dof_i][:-1]
            wi = mesh_velocity_dof[dof_i][:-1]
            vi = vi - wi
            ele_max_speed = max(
                [np.sqrt(vi[0] * vi[0] + vi[1] * vi[1]), np.sqrt(wi[0] * wi[0] + wi[1] * wi[1]), ele_max_speed])
        cfl[e, :] = ele_max_speed / elementDiameter[e]

        # compute C matrix
        c_x = np.zeros((3, 3), 'd')
        c_y = np.zeros((3, 3), 'd')

        dij_e = np.zeros((3, 3), 'd')
        for i in range(3):
            for j in range(3):
                for k in range(n_pts):
                    c_x[i, j] += v_basis[k, i] * \
                        v_grad_basis[k][j, 0] * dV[k]
                    c_y[i, j] += v_basis[k, i] * \
                        v_grad_basis[k][j, 1] * dV[k]
#                     c_x[i, j] += v_basis[k, i] * v_basis[k, j] * dV[k]
#                     c_y[i, j] += v_basis[k, j] * v_basis[k, i] * dV[k]
                    dij_e[i, j] += np.dot(v_grad_basis[k][j, :],
                                          v_grad_basis[k][i, :]) * dV[k]
                # end-of-loop-over-k

                Cx[csrRowIndeces_CellLoops[e, i]
                   +
                   csrColumnOffsets_CellLoops[e, i, j]] += c_x[i, j]
                Cy[csrRowIndeces_CellLoops[e, i]
                   +
                   csrColumnOffsets_CellLoops[e, i, j]] += c_y[i, j]

                Cx_T[csrRowIndeces_CellLoops[e, i]
                     +
                     csrColumnOffsets_CellLoops[e, i, j]] += c_x[j, i]
                Cy_T[csrRowIndeces_CellLoops[e, i]
                     +
                     csrColumnOffsets_CellLoops[e, i, j]] += c_y[j, i]

                dij[csrRowIndeces_CellLoops[e, i]
                    +
                    csrColumnOffsets_CellLoops[e, i, j]] += 0.1 * elementDiameter[e] * dij_e[i, j]
            # end-of-loop-over-j
        # end-of-loop-over-i
    # end-of-loop-over-e
    #dij = np.zeros_like(Cx, 'd')

    run_times += 1

    globalResidual[:] = 0.0

    edge_based_cfl[:] = 0.0

    for dof_i in xrange(numDOFs):  # Serious error: numDOFs=nnz-1

        # get new APPROXIMATE lumped mass matrix
        div_w_i = 0.0
        for j in xrange(csrRowIndeces_DofLoops[dof_i], csrRowIndeces_DofLoops[dof_i + 1]):
            dof_j = csrColumnOffsets_DofLoops[j]

            div_w_i += dt * (Cx[j] * mesh_velocity_dof[dof_j, 0] +
                             Cy[j] * mesh_velocity_dof[dof_j, 1])  # assume dof of function is equal to dof of node

        ML_new[dof_i] = ML[dof_i] + dt * div_w_i

        # main problem
        if div_w_i < 0.0:
            #             import pdb
            #             pdb.set_trace()
            #             print dof_i, ML[dof_i], ML_new[dof_i]
            edge_based_cfl[dof_i] = max(
                [edge_based_cfl[dof_i],  -div_w_i / ML[dof_i]])

        if ML_new[dof_i] < 0.0:
            raise Exception(`ML_new[dof_i]`, `ML[dof_i]`, `dof_i`, `div_w_i`)

        dii = 0.0
        spatial_residual_dof_i = 0.0
        dof_ii_index = 0
        for j in xrange(csrRowIndeces_DofLoops[dof_i], csrRowIndeces_DofLoops[dof_i + 1]):
            dof_j = csrColumnOffsets_DofLoops[j]
            if dof_i != dof_j:

                cij = np.array([Cx[j], Cy[j]])
                cji = np.array([Cx_T[j], Cy_T[j]])

                vj = velocity[dof_j, :-1]  # =[uj,uj]
                vi = velocity[dof_i, :-1]  # =[ui,ui]

                ui = u_dof_old[dof_i]
                uj = u_dof_old[dof_j]

                wj = mesh_velocity_dof[dof_j][:-1]
                wi = mesh_velocity_dof[dof_i][:-1]

                fi = np.array([0.5 * ui * ui, 0.5 * ui * ui]) - ui * wi
                fj = np.array([0.5 * uj * uj, 0.5 * uj * uj]) - uj * wj

                # maximum wave speed

                dij[j] = max([fabs(cij[0] + cij[1]) * maximum_wave_speed(ui, uj, np.dot(cij, wj) / (cij[0] + cij[1] + 1e-8)),  # for Burgers equation
                              fabs(cji[0] + cji[1]) * maximum_wave_speed(uj, ui, np.dot(cji, wi) / (cji[0] + cji[1] + 1e-8))])  # for symmetry

                dii -= dij[j]

                # This is because sum dij = 0.
                # It implies LED.
                # The value of GP method is to use 1d Riemann problem to get dt condition for total discrete method,
                # while Kuzmin's argument is from ODE and LED.
                # They are almost the same, the code may be exact same. Numerical viscosity larger is better.
                # The difference is small.
                # Both are low-order method.
                spatial_residual_dof_i += dij[j] * (
                    u_dof_old[dof_j] - u_dof_old[dof_i]) - np.dot(cij, fj - fi)
            else:
                dof_ii_index = j
        # end-of-loop-over-j

        quantDOFs[dof_i] = dij[dof_ii_index] = dii

        globalResidual[dof_i] = (u_dof_old[dof_i] * ML[dof_i] +
                                 dt * spatial_residual_dof_i) / ML_new[dof_i]
        ML[dof_i] = ML_new[dof_i]

        edge_based_cfl[dof_i] = max(
            [edge_based_cfl[dof_i],  2.0 * fabs(dii) / ML[dof_i]])
    # end-of-loop-over-dof_i


def P1_calculateMapping_element(nodes_coord, geo_basis, grad_geo_basis):
    n_pts = geo_basis.shape[0]
    J = np.zeros((n_pts, 2, 2), 'd')
    invJ = np.zeros((n_pts, 2, 2), 'd')
    invJT = np.zeros((n_pts, 2, 2), 'd')
    detJ = np.zeros((n_pts,), 'd')

    quad_pts = np.dot(geo_basis, nodes_coord)  # [npx2]=[np,3]x[3,2]
    # [npx2]=[np,3]x[3,2]
    x_ksi_pts = np.dot(grad_geo_basis[:, :, 0], nodes_coord)
    # [npx2]=[np,3]x[3,2]
    x_eta_pts = np.dot(grad_geo_basis[:, :, 1], nodes_coord)
    for i in xrange(n_pts):
        J[i] = np.array([x_ksi_pts[i][:2], x_eta_pts[i][:2]])
        detJ[i] = x_ksi_pts[i][0] * x_eta_pts[i][1] - \
            x_ksi_pts[i][1] * x_eta_pts[i][0]
        invJ[i] = np.array([[x_eta_pts[i][1], -x_ksi_pts[i][1]],
                            [-x_eta_pts[i][0], x_ksi_pts[i][0]]]) / detJ[i]
        invJT[i] = np.array([[x_eta_pts[i][1], -x_eta_pts[i][0]],
                             [-x_ksi_pts[i][1], x_ksi_pts[i][0]]]) / detJ[i]

    return quad_pts, J, invJ, invJT, detJ


def maximum_wave_speed(uL, uR, w=0):
    maximum_speed = 0.0
    if uL >= uR:
        maximum_speed = fabs(0.5 * (uL + uR) - w)
    else:
        maximum_speed = max([fabs(uL - w), fabs(uR - w)])

    return maximum_speed
