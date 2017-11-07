from proteus import *
from proteus.default_p import *
from math import *
from vortex2D import *
from proteus.mprans import NCLS
import numpy as np


def get_error(nodeArray,
              elementNodesArray,
              func_distance_to_0_ls,
              phi_distance,
              u,
              true_u):

    nElements_global = elementNodesArray.shape[0]

    n_pts = 3
    geo_basis = np.zeros((n_pts, 3), 'd')
    geo_grad_basis = np.zeros((n_pts, 3, 2), 'd')
    geo_basis[0, :] = np.array([1, 0, 0], 'd')
    geo_basis[1, :] = np.array([0, 1, 0], 'd')
    geo_basis[2, :] = np.array([0, 0, 1], 'd')

    geo_grad_basis[:, :, 0] = np.array([-1, 0, 1], 'd')
    geo_grad_basis[:, :, 1] = np.array([-1, 1, 0], 'd')

    dV_ref = np.zeros((n_pts, 1), 'd')
    dV_ref[:] = 1.0 / 3

    u_trial_ref = np.zeros((n_pts, 3), 'd')
    u_trial_ref[0, :] = np.array([1, 0, 0], 'd')
    u_trial_ref[1, :] = np.array([0, 1, 0], 'd')
    u_trial_ref[2, :] = np.array([0, 0, 1], 'd')

    u_grad_trial_ref = np.zeros((n_pts, 3, 2), 'd')
    u_grad_trial_ref[:, 0, :] = np.array([-1, -1], 'd')
    u_grad_trial_ref[:, 1, :] = np.array([0, 1], 'd')
    u_grad_trial_ref[:, 2, :] = np.array([1, 0], 'd')

    u_error_L_infty = 0.0
    grad_u_error_L_infty = 0.0
    ls_grad_u_error_L_infty = 0.0

    for e in xrange(nElements_global):
        qpt, J, invJ, invJT, detJ = P1_calculateMapping_element(
            nodeArray[elementNodesArray[e]], geo_basis, geo_grad_basis)

        dV = np.abs(detJ) * dV_ref

        v_basis = u_trial_ref
        v_grad_basis = [np.dot(u_grad_trial_ref[k, :, :], invJ[k, :, :])
                        for k in xrange(n_pts)]

        distance_at_node = np.zeros((n_pts,), 'd')
#         for k in range(n_pts):
#             distance_at_node[k] = func_distance_to_0_ls(qpt[k])
        distance_at_node[:] = phi_distance[elementNodesArray[e]]

        u_pt = np.zeros((n_pts,), 'd')
        true_u_pt = np.zeros((n_pts,), 'd')
        grad_u_pt = np.zeros((n_pts, 2), 'd')

#         print e, np.min(distance_at_node), np.max(distance_at_node)
#         print qpt

        for k in range(n_pts):
            for i in range(3):
                u_pt[k] += v_basis[k][i] * u[elementNodesArray[e, i]]
                grad_u_pt[k, :] += v_grad_basis[k][i, :] * \
                    u[elementNodesArray[e, i]]
                true_u_pt[k] += v_basis[k][i] * true_u[elementNodesArray[e, i]]
            # end-i
            u_error_L_infty = max(
                [u_error_L_infty, np.abs(u_pt[k] - true_u_pt[k])])

        # check 1 point is enough for P1
        grad_u_error_L_infty = max(
            [grad_u_error_L_infty, np.abs(np.linalg.norm(grad_u_pt[0]) - 1.0)])

        if np.min(distance_at_node) * np.max(distance_at_node) < 0.0:
            ls_grad_u_error_L_infty = max(
                [ls_grad_u_error_L_infty, np.abs(np.linalg.norm(grad_u_pt[0]) - 1.0)])
#             print e, np.abs(np.linalg.norm(grad_u_pt[0]) - 1.0)
        # end-k
    # end-e
    return u_error_L_infty, grad_u_error_L_infty, ls_grad_u_error_L_infty


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
        J[i] = J[i].T
        detJ[i] = x_ksi_pts[i][0] * x_eta_pts[i][1] - \
            x_ksi_pts[i][1] * x_eta_pts[i][0]
        invJT[i] = np.array([[x_eta_pts[i][1], -x_ksi_pts[i][1]],
                             [-x_eta_pts[i][0], x_ksi_pts[i][0]]]) / detJ[i]
        invJ[i] = np.array([[x_eta_pts[i][1], -x_eta_pts[i][0]],
                            [-x_ksi_pts[i][1], x_ksi_pts[i][0]]]) / detJ[i]

    return quad_pts, J, invJ, invJT, detJ
