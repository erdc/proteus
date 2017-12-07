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
node_coord = None
fes = None
lm = None #: levelmodel
q = None

quad_for_cij = {'weight':[1.0/6, 4.0/6, 1.0/6], 'point':[0.0,0.5,1.0]}#{'weight':[1.0], 'point':[1.0]}#

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

    # declare variables
    global run_times, dij, Cx_T, Cy_T, ML_new, node_coord

    Cx[:] = 0.0
    Cy[:] = 0.0

    if run_times == 0:
        Cx_T = np.zeros_like(Cx, 'd')
        Cy_T = np.zeros_like(Cy, 'd')
        dij = np.zeros_like(Cx, 'd')
        ML_new = np.zeros_like(ML, 'd')
        node_coord = np.zeros_like(mesh_dof, 'd')

    else:
        Cx_T[:] = 0.0
        Cy_T[:] = 0.0
        dij[:] = 0.0
        ML_new[:] = 0.0

    
    qpt = lm.elementQuadraturePoints
    n_pts = qpt.shape[0]
    n_eles = lm.mesh.nElements_global
    nd = lm.nSpace_global
    
    J    = np.zeros((n_eles,n_pts,nd,nd),'d')
    invJ = np.zeros((n_eles,n_pts,nd,nd),'d')
    detJ = np.zeros((n_eles,n_pts),'d')
    
    n_dofs = lm.u[0].femSpace.max_nDOF_element
    v_basis      = np.zeros((n_eles,n_pts,n_dofs),'d')
    v_grad_basis = np.zeros((n_eles,n_pts,n_dofs,nd),'d')

    geo_grad_phi = np.zeros((n_pts, n_dofs, nd),'d')
    for k in range(n_pts):
        for i in range(n_dofs):
            geo_grad_phi[k,i] = lm.u[0].femSpace.elementMaps.localFunctionSpace.basisGradients[i](qpt[k])
    
    
    cfl[:] = 1.0
    
    for it, weight_cii_at_ti in enumerate(quad_for_cij['weight']):
        # moving mesh
        node_coord[:] = mesh_dof
        node_coord[:] += quad_for_cij['point'][it]* dt * mesh_velocity_dof #: assume the velocity is the same; note dt; 
        
        #: use node_coord
        J[:]   =0.0
        detJ[:]=0.0
        invJ[:]=0.0
        cfemIntegrals.parametricMaps_getJacobianValues(geo_grad_phi, mesh_l2g, node_coord, J, detJ, invJ)#:different order of invJ2, det2
#         if  quad_for_cij['point'][it] == 0.0:
#             error_J = np.max(np.abs(q['J']-J)) 
#             assert error_J<1e-9, "J: something wrong"
#             
#             error_invJ = np.max(np.abs(q['inverse(J)']-invJ)) 
#             assert error_invJ<1e-9, "invJ: something wrong"
#             error_detJ = np.max(np.abs(q['det(J)']-detJ)) 
#             assert error_detJ<1e-9, "detJ: something wrong"
    
        
        
        lm.u[0].femSpace.getBasisValues(qpt,v_basis) #: In fact, it is independent of mesh motion
        lm.u[0].femSpace.getBasisGradientValues(qpt,invJ,v_grad_basis)
    
        w_basis = v_basis
        w_grad_basis = v_grad_basis#:depends on new mesh
    
        for e in xrange(nElements_global):
    
            dV = np.abs(detJ[e]) * dV_ref#: broadcast
    
            # Computer cell cfl number
            ele_max_speed = 1e-10
            for i in range(n_dofs):
                dof_i = u_l2g[e, i]
                vi = velocity[dof_i][:-1] #:2D
                wi = mesh_velocity_dof[dof_i][:-1]
                vi = vi - wi
                ele_max_speed = max(
                    [np.sqrt(vi[0] * vi[0] + vi[1] * vi[1]), np.sqrt(wi[0] * wi[0] + wi[1] * wi[1]), ele_max_speed])
            #end-loop-i    

            #: this is corount number; bigger is safer
            cfl[e, :] = np.maximum(cfl[e,:], ele_max_speed / get_diameter(node_coord[mesh_l2g[e]])[1])#: elementwise maximum
    
            # compute C matrix at time it
            c_x   = np.zeros((n_dofs, n_dofs), 'd')
            c_y   = np.zeros((n_dofs, n_dofs), 'd')
            ct_x   = np.zeros((n_dofs, n_dofs), 'd')
            ct_y   = np.zeros((n_dofs, n_dofs), 'd')
            
            for i in range(n_dofs):
                for j in range(n_dofs):
                    for k in range(n_pts):
                        c_x[i, j] += v_basis[e, k, i] * v_grad_basis[e, k, j, 0] * dV[k]
                        c_y[i, j] += v_basis[e, k, i] * v_grad_basis[e, k, j, 1] * dV[k]
                        ct_x[i, j] += v_basis[e, k, j] * v_grad_basis[e, k, i, 0] * dV[k]
                        ct_y[i, j] += v_basis[e, k, j] * v_grad_basis[e, k, i, 1] * dV[k]
                    # end-of-loop-over-k
    
                    Cx[csrRowIndeces_CellLoops[e, i]
                       +
                       csrColumnOffsets_CellLoops[e, i, j]] += c_x[i, j]*quad_for_cij['weight'][it]
                    Cy[csrRowIndeces_CellLoops[e, i]
                       +
                       csrColumnOffsets_CellLoops[e, i, j]] += c_y[i, j]*quad_for_cij['weight'][it]
    
                    Cx_T[csrRowIndeces_CellLoops[e, i]
                         +
                         csrColumnOffsets_CellLoops[e, i, j]] += ct_x[i, j]*quad_for_cij['weight'][it]
                    Cy_T[csrRowIndeces_CellLoops[e, i]
                         +
                         csrColumnOffsets_CellLoops[e, i, j]] += ct_y[i, j]*quad_for_cij['weight'][it]
    
                # end-loop-j
            # end-loop-i
        # end-loop-e
    # end-loop-it
    #dij = np.zeros_like(Cx, 'd')

    #: get new lumped mass matrix
    node_coord[:] = mesh_dof
    node_coord[:] += dt * mesh_velocity_dof #: assume the velocity is the same; note dt; 
    J[:]   =0.0
    detJ[:]=0.0
    invJ[:]=0.0
    cfemIntegrals.parametricMaps_getJacobianValues(geo_grad_phi, mesh_l2g, node_coord, J, detJ, invJ)#:different order of invJ2, det2
    lm.u[0].femSpace.getBasisValues(qpt,v_basis)
    get_ML(mesh_l2g, v_basis, dV_ref, detJ, ML_new)#: v_basis is independent of mesh

    assert np.all(ML_new > 0.), "some element of ML is negative"

    # update solution
    run_times += 1

    globalResidual[:] = 0.0

    edge_based_cfl[:] = 0.0

    # compute dij and compute solution and save the solution into residual
    for dof_i in xrange(numDOFs):  # Serious error: numDOFs=nnz-1

        dii = 0.0
        spatial_residual_dof_i = 0.0
        dof_ii_index = 0
        for j in xrange(csrRowIndeces_DofLoops[dof_i], csrRowIndeces_DofLoops[dof_i + 1]):
            dof_j = csrColumnOffsets_DofLoops[j]
            if dof_i != dof_j:

                cij = np.array([Cx[j], Cy[j]])
                cji = np.array([Cx_T[j], Cy_T[j]])

                cij_norm = np.linalg.norm(cij)
                cji_norm = np.linalg.norm(cji)
                nij = cij / (cij_norm+1e-10)
                nji = cji / (cji_norm+1e-10)

                vj = velocity[dof_j, :-1]  # =[uj,uj]
                vi = velocity[dof_i, :-1]  # =[ui,ui]

                ui = u_dof_old[dof_i]
                uj = u_dof_old[dof_j]

                wj = mesh_velocity_dof[dof_j][:-1]
                wi = mesh_velocity_dof[dof_i][:-1]

                fi = np.array([0.5 * ui * ui, 0.5 * ui * ui]) - ui * wi
                fj = np.array([0.5 * uj * uj, 0.5 * uj * uj]) - uj * wj

                # maximum wave speed

                dij[j] = max([cij_norm * maximum_wave_speed(ui, uj, wj, nij),  # for Burgers equation
                              cji_norm * maximum_wave_speed(uj, ui, wi, nji)])  # for symmetry

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

def maximum_wave_speed(uL, uR, w, n):

    maximum_speed = 0.0
    if uL >= uR:
        maximum_speed = fabs(0.5 * (uL + uR) * (n[0] + n[1]) - np.dot(w, n))
    else:
        maximum_speed = max(
            [fabs(uL * (n[0] + n[1]) - np.dot(w, n)), fabs(uR * (n[0] + n[1]) - np.dot(w, n))])

    return maximum_speed


def get_ML(elementNodesArray, v_basis, dV_ref, detJ, _ML):
    _ML[:] = 0.0
    
    n_ele = elementNodesArray.shape[0]
    n_dof = elementNodesArray.shape[1]
    n_pt  = dV_ref.shape[0]
    ele_ML = np.zeros((n_dof,),'d')
    
    for e in xrange(n_ele):
        ele_ML[:] = 0.0
        for i in xrange(n_dof):
            for k in xrange(n_pt):
                ele_ML[i] += v_basis[e,k,i]*dV_ref[k]*np.abs(detJ[e,k])
            #end-loop-i
            _ML[elementNodesArray[e,i]] += ele_ML[i] #: Serious error: forget += for assembling
    #end-loop-e

def get_diameter(element_nodes):
    """
        can be used for any convex polygon
    """
    
    _hmax = 0.0
    _hmin = 0.0
    N = element_nodes.shape[0]
    diff = element_nodes.reshape(N,1,3)-element_nodes
    distance=(diff**2).sum(2)

    _hmax = distance.max()

    i = np.arange(N)
    distance[i,i]=np.inf
    
    _hmin = distance.min()
    return _hmax, _hmin
    