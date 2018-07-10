cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool
from proteus.Profiling import logEvent
from proteus.mprans import MeshSmoothing as ms
from proteus.mprans cimport MeshSmoothing as ms
from proteus import Comm
from mpi4py import MPI


cdef class cCoefficients:
    cdef public:
        object pyCoefficients
        double C

    def __cinit__(self):
        self.C = 1.

    def attachPyCoefficients(self,
                             object pyCoefficients):
        self.pyCoefficients = pyCoefficients

    def preStep(self):
        pc = self.pyCoefficients
        self.cppPreStep(q_uOfX=pc.uOfXTatQuadrature,
                        q_J=pc.model.q['abs(det(J))'],
                        q_weights=pc.model.elementQuadratureWeights[('u', 0)],
                        areas_=pc.areas,
                        q_rci=pc.model.q[('r', 0)],
                        q_fci=pc.model.q[('f', 0)],
                        t=pc.t,
                        nElements_owned=pc.mesh.nElements_owned)

    cdef cppPreStep(self,
                    double[:,:] q_uOfX,
                    double[:,:] q_J,
                    double[:] q_weights,
                    double[:] areas_,
                    double[:, :] q_rci,
                    double[:, :, :] q_fci,
                    double t,
                    int nElements_owned):
        cdef double integral_1_over_f = 0.
        cdef int N_eN = q_J.shape[0]
        cdef int nE = 0
        cdef double area = 0
        cdef int nJ = len(q_J)
        cdef int nk = len(q_weights)
        for eN in xrange(nJ):
            area = 0
            for k in xrange(nk):
                area += q_J[eN, k]*q_weights[k]
                if eN < nElements_owned:
                    integral_1_over_f += q_J[eN, k]*q_weights[k]/q_uOfX[eN, k]
            areas_[eN] = area
            if eN < nElements_owned:
                nE += 1
        comm = Comm.get().comm.tompi4py()
        if comm.size > 1:
            integral_1_over_f = comm.allreduce(integral_1_over_f, op=MPI.SUM)
            nE = comm.allreduce(nE, op=MPI.SUM)
        cdef double C = integral_1_over_f/nE  # update scaling coefficient
        self.C = C
        cdef int nrci = len(q_rci)
        for eN in range(nrci):
            for k in range(nk):
                q_fci[eN, k, :] = 0.0
                q_rci[eN, k] = -(1./(q_uOfX[eN, k]*C)-1./areas_[eN])

    def postStep(self):
        pc = self.pyCoefficients
        self.cppPostStep()

    cdef cppPostStep(self):
        cdef object pc = self.pyCoefficients
        cdef object mesh = pc.mesh
        cdef int nd = pc.nd
        cdef double[:,:] elementBoundaryBarycentersArray = mesh.elementBoundaryBarycentersArray
        cdef double[:,:] elementBarycentersArray = mesh.elementBarycentersArray
        cdef double[:,:] nodeArray = mesh.nodeArray
        cdef int[:,:] elementBoundaryNodesArray = mesh.elementBoundaryNodesArray
        cdef int[:,:] elementNodesArray = mesh.elementNodesArray
        cdef int[:,:] elementBoundariesArray = mesh.elementBoundariesArray
        cdef int nElementBoundaries_global = mesh.nElementBoundaries_global
        cdef int nNodes_elementBoundary = mesh.nNodes_elementBoundary
        cdef int nElements_global = mesh.nElements_global
        cdef int nNodes_element = mesh.nNodes_element
        cdef double[:,:,:] elementBoundaryNormalsArray = pc.mesh.elementBoundaryNormalsArray
        # update element boundary barycenters
        for ebN in range(nElementBoundaries_global):
            elementBoundaryBarycentersArray[ebN, 0] = 0.
            elementBoundaryBarycentersArray[ebN, 1] = 0.
            elementBoundaryBarycentersArray[ebN, 2] = 0.
            for nN in range(nNodes_elementBoundary):
                elementBoundaryBarycentersArray[ebN, 0] += nodeArray[elementBoundaryNodesArray[ebN, nN], 0]
                elementBoundaryBarycentersArray[ebN, 1] += nodeArray[elementBoundaryNodesArray[ebN, nN], 1]
                elementBoundaryBarycentersArray[ebN, 2] += nodeArray[elementBoundaryNodesArray[ebN, nN], 2]
            elementBoundaryBarycentersArray[ebN, 0] /= nNodes_elementBoundary
            elementBoundaryBarycentersArray[ebN, 1] /= nNodes_elementBoundary
            elementBoundaryBarycentersArray[ebN, 2] /= nNodes_elementBoundary
        # update element barycenters
        for eN in range(nElements_global):
            elementBarycentersArray[eN, 0] = 0.
            elementBarycentersArray[eN, 1] = 0.
            elementBarycentersArray[eN, 2] = 0.
            for ebN in range(nNodes_element):
                elementBarycentersArray[eN, 0] += nodeArray[elementNodesArray[eN, ebN], 0]
                elementBarycentersArray[eN, 1] += nodeArray[elementNodesArray[eN, ebN], 1]
                elementBarycentersArray[eN, 2] += nodeArray[elementNodesArray[eN, ebN], 2]
            elementBarycentersArray[eN, 0] /= nNodes_element
            elementBarycentersArray[eN, 1] /= nNodes_element
            elementBarycentersArray[eN, 2] /= nNodes_element
        # update normals
        if nd == 2:
            # triangle
            ms.updateElementBoundaryNormalsTriangle(elementBoundaryNormalsArray,
                                                    nodeArray,
                                                    elementBoundariesArray,
                                                    elementBoundaryNodesArray,
                                                    elementBoundaryBarycentersArray,
                                                    elementBarycentersArray)
        elif nd == 3:
            # tetra
            ms.updateElementBoundaryNormalsTetra(elementBoundaryNormalsArray,
                                                 nodeArray,
                                                 elementBoundariesArray,
                                                 elementBoundaryNodesArray,
                                                 elementBoundaryBarycentersArray,
                                                 elementBarycentersArray)

    def pseudoTimeStepping(self,
                           xx,
                           eps=1.):
        pc = self.pyCoefficients
        return self.cppPseudoTimeSteppingParallel(xx=xx,
                                                  eps=eps)

    cdef cppPseudoTimeSteppingParallel(self,
                                       double[:,:] xx,
                                       double eps):
        logEvent("Pseudo-time stepping with dt={eps} (0<t<1)".format(eps=eps),
                 level=3)
        pc = self.pyCoefficients
        cdef double[:] t_range = np.linspace(0., 1., int(1./eps+1))[1:]
        cdef int[:] eN_phi = np.zeros(len(xx), dtype=np.int32)
        cdef double[:] normal = np.zeros(3)
        cdef double t_last = 0
        cdef double dt = 0
        cdef int flag
        cdef bool fixed
        cdef double area = 0
        cdef int eN
        cdef double ls_phi
        cdef double f
        cdef double Ccoeff = pc.C
        cdef int nd = pc.nd
        cdef double[:] dphi = np.zeros(nd)
        cdef double[:,:] grads = pc.grads
        cdef double[:] v_grad = np.zeros(nd)
        cdef double[:] areas_nodes = pc.areas_nodes,
        cdef double[:] areas = pc.areas
        cdef object u_phi = pc.u_phi

        cdef object femSpace = pc.model.u[0].femSpace
        cdef bool find_nearest_node
        # initialise mesh memoryview before loop
        cdef double[:,:] nodeArray = pc.mesh.nodeArray
        cdef int[:] nodeStarOffsets = pc.mesh.nodeStarOffsets
        cdef int[:] nodeStarArray = pc.mesh.nodeStarArray
        cdef int[:] nodeElementOffsets = pc.mesh.nodeElementOffsets
        cdef int[:] nodeElementsArray = pc.mesh.nodeElementsArray
        cdef int[:] nodeMaterialTypes = pc.mesh.nodeMaterialTypes
        cdef double[:,:] elementBarycentersArray = pc.mesh.elementBarycentersArray
        cdef int[:,:] elementBoundaryElementsArray = pc.mesh.elementBoundaryElementsArray
        cdef int[:,:] elementBoundariesArray = pc.mesh.elementBoundariesArray
        cdef int[:,:] elementNodesArray = pc.mesh.elementNodesArray
        cdef double[:,:] elementBoundaryBarycentersArray = pc.mesh.elementBoundaryBarycentersArray
        cdef int[:] exteriorElementBoundariesBoolArray = np.zeros(pc.mesh.nElementBoundaries_global, dtype=np.int32)
        for bb_i in pc.mesh.exteriorElementBoundariesArray:
            exteriorElementBoundariesBoolArray[bb_i] = 1
        cdef double[:,:,:] elementBoundaryNormalsArray = pc.mesh.elementBoundaryNormalsArray
        cdef int nNodes_owned = pc.mesh.nNodes_owned
        cdef int nNodes_global = pc.mesh.nNodes_global
        cdef int nElements_owned = pc.mesh.nElements_owned
        cdef int nElements_global = pc.mesh.nElements_global
        cdef int[:] nodeNumbering_subdomain2global = pc.mesh.globalMesh.nodeNumbering_subdomain2global
        cdef int[:] elementNumbering_subdomain2global = pc.mesh.globalMesh.elementNumbering_subdomain2global
        cdef int[:] nodeOffsets_subdomain_owned = pc.mesh.globalMesh.nodeOffsets_subdomain_owned
        cdef int[:] elementOffsets_subdomain_owned = pc.mesh.globalMesh.elementOffsets_subdomain_owned
        cdef int[:] nearest_nodes = np.array([i for i in range(len(xx))], dtype=np.int32)
        cdef int[:,:] elementNeighborsArray = pc.mesh.elementNeighborsArray
        cdef int nElementBoundaries_owned = pc.mesh.nElementBoundaries_owned
        cdef int[:] elementBoundaryNumbering_subdomain2global = pc.mesh.globalMesh.elementBoundaryNumbering_subdomain2global
        cdef int[:] elementBoundaryOffsets_subdomain_owned = pc.mesh.globalMesh.elementBoundaryOffsets_subdomain_owned
        cdef int[:,:] elementBoundaryNodesArray = pc.mesh.elementBoundaryNodesArray
        cdef int[:] fixedNodesBoolArray = pc.mesh.fixedNodesBoolArray
        cdef int nSmoothIn = pc.nSmoothIn
        cdef int nSmoothOut = pc.nSmoothOut
        cdef bool inside_eN = False
        cdef double[:] vec = np.zeros(3)
        cdef double[:] vec2 = np.zeros(3)
        cdef double vec_dist
        cdef double[:] fixed_dir = np.ones(3)
        cdef int nearestN
        cdef int typeN
        comm = Comm.get().comm.tompi4py()
        comm.barrier()
        cdef int my_rank = comm.rank
        cdef int comm_size = comm.size
        # coords to send
        cdef dict coords_2rank = {}
        # nodes to send
        cdef dict nodes0_2rank = {}
        # original rank of sending rank
        cdef dict rank0_2rank = {}
        # original node number on sending rank
        cdef dict nearestN_2rank = {}
        # elements to send
        cdef dict typeN_2rank = {}
        # direction to rank (for sliding nodes)
        cdef dict dir_2rank = {}
        cdef dict b_i_2rank = {}
        # solutions
        cdef bool solFound
        cdef bool sendBack = False
        cdef int nPending_disp = 0  # number of pending solutions from other ranks
        cdef int nPending_disp_total = 0  # number of pending solutions from other ranks
        # declare variables for speed up
        cdef int node
        cdef int nOffset
        cdef double dot
        cdef double t
        cdef int j
        cdef int i_time
        cdef int parallel_steps
        cdef double[:] coords = np.zeros(3)
        cdef int b_i, b_i_global
        cdef bool found
        if comm_size > 1:
            for rank in range(comm_size):
                # things to send to other processors to get solution
                coords_2rank[rank] = np.zeros((0, 3))
                dir_2rank[rank] = np.zeros((0, 3))
                b_i_2rank[rank] = np.zeros(0, dtype=np.int32)
                nodes0_2rank[rank] = np.zeros(0, dtype=np.int32)
                rank0_2rank[rank] = np.zeros(0, dtype=np.int32)
                nearestN_2rank[rank] = np.zeros(0, dtype=np.int32)
                typeN_2rank[rank] = np.zeros(0, dtype=np.int32)
        cdef np.ndarray coords_2doArray
        # cdef int[:] rank0_2doArray = np.array([], dtype=np.int32)
        # cdef int[:] nearestN_2doArray = np.array([], dtype=np.int32)
        # cdef int[:] typeN_2doArray = np.array([], dtype=np.int32)
        cdef double[:] nodesSentBoolArray = np.zeros(len(xx))
        cdef bool pending = False
        cdef double[:] starting_coords = np.zeros(3)
        for i_time, t in enumerate(t_range):
            logEvent("Pseudo-time stepping t={t}".format(t=t), level=3)
            nPending_disp = 0
            dt = t-t_last
            t_last = t
            for i in range(nNodes_owned):
                new_rank = my_rank
                node = i
                # coords = xx[node]
                coords = np.zeros(3)
                coords[0] = xx[node, 0]
                coords[1] = xx[node, 1]
                coords[2] = xx[node, 2]
                fixed = False
                pending = False
                area = 0.
                flag = nodeMaterialTypes[node]
                fixed_dir[0] = 1.
                fixed_dir[1] = 1.
                fixed_dir[2] = 1.
                b_i_global = -1
                for ndi in range(nd):
                    v_grad[ndi] = 0.
                    dphi[ndi] = 0.
                if nodesSentBoolArray[node] == 1:
                    pending = True
                    nPending_disp += 1
                elif nodesSentBoolArray[node] == 0:
                    if flag != 0:
                        fixed = True
                        fixed_dir[:] = 0
                        if fixedNodesBoolArray is not None:
                            if fixedNodesBoolArray[node] == 1:
                                fixed = True
                        # smooth on boundary only unless it is a corner node
                        if fixed is False:
                            # tridelat todo: works only in 2D here
                            if nd == 2:
                                vec[:] = 0.
                                vec2[:] = 0.
                                for nOffset in range(nodeStarOffsets[node],
                                                    nodeStarOffsets[node+1]):
                                    if nodeMaterialTypes[nodeStarArray[nOffset]] != 0:
                                        if vec[0] == 0. and vec[1] == 0. and vec[2] == 0.:
                                            vec[0] = nodeArray[node, 0]-nodeArray[nodeStarArray[nOffset], 0]
                                            vec[1] = nodeArray[node, 1]-nodeArray[nodeStarArray[nOffset], 1]
                                            vec[2] = nodeArray[node, 2]-nodeArray[nodeStarArray[nOffset], 2]
                                            vec_dist = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
                                            vec[0] = vec[0]/vec_dist
                                            vec[1] = vec[1]/vec_dist
                                            vec[2] = vec[2]/vec_dist
                                            fixed_dir[0] = abs(vec[0])
                                            fixed_dir[1] = abs(vec[1])
                                            fixed_dir[2] = abs(vec[2])
                                        else:
                                            vec2[0] = nodeArray[node, 0]-nodeArray[nodeStarArray[nOffset], 0]
                                            vec2[1] = nodeArray[node, 1]-nodeArray[nodeStarArray[nOffset], 1]
                                            vec2[2] = nodeArray[node, 2]-nodeArray[nodeStarArray[nOffset], 2]
                                            vec_dist = np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
                                            vec2[0] = vec2[0]/vec_dist
                                            vec2[1] = vec2[1]/vec_dist
                                            vec2[2] = vec2[2]/vec_dist
                                            dot = vec[0]*vec2[0]+vec[1]*vec2[1]+vec[2]*vec2[2]
                                            if dot == 1. or dot == -1.:
                                                dot = 1.
                                            else:
                                                fixed = True
                            elif nd == 3:
                                assert 1>2, 'works only in 2D for moving boundary nodes'
                                fixed = True
                    if not fixed:  # either flag==0 or not fixed
                        if i_time == 0:  # nodes are at original position (no search)
                            for ndi in range(nd):
                                v_grad[ndi] = grads[node, ndi]
                            area = areas_nodes[node]
                            eN = -1
                            eN_phi[node] = eN
                            if u_phi is not None:
                                ls_phi = u_phi[node]
                            else:
                                ls_phi = 1e12
                            f = pc.evaluateFunAtX(x=coords, ls_phi=ls_phi)
                            for ndi in range(nd):
                                dphi[ndi] = v_grad[ndi]/(t*1./(f*Ccoeff)+(1-t)*1./area)
                                coords[ndi] += dphi[ndi]*fixed_dir[ndi]*dt
                                xx[node, ndi] = coords[ndi]
                        else:  # find node that moved already (search)
                            # line below needs optimisation:
                            if i_time > 1:
                                find_nearest_node = True#False
                            else:
                                find_nearest_node = True
                            typeN = 0
                            nearestN = ms.getLocalNearestNode(coords=coords,
                                                              nodeArray=nodeArray,
                                                              nodeStarOffsets=nodeStarOffsets,
                                                              nodeStarArray=nodeStarArray,
                                                              node=nearest_nodes[node])
                            nearest_nodes[i] = nearestN
                            if nearestN >= nNodes_owned:
                                nearestN, new_rank = ms.cyCheckOwnedVariable(variable_nb_local=nearestN,
                                                                             rank=my_rank,
                                                                             nVariables_owned=nNodes_owned,
                                                                             variableNumbering_subdomain2global=nodeNumbering_subdomain2global,
                                                                             variableOffsets_subdomain_owned=nodeOffsets_subdomain_owned)
                            else:
                                typeN = 1
                                nearestN = ms.getLocalNearestElementAroundNode(coords=coords,
                                                                               nodeElementOffsets=nodeElementOffsets,
                                                                               nodeElementsArray=nodeElementsArray,
                                                                               elementBarycentersArray=elementBarycentersArray,
                                                                               node=nearestN)
                            if typeN == 1:
                                starting_coords = np.zeros(3)
                                starting_coords[0] = elementBarycentersArray[nearestN, 0]
                                starting_coords[1] = elementBarycentersArray[nearestN, 1]
                                starting_coords[2] = elementBarycentersArray[nearestN, 2]
                                nearestN, b_i = ms.getLocalNearestElementIntersection(coords=coords,
                                                                                      starting_coords=starting_coords,
                                                                                      elementBoundaryNormalsArray=elementBoundaryNormalsArray,
                                                                                      elementBoundariesArray=elementBoundariesArray,
                                                                                      elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                                                                      elementBoundaryElementsArray=elementBoundaryElementsArray,
                                                                                      exteriorElementBoundariesBoolArray=exteriorElementBoundariesBoolArray,
                                                                                      eN=nearestN)
                                if nearestN >= nElements_owned:
                                    nearestN, new_rank = ms.cyCheckOwnedVariable(variable_nb_local=nearestN,
                                                                                 rank=my_rank,
                                                                                 nVariables_owned=nElements_owned,
                                                                                 variableNumbering_subdomain2global=elementNumbering_subdomain2global,
                                                                                 variableOffsets_subdomain_owned=elementOffsets_subdomain_owned)
                                if nearestN == -1:
                                    b_i_global, new_rank = getGlobalVariable(variable_nb_local=b_i,
                                                                             nVariables_owned=nElementBoundaries_owned,
                                                                             variableNumbering_subdomain2global=elementBoundaryNumbering_subdomain2global,
                                                                             variableOffsets_subdomain_owned=elementBoundaryOffsets_subdomain_owned)
                                    typeN = 1
                                    pending = True
                                    for nodeEl in elementNodesArray[nearestN]:
                                        if nodeMaterialTypes[nodeEl] != 0:
                                            # tridelat hack: node is probably on exterior boundary
                                            new_rank = my_rank
                            if new_rank != my_rank:
                                pending = True
                            else:
                                pending = False
                            if pending is False:  # found an element
                                inside_eN = True
                                for j, bb_i in enumerate(elementBoundariesArray[nearestN]):
                                    normal[0] = elementBoundaryNormalsArray[nearestN, j, 0]
                                    normal[1] = elementBoundaryNormalsArray[nearestN, j, 1]
                                    normal[2] = elementBoundaryNormalsArray[nearestN, j, 2]
                                    bound_bar = elementBoundaryBarycentersArray[bb_i]
                                    dot = (bound_bar[0]-coords[0])*normal[0]+(bound_bar[1]-coords[1])*normal[1]+(bound_bar[2]-coords[2])*normal[2]
                                    if dot < 0:
                                        inside_eN = False
                                if inside_eN:
                                    xi = femSpace.elementMaps.getInverseValue(nearestN, coords)
                                    # line below needs optimisation:
                                    v_grad = pc.getGradientValue(nearestN, xi)
                                    area = pc.getAreaValue(nearestN, xi)
                                    if u_phi is not None:
                                        # line below needs optimisation:
                                        ls_phi = pc.getLevelSetValue(nearestN, coords)
                                    else:
                                        ls_phi = 1e12
                                    f = pc.evaluateFunAtX(x=coords, ls_phi=ls_phi)
                                    for ndi in range(nd):
                                        dphi[ndi] = v_grad[ndi]/(t*1./(f*Ccoeff)+(1-t)*1./area)
                                        coords[ndi] += dphi[ndi]*fixed_dir[ndi]*dt
                                        xx[node, ndi] = coords[ndi]
                                else:
                                    for nodeEl in elementNodesArray[nearestN]:
                                        if nodeEl > nNodes_owned:
                                            nearestN, new_rank = ms.cyCheckOwnedVariable(variable_nb_local=nodeEl,
                                                                                         rank=my_rank,
                                                                                         nVariables_owned=nNodes_owned,
                                                                                         variableNumbering_subdomain2global=nodeNumbering_subdomain2global,
                                                                                         variableOffsets_subdomain_owned=nodeOffsets_subdomain_owned)
                                            typeN = 0
                                            pending = True
                                if inside_eN is False and pending is False:
                                    print('outside!!', nearestN,  coords[0], coords[1], elementBarycentersArray[nearestN,0], elementBarycentersArray[nearestN,1])
                                    print('neighbours', elementNeighborsArray[nearestN, 0], elementNeighborsArray[nearestN, 1], elementNeighborsArray[nearestN, 2])
                            if pending is True:  # info to send to other processor
                                nodesSentBoolArray[node] = 1
                                coords_2rank[new_rank] = np.append(coords_2rank[new_rank], [coords], axis=0)
                                dir_2rank[new_rank] = np.append(dir_2rank[new_rank], [fixed_dir], axis=0)
                                nodes0_2rank[new_rank] = np.append(nodes0_2rank[new_rank], node)
                                nearestN_2rank[new_rank] = np.append(nearestN_2rank[new_rank], nearestN)
                                typeN_2rank[new_rank] = np.append(typeN_2rank[new_rank], typeN)
                                rank0_2rank[new_rank] = np.append(rank0_2rank[new_rank], my_rank)
                                b_i_2rank[new_rank] = np.append(b_i_2rank[new_rank], b_i_global)
                                nPending_disp += 1
            # parallel comm
            if comm_size > 1:
                # number of pending solutions to send to other processors
                nPending_disp_total = comm.allreduce(nPending_disp)
                parallel_steps = -1
                while nPending_disp_total > 0:
                    parallel_steps += 1
                    nPending_disp = 0
                    sendBack = True
                    # initialize solution from previously found solutions
                    coords_2doArray = coords_2rank[my_rank]
                    dir_2doArray = dir_2rank[my_rank]
                    nodes0_2doArray = nodes0_2rank[my_rank]
                    b_i_2doArray = b_i_2rank[my_rank]
                    rank0_2doArray = rank0_2rank[my_rank]
                    nearestN_2doArray = nearestN_2rank[my_rank]
                    typeN_2doArray = typeN_2rank[my_rank]
                    if parallel_steps > 0:
                        # if coming from this rank, solution was already found
                        solFound_2doArray = np.ones(len(coords_2doArray), dtype=np.int32)
                    else:
                        solFound_2doArray = np.zeros(len(coords_2doArray), dtype=np.int32)
                    # wipe dicts out
                    coords_2rank[my_rank] = np.zeros((0, 3))
                    dir_2rank[my_rank] = np.zeros((0, 3))
                    nodes0_2rank[my_rank] = np.zeros(0, dtype=np.int32)
                    rank0_2rank[my_rank] = np.zeros(0, dtype=np.int32)
                    nearestN_2rank[my_rank] = np.zeros(0, dtype=np.int32)
                    typeN_2rank[my_rank] = np.zeros(0, dtype=np.int32)
                    b_i_2rank[my_rank] = np.zeros(0, dtype=np.int32)
                    for rank_recv in range(comm.size):
                        for rank_send in range(comm.size):
                            if rank_send != rank_recv and rank_send == my_rank:
                                comm.send(coords_2rank[rank_recv].size, dest=rank_recv, tag=0)
                                if coords_2rank[rank_recv].size > 0:
                                    # coords
                                    comm.send(coords_2rank[rank_recv],  dest=rank_recv, tag=1)
                                    # current nearest nodes
                                    comm.send(nodes0_2rank[rank_recv], dest=rank_recv, tag=2)
                                    # rank for original nodes (to send back final solution)
                                    comm.send(rank0_2rank[rank_recv], dest=rank_recv, tag=3)
                                    # node number of nodes on original rank send back final solution)
                                    comm.send(nearestN_2rank[rank_recv], dest=rank_recv, tag=4)
                                    # current nearest elements
                                    comm.send(typeN_2rank[rank_recv], dest=rank_recv, tag=5)
                                    # node direction (if boundary nodes)
                                    comm.send(dir_2rank[rank_recv], dest=rank_recv, tag=6)
                                    comm.send(b_i_2rank[rank_recv], dest=rank_recv, tag=7)
                                    # wipe the dictionary items now that it has been sent
                                coords_2rank[rank_recv] = np.zeros((0, 3))
                                dir_2rank[rank_recv] = np.zeros((0, 3))
                                nodes0_2rank[rank_recv] = np.zeros(0, dtype=np.int32)
                                rank0_2rank[rank_recv] = np.zeros(0, dtype=np.int32)
                                nearestN_2rank[rank_recv] = np.zeros(0, dtype=np.int32)
                                typeN_2rank[rank_recv] = np.zeros(0, dtype=np.int32)
                                b_i_2rank[rank_recv] = np.zeros(0, dtype=np.int32)
                            elif rank_send!= rank_recv and rank_recv == my_rank:
                                size = comm.recv(source=rank_send, tag=0)
                                if size > 0:
                                    # coords
                                    coords_2do = comm.recv(source=rank_send, tag=1)
                                    coords_2doArray = np.append(coords_2doArray, coords_2do, axis=0)
                                    # node number on original rank
                                    nodes0_2do = comm.recv(source=rank_send, tag=2)
                                    nodes0_2doArray = np.append(nodes0_2doArray, nodes0_2do)
                                    # original rank
                                    rank0_2do = comm.recv(source=rank_send, tag=3)
                                    rank0_2doArray = np.append(rank0_2doArray, rank0_2do)
                                    # current nearest N
                                    nearestN_2do = comm.recv(source=rank_send, tag=4)
                                    nearestN_2doArray = np.append(nearestN_2doArray, nearestN_2do)
                                    # is it node (0) or element (1)
                                    typeN_2do = comm.recv(source=rank_send, tag=5)
                                    typeN_2doArray = np.append(typeN_2doArray, typeN_2do)
                                    # direction
                                    dir_2do = comm.recv(source=rank_send, tag=6)
                                    dir_2doArray = np.append(dir_2doArray, dir_2do, axis=0)
                                    b_i_2do = comm.recv(source=rank_send, tag=7)
                                    b_i_2doArray = np.append(b_i_2doArray, b_i_2do)
                                    # sent information means solution was not found
                                    solFound_2doArray = np.append(solFound_2doArray, np.zeros(len(coords_2do), dtype=np.int32))
                            comm.barrier()
                    comm.barrier()
                    for iN in range(nearestN_2doArray.size):
                        coords = np.zeros(3)
                        coords[0] = coords_2doArray[iN, 0]
                        coords[1] = coords_2doArray[iN, 1]
                        coords[2] = coords_2doArray[iN, 2]
                        nearestN = nearestN_2doArray[iN]
                        typeN = typeN_2doArray[iN]
                        new_rank = my_rank
                        fixed_dir[0] = dir_2doArray[iN, 0]
                        fixed_dir[1] = dir_2doArray[iN, 1]
                        fixed_dir[2] = dir_2doArray[iN, 2]
                        b_i_global = b_i_2doArray[iN]
                        if solFound_2doArray[iN] == 0:
                            if typeN == 0:
                                nearestN = ms.getLocalNearestNode(coords=coords,
                                                                  nodeArray=nodeArray,
                                                                  nodeStarOffsets=nodeStarOffsets,
                                                                  nodeStarArray=nodeStarArray,
                                                                  node=nearestN)
                                if nearestN >= nNodes_owned:
                                    nearestN, new_rank = ms.cyCheckOwnedVariable(variable_nb_local=nearestN,
                                                                               rank=my_rank,
                                                                               nVariables_owned=nNodes_owned,
                                                                               variableNumbering_subdomain2global=nodeNumbering_subdomain2global,
                                                                               variableOffsets_subdomain_owned=nodeOffsets_subdomain_owned)
                                else:
                                    typeN = 1
                                    nearestN = ms.getLocalNearestElementAroundNode(coords=coords,
                                                                                   nodeElementOffsets=nodeElementOffsets,
                                                                                   nodeElementsArray=nodeElementsArray,
                                                                                   elementBarycentersArray=elementBarycentersArray,
                                                                                   node=nearestN)
                            if typeN == 1:
                                if nearestN == -1:
                                    if b_i == -1:
                                        import pdb; pdb.set_trace()
                                    b_i = getLocalVariable(variable_nb_global=b_i_global,
                                                           rank=my_rank,
                                                           nVariables_owned=nElementBoundaries_owned,
                                                           variableNumbering_subdomain2global=elementBoundaryNumbering_subdomain2global,
                                                           variableOffsets_subdomain_owned=elementBoundaryOffsets_subdomain_owned)
                                    for ii, eeN in enumerate(elementBoundaryElementsArray[b_i]):
                                        if eeN == -1:
                                            nearestN = elementBoundaryElementsArray[b_i, ii-1]
                                    if nearestN == -1:
                                        nearestN = eeN
                                    assert nearestN != -1, 'did not find nearestN! {x}, {y}'.format(x=coords[0], y=coords[1])
                                starting_coords = np.zeros(3)
                                starting_coords[0] = elementBarycentersArray[nearestN, 0]
                                starting_coords[1] = elementBarycentersArray[nearestN, 1]
                                starting_coords[2] = elementBarycentersArray[nearestN, 2]
                                nearestN, b_i = ms.getLocalNearestElementIntersection(coords=coords,
                                                                                      starting_coords=starting_coords,
                                                                                      elementBoundaryNormalsArray=elementBoundaryNormalsArray,
                                                                                      elementBoundariesArray=elementBoundariesArray,
                                                                                      elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                                                                      elementBoundaryElementsArray=elementBoundaryElementsArray,
                                                                                      exteriorElementBoundariesBoolArray=exteriorElementBoundariesBoolArray,
                                                                                      eN=nearestN)
                                if nearestN >= nElements_owned:
                                    nearestN, new_rank = ms.cyCheckOwnedVariable(variable_nb_local=nearestN,
                                                                                 rank=my_rank,
                                                                                 nVariables_owned=nElements_owned,
                                                                                 variableNumbering_subdomain2global=elementNumbering_subdomain2global,
                                                                                 variableOffsets_subdomain_owned=elementOffsets_subdomain_owned)
                                if nearestN == -1:
                                    b_i_global, new_rank = getGlobalVariable(variable_nb_local=b_i,
                                                                             nVariables_owned=nElementBoundaries_owned,
                                                                             variableNumbering_subdomain2global=elementBoundaryNumbering_subdomain2global,
                                                                             variableOffsets_subdomain_owned=elementBoundaryOffsets_subdomain_owned)
                                else:
                                    # directly using nearestN for next time
                                    # (no need for b_i_global)
                                    b_i_global = -1
                            if new_rank != my_rank:
                                pending = True
                            else:
                                pending = False
                            if pending is True:
                                nPending_disp += 1
                            else:
                                solFound_2doArray[iN] += 1
                                inside_eN = True  # checking if actually true
                                for ii, bb_i in enumerate(elementBoundariesArray[nearestN]):
                                    normal[0] = elementBoundaryNormalsArray[nearestN, ii, 0]
                                    normal[1] = elementBoundaryNormalsArray[nearestN, ii, 1]
                                    normal[2] = elementBoundaryNormalsArray[nearestN, ii, 2]
                                    bound_bar = elementBoundaryBarycentersArray[bb_i]
                                    dot = (bound_bar[0]-coords[0])*normal[0]+(bound_bar[1]-coords[1])*normal[1]+(bound_bar[2]-coords[2])*normal[2]
                                    if dot < 0:
                                        inside_eN = False
                                if inside_eN:
                                    xi = femSpace.elementMaps.getInverseValue(nearestN, coords)
                                    # line below needs optimisation:
                                    v_grad = pc.getGradientValue(nearestN, xi)
                                    area = pc.getAreaValue(nearestN, xi)
                                    if u_phi is not None:
                                        # line below needs optimisation:
                                        ls_phi = pc.getLevelSetValue(nearestN, coords)
                                    else:
                                        ls_phi = 1e12
                                        f = pc.evaluateFunAtX(x=coords, ls_phi=ls_phi)
                                    for ndi in range(nd):
                                        dphi[ndi] = v_grad[ndi]/(t*1./(f*Ccoeff)+(1-t)*1./area)
                                        coords[ndi] += dphi[ndi]*fixed_dir[ndi]*dt
                                else:
                                    print('did not find coords', coords[0], coords[1])
                        coords_2rank[new_rank] = np.append(coords_2rank[new_rank], [coords], axis=0)
                        nodes0_2rank[new_rank] = np.append(nodes0_2rank[new_rank], nodes0_2doArray[iN])
                        rank0_2rank[new_rank] = np.append(rank0_2rank[new_rank], rank0_2doArray[iN])
                        nearestN_2rank[new_rank] = np.append(nearestN_2rank[new_rank], nearestN)
                        typeN_2rank[new_rank] = np.append(typeN_2rank[new_rank], typeN)
                        dir_2rank[new_rank] = np.append(dir_2rank[new_rank], [fixed_dir], axis=0)
                        b_i_2rank[new_rank] = np.append(b_i_2rank[new_rank], b_i_global)
                    nPending_disp_total = comm.allreduce(nPending_disp)


        # SEND NODES POSITION SOLUTION BACK TO ORIGINAL PROCESSORS
        if sendBack is True:
            coords_2doArray = coords_2rank[my_rank]
            nodes0_2doArray = nodes0_2rank[my_rank]
            rank0_2doArray = rank0_2rank[my_rank]
            for rank in range(comm.size):
                # things to send to other processors to get solution
                coords_2rank[rank] = np.zeros((0, 3))
                nodes0_2rank[rank] = np.zeros(0, dtype=np.int32)
                rank0_2rank[rank] = np.zeros(0, dtype=np.int32)
            for iN in range(len(coords_2doArray)):
                coords_2rank[rank0_2doArray[iN]] = np.append(coords_2rank[rank0_2doArray[iN]], [coords_2doArray[iN]], axis=0)
                nodes0_2rank[rank0_2doArray[iN]] = np.append(nodes0_2rank[rank0_2doArray[iN]], nodes0_2doArray[iN])
                rank0_2rank[rank0_2doArray[iN]] = np.append(rank0_2rank[rank0_2doArray[iN]], rank0_2doArray[iN])
            coords_2doArray = coords_2rank[my_rank]
            nodes0_2doArray = nodes0_2rank[my_rank]
            for rank_recv in range(comm.size):
                for rank_send in range(comm.size):
                    if rank_send != rank_recv and rank_send == my_rank:
                        comm.send(coords_2rank[rank_recv].size, dest=rank_recv, tag=0)
                        if coords_2rank[rank_recv].size > 0:
                            # coords
                            comm.send(coords_2rank[rank_recv],  dest=rank_recv, tag=1)
                            # original nodes
                            comm.send(nodes0_2rank[rank_recv], dest=rank_recv, tag=2)
                    elif rank_send!= rank_recv and rank_recv == my_rank:
                        size = comm.recv(source=rank_send, tag=0)
                        if size > 0:
                            # coords
                            coords_2do = comm.recv(source=rank_send, tag=1)
                            coords_2doArray = np.append(coords_2doArray, coords_2do, axis=0)
                            # original nodes
                            nodes0_2do = comm.recv(source=rank_send, tag=2)
                            nodes0_2doArray = np.append(nodes0_2doArray, nodes0_2do)
            for iN in range(len(nodes0_2doArray)):
                node0 = nodes0_2doArray[iN]
                coords = coords_2doArray[iN]
                for ind in range(nd):
                    xx[node0, ind] = coords[ind]

        # FINAL STEP: GET NON-OWNED NODES POSITION FOR CONSISTENCY
        # BUILD NON OWNED NODES ARRAY TO RETRIEVE SOLUTION
        ms.getNonOwnedNodeValues(xx,
                                 nNodes_owned,
                                 nNodes_global,
                                 nodeNumbering_subdomain2global,
                                 nodeOffsets_subdomain_owned)

        if nSmoothOut > 0:
            comm.barrier()
            logEvent('Smoothing Mesh with Laplace Smoothing - '+str(nSmoothOut))
            for iS in range(nSmoothOut):
                # elementVolumesArray = self.model.q['abs(det(J))'][:,0]
                # elementBarycentersArray = self.mesh.elementBarycentersArray
                # ms.updateElementVolumes(elementVolumesArray_=elementVolumesArray,
                #                       elementNodesArray=self.mesh.elementNodesArray,
                #                       nodeArray=self.PHI)
                # ms.updateElementBarycenters(elementBarycentersArray_=elementBarycentersArray,
                #                             elementNodesArray=self.mesh.elementNodesArray,
                #                             nodeArray=self.PHI)
                simultaneous = True
                ms.smoothNodesLaplace(nodeArray_=xx,
                                      nodeStarOffsets=nodeStarOffsets,
                                      nodeStarArray=nodeStarArray,
                                      nodeMaterialTypes=nodeMaterialTypes,
                                      nNodes_owned=nNodes_owned,
                                      simultaneous=simultaneous,
                                      smoothBoundaries=True,
                                      fixedNodesBoolArray=fixedNodesBoolArray,
                                      alpha=0.)
            comm.barrier()
            logEvent('Done smoothing')

            ms.getNonOwnedNodeValues(xx,
                                     nNodes_owned,
                                     nNodes_global,
                                     nodeNumbering_subdomain2global,
                                     nodeOffsets_subdomain_owned)

        logEvent('Done pseudo-timestepping')

    # def evaluateFunAtX(self, x, ls_phi=None):
    #     pc = self.pyCoefficients
    #     f = self.cppEvaluateFunAtX(x=x,
    #                                he_min=pc.he_min,
    #                                he_max=pc.he_max,
    #                                ls_phi=ls_phi,
    #                                fun=pc.fun,
    #                                t=pc.t)

    # cdef double cppEvaluateFunAtX(self,
    #                               double[:] x,
    #                               double he_min,
    #                               double he_max,
    #                               double ls_phi=None,
    #                               object fun=None,
    #                               double t=None):
    #     cdef double f
    #     if fun:
    #         f = fun(x, self.t)
    #     else:
    #         f = 0.
    #     if ls_phi is not None:
    #         f = min(abs(ls_phi), f)
    #     f = max(he_min, f)
    #     f = min(he_max, f)
    #     return f

    # def evaluateFunAtNodes(self):
    #     pc = self.pyCoefficients
    #     self.cppEvaluateFunAtNodes(nodeArray=pc.mesh.nodeArray,
    #                                uOfXTatNodes=pc.uOfXTatNodes,
    #                                he_min=pc.he_min,
    #                                he_max=pc.he_max,
    #                                fun=pc.func,
    #                                t=pc.t,
    #                                u_phi=pc.u_phi)

    # cdef cppEvaluateFunAtNodes(self,
    #                            double[:] nodeArray,
    #                            double[:] uOfXTatNodes,
    #                            double he_min,
    #                            double he_max,
    #                            object fun=None,
    #                            double t=None,
    #                            double[:] u_phi=None):
    #     cdef double f
    #     for i in range(len(self.mesh.nodeArray)):
    #         if fun:
    #             f = fun(nodeArray[i], t)
    #         else:
    #             f = 0.
    #         if u_phi is not None:
    #             f = min(abs(u_phi[i]), f)
    #         f = max(he_min, f)
    #         f = min(he_max, f)
    #         uOfXTatNodes[i] = f

    # def evaluateFunAtQuadraturePoints(self):
    #     pc = self.pyCoefficients
    #     self.cppEvaluateFunAtQuadraturePoints(qx=pc.model.q['x'],
    #                                           uOfXTatQuadrature=pc.uOfXTatQuadrature,
    #                                           he_min=pc.he_min,
    #                                           he_max=pc.he_max,
    #                                           q_phi=pc.q_phi,
    #                                           fun=pc.myfunc,
    #                                           t=pc.t)

    # cdef cppEvaluateFunAtQuadraturePoints(self,
    #                                       double[:,:,:] qx,
    #                                       double[:,:] uOfXTatQuadrature,
    #                                       double he_min,
    #                                       double he_max,
    #                                       double[:,:] q_phi=None,
    #                                       fun=None,
    #                                       t=None):
    #     cdef int N_k = qx.shape[1]
    #     cdef double f
    #     for e in xrange(len(qx)):
    #         for k in xrange(N_k):
    #             if fun:
    #                 f = fun(qx[e, k], self.t)
    #             if q_phi is not None:
    #                 f = min(abs(q_phi[e, k]), f)
    #             f = max(he_min, f)
    #             f = min(he_max, f)
    #             self.uOfXTatQuadrature[e, k] = f


def recoveryAtNodes(variable,
                    nodeElementsArray,
                    nodeElementOffsets):
    return cppRecoveryAtNodes(variable=variable,
                              nodeElementsArray=nodeElementsArray,
                              nodeElementOffsets=nodeElementOffsets)

cdef double[:] cppRecoveryAtNodes(double[:] variable,
                                  int[:] nodeElementsArray,
                                  int[:] nodeElementOffsets):
    """
    variable:
         Variable in element
    """
    cdef double[:] recovered_variable = np.zeros(len(nodeElementOffsets)-1)
    cdef int nb_el
    cdef double var_sum
    for node in range(len(nodeElementOffsets)-1):
        nb_el = 0
        var_sum = 0
        for eOffset in range(nodeElementOffsets[node],
                             nodeElementOffsets[node+1]):
            nb_el += 1
            var_sum += variable[nodeElementsArray[eOffset]]
        recovered_variable[node] = var_sum/nb_el
    return recovered_variable

def gradientRecoveryAtNodes(grads,
                            nodeElementsArray,
                            nodeElementOffsets,
                            nd):
    return cppGradientRecoveryAtNodes(grads=grads,
                                      nodeElementsArray=nodeElementsArray,
                                      nodeElementOffsets=nodeElementOffsets,
                                      nd=nd)

cdef double[:,:] cppGradientRecoveryAtNodes(double[:,:,:] grads,
                                            int[:] nodeElementsArray,
                                            int[:] nodeElementOffsets,
                                            int nd):
    cdef double[:, :] recovered_grads = np.zeros((len(nodeElementOffsets)-1, nd))
    cdef int nb_el
    cdef double[:] grad_sum = np.zeros(nd)
    for node in range(len(nodeElementOffsets)-1):
        nb_el = 0
        for ndi in range(nd):
            grad_sum[ndi] = 0.
        for eOffset in range(nodeElementOffsets[node],
                             nodeElementOffsets[node+1]):
            nb_el += 1
            eN = nodeElementsArray[eOffset]
            # for k in range(n_quad):
            #     grad_k = grads[eN, k]
            #     grad_eN_av += grad_k
            # grad_eN_av /= n_quad
            for ndi in range(nd):
                grad_sum[ndi] += grads[eN, 0, ndi]  # same value at all quad points
        for ndi in range(nd):
            recovered_grads[node, ndi] = grad_sum[ndi]/nb_el
    return recovered_grads

def gradientRecoveryAtNodesWeighted(double[:,:,:] grads,
                                    int[:] nodeElementsArray,
                                    int[:] nodeElementOffsets,
                                    double[:,:] detJ_array,
                                    int nd):
    return cppGradientRecoveryAtNodesWeighted(grads=grads,
                                              nodeElementsArray=nodeElementsArray,
                                              nodeElementOffsets=nodeElementOffsets,
                                              detJ_array=detJ_array,
                                              nd=nd)

cdef double[:,:] cppGradientRecoveryAtNodesWeighted(double[:,:,:] grads,
                                                    int[:] nodeElementsArray,
                                                    int[:] nodeElementOffsets,
                                                    double[:,:] detJ_array,
                                                    int nd):
    cdef double[:, :] recovered_grads = np.zeros((len(nodeElementOffsets)-1, nd))
    cdef int nb_el
    cdef double detJ_patch = 0.
    cdef double[:] grad_sum = np.zeros(nd)
    for node in range(len(nodeElementOffsets)-1):
        nb_el = 0
        detJ_patch = 0.
        for ndi in range(nd):
            grad_sum[ndi] = 0.
        for eOffset in range(nodeElementOffsets[node],
                             nodeElementOffsets[node+1]):
            nb_el += 1
            eN = nodeElementsArray[eOffset]
            # for k in range(n_quad):
            #     grad_k = grads[eN, k]
            #     grad_eN_av += grad_k
            # grad_eN_av /= n_quad
            detJ_patch += detJ_array[eN,0]
            for ndi in range(nd):
                grad_sum[ndi] += detJ_array[eN,0]*grads[eN, 0, ndi]  # same value at all quad points
        for ndi in range(nd):
            recovered_grads[node, ndi] = grad_sum[ndi]/detJ_patch
    return recovered_grads

cdef tuple pyxSearchNearestNodeElementFromMeshObject(object mesh,
                                                     double[:] x,
                                                     int[:] exteriorElementBoundariesBoolArray,
                                                     double[:,:,:] elementBoundaryNormalsArray,
                                                     int nearest_node,
                                                     int eN):
    nearest_node = ms.getLocalNearestNode(coords=x,
                                          nodeArray=mesh.nodeArray,
                                          nodeStarOffsets=mesh.nodeStarOffsets,
                                          nodeStarArray=mesh.nodeStarArray,
                                          node=nearest_node)
    eN = ms.getLocalNearestElementAroundNode(coords=x,
                                             nodeElementOffsets=mesh.nodeElementOffsets,
                                             nodeElementsArray=mesh.nodeElementsArray,
                                             elementBarycenterArray=mesh.elementBarycentersArray,
                                             node=nearest_node)
    eN = ms.getLocalNearestElementIntersection(coords=x,
                                               elementBoundaryNormalsArray=elementBoundaryNormalsArray,
                                               elementBoundariesArray=mesh.elementBoundariesArray,
                                               elementBarycentersArray=mesh.elementBoundaryBarycentersArray,
                                               elementBoundaryElementsArray=mesh.elementBoundaryElementsArray,
                                               elementBarycentersArray=mesh.elementBarycentersArray,
                                               exteriorElementBoundariesBoolArray=exteriorElementBoundariesBoolArray,
                                               eN=eN)
    # check if actually inside eN
    inside_eN = True
    for j, b_i in enumerate(mesh.elementBoundariesArray[eN]):
        normal = mesh.elementBoundaryNormalsArray[eN, j]
        bound_bar = mesh.elementBoundaryBarycentersArray[b_i]
        dot = (bound_bar[0]-x[0])*normal[0]+(bound_bar[1]-x[1])*normal[1]+(bound_bar[2]-x[2])*normal[2]
        if dot < 0:
            inside_eN = False
    return eN, nearest_node, inside_eN


cdef tuple pyxSearchNearestNodeElement(double[:] x,
                                       double[:,:] nodeArray,
                                       int[:] nodeStarOffsets,
                                       int[:] nodeStarArray,
                                       int[:] nodeElementOffsets,
                                       int[:] nodeElementsArray,
                                       double[:,:] elementBarycentersArray,
                                       int[:,:] elementBoundaryElementsArray,
                                       int[:,:] elementBoundariesArray,
                                       double[:,:] elementBoundaryBarycentersArray,
                                       int[:] exteriorElementBoundariesBoolArray,
                                       double[:,:,:] elementBoundaryNormalsArray,
                                       int nearest_node,
                                       int eN,
                                       bool find_nearest_node=True):
    if find_nearest_node is True:
        nearest_node = ms.getLocalNearestNode(coords=x,
                                              nodeArray=nodeArray,
                                              nodeStarOffsets=nodeStarOffsets,
                                              nodeStarArray=nodeStarArray,
                                              node=nearest_node)
        eN = ms.getLocalNearestElementAroundNode(coords=x,
                                                 nodeElementOffsets=nodeElementOffsets,
                                                 nodeElementsArray=nodeElementsArray,
                                                 elementBarycentersArray=elementBarycentersArray,
                                                 node=nearest_node)
    eN = ms.getLocalNearestElementIntersection(coords=x,
                                               elementBoundaryNormalsArray=elementBoundaryNormalsArray,
                                               elementBoundariesArray=elementBoundariesArray,
                                               elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                               elementBoundaryElementsArray=elementBoundaryElementsArray,
                                               elementBarycentersArray=elementBarycentersArray,
                                               exteriorElementBoundariesBoolArray=exteriorElementBoundariesBoolArray,
                                               eN=eN)
    # check if actually inside eN
    cdef bool inside_eN = True
    cdef double[:] normal = np.zeros(3)
    for j, b_i in enumerate(elementBoundariesArray[eN]):
        normal[0] = elementBoundaryNormalsArray[eN, j, 0]
        normal[1] = elementBoundaryNormalsArray[eN, j, 1]
        normal[2] = elementBoundaryNormalsArray[eN, j, 2]
        bound_bar = elementBoundaryBarycentersArray[b_i]
        dot = (bound_bar[0]-x[0])*normal[0]+(bound_bar[1]-x[1])*normal[1]+(bound_bar[2]-x[2])*normal[2]
        if dot < 0:
            inside_eN = False
    return eN, nearest_node, inside_eN


def pyCheckOwnedVariable(int variable_nb_local,
                       int rank,
                       int nVariables_owned,
                       int[:] variableNumbering_subdomain2global,
                       int[:] variableOffsets_subdomain_owned):
    return checkOwnedVariable(variable_nb_local=variable_nb_local,
                              rank=rank,
                              nVariables_owned=nVariables_owned,
                              variableNumbering_subdomain2global=variableNumbering_subdomain2global,
                              variableOffsets_subdomain_owned=variableOffsets_subdomain_owned)

cdef tuple checkOwnedVariable(int variable_nb_local,
                              int rank,
                              int nVariables_owned,
                              int[:] variableNumbering_subdomain2global,
                              int[:] variableOffsets_subdomain_owned):
    cdef int nSubdomains = len(variableOffsets_subdomain_owned)-1
    cdef int variable_nb_global
    cdef int new_variable_nb_local
    cdef int new_rank = -2  # initialised as fake rank
    if variable_nb_local >= nVariables_owned:
        # change rank ownership
        variable_nb_global = variableNumbering_subdomain2global[variable_nb_local]
        if not variableOffsets_subdomain_owned[rank] <= variable_nb_global < variableOffsets_subdomain_owned[rank+1]:
            for i in range(nSubdomains+1):
                if variableOffsets_subdomain_owned[i] > variable_nb_global:
                    # changing processor
                    if new_rank == -2:
                        new_rank = i-1
    # getting nearest variable number on new rank
    if new_rank >= 0:
        new_variable_nb_local = variable_nb_global-variableOffsets_subdomain_owned[new_rank]
    else:
        new_rank = rank
        new_variable_nb_local = variable_nb_local
    return new_variable_nb_local, new_rank

cdef tuple getGlobalVariable(int variable_nb_local,
                             int nVariables_owned,
                             int[:] variableNumbering_subdomain2global,
                             int[:] variableOffsets_subdomain_owned):
    cdef int nSubdomains = len(variableOffsets_subdomain_owned)-1
    cdef int variable_nb_global
    cdef int new_rank = -2  # initialised as fake rank
    # change rank ownership
    variable_nb_global = variableNumbering_subdomain2global[variable_nb_local]
    for i in range(nSubdomains+1):
        if variableOffsets_subdomain_owned[i] > variable_nb_global:
            # changing processor
            if new_rank == -2:
                new_rank = i-1
    variable_nb_global = variableNumbering_subdomain2global[variable_nb_local]
    return variable_nb_global, new_rank

cdef int getLocalVariable(int variable_nb_global,
                          int rank,
                          int nVariables_owned,
                          int[:] variableNumbering_subdomain2global,
                          int[:] variableOffsets_subdomain_owned):
    cdef int new_variable_nb_local
    if not variableOffsets_subdomain_owned[rank] <= variable_nb_global < variableOffsets_subdomain_owned[rank+1]:
        print('wrong global varibale nb', rank, variable_nb_global, variableOffsets_subdomain_owned[rank])
        import pdb; pdb.set_trace()
    else:
        new_variable_nb_local = variable_nb_global-variableOffsets_subdomain_owned[rank]
    return new_variable_nb_local

# cdef retrieveSolution(double[:]x,
#                       double[:,:] nodeArray,
#                       int[:] nodeStarOffsets,
#                       int[:] nodeStarArray,
#                       int nearest_node,
#                       int rank,
#                       int nNodes_owned,
# )
#                         nearest_node = ms.getLocalNearestNode(coords=x,
#                                                               nodeArray=nodeArray,
#                                                               nodeStarOffsets=nodeStarOffsets,
#                                                               nodeStarArray=nodeStarArray,
#                                                               node=nearest_node)
#                         if nearest_node >= nNodes_owned:
#                             nearest_node_new_rank, new_rank = checkOwnedVariable(variable_nb_local=nearest_node,
#                                                                                  rank=my_rank,
#                                                                                  nVariables_owned=nNodes_owned,
#                                                                                  variablesNumering_subdomain2global=nodeNumbering_subdomain2global,
#                                                                                  variableOffsets_subdomain_owned=nodeOffsets_subdomain_owned)
#                             nodes0_2rank[new_rank] = np.append(nodes0_2rank[new_rank], nearest_node_new_rank)
#                             nPending_disp += 1
#                         else:

#                             nearest_eN = ms.getLocalNearestElementAroundNode(coords=x,
#                                                                              nodeElementOffsets=nodeElementOffsets,
#                                                                              nodeElementsArray=nodeElementsArray,
#                                                                              elementBarycentersArray=elementBarycentersArray,
#                                                                              node=nearest_node)
#                             nPending_disp += 1
#                             if nearest_eN >= nElements_owned:
#                                 nearest_eN_new_rank, new_rank = checkOwnedVariable(variable_nb_local=nearest_eN,
#                                                                                    rank=my_rank,
#                                                                                    nVariables_owned=nENs_owned,
#                                                                                    variablesNumering_subdomain2global=eNNumbering_subdomain2global,
#                                                                                    variableOffsets_subdomain_owned=eNOffsets_subdomain_owned)
#                                 elements2rank[new_rank] = np.append(elements2rank[new_rank], nearest_eN_new_rank)
#                                 nPending_disp += 1
#                             else:
#                                 nearest_eN = ms.getLocalNearestElementIntersection(coords=x,
#                                                                                    elementBoundaryNormalsArray=elementBoundaryNormalsArray,
#                                                                                    elementBoundariesArray=elementBoundariesArray,
#                                                                                    elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
#                                                                                    elementNeighborsArray=elementNeighborsArray,
#                                                                                    elementBarycentersArray=elementBarycentersArray,
#                                                                                    exteriorElementBoundariesBoolArray=exteriorElementBoundariesBoolArray,
#                                                                                    eN=enearest_eN)
#                                 if nearest_eN >= nElements_owned:
#                                     nearest_eN_new_rank, new_rank = checkOwnedVariable(variable_nb_local=nearest_eN,
#                                                                                        rank=my_rank,
#                                                                                        nVariables_owned=nENs_owned,
#                                                                                        variablesNumering_subdomain2global=eNNumbering_subdomain2global,
#                                                                                        variableOffsets_subdomain_owned=eNOffsets_subdomain_owned)
#                                     elements2rank[new_rank] = np.append(elements2rank[new_rank], nearest_eN_new_rank)
#                                     nPending_disp += 1
