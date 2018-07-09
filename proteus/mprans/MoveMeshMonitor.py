import numpy as np
import copy
from proteus.Profiling import logEvent
from proteus import TransportCoefficients
from proteus.mprans import MeshSmoothing as ms
from proteus.mprans import cMoveMeshMonitor as cmm

class Coefficients(TransportCoefficients.PoissonEquationCoefficients):
    """
    Coefficients for deforming the mesh according to an area function

    Parameters
    ----------
    func: function
        function of x defining element area
    he_min: double
        minimum characteristic element size
    he_max: double
        maximum characteristic element size
    ME_MODEL: int
        index of this model
    LS_MODEL: int
        index of LS model for retrieving phi values and refine around free
        surface
    nd: int
        number of dimensions
    boundaryNormals: double[:,:]
        array of length of all domain material types, with indexes corresponding
        to the material type set as boundaryNormals.
        if e.g., boundaryNormals[2]=[0.,0.,0.], then it will be ignored
    fixedMaterialTypes: int[:]
        array of length of all domain material types, with indexes corresponding
        to the material type set as: 0 -> not fixed, 1 -> fixed.
        e.g. is nodes of material type 3 should be fixed: fixedMaterialTypes[3]=1.
        (!) must be integers
    nSmoothIn: int
        number of smoothing steps while pseudo time stepping (during each step)
    nSmoothOut: int
        number of smoothing steps after finishing pseudo time stepping
    epsFact_density: double
        epsFact*he_min where refinement is the same
    epsTimeStep: double
        pseudo time step size (0 < epsTimeStep <= 1)
    """
    def __init__(self,
                 func,
                 he_min,
                 he_max,
                 ME_MODEL,
                 LS_MODEL=None,
                 nd=2,
                 boundaryNormals=None,
                 fixedMaterialTypes=None,
                 nSmoothIn=0,
                 nSmoothOut=0,
                 epsFact_density=3,
                 epsTimeStep=1.,
                 grading=1.1):
        self.myfunc = func
        self.LS_MODEL = LS_MODEL
        self.ME_MODEL = ME_MODEL
        self.he_min = he_min
        self.he_max = he_max
        self.grading = grading
        self.C = 1.  # scaling coefficient for f (computed in preStep)
        self.integral_area = 1.
        self.nd = nd
        assert 0 < epsTimeStep <= 1, 'epsTimeStep must be between 0. and 1.'
        self.epsTimeStep = epsTimeStep
        self.epsFact_density = epsFact_density
        aOfX = [lambda x: np.array([[1., 0.], [0., 1.]])]
        fOfX = [self.uOfX]  # scaled function reciprocal
        self.t = 0
        self.boundaryNormals = boundaryNormals
        self.fixedMaterialTypes = fixedMaterialTypes
        self.nSmoothIn = nSmoothIn
        self.nSmoothOut = nSmoothOut
        self.dt_last = None
        self.u_phi = None
        #super(MyCoeff, self).__init__(aOfX, fOfX)
        TransportCoefficients.PoissonEquationCoefficients.__init__(self, aOfX, fOfX)

    def initializeMesh(self, mesh):
        self.mesh = mesh
        self.areas_nodes = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.areas = np.zeros(self.mesh.elementNodesArray.shape[0])
        # get Normals
        if self.nd == 2:
            # triangle
            self.mesh.elementBoundaryNormalsArray = np.zeros((len(self.mesh.elementBoundariesArray), 3, 3))
        elif self.nd == 3:
            # tetra
            self.mesh.elementBoundaryNormalsArray = np.zeros((len(self.mesh.elementBoundariesArray), 4, 3))
        if self.nd == 2:
            ms.updateElementBoundaryNormalsTriangle(self.mesh.elementBoundaryNormalsArray,
                                                    self.mesh.nodeArray,
                                                    self.mesh.elementBoundariesArray,
                                                    self.mesh.elementBoundaryNodesArray,
                                                    self.mesh.elementBoundaryBarycentersArray,
                                                    self.mesh.elementBarycentersArray)
        elif self.nd == 3:
            ms.updateElementBoundaryNormalsTetra(self.mesh.elementBoundaryNormalsArray,
                                                 self.mesh.nodeArray,
                                                 self.mesh.elementBoundariesArray,
                                                 self.mesh.elementBoundaryNodesArray,
                                                 self.mesh.elementBoundaryBarycentersArray,
                                                 self.mesh.elementBarycentersArray)
        # nodes to keep fixed
        self.mesh.fixedNodesBoolArray = np.zeros(len(self.mesh.nodeArray), dtype=np.int32)
        # find corner nodes
        self.mesh.cornerNodes = ms.getCornerNodesTriangle(nodeArray=self.mesh.nodeArray,
                                                          nodeStarArray=self.mesh.nodeStarArray,
                                                          nodeStarOffsets=self.mesh.nodeStarOffsets,
                                                          nodeMaterialTypes=self.mesh.nodeMaterialTypes,
                                                          nNodes_owned=self.mesh.nNodes_owned)
        for cn in self.mesh.cornerNodes:
            self.mesh.fixedNodesBoolArray[cn] = 1
        for node in range(len(self.mesh.nodeMaterialTypes)):
            if self.fixedMaterialTypes[self.mesh.nodeMaterialTypes[node]] == 1:
                self.mesh.fixedNodesBoolArray[node] = 1
        uni, count = np.unique(self.mesh.elementBoundaryElementsArray.flat, return_counts=True)

    def attachModels(self,modelList):
        self.model = modelList[self.ME_MODEL]
        self.model.bdyNullSpace = True
        self.q = self.model.q
        self.model.mesh_array0 = copy.deepcopy(self.mesh.nodeArray)
        self.PHI = np.zeros_like(self.model.mesh.nodeArray)
        self.grads = np.zeros((len(self.model.u[0].dof),2))
        self.areas_nodes = np.zeros(len(self.mesh.nodeArray))
        self.areas = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.nearest_nodes = np.array([i for i in range(len(self.mesh.nodeArray))])
        self.eN_phi = np.array([None for i in range(len(self.mesh.nodeArray))])  # element containing displaced phi (guess)
        self.nearest_nodes0 = np.array([i for i in range(len(self.mesh.nodeArray))])
        self.mesh.nodeVelocityArray = np.zeros(self.mesh.nodeArray.shape,'d')
        self.uOfXTatNodes = np.zeros(len(self.mesh.nodeArray))
        self.uOfXTatQuadrature = np.zeros((self.model.q['x'].shape[0], self.model.q['x'].shape[1]))
        if self.LS_MODEL is not None:
            self.model_ls = modelList[self.LS_MODEL]
            self.u_phi = self.model_ls.u[0].dof
            self.q_phi = self.model_ls.q[('u', 0)]
        self.cCoefficients = cmm.cCoefficients()
        self.cCoefficients.attachPyCoefficients(self)

    def preStep(self, t, firstStep=False):
        logEvent("Updating area function scaling coefficient and quadrature",
                 level=3)
        self.evaluateFunAtQuadraturePoints()
        self.cCoefficients.preStep()
        self.C = self.cCoefficients.C
        self.model.areas = self.areas
        logEvent("Finished updating area function",
                 level=3)
        if t != self.t:
            self.t = t
            self.model.mesh.nodeVelocityArray[:] = 0.

    def postStep(self, t,firstStep=False):
        # gradient recovery
        logEvent("Gradient recovery at mesh nodes", level=3)
        self.grads[:] = cmm.gradientRecoveryAtNodes(grads=self.model.q[('grad(u)', 0)],
                                                    nodeElementsArray=self.mesh.nodeElementsArray,
                                                    nodeElementOffsets=self.mesh.nodeElementOffsets,
                                                    nd=self.nd)
        logEvent("Area recovery at mesh nodes", level=3)
        self.areas_nodes[:] = cmm.recoveryAtNodes(variable=self.areas,
                                                  nodeElementsArray=self.mesh.nodeElementsArray,
                                                  nodeElementOffsets=self.mesh.nodeElementOffsets)

        nNodes_owned = self.mesh.nNodes_owned
        nNodes_global = self.mesh.nNodes_global
        from proteus import Comm
        comm = Comm.get().comm.tompi4py()
        comm_size = comm.size
        my_rank = comm.rank
        grads_2rank = {}
        nodes0_2rank = {}
        rank0_2rank = {}
        nearestN_2rank = {}
        areas_2rank = {}
        if comm_size > 1:
            logEvent("Sharing recovered value to non-owned nodes", level=3)
            comm.barrier()
            for rank in range(comm_size):
                grads_2rank[rank] = np.zeros(0, dtype=np.int32)
                nodes0_2rank[rank] = np.zeros(0, dtype=np.int32)
                nearestN_2rank[rank] = np.zeros(0, dtype=np.int32)
                rank0_2rank[rank] = np.zeros(0, dtype=np.int32)
                areas_2rank[rank] = np.zeros(0)
            for node in range(nNodes_owned, nNodes_global):
                node_new_rank, new_rank = cmm.pyCheckOwnedVariable(variable_nb_local=node,
                                                             rank=my_rank,
                                                             nVariables_owned=nNodes_owned,
                                                             variableNumbering_subdomain2global=self.mesh.globalMesh.nodeNumbering_subdomain2global,
                                                             variableOffsets_subdomain_owned=np.array([ss for ss in self.mesh.globalMesh.nodeOffsets_subdomain_owned], dtype=np.int32))
                nodes0_2rank[new_rank] = np.append(nodes0_2rank[new_rank], node)
                nearestN_2rank[new_rank] = np.append(nearestN_2rank[new_rank], node_new_rank)
                rank0_2rank[new_rank] = np.append(rank0_2rank[new_rank], my_rank)
            # SEND THOSE NODES TO RELEVANT PROCESSORS
            nodes0_2doArray = np.zeros(0, dtype=np.int32)
            rank0_2doArray = np.zeros(0, dtype=np.int32)
            nearestN_2doArray = np.zeros(0, dtype=np.int32)
            for rank_recv in range(comm.size):
                for rank_send in range(comm.size):
                    if rank_send != rank_recv and rank_send == my_rank:
                        comm.send(nodes0_2rank[rank_recv].size, dest=rank_recv, tag=0)
                        if nodes0_2rank[rank_recv].size > 0:
                            # original nodes
                            comm.send(nodes0_2rank[rank_recv],  dest=rank_recv, tag=1)
                            # nodes on other processor
                            comm.send(nearestN_2rank[rank_recv], dest=rank_recv, tag=2)
                            # rank for original nodes (to send back final solution)
                            comm.send(rank0_2rank[rank_recv], dest=rank_recv, tag=3)
                    elif rank_send!= rank_recv and rank_recv == my_rank:
                        size = comm.recv(source=rank_send, tag=0)
                        if size > 0:
                            # coords
                            nodes0_2do = comm.recv(source=rank_send, tag=1)
                            nodes0_2doArray = np.append(nodes0_2doArray, nodes0_2do, axis=0)
                            # original nodes
                            nearestN_2do = comm.recv(source=rank_send, tag=2)
                            nearestN_2doArray = np.append(nearestN_2doArray, nearestN_2do)
                            # original ranks
                            rank0_2do = comm.recv(source=rank_send, tag=3)
                            rank0_2doArray = np.append(rank0_2doArray, rank0_2do)
                    comm.barrier()
            # SEND VALUE BACK TO ORIGINAL PROCESSORS
            for rank in range(comm.size):
                grads_2rank[rank] = np.zeros((0, 3))
                nodes0_2rank[rank] = np.zeros(0, dtype=np.int32)
                areas_2rank[rank] = np.zeros(0)
            for iN in range(len(nodes0_2doArray)):
                node = nodes0_2doArray[iN]
                nearestN = nearestN_2doArray[iN]
                grad = np.zeros(3)
                for ind in range(self.nd):
                    grad[ind] = self.grads[nearestN_2doArray[iN], ind]
                grads_2rank[rank0_2doArray[iN]] = np.append(grads_2rank[rank0_2doArray[iN]], [grad], axis=0)
                nodes0_2rank[rank0_2doArray[iN]] = np.append(nodes0_2rank[rank0_2doArray[iN]], nodes0_2doArray[iN])
                areas_2rank[rank0_2doArray[iN]] = np.append(areas_2rank[rank0_2doArray[iN]], self.areas_nodes[nearestN])
            # retrieve solution
            nodes0_2doArray = np.zeros(0)
            grads_2doArray = np.zeros((0, 3))
            areas_2doArray = np.zeros(0)
            for rank_recv in range(comm.size):
                for rank_send in range(comm.size):
                    if rank_send != rank_recv and rank_send == my_rank:
                        comm.send(nodes0_2rank[rank_recv].size, dest=rank_recv, tag=0)
                        if nodes0_2rank[rank_recv].size > 0:
                            # original nodes
                            comm.send(nodes0_2rank[rank_recv],  dest=rank_recv, tag=1)
                            # grads to other processor
                            comm.send(grads_2rank[rank_recv], dest=rank_recv, tag=2)
                            # areas to other processors
                            comm.send(areas_2rank[rank_recv], dest=rank_recv, tag=3)
                    elif rank_send!= rank_recv and rank_recv == my_rank:
                        size = comm.recv(source=rank_send, tag=0)
                        if size > 0:
                            # original nodes
                            nodes0_2do = comm.recv(source=rank_send, tag=1)
                            nodes0_2doArray = np.append(nodes0_2doArray, nodes0_2do)
                            # grads
                            grads_2do = comm.recv(source=rank_send, tag=2)
                            grads_2doArray = np.append(grads_2doArray, grads_2do, axis=0)
                            # areas
                            areas_2do = comm.recv(source=rank_send, tag=3)
                            areas_2doArray = np.append(areas_2doArray, areas_2do, axis=0)
                    comm.barrier()
            # FINALLY APPLY FINAL POSITION OF NON-OWNED NODES
            for iN in range(len(nodes0_2doArray)):
                node = int(nodes0_2doArray[iN])
                grad = grads_2doArray[iN]
                area = areas_2doArray[iN]
                self.areas_nodes[node] = area
                for ndi in range(self.nd):
                    self.grads[node, ndi] = grad[ndi]
        self.model.grads = self.grads
        self.model.areas_nodes = self.areas_nodes
        # evaluate f at nodes
        self.evaluateFunAtNodes()
        # gamma parameter
        f_over_g = self.uOfXTatNodes/self.areas_nodes
        f_over_g_max = max(f_over_g)
        f_over_g_min = min(f_over_g)
        self.gamma = f_over_g_max/f_over_g_min
        self.gamma0 = 10.  # user defined parameter
        self.na = np.log(self.gamma)/np.log(self.gamma0)
        nNodes_owned = self.mesh.nNodes_owned
        if self.dt_last is not None:
            # pseudo-time step
            self.PHI[:] = self.mesh.nodeArray
            self.cCoefficients.pseudoTimeStepping(eps=self.epsTimeStep,
                                                  xx=self.PHI)
            # dt = self.model.timeIntegration.dt
            # else:
            #     dt = self.dt_last
            # dt = self.model.timeIntegration.dt
            dt = self.dt_last
            # move nodes
            self.model.mesh.nodeVelocityArray[:] = (self.PHI[:]-self.model.mesh.nodeArray[:])/dt
            #self.model.mesh.nodeVelocityArray[:] = np.zeros(self.model.mesh.nodeArray.shape)
            self.model.mesh.nodeArray[:] = self.PHI[:]
            # re-initialise nearest nodes
            self.nearest_nodes[:] = self.nearest_nodes0[:]
            # re-initialise containing element
            self.eN_phi[:] = None
        self.dt_last = self.model.timeIntegration.dt
        self.cCoefficients.postStep()
        # copyInstructions = {'clear_uList': True}
        # return copyInstructions
        IMR_nodes = ms.getInverseMeanRatioTriangleNodes(nodeArray=self.mesh.nodeArray,
                                                        elementNodesArray=self.mesh.elementNodesArray,
                                                        nodeElementOffsets=self.mesh.nodeElementOffsets,
                                                        nodeElementsArray=self.mesh.nodeElementsArray,
                                                        el_average=False)
        # self.model.u[0].dof[:] = IMR_nodes

    def uOfX(self, x, eN=None):
        """
        reaction term
        """
        # f_scaled = (self.myfunc(x)*self.C)
        # f_scaled_reciprocal = 1./f_scaled
        if eN is None:
            area = 1.
        else:
            area = self.areas[eN]
        f = 1./(self.evaluateFunAtX(x)*self.C)-1./area
        return f

    def dist_search(self, x, node, eN=-1):
        """
        search element containing containing coords x starting with a guessed
        nearest node

        Parameters
        ----------
        x: array_like
            coordinates
        node: int
            first guess of the nearest node
        """
        femSpace = self.model.u[0].femSpace
        checkedElements=[]
        patchBoundaryNodes=set()

        nearest_node = node
        min_dist = np.sqrt(np.sum((x-self.mesh.nodeArray[node])**2))
        found_node = False
        element = None
        xx = np.zeros(3)
        xx[0] = x[0]
        xx[1] = x[1]
        xx[2] = x[2]
        x = xx
        if eN >= 0:
            xi = femSpace.elementMaps.getInverseValue(eN, x)
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                element = eN
        if element < 0:
            while found_node is False:
                nearest_node0 = nearest_node
                for nOffset in range(femSpace.mesh.nodeStarOffsets[nearest_node0],
                                     femSpace.mesh.nodeStarOffsets[nearest_node0+1]):
                    node = self.mesh.nodeStarArray[nOffset]
                    node_coords = self.mesh.nodeArray[node]
                    dist = np.sqrt(np.sum((node_coords-x)**2))
                    if dist < min_dist:
                        min_dist = dist
                        nearest_node = node
                if nearest_node0 == nearest_node:
                    found_node = True
            # now find element containing x around nearest_node
            for eOffset in range(self.mesh.nodeElementOffsets[nearest_node],
                                 self.mesh.nodeElementOffsets[nearest_node+1]):
                eN = self.mesh.nodeElementsArray[eOffset]
                checkedElements.append(eN)
                # union of set
                nodes = self.mesh.nodeStarArray[femSpace.mesh.nodeStarOffsets[nearest_node]:femSpace.mesh.nodeStarOffsets[nearest_node+1]]
                patchBoundaryNodes|=set(nodes)
                # evaluate the inverse map for element eN (global to local)
                xi = femSpace.elementMaps.getInverseValue(eN, x)
                # query whether xi lies within the reference element
                if femSpace.elementMaps.referenceElement.onElement(xi):
                    element = eN
            if element < 0:
                # extra loop if case coords is in neighbour element
                for node in patchBoundaryNodes:
                    for eOffset in range(self.mesh.nodeElementOffsets[node],
                                         self.mesh.nodeElementOffsets[node+1]):
                        eN = self.mesh.nodeElementsArray[eOffset]
            #            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
                        if eN not in checkedElements:
                            checkedElements.append(eN)
                            # evaluate the inverse map for element eN
                            xi = femSpace.elementMaps.getInverseValue(eN, x)
                            # query whether xi lies within the reference element
                            if femSpace.elementMaps.referenceElement.onElement(xi):
                                element = eN
        if element is None:
            element = -1
        return element, nearest_node, xi

    def evaluateFunAtX(self, x, ls_phi=None):
        """
        reaction term
        """
        # f_scaled = (self.myfunc(x)*self.C)
        # f_scaled_reciprocal = 1./f_scaled
        if ls_phi is not None:
            f = min(abs(ls_phi), self.myfunc(x, self.t))
            # f = max(self.he_min, f)
            # f = min(self.he_max, f)
        else:
            f = self.myfunc(x, self.t)
            # f = min(self.he_max, f)
        if f < self.epsFact_density*self.he_min:
            f = 0.
        f = self.he_min*self.grading**(1.+np.log((-1./self.grading*(f-self.he_min)+f)/self.he_min)/np.log(self.grading))
        f = min(f, self.he_max)
        # f = min(self.he_max, self.he_min*self.grading**(f/self.he_min))
        return f

    def evaluateFunAtNodes(self):
        for i in range(len(self.mesh.nodeArray)):
            if self.LS_MODEL is not None:
                f = min(abs(self.u_phi[i]),
                        self.myfunc(self.mesh.nodeArray[i], self.t))
                # f = max(self.he_min, f)
                # f = min(self.he_max, f)
            else:
                f = self.myfunc(self.mesh.nodeArray[i], self.t)
                # f = min(self.he_max, f)
            if f < self.epsFact_density*self.he_min:
                f = 0.
            f = self.he_min*self.grading**(1.+np.log((-1./self.grading*(f-self.he_min)+f)/self.he_min)/np.log(self.grading))
            f = min(f, self.he_max)
            # f = min(self.he_max, self.he_min*self.grading**(f/self.he_min))
            self.uOfXTatNodes[i] = f

    def evaluateFunAtQuadraturePoints(self):
        xx = self.model.q['x']
        N_k = self.model.q['x'].shape[1]
        for e in xrange(self.mesh.elementNodesArray.shape[0]):
            for k in xrange(N_k):
                if self.LS_MODEL is not None:
                    f = min(abs(self.q_phi[e, k]),
                            self.myfunc(xx[e, k], self.t))
                    # f = max(self.he_min, f)
                    # f = min(self.he_max, f)
                else:
                    f = self.myfunc(xx[e, k], self.t)
                    # f = min(self.he_max, f)
                if f < self.epsFact_density*self.he_min:
                    f = 0.
                f = self.he_min*self.grading**(1.+np.log((-1./self.grading*(f-self.he_min)+f)/self.he_min)/np.log(self.grading))
                f = min(f, self.he_max)
                # f = min(self.he_max, self.he_min*self.grading**(f/self.he_min))
                self.uOfXTatQuadrature[e, k] = f

    def getLevelSetValue(self, eN, xi):
        value = self.model_ls.u[0].getValue(eN, xi)
        return value

    def getGradientValue(self, eN, xi):
        femSpace = self.model.u[0].femSpace
        value = 0.0
        # #  tridelat hack: the following 3 lines do not work in parallel
        # for i,psi in zip(femSpace.dofMap.l2g[eN],
        #                  femSpace.elementMaps.localFunctionSpace.basis):
        #     value+=self.grads[i]*psi(xi)
        # # reverting to grad value within element
        value = self.model.q[('grad(u)', 0)][eN, 0]
        return value

def calculate_areas(x, detJ, weights):
        N_k = len(weights)
        areas = np.zeros(len(x))
        for e in xrange(len(x)):
            area = 0
            for k in xrange(N_k):
                area += detJ[e, k]*weights[k]
            areas[e] = area
        return area

def calculate_area(x, detJ, weights):
        N_k = len(weights)
        area = 0
        for k in xrange(N_k):
            area += detJ[k]*weights[k]
        return area
