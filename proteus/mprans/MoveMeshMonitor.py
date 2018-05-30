import numpy as np
import copy
from proteus.Profiling import logEvent
from proteus import TransportCoefficients
import cMoveMeshMonitor as cmm

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
    fixedNodes: int[:]
        array of length of all domain material types, with indexes corresponding
        to the material type set as: 0 -> not fixed, 1 -> fixed.
        e.g. is nodes of material type 3 should be fixed: fixedNodes[3]=1.
        (!) must be integers
    nSmooth: int
        number of smoothing steps after solving
    """
    def __init__(self,
                 func,
                 he_min,
                 he_max,
                 ME_MODEL,
                 LS_MODEL=None,
                 nd=2,
                 boundaryNormals=None,
                 fixedNodes=None,
                 nSmooth=0,
                 epsFact=3):
        self.myfunc = func
        self.LS_MODEL = LS_MODEL
        self.ME_MODEL = ME_MODEL
        self.he_min = he_min
        self.he_max = he_max
        self.C = 1.  # scaling coefficient for f (computed in preStep)
        self.integral_area = 1.
        self.nd = nd
        aOfX = [lambda x: np.array([[1., 0.], [0., 1.]])]
        fOfX = [self.uOfX]  # scaled function reciprocal
        self.t = 0
        self.boundaryNormals = boundaryNormals
        self.fixedNodes = fixedNodes
        self.nSmooth = nSmooth
        self.dt_last = None
        #super(MyCoeff, self).__init__(aOfX, fOfX)
        TransportCoefficients.PoissonEquationCoefficients.__init__(self, aOfX, fOfX)

    def initializeMesh(self, mesh):
        self.mesh = mesh
        self.areas_nodes = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.areas = np.zeros(self.mesh.elementNodesArray.shape[0])

    def attachModels(self,modelList):
        self.model = modelList[self.ME_MODEL]
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

    def postStep(self, t,firstStep=False):
        # gradient recovery
        self.gradientRecoveryAtNodes()
        # area recovery
        self.areaRecoveryAtNodes()
        # evaluate f at nodes
        self.evaluateFunAtNodes()
        # gamma parameter
        f_over_g = self.uOfXTatNodes/self.areas_nodes
        f_over_g_max = max(f_over_g)
        f_over_g_min = min(f_over_g)
        self.gamma = f_over_g_max/f_over_g_min
        self.gamma0 = 10.  # user defined parameter
        self.na = np.log(self.gamma)/np.log(self.gamma0)
        # pseudo-time step
        eps = 1.
        self.pseudoTimeStepping(eps)
        # move nodes
        if self.dt_last is not None:
            # dt = self.model.timeIntegration.dt
            # else:
            #     dt = self.dt_last
            # dt = self.model.timeIntegration.dt
            dt = self.dt_last
            logEvent('Smoothing Mesh with Laplace Smoothing - '+str(self.nSmooth))
            from proteus.mprans import MeshSmoothing
            if self.nSmooth > 0:
                disp = MeshSmoothing.smoothNodesLaplace(nodeArray=self.PHI,
                                                 nodeStarOffsets=self.mesh.nodeStarOffsets,
                                                 nodeStarArray=self.mesh.nodeStarArray,
                                                 nodeMaterialTypes=self.mesh.nodeMaterialTypes,
                                                 nNodes_owned=self.mesh.nNodes_owned,
                                                 nSmooth=self.nSmooth,
                                                 boundaryNormals=self.boundaryNormals,
                                                 fixedNodes=self.fixedNodes,
                                                 apply_directly=False)
                self.PHI += disp
            self.model.mesh.nodeVelocityArray[:] = (self.PHI-self.model.mesh.nodeArray)/dt
            #self.model.mesh.nodeVelocityArray[:] = disp/dt
            #self.model.mesh.nodeVelocityArray[:] = np.zeros(self.model.mesh.nodeArray.shape)
            self.model.mesh.nodeArray[:] = self.PHI
            # re-initialise nearest nodes
            self.nearest_nodes[:] = self.nearest_nodes0[:]
            # re-initialise containing element
            self.eN_phi[:] = None
        self.dt_last = self.model.timeIntegration.dt
        # copyInstructions = {'clear_uList': True}
        # return copyInstructions

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

    def dist_search(self, x, node, eN=None):
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
        if eN is not None:
            xi = femSpace.elementMaps.getInverseValue(eN, x)
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                element = eN
        if element is None:
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
            if element is None:
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
            #if element is None:
            #    while element is None:
            #    # extra loop if case coords is in neighbour element
            #        print("patch", patchBoundaryNodes)
            #        for node in copy.copy(patchBoundaryNodes):
            #            for eOffset in range(self.mesh.nodeElementOffsets[node],
            #                                 self.mesh.nodeElementOffsets[node+1]):
            #                eN = self.mesh.nodeElementsArray[eOffset]
            #                patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            #                if eN not in checkedElements:
            #                    checkedElements.append(eN)
            #                    # evaluate the inverse map for element eN
            #                    xi = femSpace.elementMaps.getInverseValue(eN, coords)
            #                    # query whether xi lies within the reference element
            #                    if femSpace.elementMaps.referenceElement.onElement(xi):
            #                        element = eN
        return element, nearest_node, xi

    def evaluateFunAtX(self, x, ls_phi=None):
        """
        reaction term
        """
        # f_scaled = (self.myfunc(x)*self.C)
        # f_scaled_reciprocal = 1./f_scaled
        if ls_phi is not None:
            f = min(abs(ls_phi), self.myfunc(x, self.t))
            f = max(self.he_min, f)
            f = min(self.he_max, f)
        else:
            f = max(self.he_min, self.myfunc(x, self.t))
            f = min(self.he_max, f)
        return f

    def evaluateFunAtNodes(self):
        for i in range(len(self.mesh.nodeArray)):
            if self.LS_MODEL is not None:
                f = min(abs(self.u_phi[i]),
                        self.myfunc(self.mesh.nodeArray[i], self.t))
                f = max(self.he_min, f)
                f = min(self.he_max, f)
            else:
                f = max(self.he_min,
                        self.myfunc(self.mesh.nodeArray[i], self.t))
                f = min(self.he_max, f)
            self.uOfXTatNodes[i] = f

    def evaluateFunAtQuadraturePoints(self):
        xx = self.model.q['x']
        N_k = self.model.q['x'].shape[1]
        for e in xrange(self.mesh.elementNodesArray.shape[0]):
            for k in xrange(N_k):
                if self.LS_MODEL is not None:
                    f = min(abs(self.q_phi[e, k]),
                            self.myfunc(xx[e, k], self.t))
                    f = max(self.he_min, f)
                    f = min(self.he_max, f)
                else:
                    f = max(self.he_min,
                            self.myfunc(xx[e, k], self.t))
                    f = min(self.he_max, f)
                self.uOfXTatQuadrature[e, k] = f

    def getLevelSetValue(self, eN, xi):
        value = self.model_ls.u[0].getValue(eN, xi)
        return value

    def getGradientValue(self, eN, xi):
        femSpace = self.model.u[0].femSpace
        value = 0.0
        for i,psi in zip(
            femSpace.dofMap.l2g[eN],
            femSpace.elementMaps.localFunctionSpace.basis):
            value+=self.grads[i]*psi(xi)
        return value

    def gradientRecoveryAtNodes(self):
        logEvent("Gradient recovery at mesh nodes",
                 level=3)
        grad_v = self.model.q[('grad(u)',0)]
        # n_quad = self.mesh.nodeArray.shape[1]
        for node in range(len(self.mesh.nodeArray)):
            nb_el = 0
            grad_av = 0.
            for eOffset in range(self.mesh.nodeElementOffsets[node],
                                 self.mesh.nodeElementOffsets[node+1]):
                nb_el += 1
                eN = self.mesh.nodeElementsArray[eOffset]
                grad_eN_av = 0
                # for k in range(n_quad):
                #     grad_k = self.model.q[('grad(u)',0)][eN, k]
                #     grad_eN_av += grad_k
                # grad_eN_av /= n_quad
                grad_eN_av = self.model.q[('grad(u)',0)][eN, 0]  # same value at all quad points
                grad_av += grad_eN_av
            grad_av /= nb_el
            self.grads[node] = grad_av
        self.model.grads = self.grads

    def areaRecoveryAtNodes(self):
        logEvent("Area recovery at mesh nodes",
                 level=3)
        grad_v = self.model.q[('grad(u)',0)]
        # n_quad = self.mesh.nodeArray.shape[1]
        for node in range(len(self.mesh.nodeArray)):
            nb_el = 0
            grad_av = 0.
            for eOffset in range(self.mesh.nodeElementOffsets[node],
                                 self.mesh.nodeElementOffsets[node+1]):
                nb_el += 1
                eN = self.mesh.nodeElementsArray[eOffset]
                grad_eN_av = 0
                # for k in range(n_quad):
                #     grad_k = self.model.q[('grad(u)',0)][eN, k]
                #     grad_eN_av += grad_k
                # grad_eN_av /= n_quad
                grad_eN_av = self.areas[eN]  # same value at all quad points
                grad_av += grad_eN_av
            grad_av /= nb_el
            self.areas_nodes[node] = grad_av
        self.model.areas_nodes = self.areas_nodes

    def pseudoTimeStepping(self, eps=0.01):
        logEvent("Pseudo-time stepping with dt={eps} (0<t<1)".format(eps=eps),
                 level=3)
        t_range = np.linspace(0., 1., int(1./eps+1))[1:]
        t_last = 0
        xx = copy.deepcopy(self.mesh.nodeArray)
        for j, t in enumerate(t_range):
            logEvent("Pseudo-time stepping t={t}".format(t=t),
                    level=3)
            for i in range(len(xx)):
                phi = xx[i]
                dt = t-t_last
                flag = self.mesh.nodeMaterialTypes[i]
                bN = None
                if flag != 0 and self.boundaryNormals is not None:
                    bN = self.boundaryNormals[flag]
                    if bN[0] == 0. and bN[1] == 0. and bN[2] == 0:
                        bN = None
                fixed = False
                if self.fixedNodes is not None:
                    if self.fixedNodes[flag] == 1:
                        fixed = True
                if flag == 0 or bN is not None:
                    if fixed is False:
                        if j > 0:
                            eN, nearest_node, xi = self.dist_search(phi, self.nearest_nodes[i], self.eN_phi[i])
                            self.eN_phi[i] = eN
                            self.nearest_nodes[i] = nearest_node
                            if eN is None:
                                print("Element not found for:", i)
                                # print("node", i, self.mesh.nodeArray[i])
                                # print("new coords", phi)
                                v_grad = self.grads[i]
                                area = self.areas_nodes[i]
                            else:
                                v_grad = self.getGradientValue(eN, xi)
                                area = self.areas[eN]
                        else:
                            v_grad = self.grads[i]
                            area = self.areas_nodes[i]
                            eN = 0
                        if self.LS_MODEL is not None and eN is not None:
                            ls_phi = self.getLevelSetValue(eN, phi)
                        else:
                            ls_phi = None
                        dphi = v_grad/(t*1./(self.evaluateFunAtX(x=phi, ls_phi=ls_phi)*self.C)+(1-t)*1./area)
                        # # -----
                        # Euler
                        # # -----
                        if flag == 0:
                            phi[:self.nd] += dphi[:self.nd]*dt
                        elif bN is not None and i > 3:
                            phi[:self.nd] += dphi[:self.nd]*(1-np.abs(bN[:self.nd]))*dt
                        # # -----
                        xx[i] = phi  # not necessary but left for cythonising
            t_last = t
        self.PHI[:] = xx

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

def recoveryAtNodes(variable, grad_v, nodeElementOffsets):
    """
    variable:
         Variable in element
    """
    recovered_variable = np.zeros(len(nodeElementOffsets))
    for node in range(len(nodeElementOffsets)):
        nb_el = 0
        grad_av = 0.
        for eOffset in range(nodeElementOffsets[node],
                             nodeElementOffsets[node+1]):
            nb_el += 1
            eN = self.mesh.nodeElementsArray[eOffset]
            grad_eN_av = 0
            # for k in range(n_quad):
            #     grad_k = self.model.q[('grad(u)',0)][eN, k]
            #     grad_eN_av += grad_k
            # grad_eN_av /= n_quad
            grad_eN_av = self.areas[eN]  # same value at all quad points
            grad_av += grad_eN_av
        grad_av /= nb_el
        recovered_variable[node] = grad_av
    return recovered_variable

