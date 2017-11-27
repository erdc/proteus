import numpy as np
import copy
from proteus.Profiling import logEvent
from proteus import TransportCoefficients

class Coefficients(TransportCoefficients.PoissonEquationCoefficients):
    """
    Coefficients for deforming the mesh according to an area function

    Parameters
    ----------
    func: function
        function of x defining element area
    """
    def __init__(self, func, nd=2):
        self.myfunc=func
        self.C = 1.  # scaling coefficient for f (computed in preStep)
        self.integral_area = 1.
        self.nd = nd
        aOfX = [lambda x: np.array([[1., 0.], [0., 1.]])]
        fOfX = [self.uOfX]  # scaled function reciprocal
        #super(MyCoeff, self).__init__(aOfX, fOfX)
        TransportCoefficients.PoissonEquationCoefficients.__init__(self, aOfX, fOfX)
            
    def initializeMesh(self, mesh):
        self.mesh = mesh
        self.areas_nodes = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.areas = np.zeros(self.mesh.elementNodesArray.shape[0])
        
    def attachModels(self,modelList):
        self.model = modelList[-1]
        self.q = self.model.q
        self.model.mesh_array0 = copy.deepcopy(self.mesh.nodeArray)
        self.PHI = np.zeros_like(self.model.mesh.nodeArray) 
        self.grads = np.zeros((len(self.model.u[0].dof),2))
        self.areas_nodes = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.areas = np.zeros(self.mesh.elementNodesArray.shape[0])
        self.nearest_nodes = np.array([i for i in range(len(self.mesh.nodeArray))])
        self.nearest_nodes0 = np.array([i for i in range(len(self.mesh.nodeArray))])
        self.mesh.nodeVelocityArray = np.zeros(self.mesh.nodeArray.shape,'d')

    def preStep(self, t, firstStep=False):
        logEvent("Updating area function scaling coefficient",
                 level=3)
        integral_1_over_f = 0.
        xx = self.model.q['x']
        JJ = self.model.q['abs(det(J))']
        N_k = self.model.q['x'].shape[1]
        q_weights = self.model.elementQuadratureWeights[('u', 0)]
        nE = 0
        for e in xrange(self.mesh.elementNodesArray.shape[0]):
            area = 0
            for k in xrange(N_k):
                area += JJ[e, k]*q_weights[k]
                integral_1_over_f += JJ[e, k]*q_weights[k]/self.myfunc(xx[e,k])
            self.areas[e] = area
            nE += 1
        self.C = integral_1_over_f/(nE)  # update scaling coefficient
        self.model.areas = self.areas
        logEvent("Element Quadrature - updating uOfX",
                 level=3)
        ci = 0
        cq = self.model.q
        for ci in range(self.nc):
            cq[('f',ci)].flat[:] = 0.0
            for eN in range(len(cq[('r',ci)])):
                for k in range(len(cq[('r', ci)][eN])):
                    cq[('r',ci)][eN,k] = -self.uOfX(x=cq['x'][eN,k], eN=eN)
        logEvent("Finished updating area function",
                 level=3)

    def postStep(self, t,firstStep=False):
        # gradient recovery
        self.gradientRecoveryAtNodes()
        # area recovery
        self.areaRecoveryAtNodes()
        # pseudo-time step
        eps = 0.1
        self.pseudoTimeStepping(eps)
        # move nodes
        self.model.mesh.nodeVelocityArray[:] = (self.PHI-self.model.mesh.nodeArray)/self.model.timeIntegration.dt
        self.model.mesh.nodeArray[:] = self.PHI
        # re-initialise nearest nodes
        self.nearest_nodes[:] = self.nearest_nodes0[:]

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
        f = 1./(self.myfunc(x)*self.C)-1./area
        return f
    
    def dist_search(self, x, node):
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
        while found_node is False:
            nearest_node0 = nearest_node
            if node+1 < len(femSpace.mesh.nodeElementOffsets):
                for eOffset in range(femSpace.mesh.nodeElementOffsets[node],
                                     femSpace.mesh.nodeElementOffsets[node+1]):
                    eN = self.mesh.nodeElementsArray[eOffset]
                    checkedElements.append(eN)
                    for node in self.mesh.elementNodesArray[eN]:
                        node_coords = self.mesh.nodeArray[node]
                        dist = np.sqrt(np.sum((node_coords-x)**2))
                        if dist < min_dist:
                            min_dist = dist
                            nearest_node = node
            else:
                assert False, 'problem to find element'
            if nearest_node0 == nearest_node:
                found_node = True
        # now find element containing x around nearest_node 
        element = None
        coords = self.mesh.nodeArray[nearest_node]
        for eOffset in range(self.mesh.nodeElementOffsets[nearest_node],
                             self.mesh.nodeElementOffsets[nearest_node+1]):
            eN = self.mesh.nodeElementsArray[eOffset]
            checkedElements.append(eN)
            # union of set
            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            # evaluate the inverse map for element eN (global to local)
            xi = femSpace.elementMaps.getInverseValue(eN, coords)
            #J = femSpace.elementMaps.getJacobianValues(eN, )
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                element = eN
                #print("found1", eN)
        if element is None:
            # extra loop if case coords is in neighbour element
            for node in patchBoundaryNodes:
                for eOffset in range(self.mesh.nodeElementOffsets[node],
                                     self.mesh.nodeElementOffsets[node+1]):
                    eN = self.mesh.nodeElementsArray[eOffset]
                    if eN not in checkedElements:
                        checkedElements.append(eN)
                        # evaluate the inverse map for element eN
                        xi = femSpace.elementMaps.getInverseValue(eN, coords)
                        # query whether xi lies within the reference element
                        if femSpace.elementMaps.referenceElement.onElement(xi):
                            element = eN
                            #print("found2", eN)
        return eN, nearest_node, xi

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
                if flag == 0:
                    if j > 0:
                        eN, nearest_node, xi = self.dist_search(phi, self.nearest_nodes[i])
                        v_grad = self.getGradientValue(eN, xi)
                        area = self.areas[eN]
                        self.nearest_nodes[i] = nearest_node
                    else:
                        v_grad = self.grads[i]
                        area = self.areas_nodes[i]
                    dphi = v_grad/(t*1./(self.myfunc(phi)*self.C)+(1-t)*1./area)
                    # # RK4
                    # # -----
                    # phi1 = copy.copy(phi)
                    # phi2 = copy.copy(phi1)
                    # phi2[0] += dphi[0]*dt/2.
                    # phi2[1] += dphi[1]*dt/2.
                    # v_grad2, node2 = self.dist_search(phi2, i)
                    # dphi2 = v_grad2/(t*(self.myfunc(phi2)*self.C)+(1-t)*1)
                    # phi3 = copy.copy(phi2)
                    # phi3[0] += dphi2[0]*dt/2.
                    # phi3[1] += dphi2[1]*dt/2.
                    # v_grad3, node3 = self.dist_search(phi3, i)
                    # dphi3 = v_grad3/(t*(self.myfunc(phi3)*self.C)+(1-t)*1)
                    # phi4 = copy.copy(phi3)
                    # phi4[0] += dphi3[0]*dt/2.
                    # phi4[1] += dphi3[1]*dt/2.
                    # v_grad4, node4 = self.dist_search(phi4, i)
                    # dphi4 = v_grad4/(t*(self.myfunc(phi4)*self.C)+(1-t)*1)
                    # phii = phi1[:2]+(dt/6.)*(dphi+2*dphi2+2*dphi3+dphi4)
                    # phi[0] = phii[0]
                    # phi[1] = phii[1]
                    # # -----
                    # Euler
                    # # -----
                    phi[:self.nd] += dphi[:self.nd]*dt
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
