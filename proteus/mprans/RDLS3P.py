import proteus
import numpy
from proteus import *
from proteus.Transport import *
from proteus.Transport import OneLevelTransport
import os
from proteus import cfemIntegrals, Quadrature, Norms, Comm
from proteus.NonlinearSolvers import NonlinearEquation
from proteus.FemTools import (DOFBoundaryConditions,
                              FluxBoundaryConditions,
                              C0_AffineLinearOnSimplexWithNodalBasis)
from proteus.flcbdfWrappers import globalMax
from proteus.Profiling import memory
from proteus.Profiling import logEvent as log
from proteus.Transport import OneLevelTransport
from proteus.TransportCoefficients import TC_base
from proteus.SubgridError import SGE_base
from proteus.ShockCapturing import ShockCapturing_base
import cRDLS3P

class SubgridError(proteus.SubgridError.SGE_base):

    def __init__(self, coefficients, nd):
        proteus.SubgridError.SGE_base.__init__(self, coefficients, nd, False)

    def initializeElementQuadrature(self, mesh, t, cq):
        for ci in range(self.nc):
            cq[('dH_sge', ci, ci)] = cq[('dH', ci, ci)]

    def calculateSubgridError(self, q):
        pass

    def updateSubgridErrorHistory(self, initializationPhase=False):
        pass


class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):

    def __init__(
            self,
            coefficients,
            nd,
            shockCapturingFactor=0.25,
            lag=True,
            nStepsToDelay=None):
        proteus.ShockCapturing.ShockCapturing_base.__init__(
            self, coefficients, nd, shockCapturingFactor, lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0

    def initializeElementQuadrature(self, mesh, t, cq):
        self.mesh = mesh
        self.numDiff = []
        self.numDiff_last = []
        for ci in range(self.nc):
            self.numDiff.append(cq[('numDiff', ci, ci)])
            self.numDiff_last.append(cq[('numDiff', ci, ci)])

    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            self.lag = True
            self.numDiff_last = []
            for ci in range(self.nc):
                self.numDiff_last.append(self.numDiff[ci].copy())
        log("RDLS3P: max numDiff %e" %
            (globalMax(self.numDiff_last[0].max()),))


class PsiTC(proteus.StepControl.SC_base):

    def __init__(self, model, nOptions):
        proteus.StepControl.SC_base.__init__(self, model, nOptions)
        for ci in nOptions.atol_res.keys():
            self.atol = nOptions.atol_res[ci]
            self.rtol = nOptions.rtol_res[ci]
        self.stepExact = True
        for m in model.levelModelList:
            m.timeIntegration.isAdaptive = False
        self.nSteps = 0
        self.nStepsOsher = nOptions.psitc['nStepsForce']
        self.red_ratio = nOptions.psitc['reduceRatio']
        self.start_ratio = nOptions.psitc['startRatio']
        self.nStepsMax = nOptions.psitc['nStepsMax']

    def stepExact_model(self, tOut):
        # pseudo time step
        for m in self.model.levelModelList:
            m.timeIntegration.choose_dt()
        self.dt_model = self.start_ratio * m.timeIntegration.dt
        self.set_dt_allLevels()
        # physical time step
        self.t_model = tOut
        self.setSubsteps([tOut])

    def initialize_dt_model(self, t0, tOut):
        self.saveSolution()
        self.t_model_last = t0
        self.t_model = tOut
        # pseudo time step
        for m in self.model.levelModelList:
            m.timeIntegration.initialize_dt(t0, tOut, m.q)
        # set the starting time steps
        self.dt_model = self.start_ratio * m.timeIntegration.dt
        self.set_dt_allLevels()
        # physical time step
        self.t = tOut
        self.setSubsteps([tOut])
        self.nSteps = 0
        log("Initializing time step on model %s to dt = %12.5e" %
            (self.model.name, self.dt_model), level=1)

    def updateSubstep(self):
        # choose  a new dt and add a substep without increasing t
        # if the steady state has been reached then append the new  t to the
        # substeps
        self.solverFailures = 0
        self.errorFailures = 0
        self.saveSolution()
        # here we just need to test the error and set to tOut if steady state
        if self.nSteps == 0:
            self.res0 = self.model.solver.solverList[-1].norm_r0
            for m in self.model.levelModelList:
                m.timeIntegration.choose_dt()
            self.dt_model = self.start_ratio * \
                self.model.levelModelList[0].timeIntegration.dt
        res = self.model.solver.solverList[-1].norm_r0
        ssError = res / (self.res0 * self.rtol + self.atol)
        for m in self.model.levelModelList:
            m.updateTimeHistory(self.t_model)
            m.timeIntegration.updateTimeHistory()
        if ((self.nSteps < self.nStepsOsher or
             ssError >= 1.0) and
                self.nSteps < self.nStepsMax):
            self.nSteps += 1
            self.dt_model = m.timeIntegration.dt
            if self.nSteps >= self.nStepsOsher:  # start ramping up the time step
                self.dt_model = self.dt_model * self.red_ratio
            #log("Osher-PsiTC dt %12.5e" %(self.dt_model),level=1)
            self.set_dt_allLevels()
            # physical time step
            self.t_model = self.substeps[0]
            self.substeps.append(self.substeps[0])
            log("Osher-PsiTC iteration %d  dt = %12.5e  |res| = %12.5e %g  " %
                (self.nSteps, self.dt_model, res, (res / self.res0) * 100.0), level=1)
        elif self.nSteps >= self.nStepsMax:
            log("Osher-PsiTC DID NOT Converge |res| = %12.5e but quitting anyway" % (res,))
            log("Osher-PsiTC tolerance                %12.5e " %
                (self.res0 * self.rtol + self.atol,))
            self.nSteps = 0
        else:
            log("Osher-PsiTC converged |res| = %12.5e %12.5e" %
                (res, ssError * 100.0))
            log("Osher-PsiTC tolerance                %12.5e " %
                (self.res0 * self.rtol + self.atol,))
            self.nSteps = 0

    def choose_dt_model(self):
        # don't modify dt_model
        self.solverFailures = 0
        self.errorFailures = 0
        self.saveSolution()
        # pseudo time step
        for m in self.model.levelModelList:
            m.timeIntegration.choose_dt()
        self.dt_model = self.start_ratio * m.timeIntegration.dt
        self.set_dt_allLevels()
        # physical time step
        self.t_model = self.substeps[0]
        self.setSubsteps([self.substeps[0]])
        log("Osher-PsiTC choosing dt = %12.5e " % (self.dt_model,))


class Coefficients(proteus.TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import redistanceLevelSetCoefficientsEvaluate

    def __init__(
            self,
            applyRedistancing=True,
            epsFact=2.0,
            nModelId=None,
            u0=None,
            rdModelId=0,
            penaltyParameter=0.0,
            useMetrics=0.0,
            useConstantH=False,
            weakDirichletFactor=10.0,
            backgroundDiffusionFactor=0.01):
        self.useConstantH = useConstantH
        self.useMetrics = useMetrics
        variableNames = ['phid']
        nc = 1
        mass = {0: {0: 'linear'}}
        hamiltonian = {0: {0: 'nonlinear'}}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {0: {0: 'constant'}}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         variableNames)
        self.nModelId = nModelId
        self.rdModelId = rdModelId
        self.epsFact = epsFact
        self.q_u0 = None
        self.ebq_u0 = None
        self.ebqe_u0 = None
        self.dof_u0 = None
        self.u0 = u0
        self.applyRedistancing = applyRedistancing
        self.weakBC_on = True  # False
        self.penaltyParameter = penaltyParameter
        self.backgroundDiffusionFactor = backgroundDiffusionFactor
        self.weakDirichletFactor = weakDirichletFactor

    def attachModels(self, modelList):
        if self.nModelId is not None:
            self.nModel = modelList[self.nModelId]
            self.q_u0 = self.nModel.q[('u', 0)]
            if self.nModel.ebq.has_key(('u', 0)):
                self.ebq_u0 = self.nModel.ebq[('u', 0)]
            self.ebqe_u0 = self.nModel.ebqe[('u', 0)]
            self.dof_u0 = self.nModel.u[0].dof
        else:
            self.nModel = None
        self.rdModel = modelList[self.rdModelId]

    def initializeMesh(self, mesh):
        self.h = mesh.h
        self.eps = self.epsFact * mesh.h

    def initializeElementQuadrature(self, t, cq):
        if self.nModelId is None:
            if self.q_u0 is None:
                self.q_u0 = numpy.zeros(cq[('u', 0)].shape, 'd')
            if self.u0 is not None:
                for i in range(len(cq[('u', 0)].flat)):
                    self.q_u0.flat[i] = self.u0.uOfXT(
                        cq['x'].flat[3 * i:3 * (i + 1)], 0.)

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        if self.nModelId is None:
            if self.ebq_u0 is None:
                self.ebq_u0 = numpy.zeros(cebq[('u', 0)].shape, 'd')
            if self.u0 is not None:
                for i in range(len(cebq[('u', 0)].flat)):
                    self.ebq_u0.flat[i] = self.u0.uOfXT(
                        cebq['x'].flat[3 * i:3 * (i + 1)], 0.)

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        if self.nModelId is None:
            if self.ebqe_u0 is None:
                self.ebqe_u0 = numpy.zeros(cebqe[('u', 0)].shape, 'd')
            if self.u0 is not None:
                for i in range(len(cebqe[('u', 0)].flat)):
                    self.ebqe_u0.flat[i] = self.u0.uOfXT(
                        cebqe['x'].flat[3 * i:3 * (i + 1)], 0.)

    def preStep(self, t, firstStep=False):
        import pdb
        # pdb.set_trace()
        if self.nModel is not None:
            log("resetting signed distance level set to current level set", level=2)
            self.rdModel.u[0].dof[:] = self.nModel.u[0].dof[:]
            self.rdModel.calculateCoefficients()
            self.rdModel.calculateElementResidual()
            self.rdModel.timeIntegration.updateTimeHistory(resetFromDOF=True)
            self.rdModel.timeIntegration.resetTimeHistory(resetFromDOF=True)
            self.rdModel.updateTimeHistory(t, resetFromDOF=True)
            # now do again because of subgrid error lagging
            #\todo modify subgrid error lagging so this won't be necessary
            self.rdModel.calculateCoefficients()
            self.rdModel.calculateElementResidual()
            self.rdModel.timeIntegration.updateTimeHistory(resetFromDOF=True)
            self.rdModel.timeIntegration.resetTimeHistory(resetFromDOF=True)
            self.rdModel.updateTimeHistory(t, resetFromDOF=True)
            copyInstructions = {'copy_uList': True,
                                'uList_model': self.nModelId}
            copyInstructions = {'reset_uList': True}
            return copyInstructions
        else:
            return {}

    def postStep(self, t, firstStep=False):
        if self.nModel is not None:
            if self.applyRedistancing:
                log("resetting level set to signed distance")
                self.nModel.u[0].dof.flat[:] = self.rdModel.u[0].dof.flat[:]
                self.nModel.calculateCoefficients()
                self.nModel.calculateElementResidual()
                # save the boundary level set in the numerical flux to use for
                self.nModel.numericalFlux.ebqe[
                    ('u', 0)][:] = self.rdModel.ebqe[
                    ('u', 0)]
            copyInstructions = {}
            return copyInstructions
        else:
            return {}

    def getICDofs(self, cj):
        return self.dof_u0

    def updateToMovingDomain(self, t, c):
        # the redistancing equations is not physical so it just needs the
        # updated mesh
        pass

    def evaluate(self, t, c):
        if c[('H', 0)].shape == self.q_u0.shape:
            u0 = self.q_u0
            # print "numdiff",c[('numDiff',0,0)].min(),c[('numDiff',0,0)].max()
            # print
            # "tau",self.rdModel.stabilization.tau[0].min(),self.rdModel.stabilization.tau[0].max(),
        elif c[('H', 0)].shape == self.ebqe_u0.shape:
            u0 = self.ebqe_u0
        else:
            u0 = self.ebq_u0
        assert u0 is not None
        # \todo make redistancing epsilon depend on local element diamater instead of global max
        self.redistanceLevelSetCoefficientsEvaluate(self.eps,
                                                    u0,
                                                    c[('u', 0)],
                                                    c[('grad(u)', 0)],
                                                    c[('m', 0)],
                                                    c[('dm', 0, 0)],
                                                    c[('H', 0)],
                                                    c[('dH', 0, 0)],
                                                    c[('r', 0)])
    # weak Dirichlet conditions on level set (boundary conditions of Eikonal equation)
    # \todo clean up weak Dirichlet conditions for Eikonal equation in transport coefficents

    def setZeroLSweakDirichletBCs(vt):
        if vt.coefficients.weakBC_on:
            # print
            # "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=========================Setting
            # new weak
            # BC's=======================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            vt.dirichletNodeSetList[0] = []
            vt.dirichletGlobalNodeSet[0] = set()
            vt.dirichletValues[0] = {}
            for eN in range(vt.mesh.nElements_global):
                vt.dirichletNodeSetList[0].append(set())
                signU = 0
                j0 = 0
                eps = vt.u[0].femSpace.mesh.elementDiametersArray[eN]
                # loop over nodes looking for a node not within eps of zero
                while ((signU == 0) and
                       (j0 < vt.nDOF_trial_element[0])):
                    J0 = vt.u[0].femSpace.dofMap.l2g[eN, j0]
                    if vt.u[0].dof[J0] < -eps:
                        signU = -1
                    elif vt.u[0].dof[J0] > eps:
                        signU = 1
                    else:  # freeze this node within eps of zero
                        vt.dirichletNodeSetList[0][eN].add(j0)
                        vt.dirichletValues[0][(eN, j0)] = vt.u[0].dof[J0]
                        vt.dirichletGlobalNodeSet[0].add(J0)
                    j0 += 1
                # loop over remaining nodes to see if the zero level set cuts
                # element
                for j in range(j0, vt.nDOF_trial_element[0]):
                    J = vt.u[0].femSpace.dofMap.l2g[eN, j]
                    if (((vt.u[0].dof[J] < -eps) and
                         (signU == 1)) or
                        ((vt.u[0].dof[J] > eps) and
                         (signU == -1))):  # level set cuts element, freeze whole element
                        for jj in range(vt.nDOF_trial_element[0]):
                            JJ = vt.u[0].femSpace.dofMap.l2g[eN, jj]
                            vt.dirichletNodeSetList[0][eN].add(jj)
                            vt.dirichletValues[0][
                                (eN, jj)] = float(vt.u[0].dof[JJ])
                            vt.dirichletGlobalNodeSet[0].add(JJ)
                        break
                    # freeze this node within eps of zero
                    elif (fabs(vt.u[0].dof[J]) < eps):
                        vt.dirichletNodeSetList[0][eN].add(j)
                        vt.dirichletValues[0][(eN, j)] = float(vt.u[0].dof[J])
                        vt.dirichletGlobalNodeSet[0].add(J)
            # get all frozen dof and make sure they're frozen on each element
            for eN in range(vt.mesh.nElements_global):
                for j in range(vt.nDOF_trial_element[0]):
                    J = vt.u[0].femSpace.dofMap.l2g[eN, j]
                    if J in vt.dirichletGlobalNodeSet[0]:
                        vt.dirichletNodeSetList[0][eN].add(j)
                        vt.dirichletValues[0][(eN, j)] = float(vt.u[0].dof[J])
        else:
            # print
            # "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=========================Unsetting
            # weak
            # BC's=======================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            vt.dirichletNodeSetList[0] = []
            vt.dirichletGlobalNodeSet[0] = set()
            vt.dirichletValues[0] = {}
            for eN in range(vt.mesh.nElements_global):
                vt.dirichletNodeSetList[0].append(set())

    def setZeroLSweakDirichletBCs2(vt):
        # just look for cut edges and nodes
        vt.dirichletNodeSetList[0] = []
        vt.dirichletGlobalNodeSet[0] = set()
        vt.dirichletValues[0] = {}
        for eN in range(vt.mesh.nElements_global):
            vt.dirichletNodeSetList[0].append(set())
            signU = 0
            j0 = 0
            eps = .1 * vt.u[0].femSpace.mesh.elementDiametersArray[eN]
            eps = 0.0
            for j0 in range(vt.nDOF_trial_element[0]):
                J0 = vt.u[0].femSpace.dofMap.l2g[eN, j0]
                if vt.u[0].dof[J0] < -eps:
                    signU = -1
                elif vt.u[0].dof[J0] > eps:
                    signU = 1
                else:
                    vt.dirichletNodeSetList[0][eN].add(j0)
                    vt.dirichletValues[0][(eN, j0)] = float(vt.u[0].dof[J0])
                    vt.dirichletGlobalNodeSet[0].add(J0)
                if signU != 0:
                    for j in (
                        range(
                            0,
                            j0) +
                        range(
                            j0 +
                            1,
                            vt.nDOF_trial_element[0])):
                        J = vt.u[0].femSpace.dofMap.l2g[eN, j]
                        if (((vt.u[0].dof[J] < -eps) and
                             (signU == 1)) or
                            ((vt.u[0].dof[J] > eps) and
                             (signU == -1))):
                            vt.dirichletNodeSetList[0][eN].add(j)
                            vt.dirichletValues[0][
                                (eN, j)] = float(vt.u[0].dof[J])
                            vt.dirichletGlobalNodeSet[0].add(J)
                            vt.dirichletNodeSetList[0][eN].add(j0)
                            vt.dirichletValues[0][
                                (eN, j0)] = float(vt.u[0].dof[j0])
                            vt.dirichletGlobalNodeSet[0].add(j0)
        for eN in range(vt.mesh.nElements_global):
            for j in range(vt.nDOF_trial_element[0]):
                J = vt.u[0].femSpace.dofMap.l2g[eN, j]
                if J in vt.dirichletGlobalNodeSet[0]:
                    vt.dirichletNodeSetList[0][eN].add(j)
                    vt.dirichletValues[0][(eN, j)] = float(vt.u[0].dof[J])

    def setZeroLSweakDirichletBCs3(vt):
        if vt.coefficients.weakBC_on:
            # print
            # "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=========================Setting
            # new weak
            # BC's=======================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            vt.dirichletNodeSetList[0] = []
            vt.dirichletGlobalNodeSet[0] = set()
            vt.dirichletValues[0] = {}
            for eN in range(vt.mesh.nElements_global):
                vt.dirichletNodeSetList[0].append(set())
                eps = vt.coefficients.epsFact * \
                    vt.u[0].femSpace.mesh.elementDiametersArray[eN]
                #eps = 1.5*vt.u[0].femSpace.mesh.elementDiametersArray[eN]
                # loop over nodes looking for a node not within eps of zero
                for j in range(vt.nDOF_trial_element[0]):
                    J = vt.u[0].femSpace.dofMap.l2g[eN, j]
                    if (fabs(vt.u[0].dof[J]) <
                            eps):  # freeze this node within eps of zero
                        vt.dirichletNodeSetList[0][eN].add(j)
                        vt.dirichletValues[0][(eN, j)] = float(vt.u[0].dof[J])
                        vt.dirichletGlobalNodeSet[0].add(J)
        else:
            # print
            # "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=========================Unsetting
            # weak
            # BC's=======================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            vt.dirichletNodeSetList[0] = []
            vt.dirichletGlobalNodeSet[0] = set()
            vt.dirichletValues[0] = {}
            for eN in range(vt.mesh.nElements_global):
                vt.dirichletNodeSetList[0].append(set())

    # def
    setZeroLSweakDirichletBCs = staticmethod(setZeroLSweakDirichletBCs)
    setZeroLSweakDirichletBCs2 = staticmethod(setZeroLSweakDirichletBCs2)
    setZeroLSweakDirichletBCs3 = staticmethod(setZeroLSweakDirichletBCs3)

debugRDLS3P = False  # True


class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0

    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=True,
                 sd=True,
                 movingDomain=False):
        from proteus import Comm
        #
        # set the objects describing the method and boundary conditions
        #
        self.movingDomain = movingDomain
        self.tLast_mesh = None
        #
        self.name = name
        self.sd = sd
        self.Hess = False
        self.lowmem = True
        self.timeTerm = True  # allow turning off  the  time derivative
        # self.lowmem=False
        self.testIsTrial = True
        self.phiTrialIsTrial = True
        self.u = uDict
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # mwf try to reuse test and trial information across components if
        # spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature  # True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1, coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[
                    0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        # Simplicial Mesh
        # assume the same mesh for  all components for now
        self.mesh = self.u[0].femSpace.mesh
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.dirichletNodeSetList = None
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        # no velocity post-processing for now
        self.conservativeFlux = conservativeFluxDict
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.advection.has_key(ci):
                    for flag in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for flag in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if coefficients.potential.has_key(ci):
                    for flag in coefficients.potential[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = (
                (self.conservativeFlux is not None) or (
                    numericalFluxType is not None) or (
                    self.fluxBoundaryConditions[ci] == 'outFlow') or (
                    self.fluxBoundaryConditions[ci] == 'mixedFlow') or (
                    self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        # assume same space dim for all variables
        self.nSpace_global = self.u[0].femSpace.nSpace_global
        self.nDOF_trial_element = [
            u_j.femSpace.max_nDOF_element for u_j in self.u.values()]
        self.nDOF_phi_trial_element = [
            phi_k.femSpace.max_nDOF_element for phi_k in self.phi.values()]
        self.n_phi_ip_element = [
            phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in self.phi.values()]
        self.nDOF_test_element = [
            femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global = [
            dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self, self.nFreeVDOF_global)
        #
        # build the quadrature point dictionaries from the input (this
        # is just for convenience so that the input doesn't have to be
        # complete)
        #
        elementQuadratureDict = {}
        elemQuadIsDict = isinstance(elementQuadrature, dict)
        if elemQuadIsDict:  # set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[
                        ('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff', ci, ci)):
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            ('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            'default']
                else:
                    elementQuadratureDict[
                        ('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        # mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints, self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * \
            self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints, self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[
            0]
        self.nElementBoundaryQuadraturePoints_global = (
            self.mesh.nElements_global *
            self.mesh.nElementBoundaries_element *
            self.nElementBoundaryQuadraturePoints_elementBoundary)
#        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
#            print self.nQuadraturePoints_element
#            if self.nSpace_global == 3:
#                assert(self.nQuadraturePoints_element == 5)
#            elif self.nSpace_global == 2:
#                assert(self.nQuadraturePoints_element == 6)
#            elif self.nSpace_global == 1:
#                assert(self.nQuadraturePoints_element == 3)
#
#            print self.nElementBoundaryQuadraturePoints_elementBoundary
#            if self.nSpace_global == 3:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 2:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 1:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        #self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.q[('u', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[('m_last', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        # for time integration by VBDF and probably FLCBDF
        self.q[('m', 0)] = self.q[('u', 0)]
        # needed by PsiTCtte
        self.q[('mt', 0)] = numpy.zeros(self.q[('u', 0)].shape, 'd')
        self.q[('dH', 0, 0)] = numpy.zeros((self.mesh.nElements_global,
                                            self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[
            ('dH_sge',
             0,
             0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[('cfl', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 0, 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[
            ('u',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set(
            [('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global, i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray = numpy.zeros((self.nNodes_internal,), 'i')
        for nI, n in enumerate(self.internalNodes):
            self.internalNodesArray[nI] = n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings", 2)
        self.updateLocal2Global()
        log("Building time integration object", 2)
        log(memory("inflowBC, internalNodes,updateLocal2Global",
                   "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(
                self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration", "OneLevelTransport"), level=4)
        log("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        log(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict,
                    options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        # cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = self.numericalFlux.penalty_constant / (
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        # penalty term
        # cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = self.numericalFlux.penalty_constant / \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        # use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(
            self)
        log(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        # TODO get rid of this
        for ci, fbcObject in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(
                self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[
                        ('advectiveFlux_bc', ci)][
                        t[0], t[1]] = g(
                        self.ebqe[
                            ('x')][
                            t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
        # reduced quad
        if hasattr(self.numericalFlux, 'setDirichletValues'):
            self.numericalFlux.setDirichletValues(self.ebqe)
        if not hasattr(self.numericalFlux, 'isDOFBoundary'):
            self.numericalFlux.isDOFBoundary = {
                0: numpy.zeros(self.ebqe[('u', 0)].shape, 'i')}
        if not hasattr(self.numericalFlux, 'ebqe'):
            self.numericalFlux.ebqe = {
                ('u', 0): numpy.zeros(self.ebqe[('u', 0)].shape, 'd')}
        # TODO how to handle redistancing calls for
        # calculateCoefficients,calculateElementResidual etc
        self.globalResidualDummy = None
        # create storage for enforcing weak dirichlet boundary conditions
        # around zero level set contour
        self.freezeLevelSet = 1  # True
        self.u_dof_last = numpy.zeros(self.u[0].dof.shape, 'd')
        self.weakDirichletConditionFlags = numpy.zeros(
            self.u[0].dof.shape, 'i')
        self.dofFlag_element = numpy.zeros((self.nDOF_trial_element[0],), 'i')
        # allow Newton solves for redistancing
        if self.timeIntegration.__class__ == TimeIntegration.NoIntegration:
            self.timeIntegration.m_tmp = {
                0: numpy.zeros(self.q[('m_tmp', 0)].shape, 'd')}
            self.timeIntegration.m_last = {
                0: numpy.zeros(self.q[('m_tmp', 0)].shape, 'd')}
        compKernelFlag = 0
        if self.coefficients.useConstantH:
            self.elementDiameter = self.mesh.elementDiametersArray.copy()
            self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        else:
            self.elementDiameter = self.mesh.elementDiametersArray
        self.rdls3p = cRDLS3P.RDLS3P(
            self.nSpace_global,
            self.nQuadraturePoints_element,
            self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
            self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
            self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
            self.nElementBoundaryQuadraturePoints_elementBoundary,
            compKernelFlag)

    def calculateCoefficients(self):
        pass

    def calculateElementResidual(self):
        pass

    def getResidual(self, u, r):
        import pdb
        import copy
        # try to use 1d,2d,3d specific modules
        # mwf debug
        # pdb.set_trace()
        r.fill(0.0)
        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        # cek can put in logic to skip of BC's don't depend on t or u
        # Dirichlet boundary conditions
        if hasattr(self.numericalFlux, 'setDirichletValues'):
            self.numericalFlux.setDirichletValues(self.ebqe)
        # flux boundary conditions, SHOULDN'T HAVE
        # for  now force time integration
        useTimeIntegration = 1
        #assert self.timeIntegration.__class__ != TimeIntegration.NoIntegration
        if self.timeIntegration.__class__ == TimeIntegration.NoIntegration or not self.timeTerm:
            useTimeIntegration = 0

        if useTimeIntegration:
            alpha_bdf = self.timeIntegration.alpha_bdf
            beta_bdf = self.timeIntegration.beta_bdf
        else:
            alpha_bdf = self.timeIntegration.dt
            beta_bdf = self.timeIntegration.m_last
        # self.elementResidual[0].fill(0.0)
        # RDLS3P.calculateResidual(self.mesh.nElements_global,
        # print "beta_bdf",beta_bdf
        self.rdls3p.calculateResidual(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            # physics
            self.mesh.nElements_global,
            self.coefficients.useMetrics,
            alpha_bdf,
            self.coefficients.epsFact,
            self.coefficients.backgroundDiffusionFactor,
            self.coefficients.weakDirichletFactor,
            self.freezeLevelSet,
            useTimeIntegration,
            self.shockCapturing.lag,
            self.stabilization.lag,  # 0 nothing lagged
            # 1 dH lagged in tau calc
            # 2 dH lagged in tau and pdeResidual, Lstar*w calculations
            self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,  # self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u0,
            self.timeIntegration.m_tmp[0],
            self.q[('u', 0)],
            self.q[('grad(u)', 0)],
            self.q[('dH', 0, 0)],
            self.u_dof_last,
            beta_bdf[0],
            self.q[('dH_sge', 0, 0)],
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.weakDirichletConditionFlags,
            self.offset[0], self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_u0,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('u', 0)],
            self.ebqe[('grad(u)', 0)])
        # print "m_tmp",self.timeIntegration.m_tmp[0]
        # print "dH",self.q[('dH',0,0)]
        # print "dH_sge",self.q[('dH_sge',0,0)]
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual", level=9, data=r)
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape, 'd')

    def getJacobian(self, jacobian):
        if not debugRDLS3P:
            return self.getJacobianNew(jacobian)
        OneLevelTransport.getJacobian(self, jacobian)

    def getJacobianNew(self, jacobian):
        #import superluWrappers
        #import numpy
        import pdb
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        # for  now force time integration
        useTimeIntegration = 1
        if self.timeIntegration.__class__ == TimeIntegration.NoIntegration or not self.timeTerm:
            useTimeIntegration = 0
        if useTimeIntegration:
            alpha_bdf = self.timeIntegration.alpha_bdf
            beta_bdf = self.timeIntegration.beta_bdf
        else:
            alpha_bdf = self.timeIntegration.dt
            beta_bdf = self.timeIntegration.m_last
        self.rdls3p.calculateJacobian(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
            self.coefficients.useMetrics,
            alpha_bdf,
            self.coefficients.epsFact,
            self.coefficients.backgroundDiffusionFactor,
            self.coefficients.weakDirichletFactor,
            self.freezeLevelSet,
            useTimeIntegration,
            self.shockCapturing.lag,
            self.stabilization.lag,
            self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,  # self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.u_dof_last,
            self.coefficients.q_u0,
            beta_bdf[0],
            self.q[('dH_sge', 0, 0)],
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.weakDirichletConditionFlags,
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_u0,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.csrColumnOffsets_eb[(0, 0)])
        log("Jacobian ", level=10, data=jacobian)
        # mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    # jacobian

    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        # self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
        #                                          self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(
            self.timeIntegration.t, self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self):
        pass

    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(
            self.elementBoundaryQuadraturePoints, self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                   self.ebqe[('x')],
                                                                                   self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                   self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(
            self.timeIntegration.t, self.ebqe)

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass
        # self.q[('u',0)].flat[:]=0.0
        # for eN in range(self.u[0].femSpace.elementMaps.mesh.nElements_global):
        #     for k in range(self.nQuadraturePoints_element):
        #         for j in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
        #             J = self.u[0].femSpace.dofMap.l2g[eN,j]
        #             self.q[('u',0)][eN,k]+=self.u[0].dof[J]*self.u[0].femSpace.psi[k,j]

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def updateAfterMeshMotion(self):
        pass
# OneLevelRDLS3P
from proteus import ctransportCoefficients


def setZeroLSweakDirichletBCs(RDLS3Pvt):
    assert hasattr(RDLS3Pvt, 'freezeLevelSet')
    assert hasattr(RDLS3Pvt, 'u_dof_last')
    assert hasattr(RDLS3Pvt, 'weakDirichletConditionFlags')
    assert hasattr(RDLS3Pvt.coefficients, 'epsFact')
    assert hasattr(RDLS3Pvt, 'dofFlag_element')
    RDLS3Pvt.freezeLevelSet = 1
    useC = True  # True
    if debugRDLS3P:
        useC = False
    RDLS3Pvt.u_dof_last[:] = RDLS3Pvt.u[0].dof
    RDLS3Pvt.weakDirichletConditionFlags.fill(0)
    # to debug use original implementation and then copy over
    if not useC:
        RDLS3Pvt.coefficients.setZeroLSweakDirichletBCs(RDLS3Pvt)
        for eN in range(RDLS3Pvt.mesh.nElements_global):
            for j in range(RDLS3Pvt.nDOF_trial_element[0]):
                J = RDLS3Pvt.u[0].femSpace.dofMap.l2g[eN, j]
                if J in RDLS3Pvt.dirichletGlobalNodeSet[0]:
                    RDLS3Pvt.weakDirichletConditionFlags[J] = 1
    else:
        #
        # mwf debug
        #import pdb
        # pdb.set_trace()
        # use c directly
        ctransportCoefficients.setWeakDirichletConditionsForLevelSet(RDLS3Pvt.mesh.nElements_global,
                                                                     RDLS3Pvt.nDOF_trial_element[
                                                                         0],
                                                                     RDLS3Pvt.coefficients.epsFact,
                                                                     RDLS3Pvt.elementDiameter,  # RDLS3Pvt.mesh.elementDiametersArray,
                                                                     RDLS3Pvt.u[
                                                                         0].femSpace.dofMap.l2g,
                                                                     RDLS3Pvt.u[
                                                                         0].dof,
                                                                     RDLS3Pvt.dofFlag_element,  # temporary storage
                                                                     RDLS3Pvt.weakDirichletConditionFlags)


def setZeroLSweakDirichletBCsSimple(RDLS3Pvt):
    assert hasattr(RDLS3Pvt, 'freezeLevelSet')
    assert hasattr(RDLS3Pvt, 'u_dof_last')
    assert hasattr(RDLS3Pvt, 'weakDirichletConditionFlags')
    assert hasattr(RDLS3Pvt.coefficients, 'epsFact')
    assert hasattr(RDLS3Pvt, 'dofFlag_element')
    RDLS3Pvt.freezeLevelSet = 1
    useC = True  # True
    if debugRDLS3P:
        useC = False
    RDLS3Pvt.u_dof_last[:] = RDLS3Pvt.u[0].dof
    RDLS3Pvt.weakDirichletConditionFlags.fill(0)
    # to debug use original implementation and then copy over
    if not useC:
        RDLS3Pvt.coefficients.setZeroLSweakDirichletBCs(RDLS3Pvt)
        for eN in range(RDLS3Pvt.mesh.nElements_global):
            for j in range(RDLS3Pvt.nDOF_trial_element[0]):
                J = RDLS3Pvt.u[0].femSpace.dofMap.l2g[eN, j]
                if J in RDLS3Pvt.dirichletGlobalNodeSet[0]:
                    RDLS3Pvt.weakDirichletConditionFlags[J] = 1
    else:
        #
        # mwf debug
        #import pdb
        # pdb.set_trace()
        # use c directly
        ctransportCoefficients.setSimpleWeakDirichletConditionsForLevelSet(RDLS3Pvt.mesh.nElements_global,
                                                                           RDLS3Pvt.nDOF_trial_element[
                                                                               0],
                                                                           RDLS3Pvt.coefficients.epsFact,
                                                                           RDLS3Pvt.elementDiameter,  # RDLS3Pvt.mesh.elementDiametersArray,
                                                                           RDLS3Pvt.u[
                                                                               0].femSpace.dofMap.l2g,
                                                                           RDLS3Pvt.u[
                                                                               0].dof,
                                                                           RDLS3Pvt.dofFlag_element,  # temporary storage
                                                                           RDLS3Pvt.weakDirichletConditionFlags)
