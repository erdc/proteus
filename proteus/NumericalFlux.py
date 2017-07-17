"""
A class hierarchy for numerical flux (numerical trace) computations

.. inheritance-diagram:: proteus.NumericalFlux
   :parts: 1
"""
import cfemIntegrals,cnumericalFlux
import numpy
from .Profiling import logEvent,memory
from EGeometry import enorm
#TODO:
#  allow different flags for Dirichlet bc's, keep isDOFBoundary = 1 by default
#  ? allow flux bc's to have type flags too?
#
#
class NF_base:
    useWeakDirichletConditions=True
    useStrongDirichletConstraints=False
    hasInterior=True
    def __init__(self,vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 parallelPeriodic=False):
        import copy
        self.hasInterior=NF_base.hasInterior
        self.includeBoundaryAdjoint=False
        self.includeBoundaryAdjointInteriorOnly=False # tjp added
        self.boundaryAdjoint_sigma=0.0
        self.penalty_constant = 2.0
        self.penalty_power = 1.0
        self.vt=vt
        self.nc = vt.coefficients.nc
        self.mesh = vt.mesh
        self.DOFBoundaryConditionsDictList=[]
        self.DOFBoundaryPointDictList=[]
        self.periodicBoundaryConditionsDictList=[]
        #for now need all of the coefficient values, so don't have to distinquish
        # calls to evaluateCoefficients
        #see Transport.py elemenbt boundary quadrature initialization routine
        self.ebqeTerms = []
        self.ebqe = {}
        for k in self.vt.points_elementBoundaryQuadrature:
            self.ebqeTerms.append(k)
            self.ebqe[k] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.vt.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
        for k in self.vt.scalars_elementBoundaryQuadrature:
            self.ebqeTerms.append(k)
            self.ebqe[k]=numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.vt.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
        for k in self.vt.vectors_elementBoundaryQuadrature:
            self.ebqeTerms.append(k)
            self.ebqe[k]=numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.vt.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.vt.nSpace_global),
                'd')

        #tjp adding diffusive flux boundary conditions
        self.ebqe[('diffusiveFlux_bc',0,0)]=numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.vt.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc_flag',0,0)]=numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.vt.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)]=numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.vt.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)]=numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.vt.nElementBoundaryQuadraturePoints_elementBoundary),'i')

        for k in self.vt.tensors_elementBoundaryQuadrature:
            self.ebqeTerms.append(k)
            if (self.vt.sd and k[0] in ['a','da'] and
                self.vt.coefficients.sdInfo is not None and
                (k[1],k[2]) in self.vt.coefficients.sdInfo.keys()):
                self.ebqe[k]=numpy.zeros(
                    (self.mesh.nExteriorElementBoundaries_global,
                     self.vt.nElementBoundaryQuadraturePoints_elementBoundary,
                     self.vt.coefficients.sdInfo[(k[1],k[2])][0][self.vt.nSpace_global]),
                    'd')
            else:
                self.ebqe[k]=numpy.zeros(
                    (self.mesh.nExteriorElementBoundaries_global,
                     self.vt.nElementBoundaryQuadraturePoints_elementBoundary,
                     self.vt.nSpace_global,
                     self.vt.nSpace_global),
                    'd')
        #cek why is g getting added?
        #self.ebqeTerms.append('g')
        #self.ebqe['g'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
        #        self.vt.nElementBoundaryQuadraturePoints_elementBoundary,
        #                                      max(1,self.vt.nSpace_global-1),
        #                              max(1,self.vt.nSpace_global-1)),
                #'d')
        if self.vt.movingDomain:
            if self.vt.ebqe.has_key('xt'):
                self.ebqe['xt'] = self.vt.ebqe['xt']
                self.ebqe['n'] = self.vt.ebqe['n']
                self.ebqe['sqrt(det(g))'] = self.vt.ebqe['sqrt(det(g))']
        #copy over stuff from vt.ebq
        for term in self.ebqeTerms:
            if self.vt.ebqe.has_key(term):
                self.ebqe[term].flat[:] = self.vt.ebqe[term].flat[:]

        logEvent(memory("ebqe","NumericalFlux"),level=4)
        self.isDOFBoundary ={}
        self.isDiffusiveFluxBoundary = {}
        self.isAdvectiveFluxBoundary = {}
        self.mixedDiffusion={}
        assert(getPointwiseBoundaryConditions is not None), "must supply Dirichlet conditions"
        class ptuple:
            """
            define a dictionary key that defines points as equal if they're "close"
            """
            h=self.vt.mesh.hMin*0.001
            def __init__(self,p):
                self.p=p
            def __hash__(self):
                return hash(tuple(self.p))
            def __str__(self):
                return `self.p`
            def __eq__(self,other):
                return  enorm(self.p - other.p) < self.h
            def __lt__(self,other):
                if not self == other:
                    return tuple(self.p) < tuple(other.p)
                else:
                    return False
            def __gt__(self,other):
                if not self == other:
                    return tuple(self.p) > tuple(other.p)
                else:
                    return False
            def __ge__(self,other):
                assert(False)
            def __le__(self,other):
                assert(False)
            def __ne__(self,other):
                assert(False)
        for ci in range(self.nc):
            self.isDOFBoundary[ci] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,vt.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.isDiffusiveFluxBoundary[ci] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,vt.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.isAdvectiveFluxBoundary[ci] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,vt.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            logEvent(memory("isDOF,DiffFlux,AdvFluxBoundary[%s]" % (ci),"NumericalFlux"),level=4)
            self.mixedDiffusion[ci]=False
            self.DOFBoundaryConditionsDictList.append({})
            self.DOFBoundaryPointDictList.append({})
            self.periodicBoundaryConditionsDictList.append({})
            pset=set()
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN_global = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                materialFlag = self.mesh.elementBoundaryMaterialTypes[ebN]
                for k in range(vt.nElementBoundaryQuadraturePoints_elementBoundary):
                    x = vt.ebqe['x'][ebNE,k]
                    try:
                        #mwf now try to have flag for boundary type returned too
                        gReturn = getPointwiseBoundaryConditions[ci](x,materialFlag)
                        try:
                            #mwf debug
                            #import pdb
                            #pdb.set_trace()
                            g = gReturn[0]
                            gFlag = gReturn[1]
                        except TypeError: #didn't return a tuple?
                            g = gReturn
                            gFlag = 1
                        p = None
                        if getPeriodicBoundaryConditions is not None and not parallelPeriodic:
                            p = getPeriodicBoundaryConditions[ci](x,materialFlag)
                        if p is not None:
                            pset.add(ptuple(p))
                            #self.isDOFBoundary[ci][ebNE,k]=1
                            if ptuple(p) in self.periodicBoundaryConditionsDictList[ci].keys():#self.periodicBoundaryConditionsDictList[ci].has_key(ptuple(p)):
                                i = self.periodicBoundaryConditionsDictList[ci].keys().index(ptuple(p))
                                key = self.periodicBoundaryConditionsDictList[ci].keys()[i]
                                self.periodicBoundaryConditionsDictList[ci][key].append((ebNE,k))
                            else:
                                self.periodicBoundaryConditionsDictList[ci][ptuple(p)]=[(ebNE,k)]
                        if g is not None:
                            self.isDOFBoundary[ci][ebNE,k]=gFlag
                            self.DOFBoundaryConditionsDictList[ci][(ebNE,k)] = g
                            self.DOFBoundaryPointDictList[ci][(ebNE,k)]=x
                    except TypeError:
                        logEvent("""WARNING NumericalFlux Pointwise conditions should take arguments (x,flag) now trying without flag""")
                        #again try to get boundary condition flag too
                        gReturn = getPointwiseBoundaryConditions[ci](x)
                        try:
                            g = gReturn[0]
                            gFlag = gReturn[1]
                        except TypeError:
                            g = gReturn
                            gFlag = 1
                        p = None
                        if getPeriodicBoundaryConditions is not None:
                            p = getPeriodicBoundaryConditions[ci](x)
                        if p is not None and not parallelPeriodic:
                            if self.periodicBoundaryConditionsDictList[ci].has_key(ptuple(p)):#.hash()):
                                self.periodicBoundaryConditionsDictList[ci][ptuple(p)].append((ebNE,k))
                            else:
                                self.periodicBoundaryConditionsDictList[ci][ptuple(p)]=[(ebNE,k)]#[ptuple(p).hash()]=[(ebNE,k)]
                        elif g is not None:
                            self.isDOFBoundary[ci][ebNE,k]=gFlag
                            self.DOFBoundaryConditionsDictList[ci][(ebNE,k)] = g
                            self.DOFBoundaryPointDictList[ci][(ebNE,k)]=x
                    if getAdvectiveFluxBoundaryConditions.has_key(ci):
                        try:
                            g = getAdvectiveFluxBoundaryConditions[ci](x,materialFlag)
                            if g is not None:
                                self.isAdvectiveFluxBoundary[ci][ebNE,k]=1
                        except TypeError:
                            logEvent("""WARNING NumericalFlux Pointwise conditions should take arguments (x,flag) now trying without flag""")
                            g = getAdvectiveFluxBoundaryConditions[ci](x)
                            if g is not None:
                                self.isAdvectiveFluxBoundary[ci][ebNE,k]=1
                    if getDiffusiveFluxBoundaryConditions.has_key(ci):
                        for ck in getDiffusiveFluxBoundaryConditions[ci].keys():
                            try:
                                g = getDiffusiveFluxBoundaryConditions[ci][ck](x,materialFlag)
                                if g is not None:
                                    self.isDiffusiveFluxBoundary[ci][ebNE,k]=1
                            except TypeError:
                                logEvent("""WARNING NumericalFlux Pointwise conditions should take arguments (x,flag) now trying without flag""")
                                g = getDiffusiveFluxBoundaryConditions[ci][ck](x)
                                if g is not None:
                                    self.isDiffusiveFluxBoundary[ci][ebNE,k]=1
        #ci
        #import pdb
        #pdb.set_trace()
        #print "pcs==========================",self.periodicBoundaryConditionsDictList[0]
        #now include some notion of equation terms that the flux applies to
        self.advectiveNumericalFlux = {}; self.diffusiveNumericalFlux = {}; self.HamiltonJacobiNumericalFlux = {}
        #print "pcs=========================="
        #for t in self.periodicBoundaryConditionsDictList[1].values(): print t
        #import pdb
        #pdb.set_trace()
        for ci in range(self.nc):
            #by default advection-diffusion
            self.advectiveNumericalFlux[ci] = True
            self.diffusiveNumericalFlux[ci] = True
            self.HamiltonJacobiNumericalFlux[ci] = False

        logEvent(memory("boundary condition data structures","NumericalFlux"),level=4)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        pass
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        pass

    def setFromOptions(self,nOptions):
        """
        allow classes to set various numerical parameters
        """
        pass

class DoNothing(NF_base):
    hasInterior=False
    def __init__(self,vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        self.penalty_power=0
        self.penalty_constant=0.0
        self.hasInterior=False
        self.mixedDiffusion=[False for ci in range(vt.nc)]
        self.includeBoundaryAdjoint=False
        self.advectiveNumericalFlux = {}; self.diffusiveNumericalFlux = {}; self.HamiltonJacobiNumericalFlux = {}
    def setDirichletValues(self,ebqe):
        pass

class NoFlux(NF_base):
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            if ebq_global.has_key(('advectiveFlux',ci)):
                ebq_global[('advectiveFlux',ci)].flat[:]=0.0
                ebq_global[('dadvectiveFlux_left',ci,ci)].flat[:]=0.0
                ebq_global[('dadvectiveFlux_right',ci,ci)].flat[:]=0.0
            for ck in range(self.nc):
                if ebq_global.has_key(('diffusiveFlux',ck,ci)):
                    ebq_global[('diffusiveFlux',ck,ci)].flat[:]=0.0
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            for cj  in range(self.nc):
                fluxJacobian[ci][cj].flat[:]=0.0

class StrongDirichlet(NF_base):
    hasInterior=False
    useWeakDirichletConditions=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions,
                         getPeriodicBoundaryConditions)
        self.fluxBoundaryConditions = {}
        self.setFluxBoundaryConditions()
        self.hasInterior=StrongDirichlet.hasInterior
    def setFluxBoundaryConditions(self):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        import pdb
#        pdb.set_trace()
        for ci,cjDict in self.vt.coefficients.advection.iteritems():
            if (self.fluxBoundaryConditions[ci] == 'outFlow' or
                self.fluxBoundaryConditions[ci] == 'mixedFlow'):
                for cj in cjDict:
                    cnumericalFlux.calculateExteriorNumericalAdvectiveFlux_NoBC(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[ci],
                                                                                ebqe['n'],
                                                                                ebqe[('f',ci)],
                                                                                ebqe[('df',ci,cj)],
                                                                                ebqe[('advectiveFlux',ci)],
                                                                                ebqe[('dadvectiveFlux_left',ci,cj)])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci,cjDict in self.vt.coefficients.advection.iteritems():
            if ((self.fluxBoundaryConditions[ci] == 'outFlow' or
                 self.fluxBoundaryConditions[ci] == 'mixedFlow') and
                self.vt.timeIntegration.advectionIsImplicit[ci]):
                for cj in cjDict:
                    cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[ci],
                                                                                ebqe[('dadvectiveFlux_left',ci,cj)],
                                                                                ebqe[('v',cj)],
                                                                                fluxJacobian_exterior[ci][cj])

def StrongDirichletFactory(fluxBoundaryConditionsDict={}):
    class MyStrongDirichlet(StrongDirichlet):
        hasInterior=False
        def setFluxBoundaryConditions(self):
            self.fluxBoundaryConditions=fluxBoundaryConditionsDict
    return MyStrongDirichlet

class Advection_DiagonalUpwind(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.outFlowOnly=True
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorNumericalAdvectiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                   self.mesh.elementBoundaryElementsArray,
                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                   ebq['n'],
                                                                   ebq[('u',ci)],
                                                                   ebq[('f',ci)],
                                                                   ebq[('df',ci,ci)],
                                                                   ebq_global[('advectiveFlux',ci)],
                                                                   ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                   ebq_global[('dadvectiveFlux_right',ci,ci)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        #mwf hack need this?
        #for ci in range(self.nc):
        #    ebqe[('advectiveFlux',ci)].flat[:] = 0.0
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs ebNE=%d k=%d g=%s" % (ebNE,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                #print "bc--------------------",bci
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
#                self.ebqe[('u',ci)][bci[0][0],bci[0][1]]=ebqe[('u',ci)][bci[1][0],bci[1][1]]
#                self.ebqe[('u',ci)][bci[1][0],bci[1][1]]=ebqe[('u',ci)][bci[0][0],bci[0][1]]
#                self.ebqe[('u',ci)][ebNE_L,k_L]=ebqe[('u',ci)][ebNE_R,k_R]
#                self.ebqe[('u',ci)][ebNE_R,k_R]=ebqe[('u',ci)][ebNE_L,k_L]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        #import pdb
        #pdb.set_trace()
        for ci in range(self.nc):
            if self.outFlowOnly:
                cnumericalFlux.calculateExteriorNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       self.isDOFBoundary[ci],
                                                                       inflowFlag[ci],
                                                                       ebqe['n'],
                                                                       self.ebqe[('u',ci)],
                                                                       self.ebqe[('f',ci)],
                                                                       self.ebqe[('df',ci,ci)],
                                                                       ebqe[('u',ci)],
                                                                       ebqe[('f',ci)],
                                                                       ebqe[('df',ci,ci)],
                                                                       ebqe[('advectiveFlux',ci)],
                                                                       ebqe[('dadvectiveFlux_left',ci,ci)])
            else:
                cnumericalFlux.calculateExteriorNumericalAdvectiveFlux_free(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[ci],
                                                                            inflowFlag[ci],
                                                                            ebqe['n'],
                                                                            self.ebqe[('u',ci)],
                                                                            self.ebqe[('f',ci)],
                                                                            self.ebqe[('df',ci,ci)],
                                                                            ebqe[('u',ci)],
                                                                            ebqe[('f',ci)],
                                                                            ebqe[('df',ci,ci)],
                                                                            ebqe[('advectiveFlux',ci)],
                                                                            ebqe[('dadvectiveFlux_left',ci,ci)])
        #mwf add for inflow flux?
        for ci in range(self.nc):
            cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               inflowFlag[ci],
                                                                               ebqe[('inflowFlux',ci)],
                                                                               ebqe['n'],
                                                                               ebqe[('f',ci)],
                                                                               ebqe[('df',ci,ci)],
                                                                               ebqe[('advectiveFlux',ci)],
                                                                               ebqe[('dadvectiveFlux_left',ci,ci)])
        #mwf debug
        #import pdb
        #pdb.set_trace()
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateInteriorNumericalAdvectiveFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                           ebq_global[('dadvectiveFlux_right',ci,ci)],
                                                                           ebq[('v',ci)],
                                                                           fluxJacobian[ci][ci])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])
class Advection_Diagonal_average(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorNumericalAdvectiveFlux_average(self.mesh.interiorElementBoundariesArray,
                                                                  self.mesh.elementBoundaryElementsArray,
                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                  ebq['n'],
                                                                  ebq[('u',ci)],
                                                                  ebq[('f',ci)],
                                                                  ebq[('df',ci,ci)],
                                                                  ebq_global[('advectiveFlux',ci)],
                                                                  ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                  ebq_global[('dadvectiveFlux_right',ci,ci)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs eN=%d ebN=%d k=%d g=%s" % (eN,ebN,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFlux_average(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[ci],
                                                                           inflowFlag[ci],
                                                                           ebqe['n'],
                                                                           self.ebqe[('u',ci)],
                                                                           self.ebqe[('f',ci)],
                                                                           self.ebqe[('df',ci,ci)],
                                                                           ebqe[('u',ci)],
                                                                           ebqe[('f',ci)],
                                                                           ebqe[('df',ci,ci)],
                                                                           ebqe[('advectiveFlux',ci)],
                                                                           ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateInteriorNumericalAdvectiveFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                           ebq_global[('dadvectiveFlux_right',ci,ci)],
                                                                           ebq[('v',ci)],
                                                                           fluxJacobian[ci][ci])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])
class Advection_DiagonalUpwind_Diffusion_IIPG(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateInteriorNumericalAdvectiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                      self.mesh.elementBoundaryElementsArray,
                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                      ebq['n'],
                                                                      ebq[('u',ci)],
                                                                      ebq[('f',ci)],
                                                                      ebq[('df',ci,ci)],
                                                                      ebq_global[('advectiveFlux',ci)],
                                                                      ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                      ebq_global[('dadvectiveFlux_right',ci,ci)])
                #print "af",ebq_global[('advectiveFlux',ci)]
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.interiorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  ebq['n'],
                                                                                  ebq[('a',ci,ck)],
                                                                                  ebq[('grad(phi)',ck)],
                                                                                  ebq[('u',ck)],
                                                                                  ebq_global[('penalty')],
                                                                                  ebq_global[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               ebq[('a',ci,ck)],
                                                                               ebq[('grad(phi)',ck)],
                                                                               ebq[('u',ck)],
                                                                               ebq_global[('penalty')],
                                                                               ebq_global[('diffusiveFlux',ck,ci)])
                    #print "df",ebq_global[('diffusiveFlux',ck,ci)]
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs eN=%d ebN=%d k=%d g=%s" % (eN,ebN,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                cnumericalFlux.calculateExteriorNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       self.isDOFBoundary[ci],
                                                                       inflowFlag[ci],
                                                                       ebqe['n'],
                                                                       self.ebqe[('u',ci)],
                                                                       self.ebqe[('f',ci)],
                                                                       self.ebqe[('df',ci,ci)],
                                                                       ebqe[('u',ci)],
                                                                       ebqe[('f',ci)],
                                                                       ebqe[('df',ci,ci)],
                                                                       ebqe[('advectiveFlux',ci)],
                                                                       ebqe[('dadvectiveFlux_left',ci,ci)])
                #print "afe",ebq_global[('advectiveFlux',ci)]
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)])
                    #print "dfe",ebq_global[('diffusiveFlux',ck,ci)]
 #        for ci in range(self.nc):
#             cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
#                                                                                self.mesh.elementBoundaryElementsArray,
#                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                inflowFlag[ci],
#                                                                                ebqe[('inflowFlux',ci)],
#                                                                                ebqe['n'],
#                                                                                ebqe[('f',ci)],
#                                                                                ebqe[('df',ci,ci)],
#                                                                                ebqe[('advectiveFlux',ci)],
#                                                                                ebqe[('dadvectiveFlux_left',ci,ci)])
#         #mwf add for inflow flux?
#         for ci in range(self.nc):
#             if ebqe.has_key(('f',ci)):
#                 cnumericalFlux.calculateExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
#                                                                              self.mesh.elementBoundaryElementsArray,
#                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                              inflowFlag[ci],
#                                                                              ebqe[('inflowFlux',ci)],
#                                                                              ebqe['n'],
#                                                                              ebqe[('f',ci)],
#                                                                              ebqe[('df',ci,ci)],
#                                                                              ebqe[('advectiveFlux',ci)],
#                                                                              ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                if ebq.has_key(('f',ci)):
                    cnumericalFlux.updateInteriorNumericalAdvectiveFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                               ebq_global[('dadvectiveFlux_right',ci,ci)],
                                                                               ebq[('v',ci)],
                                                                               fluxJacobian[ci][ci])
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    for cj in range(self.nc):
                        if dphi.has_key((ck,cj)):
                            if self.vt.sd:
                                cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                               dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                               self.mesh.interiorElementBoundariesArray,
                                                                                               self.mesh.elementBoundaryElementsArray,
                                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                               ebq['n'],
                                                                                               ebq[('a',ci,ck)],
                                                                                               ebq[('da',ci,ck,cj)],
                                                                                               ebq[('grad(phi)',ck)],
                                                                                               dphi[(ck,cj)].dof,
                                                                                               ebq[('v',cj)],
                                                                                               ebq[('grad(v)',cj)],
                                                                                               ebq_global['penalty'],
                                                                                               fluxJacobian[ci][cj])
                            else:
                                cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                            self.mesh.interiorElementBoundariesArray,
                                                                                            self.mesh.elementBoundaryElementsArray,
                                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                            ebq['n'],
                                                                                            ebq[('a',ci,ck)],
                                                                                            ebq[('da',ci,ck,cj)],
                                                                                            ebq[('grad(phi)',ck)],
                                                                                            dphi[(ck,cj)].dof,
                                                                                            ebq[('v',cj)],
                                                                                            ebq[('grad(v)',cj)],
                                                                                            ebq_global['penalty'],
                                                                                            fluxJacobian[ci][cj])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        #cek hack
        #pass
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                if ebqe.has_key(('f',ci)):
                    cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[ci],
                                                                                ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                                ebqe[('v',ci)],
                                                                                fluxJacobian_exterior[ci][ci])
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj])
class Advection_DiagonalUpwind_Diffusion_SIPG(Advection_DiagonalUpwind_Diffusion_IIPG):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG.__init__(self,vt,getPointwiseBoundaryConditions,
                                                         getAdvectiveFluxBoundaryConditions,
                                                         getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0

class Advection_DiagonalUpwind_Diffusion_NIPG(Advection_DiagonalUpwind_Diffusion_IIPG):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG.__init__(self,vt,getPointwiseBoundaryConditions,
                                                         getAdvectiveFluxBoundaryConditions,
                                                         getDiffusiveFluxBoundaryConditions)
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=-1.0

class Advection_DiagonalUpwind_exterior(Advection_DiagonalUpwind):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        Advection_DiagonalUpwind.__init__(self,vt,getPointwiseBoundaryConditions,
                                          getAdvectiveFluxBoundaryConditions,
                                          getDiffusiveFluxBoundaryConditions,
                                          getPeriodicBoundaryConditions)
        self.hasInterior=False
        self.outFlowOnly=True
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass

class Advection_DiagonalUpwind_Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                   self.mesh.elementBoundaryElementsArray,
                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                   self.isDOFBoundary[ci],
                                                                   inflowFlag[ci],
                                                                   ebqe['n'],
                                                                   self.ebqe[('u',ci)],
                                                                   self.ebqe[('f',ci)],
                                                                   self.ebqe[('df',ci,ci)],
                                                                   ebqe[('u',ci)],
                                                                   ebqe[('f',ci)],
                                                                   ebqe[('df',ci,ci)],
                                                                   ebqe[('advectiveFlux',ci)],
                                                                   ebqe[('dadvectiveFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ci],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)])
        #mwf add for inflow flux?
        #for ci in range(self.nc):
        #    cnumericalFlux.calculateExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
        #                                                                 self.mesh.elementBoundaryElementsArray,
        #                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                                                 inflowFlag[ci],
        #                                                                 ebqe[('inflowFlux',ci)],
        #                                                                 ebqe['n'],
        #                                                                 ebqe[('f',ci)],
        #                                                                 ebqe[('df',ci,ci)],
        #                                                                 ebqe[('advectiveFlux',ci)],
        #                                                                 ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj])
class ConstantAdvection_Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
        self.scale_penalty = 1
        self.penalty_floor = 0.0
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                ebqe[('advectiveFlux',ci)][:] = (ebqe[('f',ci)]*ebqe['n']).sum(-1)
            ebqe[('dadvectiveFlux_left',ci,ci)][:] = 0.0
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ci],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)],
                                                                                  self.scale_penalty,
                                                                                  self.penalty_floor)
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                fluxJacobian_exterior[ci][cj][:] = 0.0
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj],
                                                                                                   self.scale_penalty,
                                                                                                   self.penalty_floor)
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj],
                                                                                                self.scale_penalty,
                                                                                                self.penalty_floor)
class ConstantAdvection_exterior(NF_base):
    useStrongDirichletConstraints=True
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            ebqe[('advectiveFlux',ci)][:] = (ebqe[('f',ci)]*ebqe['n']).sum(-1)
            ebqe[('dadvectiveFlux_left',ci,ci)][:] = 0.0
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        pass

class MixedDarcy_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        cnumericalFlux.calculateExteriorNumericalAdvectiveFlux_NoBC(self.mesh.exteriorElementBoundariesArray,
                                                                    self.mesh.elementBoundaryElementsArray,
                                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                    inflowFlag[1],
                                                                    ebqe['n'],
                                                                    ebqe[('f',1)],
                                                                    ebqe[('df',1,1)],
                                                                    ebqe[('advectiveFlux',1)],
                                                                    ebqe[('dadvectiveFlux_left',1,1)])
        for ci in range(self.nc):
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ci],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)])
        #mwf add for inflow flux?
        #for ci in range(self.nc):
        #    cnumericalFlux.calculateExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
        #                                                                 self.mesh.elementBoundaryElementsArray,
        #                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                                                 inflowFlag[ci],
        #                                                                 ebqe[('inflowFlux',ci)],
        #                                                                 ebqe['n'],
        #                                                                 ebqe[('f',ci)],
        #                                                                 ebqe[('df',ci,ci)],
        #                                                                 ebqe[('advectiveFlux',ci)],
        #                                                                 ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        if self.vt.timeIntegration.advectionIsImplicit[ci]:
            cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[1],
                                                                            ebqe[('dadvectiveFlux_left',1,1)],
                                                                            ebqe[('v',1)],
                                                                            fluxJacobian_exterior[1][1])
        for ci in range(self.nc):
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj])

class Advection_DiagonalUpwind_Diffusion_NIPG_exterior(Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,
                                                                  vt,
                                                                  getPointwiseBoundaryConditions,
                                                                  getAdvectiveFluxBoundaryConditions,
                                                                  getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=-1.0

class Advection_DiagonalUpwind_Diffusion_SIPG_exterior(Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,
                                                                  vt,
                                                                  getPointwiseBoundaryConditions,
                                                                  getAdvectiveFluxBoundaryConditions,
                                                                  getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0

class ConstantAdvection_Diffusion_NIPG_exterior(ConstantAdvection_Diffusion_IIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        ConstantAdvection_Diffusion_IIPG_exterior.__init__(self,
                                                           vt,
                                                           getPointwiseBoundaryConditions,
                                                           getAdvectiveFluxBoundaryConditions,
                                                           getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=-1.0

class ConstantAdvection_Diffusion_SIPG_exterior(ConstantAdvection_Diffusion_IIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        ConstantAdvection_Diffusion_IIPG_exterior.__init__(self,
                                                           vt,
                                                           getPointwiseBoundaryConditions,
                                                           getAdvectiveFluxBoundaryConditions,
                                                           getDiffusiveFluxBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.scale_penalty = 1; self.penalty_floor = 0.0

class Advection_DiagonalUpwind_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                   self.mesh.elementBoundaryElementsArray,
                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                   self.isDOFBoundary[ci],
                                                                   inflowFlag[ci],
                                                                   ebqe['n'],
                                                                   self.ebqe[('u',ci)],
                                                                   self.ebqe[('f',ci)],
                                                                   self.ebqe[('df',ci,ci)],
                                                                   ebqe[('u',ci)],
                                                                   ebqe[('f',ci)],
                                                                   ebqe[('df',ci,ci)],
                                                                   ebqe[('advectiveFlux',ci)],
                                                                   ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])
class Curvature_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFlux_free(self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('f',ci)],
                                                                        self.ebqe[('df',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('f',ci)],
                                                                        ebqe[('df',ci,ci)],
                                                                        ebqe[('advectiveFlux',ci)],
                                                                        ebqe[('dadvectiveFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_free_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                       self.mesh.exteriorElementBoundariesArray,
                                                                                       self.mesh.elementBoundaryElementsArray,
                                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                       self.isDOFBoundary[ci],
                                                                                       ebqe['n'],
                                                                                       self.ebqe[('a',ci,ck)],
                                                                                       self.ebqe[('grad(phi)',ck)],
                                                                                       self.ebqe[('u',ck)],
                                                                                       ebqe[('a',ci,ck)],
                                                                                       ebqe[('grad(phi)',ck)],
                                                                                       ebqe[('u',ck)],
                                                                                       ebqe[('penalty')],
                                                                                       ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_free(self.mesh.exteriorElementBoundariesArray,
                                                                                    self.mesh.elementBoundaryElementsArray,
                                                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                    self.isDOFBoundary[ci],
                                                                                    ebqe['n'],
                                                                                    self.ebqe[('a',ci,ck)],
                                                                                    self.ebqe[('grad(phi)',ck)],
                                                                                    self.ebqe[('u',ck)],
                                                                                    ebqe[('a',ci,ck)],
                                                                                    ebqe[('grad(phi)',ck)],
                                                                                    ebqe[('u',ck)],
                                                                                    ebqe[('penalty')],
                                                                                    ebqe[('diffusiveFlux',ck,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian_free(self.mesh.exteriorElementBoundariesArray,
                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                 inflowFlag[ci],
                                                                                 ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                                 ebqe[('v',ci)],
                                                                                 fluxJacobian_exterior[ci][ci])
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_free_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                        dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                                                        self.mesh.elementBoundaryElementsArray,
                                                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                        self.isDOFBoundary[ck],
                                                                                                        ebqe['n'],
                                                                                                        ebqe[('a',ci,ck)],
                                                                                                        ebqe[('da',ci,ck,cj)],
                                                                                                        ebqe[('grad(phi)',ck)],
                                                                                                        dphi[(ck,cj)].dof,
                                                                                                        ebqe[('v',cj)],
                                                                                                        ebqe[('grad(v)',cj)],
                                                                                                        ebqe['penalty'],
                                                                                                        fluxJacobian_exterior[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_free(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                     self.mesh.exteriorElementBoundariesArray,
                                                                                                     self.mesh.elementBoundaryElementsArray,
                                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                     self.isDOFBoundary[ck],
                                                                                                     ebqe['n'],
                                                                                                     ebqe[('a',ci,ck)],
                                                                                                     ebqe[('da',ci,ck,cj)],
                                                                                                     ebqe[('grad(phi)',ck)],
                                                                                                     dphi[(ck,cj)].dof,
                                                                                                     ebqe[('v',cj)],
                                                                                                     ebqe[('grad(v)',cj)],
                                                                                                     ebqe['penalty'],
                                                                                                     fluxJacobian_exterior[ci][cj])
class Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    """
    To use with regular Stokes, takes advantage of existence of 'advectiveFlux' flag
    even when there is no advective term
    """
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,
                         vt,
                         getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions,
                         getPeriodicBoundaryConditions)
        self.hasInterior=False
        self.scale_penalty = 1
        self.penalty_floor = 0.0
        self.penalty_constant = 100.0
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.nSpace_global == 2:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxStokes2D(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[0],
                                                                           self.isDOFBoundary[1],
                                                                           self.isDOFBoundary[2],
                                                                           ebqe['n'],
                                                                           self.ebqe[('u',0)],
                                                                           self.ebqe[('f',0)],
                                                                           ebqe[('u',0)],
                                                                           ebqe[('f',0)],
                                                                           ebqe[('df',0,1)],
                                                                           ebqe[('df',0,2)],
                                                                           ebqe[('advectiveFlux',0)],
                                                                           ebqe[('advectiveFlux',1)],
                                                                           ebqe[('advectiveFlux',2)],
                                                                           ebqe[('dadvectiveFlux_left',0,1)],
                                                                           ebqe[('dadvectiveFlux_left',0,2)],
                                                                           ebqe[('dadvectiveFlux_left',1,0)],
                                                                           ebqe[('dadvectiveFlux_left',2,0)],
                                                                           ebqe[('velocity',0)])
        elif self.vt.nSpace_global == 3:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxStokes3D(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[0],
                                                                           self.isDOFBoundary[1],
                                                                           self.isDOFBoundary[2],
                                                                           self.isDOFBoundary[3],
                                                                           ebqe['n'],
                                                                           self.ebqe[('u',0)],
                                                                           self.ebqe[('f',0)],
                                                                           ebqe[('u',0)],
                                                                           ebqe[('f',0)],
                                                                           ebqe[('df',0,1)],
                                                                           ebqe[('df',0,2)],
                                                                           ebqe[('df',0,3)],
                                                                           ebqe[('advectiveFlux',0)],
                                                                           ebqe[('advectiveFlux',1)],
                                                                           ebqe[('advectiveFlux',2)],
                                                                           ebqe[('advectiveFlux',3)],
                                                                           ebqe[('dadvectiveFlux_left',0,1)],
                                                                           ebqe[('dadvectiveFlux_left',0,2)],
                                                                           ebqe[('dadvectiveFlux_left',0,3)],
                                                                           ebqe[('dadvectiveFlux_left',1,0)],
                                                                           ebqe[('dadvectiveFlux_left',2,0)],
                                                                           ebqe[('dadvectiveFlux_left',3,0)],
                                                                           ebqe[('velocity',0)])
        for ci in range(1,self.nc):
            if ebqe.has_key(('a',ci,ci)):
                if self.vt.sd:
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],
                                                                              self.vt.coefficients.sdInfo[(ci,ci)][1],
                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              self.isDOFBoundary[ci],
                                                                              ebqe['n'],
                                                                              self.ebqe[('a',ci,ci)],
                                                                              self.ebqe[('grad(phi)',ci)],
                                                                              self.ebqe[('u',ci)],
                                                                              ebqe[('a',ci,ci)],
                                                                              ebqe[('grad(phi)',ci)],
                                                                              ebqe[('u',ci)],
                                                                              ebqe[('penalty')],
                                                                              ebqe[('diffusiveFlux',ci,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
                else:
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[ci],
                                                                           ebqe['n'],
                                                                           self.ebqe[('a',ci,ci)],
                                                                           self.ebqe[('grad(phi)',ci)],
                                                                           self.ebqe[('u',ci)],
                                                                           ebqe[('a',ci,ci)],
                                                                           ebqe[('grad(phi)',ci)],
                                                                           ebqe[('u',ci)],
                                                                           ebqe[('penalty')],
                                                                           ebqe[('diffusiveFlux',ci,ci)],
                                                                           self.scale_penalty,
                                                                           self.penalty_floor)
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb, fluxJacobian_hj):
        for ci in range(self.nc):
            for cj in range(self.nc):
                if ebqe.has_key(('dadvectiveFlux_left',ci,cj)):
                    cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[ci],
                                                                                ebqe[('dadvectiveFlux_left',ci,cj)],
                                                                                ebqe[('v',cj)],
                                                                                fluxJacobian_exterior[ci][cj])
        for ci in range(1,self.nc):
            if self.vt.sd:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],self.vt.coefficients.sdInfo[(ci,ci)][1],
                                                                               dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                               self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               ebqe[('a',ci,ci)],
                                                                               ebqe[('da',ci,ci,ci)],
                                                                               ebqe[('grad(phi)',ci)],
                                                                               dphi[(ci,ci)].dof,
                                                                               ebqe[('v',ci)],
                                                                               ebqe[('grad(v)',ci)],
                                                                               ebqe['penalty'],
                                                                               fluxJacobian_exterior[ci][ci],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
            else:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                            self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[ci],
                                                                            ebqe['n'],
                                                                            ebqe[('a',ci,ci)],
                                                                            ebqe[('da',ci,ci,ci)],
                                                                            ebqe[('grad(phi)',ci)],
                                                                            dphi[(ci,ci)].dof,
                                                                            ebqe[('v',ci)],
                                                                            ebqe[('grad(v)',ci)],
                                                                            ebqe['penalty'],
                                                                            fluxJacobian_exterior[ci][ci],
                                                                            self.scale_penalty,
                                                                            self.penalty_floor)

class Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior(Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,
                                                                         vt,
                                                                         getPointwiseBoundaryConditions,
                                                                         getAdvectiveFluxBoundaryConditions,
                                                                         getDiffusiveFluxBoundaryConditions,
                                                                         getPeriodicBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0

class StokesP_Advection_DiagonalUpwind_Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    """
    Weak boundary fluxes for StokesP
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs eN=%d ebN=%d k=%d g=%s" % (eN,ebN,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.nSpace_global == 2:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxStokesP2D(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[0],
                                                                           self.isDOFBoundary[1],
                                                                           self.isDOFBoundary[2],
                                                                           ebqe['n'],
                                                                           self.ebqe[('f',0)],
                                                                           self.ebqe[('f',1)],
                                                                           self.ebqe[('f',2)],
                                                                           ebqe[('f',0)],
                                                                           ebqe[('f',1)],
                                                                           ebqe[('f',2)],
                                                                           ebqe[('df',0,1)],
                                                                           ebqe[('df',0,2)],
                                                                           ebqe[('df',1,0)],
                                                                           ebqe[('df',2,0)],
                                                                           ebqe[('advectiveFlux',0)],
                                                                           ebqe[('advectiveFlux',1)],
                                                                           ebqe[('advectiveFlux',2)],
                                                                           ebqe[('dadvectiveFlux_left',0,1)],
                                                                           ebqe[('dadvectiveFlux_left',0,2)],
                                                                           ebqe[('dadvectiveFlux_left',1,0)],
                                                                           ebqe[('dadvectiveFlux_left',2,0)])
        elif self.vt.nSpace_global == 3:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxStokesP3D(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[0],
                                                                           self.isDOFBoundary[1],
                                                                           self.isDOFBoundary[2],
                                                                           self.isDOFBoundary[3],
                                                                           ebqe['n'],
                                                                           self.ebqe[('f',0)],
                                                                           self.ebqe[('f',1)],
                                                                           self.ebqe[('f',2)],
                                                                           self.ebqe[('f',3)],
                                                                           ebqe[('f',0)],
                                                                           ebqe[('f',1)],
                                                                           ebqe[('f',2)],
                                                                           ebqe[('f',3)],
                                                                           ebqe[('df',0,1)],
                                                                           ebqe[('df',0,2)],
                                                                           ebqe[('df',0,3)],
                                                                           ebqe[('df',1,0)],
                                                                           ebqe[('df',2,0)],
                                                                           ebqe[('df',3,0)],
                                                                           ebqe[('advectiveFlux',0)],
                                                                           ebqe[('advectiveFlux',1)],
                                                                           ebqe[('advectiveFlux',2)],
                                                                           ebqe[('advectiveFlux',3)],
                                                                           ebqe[('dadvectiveFlux_left',0,1)],
                                                                           ebqe[('dadvectiveFlux_left',0,2)],
                                                                           ebqe[('dadvectiveFlux_left',0,3)],
                                                                           ebqe[('dadvectiveFlux_left',1,0)],
                                                                           ebqe[('dadvectiveFlux_left',2,0)],
                                                                           ebqe[('dadvectiveFlux_left',3,0)])

        for ci in range(1,self.nc):
            if ebqe.has_key(('a',ci,ci)):
                if self.vt.sd:
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],self.vt.coefficients.sdInfo[(ci,ci)][1],
                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              self.isDOFBoundary[ci],
                                                                              ebqe['n'],
                                                                              self.ebqe[('a',ci,ci)],
                                                                              self.ebqe[('grad(phi)',ci)],
                                                                              self.ebqe[('u',ci)],
                                                                              ebqe[('a',ci,ci)],
                                                                              ebqe[('grad(phi)',ci)],
                                                                              ebqe[('u',ci)],
                                                                              ebqe[('penalty')],
                                                                              ebqe[('diffusiveFlux',ci,ci)])
                else:
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[ci],
                                                                           ebqe['n'],
                                                                           self.ebqe[('a',ci,ci)],
                                                                           self.ebqe[('grad(phi)',ci)],
                                                                           self.ebqe[('u',ci)],
                                                                           ebqe[('a',ci,ci)],
                                                                           ebqe[('grad(phi)',ci)],
                                                                           ebqe[('u',ci)],
                                                                           ebqe[('penalty')],
                                                                           ebqe[('diffusiveFlux',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for cj in range(1,self.nc):
            cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        inflowFlag[0],#mwf should this be [cj]
                                                                        ebqe[('dadvectiveFlux_left',0,cj)],
                                                                        ebqe[('v',cj)],
                                                                        fluxJacobian_exterior[0][cj])
            cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        inflowFlag[0],
                                                                        ebqe[('dadvectiveFlux_left',cj,0)],
                                                                        ebqe[('v',0)],
                                                                        fluxJacobian_exterior[cj][0])
        for ci in range(1,self.nc):
            if self.vt.sd:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],self.vt.coefficients.sdInfo[(ci,ci)][0],
                                                                               dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                               self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               ebqe[('a',ci,ci)],
                                                                               ebqe[('da',ci,ci,ci)],
                                                                               ebqe[('grad(phi)',ci)],
                                                                               dphi[(ci,ci)].dof,
                                                                               ebqe[('v',ci)],
                                                                               ebqe[('grad(v)',ci)],
                                                                               ebqe['penalty'],
                                                                               fluxJacobian_exterior[ci][ci])
            else:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                            self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[ci],
                                                                            ebqe['n'],
                                                                            ebqe[('a',ci,ci)],
                                                                            ebqe[('da',ci,ci,ci)],
                                                                            ebqe[('grad(phi)',ci)],
                                                                            dphi[(ci,ci)].dof,
                                                                            ebqe[('v',ci)],
                                                                            ebqe[('grad(v)',ci)],
                                                                            ebqe['penalty'],
                                                                            fluxJacobian_exterior[ci][ci])
class NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions,
                         getPeriodicBoundaryConditions,
                         parallelPeriodic=True)
        self.penalty_constant=100.0
        self.hasInterior=False
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        import ctransportCoefficients
        self.setDirichletValues(ebqe)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.nSpace_global == 2:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxNavierStokes2D(self.mesh.exteriorElementBoundariesArray,
                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                 self.isDOFBoundary[0],
                                                                                 self.isDOFBoundary[1],
                                                                                 self.isDOFBoundary[2],
                                                                                 ebqe['n'],
                                                                                 self.ebqe[('u',0)],
                                                                                 self.ebqe[('f',0)],
                                                                                 self.ebqe[('f',1)],
                                                                                 self.ebqe[('f',2)],
                                                                                 ebqe[('u',0)],
                                                                                 ebqe[('dm',1,1)],
                                                                                 ebqe[('f',0)],
                                                                                 ebqe[('f',1)],
                                                                                 ebqe[('f',2)],
                                                                                 ebqe[('df',0,1)],
                                                                                 ebqe[('df',0,2)],
                                                                                 ebqe[('df',1,0)],
                                                                                 ebqe[('df',1,1)],
                                                                                 ebqe[('df',1,2)],
                                                                                 ebqe[('df',2,0)],
                                                                                 ebqe[('df',2,1)],
                                                                                 ebqe[('df',2,2)],
                                                                                 ebqe[('advectiveFlux',0)],
                                                                                 ebqe[('advectiveFlux',1)],
                                                                                 ebqe[('advectiveFlux',2)],
                                                                                 ebqe[('dadvectiveFlux_left',0,0)],
                                                                                 ebqe[('dadvectiveFlux_left',0,1)],
                                                                                 ebqe[('dadvectiveFlux_left',0,2)],
                                                                                 ebqe[('dadvectiveFlux_left',1,0)],
                                                                                 ebqe[('dadvectiveFlux_left',1,1)],
                                                                                 ebqe[('dadvectiveFlux_left',1,2)],
                                                                                 ebqe[('dadvectiveFlux_left',2,0)],
                                                                                 ebqe[('dadvectiveFlux_left',2,1)],
                                                                                 ebqe[('dadvectiveFlux_left',2,2)],
                                                                                 ebqe[('velocity',0)])
        elif self.vt.nSpace_global == 3:
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxNavierStokes3D(self.mesh.exteriorElementBoundariesArray,
                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                 self.isDOFBoundary[0],
                                                                                 self.isDOFBoundary[1],
                                                                                 self.isDOFBoundary[2],
                                                                                 self.isDOFBoundary[3],
                                                                                 ebqe['n'],
                                                                                 self.ebqe[('u',0)],
                                                                                 self.ebqe[('f',0)],
                                                                                 self.ebqe[('f',1)],
                                                                                 self.ebqe[('f',2)],
                                                                                 self.ebqe[('f',3)],
                                                                                 ebqe[('u',0)],
                                                                                 ebqe[('f',0)],
                                                                                 ebqe[('f',1)],
                                                                                 ebqe[('f',2)],
                                                                                 ebqe[('f',3)],
                                                                                 ebqe[('df',0,1)],
                                                                                 ebqe[('df',0,2)],
                                                                                 ebqe[('df',0,3)],
                                                                                 ebqe[('df',1,0)],
                                                                                 ebqe[('df',1,1)],
                                                                                 ebqe[('df',1,2)],
                                                                                 ebqe[('df',1,3)],
                                                                                 ebqe[('df',2,0)],
                                                                                 ebqe[('df',2,1)],
                                                                                 ebqe[('df',2,2)],
                                                                                 ebqe[('df',2,3)],
                                                                                 ebqe[('df',3,0)],
                                                                                 ebqe[('df',3,1)],
                                                                                 ebqe[('df',3,2)],
                                                                                 ebqe[('df',3,3)],
                                                                                 ebqe[('advectiveFlux',0)],
                                                                                 ebqe[('advectiveFlux',1)],
                                                                                 ebqe[('advectiveFlux',2)],
                                                                                 ebqe[('advectiveFlux',3)],
                                                                                 ebqe[('dadvectiveFlux_left',0,1)],
                                                                                 ebqe[('dadvectiveFlux_left',0,2)],
                                                                                 ebqe[('dadvectiveFlux_left',0,3)],
                                                                                 ebqe[('dadvectiveFlux_left',1,0)],
                                                                                 ebqe[('dadvectiveFlux_left',1,1)],
                                                                                 ebqe[('dadvectiveFlux_left',1,2)],
                                                                                 ebqe[('dadvectiveFlux_left',1,3)],
                                                                                 ebqe[('dadvectiveFlux_left',2,0)],
                                                                                 ebqe[('dadvectiveFlux_left',2,1)],
                                                                                 ebqe[('dadvectiveFlux_left',2,2)],
                                                                                 ebqe[('dadvectiveFlux_left',2,3)],
                                                                                 ebqe[('dadvectiveFlux_left',3,0)],
                                                                                 ebqe[('dadvectiveFlux_left',3,1)],
                                                                                 ebqe[('dadvectiveFlux_left',3,2)],
                                                                                 ebqe[('dadvectiveFlux_left',3,3)],
                                                                                 ebqe[('velocity',0)])
        for ci in range(1,self.nc):
            if ebqe.has_key(('a',ci,ci)):
                if self.vt.sd:
                    ebqe[('diffusiveFlux',ci,ci)].fill(0.0)
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],self.vt.coefficients.sdInfo[(ci,ci)][1],
                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              self.isDOFBoundary[ci],
                                                                              ebqe['n'],
                                                                              self.ebqe[('a',ci,ci)],
                                                                              self.ebqe[('grad(phi)',ci)],
                                                                              self.ebqe[('u',ci)],
                                                                              ebqe[('a',ci,ci)],
                                                                              ebqe[('grad(phi)',ci)],
                                                                              ebqe[('u',ci)],
                                                                              ebqe[('penalty')],
                                                                              ebqe[('diffusiveFlux',ci,ci)])
                else:
                    cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           self.isDOFBoundary[ci],
                                                                           ebqe['n'],
                                                                           self.ebqe[('a',ci,ci)],
                                                                           self.ebqe[('grad(phi)',ci)],
                                                                           self.ebqe[('u',ci)],
                                                                           ebqe[('a',ci,ci)],
                                                                           ebqe[('grad(phi)',ci)],
                                                                           ebqe[('u',ci)],
                                                                           ebqe[('penalty')],
                                                                           ebqe[('diffusiveFlux',ci,ci)])
#                 ctransportCoefficients.applyContactLineSlip(2.0*self.vt.coefficients.eps_viscosity,
#                                                             self.isDOFBoundary[ci],
#                                                             self.vt.coefficients.ebqe_phi,
#                                                             ebqe[('advectiveFlux',ci)],
#                                                             ebqe[('diffusiveFlux',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        import ctransportCoefficients
        for ci in range(self.nc):
            for cj in range(self.nc):
                if ebqe.has_key(('dadvectiveFlux_left',ci,cj)):
                    cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[0],#mwf should this be [cj]
                                                                                ebqe[('dadvectiveFlux_left',ci,cj)],
                                                                                ebqe[('v',cj)],
                                                                                fluxJacobian_exterior[ci][cj])
        for ci in range(1,self.nc):
            if self.vt.sd:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ci)][0],self.vt.coefficients.sdInfo[(ci,ci)][1],
                                                                               dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                               self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ci],
                                                                               ebqe['n'],
                                                                               ebqe[('a',ci,ci)],
                                                                               ebqe[('da',ci,ci,ci)],
                                                                               ebqe[('grad(phi)',ci)],
                                                                               dphi[(ci,ci)].dof,
                                                                               ebqe[('v',ci)],
                                                                               ebqe[('grad(v)',ci)],
                                                                               ebqe['penalty'],
                                                                               fluxJacobian_exterior[ci][ci])
            else:
                cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ci,ci)].femSpace.dofMap.l2g,
                                                                            self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[ci],
                                                                            ebqe['n'],
                                                                            ebqe[('a',ci,ci)],
                                                                            ebqe[('da',ci,ci,ci)],
                                                                            ebqe[('grad(phi)',ci)],
                                                                            dphi[(ci,ci)].dof,
                                                                            ebqe[('v',ci)],
                                                                            ebqe[('grad(v)',ci)],
                                                                            ebqe['penalty'],
                                                                            fluxJacobian_exterior[ci][ci])
#             ctransportCoefficients.applyContactLineSlipJacobian(2.0*self.vt.coefficients.eps_viscosity,
#                                                                 self.isDOFBoundary[ci],
#                                                                 self.vt.coefficients.ebqe_phi,
#                                                                 fluxJacobian_exterior[ci][ci])
class NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior(NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                               getAdvectiveFluxBoundaryConditions,
                                                                               getDiffusiveFluxBoundaryConditions,getPeriodicBoundaryConditions)
        self.penalty_constant = 10.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.hasInterior=False

class Diffusion_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.hasInterior=False
        self.scale_penalty = 1; self.penalty_floor = 0.0
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs eN=%d ebN=%d k=%d g=%s" % (eN,ebN,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    #print self.ebqe[('u',ck)]
                    ebqe[('advectiveFlux',ci)].flat[:]=0.0
                    ebqe[('diffusiveFlux',ck,ci)].flat[:]=0.0
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)],
                                                                                  self.scale_penalty,
                                                                                  self.penalty_floor)
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj],
                                                                                                   self.scale_penalty,
                                                                                                   self.penalty_floor)
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,

                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj],
                                                                                                self.scale_penalty,
                                                                                                self.penalty_floor)

class Diffusion_SIPG_exterior(Diffusion_IIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        Diffusion_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                         getAdvectiveFluxBoundaryConditions,
                                         getDiffusiveFluxBoundaryConditions,
                                         getPeriodicBoundaryConditions)
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0

class DarcySplitPressure_IIPG_exterior(NF_base):
    """
    weak dirichlet boundary conditions for Twophase_split_pressure class

    Diffusion_IIPG_exterior is ok except for need to switch between psi_w and
    psi_n bc types

    .. todo::

       put in bc that switches flux and dirichlet types
    """
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        #hold Dirichlet values for \psi_n (non-wetting phase head)
        #also need extra psi_n entries that aren't part of default quadrature
        for term in ['psi_n_bc','psi_n',('dpsi_n',0)]:
            self.ebqe[term] = numpy.zeros(self.ebqe[('u',0)].shape,'d')
        self.hasInterior=False
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                if self.isDOFBoundary[ci][ebNE,k] == 2:
                    self.ebqe['psi_n_bc'][ebNE,k] = g(x,self.vt.timeIntegration.t)
                else:
                    self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        ebqe[('diffusiveFlux',0,0)].fill(0.0)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if self.vt.sd:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcySplitPressure_sd(self.vt.coefficients.sdInfo[(0,0)][0],
                                                                                     self.vt.coefficients.sdInfo[(0,0)][1],
                                                                                     self.mesh.exteriorElementBoundariesArray,
                                                                                     self.mesh.elementBoundaryElementsArray,
                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                     self.isDOFBoundary[0],
                                                                                     ebqe['n'],
                                                                                     self.ebqe[('a',0,0)],
                                                                                     self.ebqe[('grad(phi)',0)],
                                                                                     self.ebqe[('u',0)],
                                                                                     self.ebqe['psi_n_bc'],
                                                                                     ebqe[('a',0,0)],
                                                                                     ebqe[('grad(phi)',0)],
                                                                                     ebqe[('u',0)],
                                                                                     ebqe['psi_n'],
                                                                                     ebqe[('penalty')],
                                                                                     ebqe[('diffusiveFlux',0,0)])

        else:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcySplitPressure(self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[0],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',0,0)],
                                                                                  self.ebqe[('grad(phi)',0)],
                                                                                  self.ebqe[('u',0)],
                                                                                  self.ebqe['psi_n_bc'],
                                                                                  ebqe[('a',0,0)],
                                                                                  ebqe[('grad(phi)',0)],
                                                                                  ebqe[('u',0)],
                                                                                  ebqe['psi_n'],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',0,0)])


    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        #dependence of psi_n and psi_w same so don't need special evaluation routine
        for ci in range(self.nc):
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,

                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj])


from math import *
class Diffusion_LDG(NF_base):
    """
    initial LDG numerical flux  for pure diffusion
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):

        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions)
        self.mixedDiffusion={} #mwf was True
        self.phi_trace = {}
        self.dphi_trace_left = {}
        self.dphi_trace_right = {}
        self.b={}
        self.db={}
        self.db_eb={}
        self.qV = {}
        self.ebqV = {}
        self.qDV = {}
        self.qDV_eb = {}
        self.ebqDV = {}
        self.ebqDV_eb = {}
        self.eb_aTilde={}
        self.aTilde={}
        self.eb_aHat={}
        self.aHat={}
        assert self.vt.sd, "LDG factor code currently assumes sparse diffusion rep"
        self.diagEntry_a = {} #only need if not using c for splitting calculations
        self.aSplit = 1 #how LDG diffusion matrix is factored
        self.rho_split=0 # for variable density flow
        self.aSplittingsHaveBeenAllocated= False
        for ci in range(self.nc):
            for ck in range(self.nc):
                if self.vt.coefficients.sdInfo.has_key((ci,ck)):
                    self.diagEntry_a[(ci,ck)] = numpy.zeros((self.vt.nSpace_global,),'i')
                    rowptr = self.vt.coefficients.sdInfo[(ci,ck)][0]
                    colind = self.vt.coefficients.sdInfo[(ci,ck)][1]
                    nnz    = self.vt.coefficients.sdInfo[(ci,ck)][0][self.vt.nSpace_global]
                    for I in range(self.vt.nSpace_global):
                        for m in range(rowptr[I],
                                       rowptr[I+1]):
                            if colind.flat[m] == I:
                                self.diagEntry_a[(ci,ck)][I] = m
        for ci in range(self.nc):
            self.mixedDiffusion[ci] = False
        for ck in vt.coefficients.potential.keys():
            self.mixedDiffusion[ck] = True #mwf inserted here
            self.phi_trace[ck]=numpy.zeros((self.mesh.nElementBoundaries_global,
                                            vt.nElementBoundaryQuadraturePoints_elementBoundary),
                                           'd')
            self.dphi_trace_left[ck]=numpy.zeros((self.mesh.nElementBoundaries_global,
                                                  vt.nElementBoundaryQuadraturePoints_elementBoundary),
                                                 'd')
            self.dphi_trace_right[ck]=numpy.zeros((self.mesh.nElementBoundaries_global,
                                                   vt.nElementBoundaryQuadraturePoints_elementBoundary),
                                                  'd')
            self.A_inv = numpy.zeros((self.mesh.nElements_global,
                                      vt.nDOF_trial_element[ck],
                                      vt.nDOF_trial_element[ck]),
                                     'd')
            self.vXw =numpy.zeros(
                (vt.mesh.nElements_global,
                 vt.nQuadraturePoints_element,
                 vt.nDOF_trial_element[ck],
                 vt.nDOF_test_element[ck]),
                'd')
            cfemIntegrals.calculateShape_X_weightedShape(vt.q[('v',ck)],
                                                         vt.q[('w*dV_a',ck,ck)],
                                                         self.vXw)
            cfemIntegrals.calculateVelocityProjectionMatrixLDG(self.vXw,
                                                               self.A_inv)
    #             cfemIntegrals.calculateVelocityProjectionMatrixLDG(vt.q[('vXw*dV_a',ck,ck,ck)],
    #                                                                self.A_inv)
            self.b[ck]=numpy.zeros((self.mesh.nElements_global,
                                    vt.q[('grad(phi)',ck)].shape[-1],
                                    vt.nDOF_trial_element[ck]),
                                   'd')
            self.db[ck]=numpy.zeros((self.mesh.nElements_global,
                                     vt.q[('grad(phi)',ck)].shape[-1],
                                     vt.nDOF_trial_element[ck],
                                     vt.nDOF_trial_element[ck]),
                                    'd')
            self.db_eb[ck] = numpy.zeros((self.mesh.nElements_global,
                                          self.mesh.nElementBoundaries_element,
                                          vt.q[('grad(phi)',ck)].shape[-1],
                                          vt.nDOF_trial_element[ck],
                                          vt.nDOF_trial_element[ck]),
                                         'd')
            self.qV[ck] = numpy.zeros((self.mesh.nElements_global,
                                       vt.nQuadraturePoints_element,
                                       vt.q[('grad(phi)',ck)].shape[-1]),
                                      'd')
            self.ebqV[ck] = numpy.zeros((self.mesh.nElements_global,
                                         self.mesh.nElementBoundaries_element,
                                         vt.nElementBoundaryQuadraturePoints_elementBoundary,
                                         vt.q[('grad(phi)',ck)].shape[-1]),
                                        'd')
            self.qDV[ck] = numpy.zeros((self.mesh.nElements_global,
                                        vt.nQuadraturePoints_element,
                                        vt.nDOF_trial_element[ck],
                                        vt.q[('grad(phi)',ck)].shape[-1]),
                                       'd')
            self.qDV_eb[ck] = numpy.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           vt.nQuadraturePoints_element,
                                           vt.nDOF_trial_element[ck],
                                           vt.q[('grad(phi)',ck)].shape[-1]),
                                          'd')
            self.ebqDV[ck] = numpy.zeros((self.mesh.nElements_global,
                                          self.mesh.nElementBoundaries_element,
                                          vt.nElementBoundaryQuadraturePoints_elementBoundary,
                                          vt.nDOF_trial_element[ck],
                                          vt.q[('grad(phi)',ck)].shape[-1]),
                                         'd')
            self.ebqDV_eb[ck] = numpy.zeros((self.mesh.nElements_global,
                                             self.mesh.nElementBoundaries_element,
                                             self.mesh.nElementBoundaries_element,
                                             vt.nElementBoundaryQuadraturePoints_elementBoundary,
                                             vt.nDOF_trial_element[ck],
                                             vt.q[('grad(phi)',ck)].shape[-1]),
                                            'd')

    def calculateDiffusionSplittings(self,q,ebq,ebq_global):
        if not self.aSplittingsHaveBeenAllocated:
            for ci in range(self.nc):
                for ck in range(self.nc):
                    if self.vt.coefficients.sdInfo.has_key((ci,ck)):
                        self.eb_aTilde[(ci,ck)] = numpy.zeros(ebq[('a',ci,ck)].shape,dtype='d')
                        self.eb_aHat[(ci,ck)]   = numpy.zeros(ebq[('a',ci,ck)].shape,dtype='d')
                        self.aTilde[(ci,ck)]    = numpy.zeros(q[('a',ci,ck)].shape,dtype='d')
                        self.aHat[(ci,ck)]      = numpy.zeros(q[('a',ci,ck)].shape,dtype='d')
            self.aSplittingsHaveBeenAllocated = True
        for ci in range(self.nc):
            for ck in range(self.nc):
                if self.vt.coefficients.sdInfo.has_key((ci,ck)):
                    if self.aSplit == 0:
                        self.eb_aTilde[(ci,ck)]=ebq[('a',ci,ck)]
                        self.aTilde[(ci,ck)]=q[('a',ci,ck)]
                    elif self.aSplit == 1:
                        self.eb_aHat[(ci,ck)]=ebq[('a',ci,ck)]
                        self.aHat[(ci,ck)]=q[('a',ci,ck)]
                    if self.vt.sd:
                        cnumericalFlux.calculateDiffusionMatrixSplittings_LDG_sd(self.aSplit,
                                                                                 self.vt.nSpace_global,
                                                                                 self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                 self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                 ebq[('a',ci,ck)],
                                                                                 q[('a',ci,ck)],
                                                                                 self.eb_aHat[(ci,ck)],
                                                                                 self.eb_aTilde[(ci,ck)],
                                                                                 self.aHat[(ci,ck)],
                                                                                 self.aTilde[(ci,ck)])
                    else:
                        for eN in range(self.mesh.nElements_global):
                            for ebN in range(self.mesh.nElementBoundaries_element):
                                for k in range(self.vt.nElementBoundaryQuadraturePoints_elementBoundary):
                                    for I in range(self.vt.nSpace_global):
                                        if self.aSplit==0:
                                            self.eb_aHat[(ci,ck)][eN,ebN,k,self.diagEntry_a[(ci,ck)][I]]=1.0
                                        elif self.aSplit==1:
                                            self.eb_aTilde[(ci,ck)][eN,ebN,k,self.diagEntry_a[(ci,ck)][I]]=1.0
                                        elif self.aSplit==2:
                                            self.eb_aTilde[(ci,ck)][eN,ebN,k,self.diagEntry_a[(ci,ck)][I]]=sqrt(ebq[('a',ck,ck)][eN,ebN,k,self.diagEntry_a[(ci,ck)][I]])
                        for eN in range(self.mesh.nElements_global):
                            for k in range(self.vt.nQuadraturePoints_element):
                                for I in range(self.vt.nSpace_global):
                                    if self.aSplit==0:
                                        self.aHat[(ci,ck)][eN,k,self.diagEntry_a[(ci,ck)][I]]=1.0
                                    elif self.aSplit==1:
                                        self.aTilde[(ci,ck)][eN,k,self.diagEntry_a[(ci,ck)][I]]=1.0
                                    elif self.aSplit==2:
                                        self.aTilde[(ci,ck)][eN,k,self.diagEntry_a[(ci,ck)][I]]=sqrt(q[('a',ck,ck)][eN,k,self.diagEntry_a[(ci,ck)][I]])



    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cfemIntegrals.copyExteriorElementBoundaryValuesFromElementBoundaryValues(self.mesh.exteriorElementBoundariesArray,
                                                                                     self.mesh.elementBoundaryElementsArray,
                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                     ebq[('u',ci)],
                                                                                     self.ebqe[('u',ci)])
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind computing bcs ebNE=%d k=%d g=%s" % (ebNE,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        #calculate diffusive velocity projection--requires interior and exterior loops
        self.calculateDiffusionSplittings(q,ebq,ebq_global)
        for ck in self.vt.coefficients.potential.keys():
            assert (len(self.vt.coefficients.potential[ck].keys()) == 1 and
                    self.vt.coefficients.potential[ck].keys()[0] == ck), "No off-diagonal dependence in phi currently allowed in LDG"
            if self.vt.coefficients.potential[ck][ck] == 'u':
                self.ebqe[('phi',ck)].flat[:]=self.ebqe[('u',ck)].flat
                self.ebqe[('dphi',ck,ck)].flat[:]=1.0
            cfemIntegrals.calculateInteriorNumericalTrace_Potential(self.mesh.interiorElementBoundariesArray,
                                                                    self.mesh.elementBoundaryElementsArray,
                                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                    ebq[('phi',ck)],
                                                                    ebq[('dphi',ck,ck)],
                                                                    self.phi_trace[ck],
                                                                    self.dphi_trace_left[ck],
                                                                    self.dphi_trace_right[ck])
            cfemIntegrals.calculateExteriorNumericalTrace_Potential(self.isDOFBoundary[ck],
                                                                    self.mesh.exteriorElementBoundariesArray,
                                                                    self.mesh.elementBoundaryElementsArray,
                                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                    self.ebqe[('phi',ck)],
                                                                    ebq[('phi',ck)],
                                                                    ebq[('dphi',ck,ck)],
                                                                    self.phi_trace[ck],
                                                                    self.dphi_trace_left[ck])
            #initialization
            self.b[ck].fill(0.0)
            #updates
            cfemIntegrals.updateInteriorElementBoundary_MixedForm_weak(self.mesh.interiorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       ebq['n'],
                                                                       self.phi_trace[ck],
                                                                       ebq[('w*dS_a',ck,ck)],
                                                                       self.b[ck])
            cfemIntegrals.updateExteriorElementBoundary_MixedForm_weak(self.mesh.exteriorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       ebq['n'],
                                                                       self.phi_trace[ck],
                                                                       ebq[('w*dS_a',ck,ck)],
                                                                       self.b[ck])
            cfemIntegrals.updatePotential_MixedForm_weak(q[('phi',ck)],
                                                         q[('grad(w)*dV_a',ck,ck)],
                                                         self.b[ck])
            self.useMixedForm1=False
            if self.useMixedForm1:
                cfemIntegrals.calculateVelocityQuadrature_MixedForm(self.A_inv,
                                                                    self.b[ck],
                                                                    ebq[('v',ck)],
                                                                    self.ebqV[ck],
                                                                    q[('v',ck)],
                                                                    self.qV[ck])
                self.eb_aTilde[(ck,ck)]=ebq[('a',ck,ck)]
                self.aTilde[(ck,ck)]=q[('a',ck,ck)]
            else:
                if self.vt.sd == True:
                    cfemIntegrals.calculateVelocityQuadrature_MixedForm2_sd(self.vt.coefficients.sdInfo[(ck,ck)][0],
                                                                            self.vt.coefficients.sdInfo[(ck,ck)][1],
                                                                            self.aHat[(ck,ck)],#q[('a',ck,ck)],
                                                                            q[('w*dV_a',ck,ck)],
                                                                            self.b[ck],
                                                                            ebq[('v',ck)],
                                                                            self.ebqV[ck],
                                                                            q[('v',ck)],
                                                                            self.qV[ck])
                else:
                    cfemIntegrals.calculateVelocityQuadrature_MixedForm2(self.aHat[(ck,ck)],#q[('a',ck,ck)],
                                                                         q[('w*dV_a',ck,ck)],
                                                                         self.b[ck],
                                                                         ebq[('v',ck)],
                                                                         self.ebqV[ck],
                                                                         q[('v',ck)],
                                                                         self.qV[ck])

        for ci in range(self.nc):
            ##\todo fix numerical diffusive flux for multiple diffusion terms
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    if self.vt.sd == True:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux_LDG_upwind_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                             self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                             self.mesh.interiorElementBoundariesArray,
                                                                                             self.mesh.elementBoundaryElementsArray,
                                                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                             ebq['n'],
                                                                                             ebq[('u',ck)],
                                                                                             self.eb_aTilde[(ci,ck)],
                                                                                             #ebq[('a',ci,ck)],
                                                                                             ebq[('phi',ck)],
                                                                                             self.ebqV[ck],
                                                                                             ebq_global[('penalty')],
                                                                                             ebq_global[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux_LDG_upwind(self.mesh.interiorElementBoundariesArray,
                                                                                          self.mesh.elementBoundaryElementsArray,
                                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                          ebq['n'],
                                                                                          ebq[('u',ck)],
                                                                                          self.eb_aTilde[(ci,ck)],
                                                                                          #ebq[('a',ci,ck)],
                                                                                          ebq[('phi',ck)],
                                                                                          self.ebqV[ck],
                                                                                          ebq_global[('penalty')],
                                                                                          ebq_global[('diffusiveFlux',ck,ci)])
        ####
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_LDG_upwind_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                             self.mesh.exteriorElementBoundariesArray,
                                                                                             self.mesh.elementBoundaryElementsArray,
                                                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                             ebqe['n'],
                                                                                             ebqe[('u',ck)],
                                                                                             self.eb_aTilde[(ci,ck)],
                                                                                             #ebqe[('a',ci,ck)],
                                                                                             self.ebqe[('phi',ck)],
                                                                                             ebqe[('phi',ck)],
                                                                                             self.ebqV[ck],
                                                                                             ebqe[('penalty')],
                                                                                             ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_LDG_upwind(self.mesh.exteriorElementBoundariesArray,
                                                                                          self.mesh.elementBoundaryElementsArray,
                                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                          ebqe['n'],
                                                                                          ebqe[('u',ck)],
                                                                                          self.eb_aTilde[(ci,ck)],
                                                                                          #ebqe[('a',ci,ck)],
                                                                                          self.ebqe[('phi',ck)],
                                                                                          ebqe[('phi',ck)],
                                                                                          self.ebqV[ck],
                                                                                          ebqe[('penalty')],
                                                                                          ebqe[('diffusiveFlux',ck,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        import pdb
        for ck in self.vt.coefficients.potential.keys():
            #initializations
            self.db[ck].fill(0.0)
            self.db_eb[ck].fill(0.0)
            #updates
            cfemIntegrals.updateInteriorElementBoundary_MixedForm_weakJacobian(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               self.dphi_trace_left[ck],
                                                                               self.dphi_trace_right[ck],
                                                                               ebq[('v',ck)],
                                                                               ebq[('w*dS_a',ck,ck)],
                                                                               self.db[ck],
                                                                               self.db_eb[ck])
            cfemIntegrals.updateExteriorElementBoundary_MixedForm_weakJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               self.dphi_trace_left[ck],
                                                                               ebq[('v',ck)],
                                                                               ebq[('w*dS_a',ck,ck)],
                                                                               self.db[ck],
                                                                               self.db_eb[ck])
            cfemIntegrals.updatePotential_MixedForm_weakJacobian(q[('dphi',ck,ck)],
                                                                 q[('v',ck)],
                                                                 q[('grad(w)*dV_a',ck,ck)],
                                                                 self.db[ck])
            if self.useMixedForm1:
                cfemIntegrals.calculateVelocityQuadrature_MixedForm_Jacobian(self.A_inv,
                                                                             self.db[ck],
                                                                             self.db_eb[ck],
                                                                             ebq[('v',ck)],
                                                                             self.ebqDV[ck],
                                                                             self.ebqDV_eb[ck],
                                                                             q[('v',ck)],
                                                                             self.qDV[ck],
                                                                             self.qDV_eb[ck])
            else:
                if self.vt.sd == True:
                    cfemIntegrals.calculateVelocityQuadrature_MixedForm2_Jacobian_sd(self.vt.coefficients.sdInfo[(ck,ck)][0],
                                                                                     self.vt.coefficients.sdInfo[(ck,ck)][1],
                                                                                     self.aHat[(ck,ck)],
                                                                                     q[('w*dV_a',ck,ck)],
                                                                                     self.db[ck],
                                                                                     self.db_eb[ck],
                                                                                     ebq[('v',ck)],
                                                                                     self.ebqDV[ck],
                                                                                     self.ebqDV_eb[ck],
                                                                                     q[('v',ck)],
                                                                                     self.qDV[ck],
                                                                                     self.qDV_eb[ck])

                else:
                    cfemIntegrals.calculateVelocityQuadrature_MixedForm2_Jacobian(self.aHat[(ck,ck)],
                                                                                  q[('w*dV_a',ck,ck)],
                                                                                  self.db[ck],
                                                                                  self.db_eb[ck],
                                                                                  ebq[('v',ck)],
                                                                                  self.ebqDV[ck],
                                                                                  self.ebqDV_eb[ck],
                                                                                  q[('v',ck)],
                                                                                  self.qDV[ck],
                                                                                  self.qDV_eb[ck])

        for ci in range(self.nc):
            #ignore advection contribution for right now
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    for cj in range(self.nc):
                        if dphi.has_key((ck,cj)):
                            if self.vt.sd:
                                cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                          self.mesh.interiorElementBoundariesArray,
                                                                                                          self.mesh.elementBoundaryElementsArray,
                                                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                          ebq['n'],
                                                                                                          self.eb_aTilde[(ci,ck)],
                                                                                                          #ebq[('a',ci,ck)],
                                                                                                          ebq[('da',ci,ck,cj)],
                                                                                                          ebq[('dphi',ck,cj)],
                                                                                                          self.ebqV[ck],
                                                                                                          self.ebqDV[ck],
                                                                                                          self.ebqDV_eb[ck],
                                                                                                          ebq[('v',cj)],
                                                                                                          ebq_global[('penalty')],
                                                                                                          fluxJacobian[ci][cj],
                                                                                                          fluxJacobian_eb[ci][cj])
                            else:
                                cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind(self.mesh.interiorElementBoundariesArray,
                                                                                                       self.mesh.elementBoundaryElementsArray,
                                                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                       ebq['n'],
                                                                                                       self.eb_aTilde[(ci,ck)],
                                                                                                       #ebq[('a',ci,ck)],
                                                                                                       ebq[('da',ci,ck,cj)],
                                                                                                       ebq[('dphi',ck,cj)],
                                                                                                       self.ebqV[ck],
                                                                                                       self.ebqDV[ck],
                                                                                                       self.ebqDV_eb[ck],
                                                                                                       ebq[('v',cj)],
                                                                                                       ebq_global[('penalty')],
                                                                                                       fluxJacobian[ci][cj],
                                                                                                       fluxJacobian_eb[ci][cj])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            #ignore advection contribution for now
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                                              self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                              self.isDiffusiveFluxBoundary[ci],
                                                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                                                              self.mesh.elementBoundaryElementsArray,
                                                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                              ebqe['n'],
                                                                                                              self.eb_aTilde[(ci,ck)],
                                                                                                              #ebqe[('a',ci,ck)],
                                                                                                              ebqe[('da',ci,ck,cj)],
                                                                                                              ebqe[('dphi',ck,cj)],
                                                                                                              self.ebqV[ck],
                                                                                                              self.ebqDV[ck],
                                                                                                              self.ebqDV_eb[ck],
                                                                                                              ebqe[('v',cj)],
                                                                                                              ebqe[('penalty')],
                                                                                                              fluxJacobian_exterior[ci][cj],
                                                                                                              fluxJacobian_eb[ci][cj])
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind(self.isDiffusiveFluxBoundary[ci],
                                                                                                           self.mesh.exteriorElementBoundariesArray,
                                                                                                           self.mesh.elementBoundaryElementsArray,
                                                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                           ebqe['n'],
                                                                                                           self.eb_aTilde[(ci,ck)],
                                                                                                           #ebqe[('a',ci,ck)],
                                                                                                           ebqe[('da',ci,ck,cj)],
                                                                                                           ebqe[('dphi',ck,cj)],
                                                                                                           self.ebqV[ck],
                                                                                                           self.ebqDV[ck],
                                                                                                           self.ebqDV_eb[ck],
                                                                                                           ebqe[('v',cj)],
                                                                                                           ebqe[('penalty')],
                                                                                                           fluxJacobian_exterior[ci][cj],
                                                                                                           fluxJacobian_eb[ci][cj])



class Advection_DiagonalUpwind_Diffusion_LDG(Diffusion_LDG):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Diffusion_LDG.__init__(self,vt,getPointwiseBoundaryConditions,
                               getAdvectiveFluxBoundaryConditions,
                               getDiffusiveFluxBoundaryConditions)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        Diffusion_LDG.calculateInteriorNumericalFlux(self,q,ebq,ebq_global)
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateInteriorNumericalAdvectiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                      self.mesh.elementBoundaryElementsArray,
                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                      ebq['n'],
                                                                      ebq[('u',ci)],
                                                                      ebq[('f',ci)],
                                                                      ebq[('df',ci,ci)],
                                                                      ebq_global[('advectiveFlux',ci)],
                                                                      ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                      ebq_global[('dadvectiveFlux_right',ci,ci)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        Diffusion_LDG.calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe)
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                cnumericalFlux.calculateExteriorNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       self.isDOFBoundary[ci],
                                                                       inflowFlag[ci],
                                                                       ebqe['n'],
                                                                       self.ebqe[('u',ci)],
                                                                       self.ebqe[('f',ci)],
                                                                       self.ebqe[('df',ci,ci)],
                                                                       ebqe[('u',ci)],
                                                                       ebqe[('f',ci)],
                                                                       ebqe[('df',ci,ci)],
                                                                       ebqe[('advectiveFlux',ci)],
                                                                       ebqe[('dadvectiveFlux_left',ci,ci)])
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                   inflowFlag[ci],
                                                                                   ebqe[('inflowFlux',ci)],
                                                                                   ebqe['n'],
                                                                                   ebqe[('f',ci)],
                                                                                   ebqe[('df',ci,ci)],
                                                                                   ebqe[('advectiveFlux',ci)],
                                                                                   ebqe[('dadvectiveFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        import pdb
        Diffusion_LDG.updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj)
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                if ebq.has_key(('f',ci)):
                    cnumericalFlux.updateInteriorNumericalAdvectiveFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                               ebq_global[('dadvectiveFlux_right',ci,ci)],
                                                                               ebq[('v',ci)],
                                                                               fluxJacobian[ci][ci])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        Diffusion_LDG.updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj)
        for ci in range(self.nc):
            if self.vt.timeIntegration.advectionIsImplicit[ci]:
                if ebqe.has_key(('f',ci)):
                    cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                inflowFlag[ci],
                                                                                ebqe[('dadvectiveFlux_left',ci,ci)],
                                                                                ebqe[('v',ci)],
                                                                                fluxJacobian_exterior[ci][ci])
#mwf add some scalar numerical fluxes
class RusanovNumericalFlux_Diagonal(Advection_DiagonalUpwind):
    r"""
    apply numerical flux :math:`f_{num}(a,b) = 1/2(f(a)+f(b)-\bar{\lambda}(b-a)` where
    :math:`\lambda >= max |f^{\prime}| for a<= u <= b`

    this one applies flux to each component of flux separately
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicConditions):
        Advection_DiagonalUpwind.__init__(self,vt,
                                          getPointwiseBoundaryConditions,
                                          getAdvectiveFluxBoundaryConditions,
                                          getDiffusiveFluxBoundaryConditions,
                                          getPeriodicConditions)
        self.safetyFactor=1.1
        #add extra terms that can be lagged specifically for advective flux. Time integrator has to do this though
        for ci in range(self.nc):
            #ebq
            vt.ebq[('u_advectiveNumericalFlux',ci)]= vt.ebq[('u',ci)]
            if vt.ebq.has_key(('f',ci)):
                vt.ebq[('f_advectiveNumericalFlux',ci)]= vt.ebq[('f',ci)]
            for cj in range(self.nc):
                if vt.ebq.has_key(('df',ci,cj)):
                    vt.ebq[('df_advectiveNumericalFlux',ci,cj)]= vt.ebq[('df',ci,cj)]
                if vt.q.has_key(('df',ci,cj)):
                    vt.q[('df_advectiveNumericalFlux',ci,cj)] = vt.q[('df',ci,cj)]
            #ebqe
            vt.ebqe[('u_advectiveNumericalFlux',ci)]= vt.ebqe[('u',ci)]
            if vt.ebqe.has_key(('f',ci)):
                vt.ebqe[('f_advectiveNumericalFlux',ci)]= vt.ebqe[('f',ci)]
            for cj in range(self.nc):
                if vt.ebqe.has_key(('df',ci,cj)):
                    vt.ebqe[('df_advectiveNumericalFlux',ci,cj)]= vt.ebqe[('df',ci,cj)]

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                          self.mesh.interiorElementBoundariesArray,
                                                                          self.mesh.elementBoundaryElementsArray,
                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                          ebq['n'],
                                                                          ebq[('u_advectiveNumericalFlux',ci)],
                                                                          ebq[('f_advectiveNumericalFlux',ci)],
                                                                          ebq[('df_advectiveNumericalFlux',ci,ci)],
                                                                          q[('df_advectiveNumericalFlux',ci,ci)],
                                                                          ebq_global[('advectiveFlux',ci)],
                                                                          ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                          ebq_global[('dadvectiveFlux_right',ci,ci)])

    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Rusanove_DiagonalUpwind computing bcs ebNE=%d k=%d x=%s g=%s " % (ebNE,k,x,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                          self.mesh.exteriorElementBoundariesArray,
                                                                          self.mesh.elementBoundaryElementsArray,
                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                          self.isDOFBoundary[ci],
                                                                          inflowFlag[ci],
                                                                          ebqe['n'],
                                                                          self.ebqe[('u',ci)],
                                                                          self.ebqe[('f',ci)],
                                                                          self.ebqe[('df',ci,ci)],
                                                                          ebqe[('u_advectiveNumericalFlux',ci)],
                                                                          ebqe[('f_advectiveNumericalFlux',ci)],
                                                                          ebqe[('df_advectiveNumericalFlux',ci,ci)],
                                                                          q[('df_advectiveNumericalFlux',ci,ci)],
                                                                          ebqe[('advectiveFlux',ci)],
                                                                          ebqe[('dadvectiveFlux_left',ci,ci)])
            #mwf debug
            #import pdb
            #pdb.set_trace()
        #mwf add for inflow flux?
        for ci in range(self.nc):
            cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               inflowFlag[ci],
                                                                               ebqe[('inflowFlux',ci)],
                                                                               ebqe['n'],
                                                                               ebqe[('f_advectiveNumericalFlux',ci)],
                                                                               ebqe[('df_advectiveNumericalFlux',ci,ci)],
                                                                               ebqe[('advectiveFlux',ci)],
                                                                               ebqe[('dadvectiveFlux_left',ci,ci)])
class RusanovNumericalFlux_Diagonal_Diffusion_IIPG(Advection_DiagonalUpwind_Diffusion_IIPG):
    r"""
    apply numerical flus :math:`f_{num}(a,b) = 1/2(f(a)+f(b)-\bar{\lambda}(b-a)` where
    :math:`\lambda >= max |f^{\prime}|` for :math:`a<= u <= b`
    this one applies flux to each component of flux separately
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Advection_DiagonalUpwind_Diffusion_IIPG.__init__(self,vt,getPointwiseBoundaryConditions,
                                                         getAdvectiveFluxBoundaryConditions,
                                                         getDiffusiveFluxBoundaryConditions)
        self.safetyFactor=1.1
        #add extra terms that can be lagged specifically for advective flux. Time integrator has to do this though
        for ci in range(self.nc):
            #ebq
            vt.ebq[('u_advectiveNumericalFlux',ci)]= vt.ebq[('u',ci)]
            if vt.ebq.has_key(('f',ci)):
                vt.ebq[('f_advectiveNumericalFlux',ci)]= vt.ebq[('f',ci)]
            for cj in range(self.nc):
                if vt.ebq.has_key(('df',ci,cj)):
                    vt.ebq[('df_advectiveNumericalFlux',ci,cj)]= vt.ebq[('df',ci,cj)]
                if vt.q.has_key(('df',ci,cj)):
                    vt.q[('df_advectiveNumericalFlux',ci,cj)] = vt.q[('df',ci,cj)]
            #ebqe
            vt.ebqe[('u_advectiveNumericalFlux',ci)]= vt.ebqe[('u',ci)]
            if vt.ebqe.has_key(('f',ci)):
                vt.ebqe[('f_advectiveNumericalFlux',ci)]= vt.ebqe[('f',ci)]
            for cj in range(self.nc):
                if vt.ebqe.has_key(('df',ci,cj)):
                    vt.ebqe[('df_advectiveNumericalFlux',ci,cj)]= vt.ebqe[('df',ci,cj)]

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateInteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                              self.mesh.interiorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              ebq['n'],
                                                                              ebq[('u_advectiveNumericalFlux',ci)],
                                                                              ebq[('f_advectiveNumericalFlux',ci)],
                                                                              ebq[('df_advectiveNumericalFlux',ci,ci)],
                                                                              q[('df_advectiveNumericalFlux',ci,ci)],
                                                                              ebq_global[('advectiveFlux',ci)],
                                                                              ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                              ebq_global[('dadvectiveFlux_right',ci,ci)])
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.interiorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  ebq['n'],
                                                                                  ebq[('a',ci,ck)],
                                                                                  ebq[('grad(phi)',ck)],
                                                                                  ebq[('u',ck)],
                                                                                  ebq_global[('penalty')],
                                                                                  ebq_global[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               ebq[('a',ci,ck)],
                                                                               ebq[('grad(phi)',ck)],
                                                                               ebq[('u',ck)],
                                                                               ebq_global[('penalty')],
                                                                               ebq_global[('diffusiveFlux',ck,ci)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Rusanove_DiagonalUpwind computing bcs ebNE=%d k=%d x=%s g=%s " % (ebNE,k,x,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                cnumericalFlux.calculateExteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              self.isDOFBoundary[ci],
                                                                              inflowFlag[ci],
                                                                              ebqe['n'],
                                                                              self.ebqe[('u',ci)],
                                                                              self.ebqe[('f',ci)],
                                                                              self.ebqe[('df',ci,ci)],
                                                                              ebqe[('u_advectiveNumericalFlux',ci)],
                                                                              ebqe[('f_advectiveNumericalFlux',ci)],
                                                                              ebqe[('df_advectiveNumericalFlux',ci,ci)],
                                                                              q[('df_advectiveNumericalFlux',ci,ci)],
                                                                              ebqe[('advectiveFlux',ci)],
                                                                              ebqe[('dadvectiveFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)])
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)])
        #mwf add for inflow flux? needs to be removed
        for ci in range(self.nc):
            cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               inflowFlag[ci],
                                                                               ebqe[('inflowFlux',ci)],
                                                                               ebqe['n'],
                                                                               ebqe[('f_advectiveNumericalFlux',ci)],
                                                                               ebqe[('df_advectiveNumericalFlux',ci,ci)],
                                                                               ebqe[('advectiveFlux',ci)],
                                                                               ebqe[('dadvectiveFlux_left',ci,ci)])
class ConvexOneSonicPointNumericalFlux(Advection_DiagonalUpwind):
    r"""
    basic Godunov flux :math:`f_{num}(a,b) = max_{b<= u <= a} f(u)` if :math:`a >= b = min_{a<= u <= b} f(u)` otherwise
    where there is only one sonic point, :math:`u_s` with :math:`f^{\prime}(u_s) = 0` and
    :math:`f` is convex so :math:`f(u_s)` is a minimum

    This class typically has to be "wrapped" for a given problem to specify the
    correct sonic point and sonic flux
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions,
                 sonicPoint=0.0,sonicFlux=0.0):
        Advection_DiagonalUpwind.__init__(self,vt,getPointwiseBoundaryConditions,
                                          getAdvectiveFluxBoundaryConditions,
                                          getDiffusiveFluxBoundaryConditions,
                                          getPeriodicBoundaryConditions)
        self.sonicPoint = sonicPoint
        self.sonicFlux  = sonicFlux
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorNumericalAdvectiveFluxConvexOneSonicPoint(self.sonicPoint,
                                                                                      self.sonicFlux,
                                                                                      self.mesh.interiorElementBoundariesArray,
                                                                                      self.mesh.elementBoundaryElementsArray,
                                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                      ebq['n'],
                                                                                      ebq[('u',ci)],
                                                                                      ebq[('f',ci)],
                                                                                      ebq[('df',ci,ci)],
                                                                                      ebq_global[('advectiveFlux',ci)],
                                                                                      ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                                      ebq_global[('dadvectiveFlux_right',ci,ci)])




class HamiltonJacobi_DiagonalLesaintRaviart(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 speedEvaluationType=2):#mwf hack
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.speedEvaluationType = int(speedEvaluationType)#1 -- use max, otherwise allow disc.
        #self.outFlowOnly=True
        for ci in range(self.nc):
            #by default advection-diffusion
            self.advectiveNumericalFlux[ci] = False
            self.diffusiveNumericalFlux[ci] = False
            self.HamiltonJacobiNumericalFlux[ci] = True
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.interiorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        ebq['n'],
                                                                        ebq[('u',ci)],
                                                                        ebq[('H',ci)],
                                                                        ebq[('dH',ci,ci)],
                                                                        ebq[('HamiltonJacobiFlux',ci)],
                                                                        ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                        ebq[('dHamiltonJacobiFlux_right',ci,ci)])

        #mwf debug
        #import pdb
        #pdb.set_trace()
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                #print "bc--------------------",bci
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('H',ci)],
                                                                        self.ebqe[('dH',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('H',ci)],
                                                                        ebqe[('dH',ci,ci)],
                                                                        ebqe[('HamiltonJacobiFlux',ci)],
                                                                        ebqe[('dHamiltonJacobiFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateInteriorTwoSidedNumericalFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                           ebq[('dHamiltonJacobiFlux_right',ci,ci)],
                                                                           ebq[('v',ci)],
                                                                           fluxJacobian_hj[ci][ci])
                #mwf debug
                #import pdb
                #pdb.set_trace()

    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])




class HamiltonJacobi_DiagonalLesaintRaviart_Diffusion_IIPG(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 speedEvaluationType=1):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.speedEvaluationType = int(speedEvaluationType)#1 -- use max, otherwise allow disc.
        #self.outFlowOnly=True
        for ci in range(self.nc):
            self.advectiveNumericalFlux[ci] = False
            self.diffusiveNumericalFlux[ci] = True
            self.HamiltonJacobiNumericalFlux[ci] = True
         #
        self.scale_penalty = 1; self.penalty_floor = 0.0
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.interiorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        ebq['n'],
                                                                        ebq[('u',ci)],
                                                                        ebq[('H',ci)],
                                                                        ebq[('dH',ci,ci)],
                                                                        ebq[('HamiltonJacobiFlux',ci)],
                                                                        ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                        ebq[('dHamiltonJacobiFlux_right',ci,ci)])

            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                               self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               ebq[('a',ci,ck)],
                                                                               ebq[('grad(phi)',ck)],
                                                                               ebq[('u',ck)],
                                                                               ebq_global[('penalty')],
                                                                               ebq_global[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
                    else:
                        cnumericalFlux.calculateInteriorNumericalDiffusiveFlux(self.mesh.interiorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               ebq['n'],
                                                                               ebq[('a',ci,ck)],
                                                                               ebq[('grad(phi)',ck)],
                                                                               ebq[('u',ck)],
                                                                               ebq_global[('penalty')],
                                                                               ebq_global[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)
        #mwf debug
        #import pdb
        #pdb.set_trace()
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                #print "bc--------------------",bci
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('H',ci)],
                                                                        self.ebqe[('dH',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('H',ci)],
                                                                        ebqe[('dH',ci,ci)],
                                                                        ebqe[('HamiltonJacobiFlux',ci)],
                                                                        ebqe[('dHamiltonJacobiFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)],
                                                                                  self.scale_penalty,
                                                                                  self.penalty_floor)
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)

    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateInteriorTwoSidedNumericalFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                           ebq[('dHamiltonJacobiFlux_right',ci,ci)],
                                                                           ebq[('v',ci)],
                                                                           fluxJacobian_hj[ci][ci])
            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebq.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.interiorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   ebq['n'],
                                                                                                   ebq[('a',ci,ck)],
                                                                                                   ebq[('da',ci,ck,cj)],
                                                                                                   ebq[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebq[('v',cj)],
                                                                                                   ebq[('grad(v)',cj)],
                                                                                                   ebq_global['penalty'],
                                                                                                   fluxJacobian[ci][cj],
                                                                                                   self.scale_penalty,
                                                                                                   self.penalty_floor)
                                else:
                                    cnumericalFlux.updateInteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.interiorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                ebq['n'],
                                                                                                ebq[('a',ci,ck)],
                                                                                                ebq[('da',ci,ck,cj)],
                                                                                                ebq[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebq[('v',cj)],
                                                                                                ebq[('grad(v)',cj)],
                                                                                                ebq_global['penalty'],
                                                                                                fluxJacobian[ci][cj],
                                                                                                self.scale_penalty,
                                                                                                self.penalty_floor)
                #mwf debug
                #import pdb
                #pdb.set_trace()

    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])

            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj],
                                                                                                   self.scale_penalty,
                                                                                                   self.penalty_floor)
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj],
                                                                                                self.scale_penalty,
                                                                                                self.penalty_floor)

    #        if self.vt.timeIntegration.diffusionIsImplicit[ci]:
     #           for ck in range(self.nc):
      #              if ebqe.has_key(('a',ci,ck)):
       #                 for cj in range(self.nc):
        #                    if dphi.has_key((ck,cj)):
         #                       cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
          #                                                                                  self.mesh.exteriorElementBoundariesArray,
            #                                                                                 self.mesh.elementBoundaryElementsArray,
            #                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
             #                                                                               self.isDOFBoundary[ck],
              #                                                                              ebqe['n'],
                #                                                                             ebqe[('a',ci,ck)],
                #                                                                            ebqe[('da',ci,ck,cj)],
                 #                                                                           ebqe[('grad(phi)',ck)],
                  #                                                                          dphi[(ck,cj)].dof,
                    #                                                                         ebqe[('v',cj)],
                    #                                                                        ebqe[('grad(v)',cj)],
                     #                                                                       ebqe['penalty'],
                      #                                                                      fluxJacobian_exterior[ci][cj],
                       #                                                                     self.scale_penalty,
                        #                                                                    self.penalty_floor)

class HamiltonJacobi_DiagonalLesaintRaviart_Diffusion_SIPG_exterior(Diffusion_SIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 speedEvaluationType=1):
        Diffusion_SIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.speedEvaluationType = int(speedEvaluationType)#1 -- use max, otherwise allow disc.
        for ci in range(self.nc):
            self.advectiveNumericalFlux[ci] = False
            self.diffusiveNumericalFlux[ci] = True
            self.HamiltonJacobiNumericalFlux[ci] = True
        self.scale_penalty = 1
        self.penalty_floor = 0.0
        self.penalty_constant = 10.0
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)

    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('H',ci)],
                                                                        self.ebqe[('dH',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('H',ci)],
                                                                        ebqe[('dH',ci,ci)],
                                                                        ebqe[('HamiltonJacobiFlux',ci)],
                                                                        ebqe[('dHamiltonJacobiFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)],
                                                                                  self.scale_penalty,
                                                                                  self.penalty_floor)
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)

    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])

            if self.vt.timeIntegration.diffusionIsImplicit[ci]:
                for ck in range(self.nc):
                    if ebqe.has_key(('a',ci,ck)):
                        for cj in range(self.nc):
                            if dphi.has_key((ck,cj)):
                                if self.vt.sd:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][1],
                                                                                                   dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[ck],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',ci,ck)],
                                                                                                   ebqe[('da',ci,ck,cj)],
                                                                                                   ebqe[('grad(phi)',ck)],
                                                                                                   dphi[(ck,cj)].dof,
                                                                                                   ebqe[('v',cj)],
                                                                                                   ebqe[('grad(v)',cj)],
                                                                                                   ebqe['penalty'],
                                                                                                   fluxJacobian_exterior[ci][cj],
                                                                                                   self.scale_penalty,
                                                                                                   self.penalty_floor)
                                else:
                                    cnumericalFlux.updateExteriorNumericalDiffusiveFluxJacobian(dphi[(ck,cj)].femSpace.dofMap.l2g,
                                                                                                self.mesh.exteriorElementBoundariesArray,
                                                                                                self.mesh.elementBoundaryElementsArray,
                                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                self.isDOFBoundary[ck],
                                                                                                ebqe['n'],
                                                                                                ebqe[('a',ci,ck)],
                                                                                                ebqe[('da',ci,ck,cj)],
                                                                                                ebqe[('grad(phi)',ck)],
                                                                                                dphi[(ck,cj)].dof,
                                                                                                ebqe[('v',cj)],
                                                                                                ebqe[('grad(v)',cj)],
                                                                                                ebqe['penalty'],
                                                                                                fluxJacobian_exterior[ci][cj],
                                                                                                self.scale_penalty,
                                                                                                self.penalty_floor)

# If using this flux, it adds the term  + p*\n to the existing HamiltonJacobiFlux.  Because this is
# in a sense treating pressure (or pI) as a conservative advection velocity, but pressure is
# not really being advected, it is important not to use this flux in combination with
# an upwinding numericalFlux for the pressure variable.  Recommended ones are
# ConstantAdvection_exterior or others of that type for pressure.
#
class HamiltonJacobi_Pressure_DiagonalLesaintRaviart_Diffusion_SIPG_exterior(HamiltonJacobi_DiagonalLesaintRaviart_Diffusion_SIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 speedEvaluationType=1):
        HamiltonJacobi_DiagonalLesaintRaviart_Diffusion_SIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions,
                 speedEvaluationType)
        self.speedEvaluationType = int(speedEvaluationType)#1 -- use max, otherwise allow disc.
        for ci in range(self.nc):
            self.advectiveNumericalFlux[ci] = False
            self.diffusiveNumericalFlux[ci] = True
            self.HamiltonJacobiNumericalFlux[ci] = True
        self.scale_penalty = 1
        self.penalty_floor = 0.0
        self.penalty_constant = 10.0
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('H',ci)],
                                                                        self.ebqe[('dH',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('H',ci)],
                                                                        ebqe[('dH',ci,ci)],
                                                                        ebqe[('HamiltonJacobiFlux',ci)],
                                                                        ebqe[('dHamiltonJacobiFlux_left',ci,ci)])

            # we add the term that technically should be an advective flux but since
            # it is constant wrt variables here, we just add it to the HamiltonJacobiFlux
            #  < div(p I), w> = < -p I, grad w > + <pI n, grad w>_{\partial\Lambda}
            #
            # so our term is  pIn = p \n  and since 'f' = pI, we just innerproduct
            # with 'n' to get our term.
            ebqe[('HamiltonJacobiFlux',ci)][:] += (ebqe[('f',ci)]*ebqe['n']).sum(-1)   #  + p N from advection terms
            ebqe[('dHamiltonJacobiFlux_left',ci,ci)][:] += 0.0

            for ck in range(self.nc):
                if ebqe.has_key(('a',ci,ck)):
                    if self.vt.sd:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux_sd(self.vt.coefficients.sdInfo[(ci,ck)][0],self.vt.coefficients.sdInfo[(ci,ck)][0],
                                                                                  self.mesh.exteriorElementBoundariesArray,
                                                                                  self.mesh.elementBoundaryElementsArray,
                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                  self.isDOFBoundary[ck],
                                                                                  ebqe['n'],
                                                                                  self.ebqe[('a',ci,ck)],
                                                                                  self.ebqe[('grad(phi)',ck)],
                                                                                  self.ebqe[('u',ck)],
                                                                                  ebqe[('a',ci,ck)],
                                                                                  ebqe[('grad(phi)',ck)],
                                                                                  ebqe[('u',ck)],
                                                                                  ebqe[('penalty')],
                                                                                  ebqe[('diffusiveFlux',ck,ci)],
                                                                                  self.scale_penalty,
                                                                                  self.penalty_floor)
                    else:
                        cnumericalFlux.calculateExteriorNumericalDiffusiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                               self.mesh.elementBoundaryElementsArray,
                                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                               self.isDOFBoundary[ck],
                                                                               ebqe['n'],
                                                                               self.ebqe[('a',ci,ck)],
                                                                               self.ebqe[('grad(phi)',ck)],
                                                                               self.ebqe[('u',ck)],
                                                                               ebqe[('a',ci,ck)],
                                                                               ebqe[('grad(phi)',ck)],
                                                                               ebqe[('u',ck)],
                                                                               ebqe[('penalty')],
                                                                               ebqe[('diffusiveFlux',ck,ci)],
                                                                               self.scale_penalty,
                                                                               self.penalty_floor)


class DarcyFCFF_IIPG_exterior(NF_base):
    """
    weak dirichlet boundary conditions for Twophase_fc_ff class
    """
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Advection_DiagonalUpwind_ext computing bcs ebNE=%d k=%d g=%s" % (ebNE,k,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.sd:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCFF_sd(self.vt.coefficients.sdInfo[(0,1)][0],self.vt.coefficients.sdInfo[(0,1)][1],
                                                                            self.vt.coefficients.sdInfo[(0,1)][0],self.vt.coefficients.sdInfo[(1,0)][1],
                                                                            self.vt.coefficients.sdInfo[(0,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                            self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[0],
                                                                            self.isDOFBoundary[1],
                                                                            ebqe['n'],
                                                                            self.ebqe[('f',1)],
                                                                            self.ebqe[('a',0,1)],
                                                                            self.ebqe[('a',1,0)],
                                                                            self.ebqe[('a',1,1)],
                                                                            self.ebqe[('grad(phi)',0)],
                                                                            self.ebqe[('grad(phi)',1)],
                                                                            self.ebqe[('u',0)],
                                                                            self.ebqe[('u',1)],
                                                                            ebqe[('f',1)],
                                                                            ebqe[('df',1,0)],#no dependence on u_m unless compressible
                                                                            ebqe[('a',0,1)],
                                                                            ebqe[('a',1,0)],
                                                                            ebqe[('a',1,1)],
                                                                            ebqe[('grad(phi)',0)],
                                                                            ebqe[('grad(phi)',1)],
                                                                            ebqe[('u',0)],
                                                                            ebqe[('u',1)],
                                                                            ebqe[('penalty')],
                                                                            ebqe[('penalty')],
                                                                            ebqe[('advectiveFlux',1)],
                                                                            ebqe[('dadvectiveFlux_left',1,0)],#(ci,cj)
                                                                            ebqe[('diffusiveFlux',1,0)], #(ck,ci)
                                                                            ebqe[('diffusiveFlux',0,1)],
                                                                            ebqe[('diffusiveFlux',1,1)])
        else:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCFF(self.mesh.exteriorElementBoundariesArray,
                                                                         self.mesh.elementBoundaryElementsArray,
                                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                         self.isDOFBoundary[0],
                                                                         self.isDOFBoundary[1],
                                                                         ebqe['n'],
                                                                         self.ebqe[('f',1)],
                                                                         self.ebqe[('a',0,1)],
                                                                         self.ebqe[('a',1,0)],
                                                                         self.ebqe[('a',1,1)],
                                                                         self.ebqe[('grad(phi)',0)],
                                                                         self.ebqe[('grad(phi)',1)],
                                                                         self.ebqe[('u',0)],
                                                                         self.ebqe[('u',1)],
                                                                         ebqe[('f',1)],
                                                                         ebqe[('df',1,0)],#no dependence on u_m unless compressible
                                                                         ebqe[('a',0,1)],
                                                                         ebqe[('a',1,0)],
                                                                         ebqe[('a',1,1)],
                                                                         ebqe[('grad(phi)',0)],
                                                                         ebqe[('grad(phi)',1)],
                                                                         ebqe[('u',0)],
                                                                         ebqe[('u',1)],
                                                                         ebqe[('penalty')],
                                                                         ebqe[('penalty')],
                                                                         ebqe[('advectiveFlux',1)],
                                                                         ebqe[('dadvectiveFlux_left',1,0)],#(ci,cj)
                                                                         ebqe[('diffusiveFlux',1,0)], #(ck,ci)
                                                                         ebqe[('diffusiveFlux',0,1)],
                                                                         ebqe[('diffusiveFlux',1,1)])


    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        ci = 1;
        if self.vt.timeIntegration.advectionIsImplicit[ci]:
            for cj in [0]:#no dependence on u_m unless compressible
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],#not evaluated right now
                                                                            ebqe[('dadvectiveFlux_left',ci,cj)],
                                                                            ebqe[('v',cj)],
                                                                            fluxJacobian_exterior[ci][cj])
        if (self.vt.timeIntegration.diffusionIsImplicit[0] or self.vt.timeIntegration.diffusionIsImplicit[1]):
            #assume the u_w and u_m spaces are the same right now
            l2g = dphi[(1,1)].femSpace.dofMap.l2g; v   = ebqe[('v',1)]; grad_v = ebqe[('grad(v)',1)];
            if self.vt.sd:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(0,1)][0],self.vt.coefficients.sdInfo[(0,1)][1],
                                                                                                      self.vt.coefficients.sdInfo[(1,0)][0],self.vt.coefficients.sdInfo[(1,0)][1],
                                                                                                      self.vt.coefficients.sdInfo[(1,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                                                      l2g,
                                                                                                      self.mesh.exteriorElementBoundariesArray,
                                                                                                      self.mesh.elementBoundaryElementsArray,
                                                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                      self.isDOFBoundary[0],
                                                                                                      self.isDOFBoundary[1],
                                                                                                      ebqe['n'],
                                                                                                      ebqe[('f',1)],
                                                                                                      ebqe[('df',1,0)],
                                                                                                      ebqe[('a',0,1)],
                                                                                                      ebqe[('da',0,1,0)],
                                                                                                      ebqe[('da',0,1,1)],
                                                                                                      ebqe[('a',1,0)],
                                                                                                      ebqe[('da',1,0,0)],
                                                                                                      ebqe[('da',1,0,1)],
                                                                                                      ebqe[('a',1,1)],
                                                                                                      ebqe[('da',1,1,0)],
                                                                                                      ebqe[('da',1,1,1)],
                                                                                                      ebqe[('grad(phi)',0)],
                                                                                                      ebqe[('grad(phi)',1)],
                                                                                                      dphi[(0,0)].dof,
                                                                                                      dphi[(0,1)].dof,
                                                                                                      dphi[(1,0)].dof,
                                                                                                      dphi[(1,1)].dof,
                                                                                                      ebqe[('u',0)],
                                                                                                      ebqe[('u',1)],
                                                                                                      v,
                                                                                                      grad_v,
                                                                                                      ebqe[('penalty')],
                                                                                                      ebqe[('penalty')],
                                                                                                      fluxJacobian_exterior[0][0],
                                                                                                      fluxJacobian_exterior[0][1],
                                                                                                      fluxJacobian_exterior[1][0],
                                                                                                      fluxJacobian_exterior[1][1])
            else:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian(l2g,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[0],
                                                                                                   self.isDOFBoundary[1],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('f',1)],
                                                                                                   ebqe[('df',1,0)],
                                                                                                   ebqe[('a',0,1)],
                                                                                                   ebqe[('da',0,1,0)],
                                                                                                   ebqe[('da',0,1,1)],
                                                                                                   ebqe[('a',1,0)],
                                                                                                   ebqe[('da',1,0,0)],
                                                                                                   ebqe[('da',1,0,1)],
                                                                                                   ebqe[('a',1,1)],
                                                                                                   ebqe[('da',1,1,0)],
                                                                                                   ebqe[('da',1,1,1)],
                                                                                                   ebqe[('grad(phi)',0)],
                                                                                                   ebqe[('grad(phi)',1)],
                                                                                                   dphi[(0,0)].dof,
                                                                                                   dphi[(0,1)].dof,
                                                                                                   dphi[(1,0)].dof,
                                                                                                   dphi[(1,1)].dof,
                                                                                                   ebqe[('u',0)],
                                                                                                   ebqe[('u',1)],
                                                                                                   v,
                                                                                                   grad_v,
                                                                                                   ebqe[('penalty')],
                                                                                                   ebqe[('penalty')],
                                                                                                   fluxJacobian_exterior[0][0],
                                                                                                   fluxJacobian_exterior[0][1],
                                                                                                   fluxJacobian_exterior[1][0],
                                                                                                   fluxJacobian_exterior[1][1])


class DarcyFC_IIPG_exterior(NF_base):
    r"""
    weak dirichlet boundary conditions for Twophase_fc class

    .. todo::
       - put in nonlinear bc for setting :math:`\psi_n = \psi_n^b`
       - put in bc that switches flux and dirichlet types
    """
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        #hold Dirichlet values for \psi_n (non-wetting phase head)
        #also need extra psi_n entries that aren't part of default quadrature
        for term in ['psi_n_bc','psi_n',('dpsi_n',0),('dpsi_n',1),'sw']:
            self.ebqe[term] = numpy.zeros(self.ebqe[('u',1)].shape,'d')
        self.hasInterior=False
        self.penalty_constant = 2.0
        self.penalty_power = 1.0
        self.fluxBoundaryFlags = {0:0,1:0}#default is no-flow
        for ci in range(2):
            if vt.fluxBoundaryConditions[ci] == 'outFlow':
                self.fluxBoundaryFlags[ci] = 1
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                if self.isDOFBoundary[ci][ebNE,k] == 2:
                    self.ebqe['psi_n_bc'][ebNE,k] = g(x,self.vt.timeIntegration.t)
                else:
                    self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
#         if ebqe.has_key(('f',0)) and ebqe.has_key(('f',1)):
#             cnumericalFlux.calculateGlobalExteriorNumericalAdvectiveFlux_DarcyFC(self.mesh.exteriorElementBoundariesArray,
#                                                                                  self.mesh.elementBoundaryElementsArray,
#                                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                  self.isDOFBoundary[0],
#                                                                                  self.isDOFBoundary[1],
#                                                                                  ebqe['n'],
#                                                                                  self.ebqe[('u',0)],
#                                                                                  self.ebqe[('u',1)],
#                                                                                  self.ebqe[('f',0)],
#                                                                                  self.ebqe[('df',0,0)],
#                                                                                  self.ebqe[('df',0,1)],
#                                                                                  self.ebqe[('f',1)],
#                                                                                  self.ebqe[('df',1,0)],
#                                                                                  self.ebqe[('df',1,1)],
#                                                                                  ebqe[('u',0)],
#                                                                                  ebqe[('u',1)],
#                                                                                  ebqe[('f',0)],
#                                                                                  ebqe[('df',0,0)],
#                                                                                  ebqe[('df',0,1)],
#                                                                                  ebqe[('f',1)],
#                                                                                  ebqe[('df',1,0)],
#                                                                                  ebqe[('df',1,1)],
#                                                                                  ebqe[('advectiveFlux',0)],
#                                                                                  ebqe[('dadvectiveFlux_left',0,0)],
#                                                                                  ebqe[('dadvectiveFlux_left',0,1)],
#                                                                                  ebqe[('advectiveFlux',1)],
#                                                                                  ebqe[('dadvectiveFlux_left',1,0)],
#                                                                                  ebqe[('dadvectiveFlux_left',1,1)])
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if self.vt.sd:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFC_sd(self.vt.coefficients.sdInfo[(0,0)][0],self.vt.coefficients.sdInfo[(0,0)][1],
                                                                          self.vt.coefficients.sdInfo[(1,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                          self.mesh.exteriorElementBoundariesArray,
                                                                          self.mesh.elementBoundaryElementsArray,
                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                          self.isDOFBoundary[0],
                                                                          self.isDOFBoundary[1],
                                                                          self.fluxBoundaryFlags[0],
                                                                          self.fluxBoundaryFlags[1],
                                                                          ebqe['n'],
                                                                          self.ebqe[('a',0,0)],
                                                                          self.ebqe[('a',1,1)],
                                                                          self.ebqe[('grad(phi)',0)],
                                                                          self.ebqe[('grad(phi)',1)],
                                                                          self.ebqe[('u',0)],
                                                                          self.ebqe[('u',1)],
                                                                          self.ebqe['psi_n_bc'],
                                                                          ebqe[('a',0,0)],
                                                                          ebqe[('a',1,1)],
                                                                          ebqe[('grad(phi)',0)],
                                                                          ebqe[('grad(phi)',1)],
                                                                          ebqe[('u',0)],
                                                                          ebqe[('u',1)],
                                                                          ebqe['psi_n'],
                                                                          ebqe[('penalty')],
                                                                          ebqe[('penalty')],
                                                                          ebqe[('diffusiveFlux',0,0)], #(ck,ci)
                                                                          ebqe[('diffusiveFlux',1,1)])
        else:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFC(self.mesh.exteriorElementBoundariesArray,
                                                                       self.mesh.elementBoundaryElementsArray,
                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                       self.isDOFBoundary[0],
                                                                       self.isDOFBoundary[1],
                                                                       self.fluxBoundaryFlags[0],
                                                                       self.fluxBoundaryFlags[1],
                                                                       ebqe['n'],
                                                                       self.ebqe[('a',0,0)],
                                                                       self.ebqe[('a',1,1)],
                                                                       self.ebqe[('grad(phi)',0)],
                                                                       self.ebqe[('grad(phi)',1)],
                                                                       self.ebqe[('u',0)],
                                                                       self.ebqe[('u',1)],
                                                                       self.ebqe['psi_n_bc'],
                                                                       ebqe[('a',0,0)],
                                                                       ebqe[('a',1,1)],
                                                                       ebqe[('grad(phi)',0)],
                                                                       ebqe[('grad(phi)',1)],
                                                                       ebqe[('u',0)],
                                                                       ebqe[('u',1)],
                                                                       ebqe['psi_n'],
                                                                       ebqe[('penalty')],
                                                                       ebqe[('penalty')],
                                                                       ebqe[('diffusiveFlux',0,0)], #(ck,ci)
                                                                       ebqe[('diffusiveFlux',1,1)])
        #print "DarcyFC exterior advectiveFlux[0]=%s \n advectiveFlux[1]=%s " % (ebqe[('advectiveFlux',0)],ebqe[('advectiveFlux',1)])
        #print "DarcyFC exterior diffusiveFlux[0]=%s \n diffusiveFlux[1]=%s " % (ebqe[('diffusiveFlux',0,0)],ebqe[('diffusiveFlux',1,1)])

    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        #mwf debug
        #import pdb
        #pdb.set_trace()
#         for ci in range(self.nc):
#              if self.vt.timeIntegration.advectionIsImplicit[ci]:
#                 for cj in range(self.nc):
#                     if ebqe.has_key(('df',ci,cj)):
#                         cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
#                                                                                     self.mesh.elementBoundaryElementsArray,
#                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                     inflowFlag[ci],#not used
#                                                                                     ebqe[('dadvectiveFlux_left',ci,cj)],
#                                                                                     ebqe[('v',cj)],
#                                                                                     fluxJacobian_exterior[ci][cj])
        #mwf debug
        #import pdb
        #pdb.set_trace()

        if (self.vt.timeIntegration.diffusionIsImplicit[0] or self.vt.timeIntegration.diffusionIsImplicit[1]):
            #assume the u_w and u_m spaces are the same right now
            l2g_local = dphi[(0,0)].femSpace.dofMap.l2g; v   = ebqe[('v',0)]; grad_v = ebqe[('grad(v)',0)];
            if self.vt.sd:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(0,0)][0],self.vt.coefficients.sdInfo[(0,0)][1],
                                                                                                    self.vt.coefficients.sdInfo[(1,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                                                    l2g_local,
                                                                                                    self.mesh.exteriorElementBoundariesArray,
                                                                                                    self.mesh.elementBoundaryElementsArray,
                                                                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                    self.isDOFBoundary[0],
                                                                                                    self.isDOFBoundary[1],
                                                                                                    self.fluxBoundaryFlags[0],
                                                                                                    self.fluxBoundaryFlags[1],
                                                                                                    ebqe['n'],
                                                                                                    ebqe[('a',0,0)],
                                                                                                    ebqe[('da',0,0,0)],
                                                                                                    ebqe[('da',0,0,1)],
                                                                                                    ebqe[('a',1,1)],
                                                                                                    ebqe[('da',1,1,0)],
                                                                                                    ebqe[('da',1,1,1)],
                                                                                                    ebqe[('grad(phi)',0)],
                                                                                                    ebqe[('grad(phi)',1)],
                                                                                                    dphi[(0,0)].dof,
                                                                                                    dphi[(0,1)].dof,
                                                                                                    dphi[(1,0)].dof,
                                                                                                    dphi[(1,1)].dof,
                                                                                                    ebqe[('u',0)],
                                                                                                    ebqe[('u',1)],
                                                                                                    ebqe['psi_n'],
                                                                                                    ebqe[('dpsi_n',0)],
                                                                                                    ebqe[('dpsi_n',1)],
                                                                                                    v,
                                                                                                    grad_v,
                                                                                                    ebqe[('penalty')],
                                                                                                    ebqe[('penalty')],
                                                                                                    fluxJacobian_exterior[0][0],
                                                                                                    fluxJacobian_exterior[0][1],
                                                                                                    fluxJacobian_exterior[1][0],
                                                                                                    fluxJacobian_exterior[1][1])
            else:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian(l2g_local,
                                                                                                 self.mesh.exteriorElementBoundariesArray,
                                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                 self.isDOFBoundary[0],
                                                                                                 self.isDOFBoundary[1],
                                                                                                 self.fluxBoundaryFlags[0],
                                                                                                 self.fluxBoundaryFlags[1],
                                                                                                 ebqe['n'],
                                                                                                 ebqe[('a',0,0)],
                                                                                                 ebqe[('da',0,0,0)],
                                                                                                 ebqe[('da',0,0,1)],
                                                                                                 ebqe[('a',1,1)],
                                                                                                 ebqe[('da',1,1,0)],
                                                                                                 ebqe[('da',1,1,1)],
                                                                                                 ebqe[('grad(phi)',0)],
                                                                                                 ebqe[('grad(phi)',1)],
                                                                                                 dphi[(0,0)].dof,
                                                                                                 dphi[(0,1)].dof,
                                                                                                 dphi[(1,0)].dof,
                                                                                                 dphi[(1,1)].dof,
                                                                                                 ebqe[('u',0)],
                                                                                                 ebqe[('u',1)],
                                                                                                 ebqe['psi_n'],
                                                                                                 ebqe[('dpsi_n',0)],
                                                                                                 ebqe[('dpsi_n',1)],
                                                                                                 v,
                                                                                                 grad_v,
                                                                                                 ebqe[('penalty')],
                                                                                                 ebqe[('penalty')],
                                                                                                 fluxJacobian_exterior[0][0],
                                                                                                 fluxJacobian_exterior[0][1],
                                                                                                 fluxJacobian_exterior[1][0],
                                                                                                 fluxJacobian_exterior[1][1])

            #
            #print "DarcyFC exterior fluxJacobian[0][0]=%s \n fluxJacobian[0][1]=%s " % (fluxJacobian_exterior[0][0],fluxJacobian_exterior[0][1])
            #mwf debug
            #import pdb
            #pdb.set_trace()

class DarcyFCPP_IIPG_exterior(NF_base):
    r"""
    weak dirichlet boundary conditions for Twophase_fc class

    .. todo::
       - put in nonlinear bc for setting :math:`\psi_n = \psi_n^b`
       - put in bc that switches flux and dirichlet types
    """
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        #hold Dirichlet values for \psi_n (non-wetting phase head)
        #also need extra psi_n entries that aren't part of default quadrature
        for term in ['psi_n_bc','psi_n',('dpsi_n',0),('dpsi_n',1),'sw']:
            self.ebqe[term] = numpy.zeros(self.ebqe[('u',1)].shape,'d')
        self.hasInterior=False
        self.penalty_constant = 2.0
        self.penalty_power = 1.0
        self.fluxBoundaryFlags = {0:0,1:0}#default is no-flow
        for ci in range(2):
            if vt.fluxBoundaryConditions[ci] == 'outFlow':
                self.fluxBoundaryFlags[ci] = 1
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                if self.isDOFBoundary[ci][ebNE,k] == 2:
                    self.ebqe['psi_n_bc'][ebNE,k] = g(x,self.vt.timeIntegration.t)
                else:
                    self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
                #mwf debug
                #import pdb
                #pdb.set_trace()
                #print "DarcyFCPP g(%s,%s)= %s " % (x,self.vt.timeIntegration.t,g(x,self.vt.timeIntegration.t))
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #how to handle flux boundary conditions, since usual advective/diffusive split doesn't apply here

        if self.vt.sd:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCPP_sd(self.vt.coefficients.sdInfo[(0,0)][0],self.vt.coefficients.sdInfo[(0,0)][1],
                                                                            self.vt.coefficients.sdInfo[(1,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                            self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            self.isDOFBoundary[0],
                                                                            self.isDOFBoundary[1],
                                                                            self.fluxBoundaryFlags[0],
                                                                            self.fluxBoundaryFlags[1],
                                                                            ebqe['n'],
                                                                            self.ebqe[('a',0,0)],
                                                                            self.ebqe[('a',1,1)],
                                                                            self.ebqe[('grad(phi)',0)],
                                                                            self.ebqe[('grad(phi)',1)],
                                                                            self.ebqe[('u',0)],
                                                                            self.ebqe[('u',1)],
                                                                            self.ebqe['psi_n_bc'],
                                                                            ebqe[('a',0,0)],
                                                                            ebqe[('a',1,1)],
                                                                            ebqe[('grad(phi)',0)],
                                                                            ebqe[('grad(phi)',1)],
                                                                            ebqe[('u',0)],
                                                                            ebqe[('u',1)],
                                                                            ebqe['psi_n'],
                                                                            ebqe[('penalty')],
                                                                            ebqe[('penalty')],
                                                                            ebqe[('diffusiveFlux',0,0)], #(ck,ci)
                                                                            ebqe[('diffusiveFlux',1,1)])
        else:
            cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCPP(self.mesh.exteriorElementBoundariesArray,
                                                                         self.mesh.elementBoundaryElementsArray,
                                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                         self.isDOFBoundary[0],
                                                                         self.isDOFBoundary[1],
                                                                         self.fluxBoundaryFlags[0],
                                                                         self.fluxBoundaryFlags[1],
                                                                         ebqe['n'],
                                                                         self.ebqe[('a',0,0)],
                                                                         self.ebqe[('a',1,1)],
                                                                         self.ebqe[('grad(phi)',0)],
                                                                         self.ebqe[('grad(phi)',1)],
                                                                         self.ebqe[('u',0)],
                                                                         self.ebqe[('u',1)],
                                                                         self.ebqe['psi_n_bc'],
                                                                         ebqe[('a',0,0)],
                                                                         ebqe[('a',1,1)],
                                                                         ebqe[('grad(phi)',0)],
                                                                         ebqe[('grad(phi)',1)],
                                                                         ebqe[('u',0)],
                                                                         ebqe[('u',1)],
                                                                         ebqe['psi_n'],
                                                                         ebqe[('penalty')],
                                                                         ebqe[('penalty')],
                                                                         ebqe[('diffusiveFlux',0,0)], #(ck,ci)
                                                                         ebqe[('diffusiveFlux',1,1)])
        #print "DarcyFCPP exterior diffusiveFlux[0]=%s \n diffusiveFlux[1]=%s " % (ebqe[('diffusiveFlux',0,0)],ebqe[('diffusiveFlux',1,1)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        if (self.vt.timeIntegration.diffusionIsImplicit[0] or self.vt.timeIntegration.diffusionIsImplicit[1]):
            #assume the u_w and u_m spaces are the same right now
            l2g_local = dphi[(0,0)].femSpace.dofMap.l2g; v   = ebqe[('v',0)]; grad_v = ebqe[('grad(v)',0)];
            if self.vt.sd:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian_sd(self.vt.coefficients.sdInfo[(0,0)][0],self.vt.coefficients.sdInfo[(0,0)][1],
                                                                                                      self.vt.coefficients.sdInfo[(1,1)][0],self.vt.coefficients.sdInfo[(1,1)][1],
                                                                                                      l2g_local,
                                                                                                      self.mesh.exteriorElementBoundariesArray,
                                                                                                      self.mesh.elementBoundaryElementsArray,
                                                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                      self.isDOFBoundary[0],
                                                                                                      self.isDOFBoundary[1],
                                                                                                      self.fluxBoundaryFlags[0],
                                                                                                      self.fluxBoundaryFlags[1],
                                                                                                      ebqe['n'],
                                                                                                      ebqe[('a',0,0)],
                                                                                                      ebqe[('da',0,0,0)],
                                                                                                      ebqe[('da',0,0,1)],
                                                                                                      ebqe[('a',1,1)],
                                                                                                      ebqe[('da',1,1,0)],
                                                                                                      ebqe[('da',1,1,1)],
                                                                                                      ebqe[('grad(phi)',0)],
                                                                                                      ebqe[('grad(phi)',1)],
                                                                                                      dphi[(0,0)].dof,
                                                                                                      dphi[(0,1)].dof,
                                                                                                      dphi[(1,0)].dof,
                                                                                                      dphi[(1,1)].dof,
                                                                                                      ebqe[('u',0)],
                                                                                                      ebqe[('u',1)],
                                                                                                      ebqe['psi_n'],
                                                                                                      ebqe[('dpsi_n',0)],
                                                                                                      ebqe[('dpsi_n',1)],
                                                                                                      v,
                                                                                                      grad_v,
                                                                                                      ebqe[('penalty')],
                                                                                                      ebqe[('penalty')],
                                                                                                      fluxJacobian_exterior[0][0],
                                                                                                      fluxJacobian_exterior[0][1],
                                                                                                      fluxJacobian_exterior[1][0],
                                                                                                      fluxJacobian_exterior[1][1])
            else:
                cnumericalFlux.calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian(l2g_local,
                                                                                                   self.mesh.exteriorElementBoundariesArray,
                                                                                                   self.mesh.elementBoundaryElementsArray,
                                                                                                   self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                   self.isDOFBoundary[0],
                                                                                                   self.isDOFBoundary[1],
                                                                                                   self.fluxBoundaryFlags[0],
                                                                                                   self.fluxBoundaryFlags[1],
                                                                                                   ebqe['n'],
                                                                                                   ebqe[('a',0,0)],
                                                                                                   ebqe[('da',0,0,0)],
                                                                                                   ebqe[('da',0,0,1)],
                                                                                                   ebqe[('a',1,1)],
                                                                                                   ebqe[('da',1,1,0)],
                                                                                                   ebqe[('da',1,1,1)],
                                                                                                   ebqe[('grad(phi)',0)],
                                                                                                   ebqe[('grad(phi)',1)],
                                                                                                   dphi[(0,0)].dof,
                                                                                                   dphi[(0,1)].dof,
                                                                                                   dphi[(1,0)].dof,
                                                                                                   dphi[(1,1)].dof,
                                                                                                   ebqe[('u',0)],
                                                                                                   ebqe[('u',1)],
                                                                                                   ebqe['psi_n'],
                                                                                                   ebqe[('dpsi_n',0)],
                                                                                                   ebqe[('dpsi_n',1)],
                                                                                                   v,
                                                                                                   grad_v,
                                                                                                   ebqe[('penalty')],
                                                                                                   ebqe[('penalty')],
                                                                                                   fluxJacobian_exterior[0][0],
                                                                                                   fluxJacobian_exterior[0][1],
                                                                                                   fluxJacobian_exterior[1][0],
                                                                                                   fluxJacobian_exterior[1][1])

class ShallowWater_1D(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-8,
                 tol_u=1.0e-8):
        self.h_eps=h_eps
        self.tol_u=tol_u
        self.g = vt.coefficients.g
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        cnumericalFlux.calculateInteriorNumericalFluxShallowWater_1D(self.h_eps,
                                                                     self.tol_u,
                                                                     self.g,
                                                                     self.mesh.interiorElementBoundariesArray,
                                                                     self.mesh.elementBoundaryElementsArray,
                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                     ebq['n'],
                                                                     ebq[('u',0)],
                                                                     ebq[('u',1)],
                                                                     ebq_global[('advectiveFlux',0)],
                                                                     ebq_global[('advectiveFlux',1)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)

        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        cnumericalFlux.calculateExteriorNumericalFluxShallowWater_1D(self.mesh.nExteriorElementBoundaries_global,
                                                                     self.h_eps,
                                                                     self.tol_u,
                                                                     self.g,
                                                                     ebqe['n'],
                                                                     ebqe[('u',0)],
                                                                     ebqe[('u',1)],
                                                                     self.ebqe[('u',0)],
                                                                     self.ebqe[('u',1)],
                                                                     ebqe[('advectiveFlux',0)],
                                                                     ebqe[('advectiveFlux',1)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        pass
class ShallowWaterHLL_1D(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-8,#1.0e-8
                 tol_u=1.0e-8):
        self.h_eps=h_eps
        self.tol_u=tol_u
        self.g = vt.coefficients.g
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        cnumericalFlux.calculateInteriorNumericalFluxShallowWaterHLL_1D(self.h_eps,
                                                                        self.tol_u,
                                                                        self.g,
                                                                        self.mesh.interiorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        ebq['n'],
                                                                        ebq[('u',0)],
                                                                        ebq[('u',1)],
                                                                        ebq_global[('advectiveFlux',0)],
                                                                        ebq_global[('advectiveFlux',1)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)

        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        cnumericalFlux.calculateExteriorNumericalFluxShallowWaterHLL_1D(self.mesh.nExteriorElementBoundaries_global,
                                                                        self.h_eps,
                                                                        self.tol_u,
                                                                        self.g,
                                                                        ebqe['n'],
                                                                        ebqe[('u',0)],
                                                                        ebqe[('u',1)],
                                                                        self.ebqe[('u',0)],
                                                                        self.ebqe[('u',1)],
                                                                        ebqe[('advectiveFlux',0)],
                                                                        ebqe[('advectiveFlux',1)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        pass
class ShallowWater_2D(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-8,
                 tol_u=1.0e-8):
        self.h_eps=h_eps
        self.tol_u=tol_u
        self.g = vt.coefficients.g
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        cnumericalFlux.calculateInteriorNumericalFluxShallowWater_2D(self.h_eps,
                                                                     self.tol_u,
                                                                     self.g,
                                                                     self.mesh.interiorElementBoundariesArray,
                                                                     self.mesh.elementBoundaryElementsArray,
                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                     ebq['n'],
                                                                     ebq[('u',0)],
                                                                     ebq[('u',1)],
                                                                     ebq[('u',2)],
                                                                     ebq_global[('advectiveFlux',0)],
                                                                     ebq_global[('advectiveFlux',1)],
                                                                     ebq_global[('advectiveFlux',2)])
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        cnumericalFlux.calculateExteriorNumericalFluxShallowWater_2D(self.mesh.nExteriorElementBoundaries_global,
                                                                     self.h_eps,
                                                                     self.tol_u,
                                                                     self.g,
                                                                     ebqe['n'],
                                                                     ebqe[('u',0)],
                                                                     ebqe[('u',1)],
                                                                     ebqe[('u',2)],
                                                                     self.ebqe[('u',0)],
                                                                     self.ebqe[('u',1)],
                                                                     self.ebqe[('u',2)],
                                                                     ebqe[('advectiveFlux',0)],
                                                                     ebqe[('advectiveFlux',1)],
                                                                     ebqe[('advectiveFlux',2)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        pass


class RusanovNumericalFlux(RusanovNumericalFlux_Diagonal):
    r"""Base class for Rusanov scheme for generic systems, relies on a
    user-specified estimate for the maximum (magnitude) eigenvalue for the
    system.

    Default is to take

    .. math::

       \bar{\lamda}= max |\vec f^{\prime} \cdot n\|

    Then we just apply apply numerical flux

    .. math::

       f_{num}(a,b) = 1/2(f(a)+f(b)-\bar{\lambda}(b-a)

    For now, we will try to piggy back on the cfl calculation, assuming the
    user has implemented this correctly for the system in question
    """

    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        RusanovNumericalFlux_Diagonal.__init__(self,vt,getPointwiseBoundaryConditions,
                                               getAdvectiveFluxBoundaryConditions,
                                               getDiffusiveFluxBoundaryConditions)

        self.ebqe['eigen_bound']=numpy.copy(self.ebqe[('u',0)])
        if not self.vt.q.has_key('eigen_bound'):
            self.vt.q['eigen_bound']=numpy.copy(self.vt.q[('u',0)])
        #mwf debug
        #self.safetyFactor=1.5
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        self.vt.q['eigen_bound'].fill(0.0)
        for ci in range(self.nc):
            self.vt.q['eigen_bound']=numpy.maximum(numpy.absolute(q[('cfl',ci)]),self.vt.q['eigen_bound'],self.vt.q['eigen_bound'])
        #cfl has h^-1 included
        for k in range(self.vt.q['eigen_bound'].shape[-1]):
            self.vt.q['eigen_bound'][:,k] *= self.vt.mesh.elementDiametersArray
        #mwf debug
        #print "Rusanov System \nh=%s \nu=%s \ncfl_0=%s \ncfl_1=%s \neigen_bound=%s " % (q[('u',0)],q[('u',1)],q[('cfl',0)],q[('cfl',1)],self.vt.q['eigen_bound'])
        #TODO need an additional loop through j to get nondiagonal flux derivatives for Jacobian
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(self.safetyFactor,
                                                                                             self.mesh.interiorElementBoundariesArray,
                                                                                             self.mesh.elementBoundaryElementsArray,
                                                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                             ebq['n'],
                                                                                             ebq[('u_advectiveNumericalFlux',ci)],
                                                                                             ebq[('f_advectiveNumericalFlux',ci)],
                                                                                             self.vt.q['eigen_bound'],
                                                                                             ebq_global[('advectiveFlux',ci)])

    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Rusanove_DiagonalUpwind computing bcs ebNE=%d k=%d x=%s g=%s " % (ebNE,k,x,g(x,self.vt.timeIntegration.t))
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        #Assume q['eigen_bound'] is up to date from interior numerical flux calc?
        self.vt.q['eigen_bound'].fill(0.0)
        for ci in range(self.nc):
            self.vt.q['eigen_bound']=numpy.maximum(numpy.absolute(q[('cfl',ci)]),self.vt.q['eigen_bound'],self.vt.q['eigen_bound'])
        #cfl has h^-1 included
        for k in range(self.vt.q['eigen_bound'].shape[-1]):
            self.vt.q['eigen_bound'][:,k] *= self.vt.mesh.elementDiametersArray
        #TODO need an additional loop through j to get nondiagonal flux derivatives for Jacobian
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(self.safetyFactor,
                                                                                             self.mesh.exteriorElementBoundariesArray,
                                                                                             self.mesh.elementBoundaryElementsArray,
                                                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                             self.isDOFBoundary[ci],
                                                                                             inflowFlag[ci],
                                                                                             ebqe['n'],
                                                                                             self.ebqe[('u',ci)],
                                                                                             self.ebqe[('f',ci)],
                                                                                             ebqe[('u_advectiveNumericalFlux',ci)],
                                                                                             ebqe[('f_advectiveNumericalFlux',ci)],
                                                                                             self.vt.q['eigen_bound'],
                                                                                             ebqe[('advectiveFlux',ci)])
            #mwf debug
            #if ci == self.nc-1:
            #    import pdb
            #    pdb.set_trace()
        #mwf add for inflow flux?
#         for ci in range(self.nc):
#             cnumericalFlux.calculateGlobalExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
#                                                                                self.mesh.elementBoundaryElementsArray,
#                                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                inflowFlag[ci],
#                                                                                ebqe[('inflowFlux',ci)],
#                                                                                ebqe['n'],
#                                                                                ebqe[('f',ci)],
#                                                                                ebqe[('df',ci,ci)],
#                                                                                ebqe[('advectiveFlux',ci)],
#                                                                                ebqe[('dadvectiveFlux_left',ci,ci)])



class RusanovLDG(Diffusion_LDG):
    """
    combine Rusanov flux for advection and LDG for diffusion
    right now advection jacobian is not calculated
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):

        Diffusion_LDG.__init__(self,vt,getPointwiseBoundaryConditions,
                              getAdvectiveFluxBoundaryConditions,
                              getDiffusiveFluxBoundaryConditions)
        self.safetyFactor=1.1
        self.ebqe['eigen_bound']=numpy.copy(self.ebqe[('u',0)])
        if not self.vt.q.has_key('eigen_bound'):
            self.vt.q['eigen_bound']=numpy.copy(self.vt.q[('u',0)])
        #add extra terms that can be lagged specifically for advective flux. Time integrator has to do this though
        for ci in range(self.nc):
            #ebq
            vt.ebq[('u_advectiveNumericalFlux',ci)]= vt.ebq[('u',ci)]
            if vt.ebq.has_key(('f',ci)):
                vt.ebq[('f_advectiveNumericalFlux',ci)]= vt.ebq[('f',ci)]
            for cj in range(self.nc):
                if vt.ebq.has_key(('df',ci,cj)):
                    vt.ebq[('df_advectiveNumericalFlux',ci,cj)]= vt.ebq[('df',ci,cj)]
                if vt.q.has_key(('df',ci,cj)):
                    vt.q[('df_advectiveNumericalFlux',ci,cj)] = vt.q[('df',ci,cj)]
            #ebqe
            vt.ebqe[('u_advectiveNumericalFlux',ci)]= vt.ebqe[('u',ci)]
            if vt.ebqe.has_key(('f',ci)):
                vt.ebqe[('f_advectiveNumericalFlux',ci)]= vt.ebqe[('f',ci)]
            for cj in range(self.nc):
                if vt.ebqe.has_key(('df',ci,cj)):
                    vt.ebqe[('df_advectiveNumericalFlux',ci,cj)]= vt.ebqe[('df',ci,cj)]

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        Diffusion_LDG.calculateInteriorNumericalFlux(self,q,ebq,ebq_global)
        #now Rusanov for advection
        self.vt.q['eigen_bound'].fill(0.0)
        for ci in range(self.nc):
            self.vt.q['eigen_bound']=numpy.maximum(numpy.absolute(q[('cfl',ci)]),self.vt.q['eigen_bound'],self.vt.q['eigen_bound'])
        #cfl has h^-1 included
        for k in range(self.vt.q['eigen_bound'].shape[-1]):
            self.vt.q['eigen_bound'][:,k] *= self.vt.mesh.elementDiametersArray
        #mwf debug
        #print "Rusanov System \nh=%s \nu=%s \ncfl_0=%s \ncfl_1=%s \neigen_bound=%s " % (q[('u',0)],q[('u',1)],q[('cfl',0)],q[('cfl',1)],self.vt.q['eigen_bound'])
        #TODO need an additional loop through j to get nondiagonal flux derivatives for Jacobian
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(self.safetyFactor,
                                                                                                 self.mesh.interiorElementBoundariesArray,
                                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                 ebq['n'],
                                                                                                 ebq[('u_advectiveNumericalFlux',ci)],
                                                                                                 ebq[('f_advectiveNumericalFlux',ci)],
                                                                                                 self.vt.q['eigen_bound'],
                                                                                                 ebq_global[('advectiveFlux',ci)])


    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        Diffusion_LDG.calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe)
        ####
        #now Rusanov contr.
        ####
        self.vt.q['eigen_bound'].fill(0.0)
        for ci in range(self.nc):
            self.vt.q['eigen_bound']=numpy.maximum(numpy.absolute(q[('cfl',ci)]),self.vt.q['eigen_bound'],self.vt.q['eigen_bound'])
        #cfl has h^-1 included
        for k in range(self.vt.q['eigen_bound'].shape[-1]):
            self.vt.q['eigen_bound'][:,k] *= self.vt.mesh.elementDiametersArray
        #TODO need an additional loop through j to get nondiagonal flux derivatives for Jacobian
        for ci in range(self.nc):
            if ebqe.has_key(('f',ci)):
                cnumericalFlux.calculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(self.safetyFactor,
                                                                                                 self.mesh.exteriorElementBoundariesArray,
                                                                                                 self.mesh.elementBoundaryElementsArray,
                                                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                                 self.isDOFBoundary[ci],
                                                                                                 inflowFlag[ci],
                                                                                                 ebqe['n'],
                                                                                                 self.ebqe[('u',ci)],
                                                                                                 self.ebqe[('f',ci)],
                                                                                                 ebqe[('u',ci)],
                                                                                                 ebqe[('f',ci)],
                                                                                                 self.vt.q['eigen_bound'],
                                                                                                 ebqe[('advectiveFlux',ci)])
class HamiltonJacobi_DiagonalChengShu(NF_base):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 speedEvaluationType=0):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions)
        self.speedEvaluationType = int(speedEvaluationType)#1 -- use max, otherwise allow disc.
        #self.outFlowOnly=True
        for ci in range(self.nc):
            #by default advection-diffusion
            self.advectiveNumericalFlux[ci] = False
            self.diffusiveNumericalFlux[ci] = False
            self.HamiltonJacobiNumericalFlux[ci] = True

    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            cnumericalFlux.calculateInteriorChengShuNumericalFlux(self.speedEvaluationType,
                                                                  self.mesh.interiorElementBoundariesArray,
                                                                  self.mesh.elementBoundaryElementsArray,
                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                  ebq['n'],
                                                                  ebq[('u',ci)],
                                                                  ebq[('H',ci)],
                                                                  ebq[('dH',ci,ci)],
                                                                  q[('H',ci)],
                                                                  q[('dH',ci,ci)],
                                                                  ebq[('HamiltonJacobiFlux',ci)],
                                                                  ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                  ebq[('dHamiltonJacobiFlux_right',ci,ci)])

        #mwf debug
        #import pdb
        #pdb.set_trace()
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                #print "bc--------------------",bci
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
        for ci in range(self.nc):
            cnumericalFlux.calculateExteriorLesaintRaviartNumericalFlux(self.speedEvaluationType,
                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                        self.mesh.elementBoundaryElementsArray,
                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                        self.isDOFBoundary[ci],
                                                                        inflowFlag[ci],
                                                                        ebqe['n'],
                                                                        self.ebqe[('u',ci)],
                                                                        self.ebqe[('H',ci)],
                                                                        self.ebqe[('dH',ci,ci)],
                                                                        ebqe[('u',ci)],
                                                                        ebqe[('H',ci)],
                                                                        ebqe[('dH',ci,ci)],
                                                                        ebqe[('HamiltonJacobiFlux',ci)],
                                                                        ebqe[('dHamiltonJacobiFlux_left',ci,ci)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateInteriorTwoSidedNumericalFluxJacobian(self.mesh.interiorElementBoundariesArray,
                                                                           self.mesh.elementBoundaryElementsArray,
                                                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                           ebq[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                           ebq[('dHamiltonJacobiFlux_right',ci,ci)],
                                                                           ebq[('v',ci)],
                                                                           fluxJacobian_hj[ci][ci])
                #mwf debug
                #import pdb
                #pdb.set_trace()

    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        for ci in range(self.nc):
            if self.vt.timeIntegration.hamiltonianIsImplicit[ci]:
                cnumericalFlux.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebqe[('dHamiltonJacobiFlux_left',ci,ci)],
                                                                            ebqe[('v',ci)],
                                                                            fluxJacobian_exterior[ci][ci])

class Stress_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                         getAdvectiveFluxBoundaryConditions,
                         getDiffusiveFluxBoundaryConditions,
                         getPeriodicBoundaryConditions,
                         parallelPeriodic=True)
        self.hasInterior=False
        self.penalty_constant = 100.0e5
        self.penalty_power = 1.0
        self.isStressBoundary={}
        for ci,sbcObject  in vt.stressFluxBoundaryConditionsObjectsDict.iteritems():
            self.isStressBoundary[ci] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,vt.ebqe['x'].shape[-2]),'i')
            for t,g in sbcObject.stressFluxBoundaryConditionsDict.iteritems():
                self.isStressBoundary[ci][t[0],t[1]] = 1
    def setDirichletValues(self,ebqe):
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                #cek todo: this needs to be generlzied for  all the numerical  fluxes when domain moves
                self.ebqe[('u',ci)][ebNE,k]=g(ebqe['x'][ebNE,k],self.vt.timeIntegration.t)
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
        if self.vt.movingDomain:
            self.vt.coefficients.updateToMovingDomain(self.vt.timeIntegration.t,self.ebqe)
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        pass
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        import ctransportCoefficients
        self.setDirichletValues(ebqe)
        #not using any coefficients on boundary
        #self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        cnumericalFlux.calculateGlobalExteriorNumericalStressFlux(self.mesh.exteriorElementBoundariesArray,
                                                                  self.mesh.elementBoundaryElementsArray,
                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                  self.isDOFBoundary[0],
                                                                  self.isDOFBoundary[1],
                                                                  self.isDOFBoundary[2],
                                                                  ebqe['n'],
                                                                  self.ebqe[('u',0)],
                                                                  self.ebqe[('u',1)],
                                                                  self.ebqe[('u',2)],
                                                                  ebqe['sigma'],
                                                                  ebqe[('u',0)],
                                                                  ebqe[('u',1)],
                                                                  ebqe[('u',2)],
                                                                  ebqe[('penalty')],
                                                                  ebqe[('stressFlux',0)],
                                                                  ebqe[('stressFlux',1)],
                                                                  ebqe[('stressFlux',2)])
    def updateInteriorNumericalFluxJacobian(self,l2g,q,ebq,ebq_global,dphi,fluxJacobian,fluxJacobian_eb,fluxJacobian_hj):
        pass
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        cnumericalFlux.updateExteriorNumericalStressFluxJacobian(self.mesh.exteriorElementBoundariesArray,
                                                                 self.mesh.elementBoundaryElementsArray,
                                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                 self.isDOFBoundary[0],
                                                                 self.isDOFBoundary[1],
                                                                 self.isDOFBoundary[2],
                                                                 self.isStressBoundary[0],
                                                                 self.isStressBoundary[1],
                                                                 self.isStressBoundary[2],
                                                                 ebqe['n'],
                                                                 ebqe[('dsigma',0,0)],
                                                                 ebqe[('dsigma',0,1)],
                                                                 ebqe[('dsigma',0,2)],
                                                                 ebqe[('dsigma',1,0)],
                                                                 ebqe[('dsigma',1,1)],
                                                                 ebqe[('dsigma',1,2)],
                                                                 ebqe[('dsigma',2,0)],
                                                                 ebqe[('dsigma',2,1)],
                                                                 ebqe[('dsigma',2,2)],
                                                                 ebqe[('v',0)],
                                                                 ebqe[('grad(v)',0)],
                                                                 ebqe['penalty'],
                                                                 fluxJacobian_exterior[0][0],
                                                                 fluxJacobian_exterior[0][1],
                                                                 fluxJacobian_exterior[0][2],
                                                                 fluxJacobian_exterior[1][0],
                                                                 fluxJacobian_exterior[1][1],
                                                                 fluxJacobian_exterior[1][2],
                                                                 fluxJacobian_exterior[2][0],
                                                                 fluxJacobian_exterior[2][1],
                                                                 fluxJacobian_exterior[2][2])

class Stress_SIPG_exterior(Stress_IIPG_exterior):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        Stress_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                      getAdvectiveFluxBoundaryConditions,
                                      getDiffusiveFluxBoundaryConditions,getPeriodicBoundaryConditions)
        self.includeBoundaryAdjoint=False
        self.boundaryAdjoint_sigma=0.0
        self.hasInterior=False

class Richards_IIPG_exterior(NF_base):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NF_base.__init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions)
        self.epsSeepage = 3.0
        self.hasInterior=False
        self.penalty_constant = 100.0
        self.penalty_power = 1.0

    def setDirichletValues(self,ebqe):
        self.isDOFBoundary[0][:]=0
        for ci in range(self.nc):
            self.ebqe[('u',ci)].flat[:] = ebqe[('u',ci)].flat[:]
            for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                    self.DOFBoundaryConditionsDictList[ci].values(),
                                    self.DOFBoundaryPointDictList[ci].values()):
                self.ebqe[('u',ci)][ebNE,k]=g(x,self.vt.timeIntegration.t)
                self.isDOFBoundary[0][ebNE,k] = 1#this will get turned off if on the seepage boundary and flow is inward
        for ci in range(self.nc):
            for bci in self.periodicBoundaryConditionsDictList[ci].values():
                self.ebqe[('u',ci)][bci[0]]=ebqe[('u',ci)][bci[1]]
                self.ebqe[('u',ci)][bci[1]]=ebqe[('u',ci)][bci[0]]
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebqe):
        ebqe[('advectiveFlux',0)].flat[:]=0.0#just put everthing in the diffusive flux
        self.ebqe[('u',0)].flat[:] = ebqe[('u',0)].flat[:]
        for (ebNE,k),g,x in zip(self.DOFBoundaryConditionsDictList[0].keys(),
                                self.DOFBoundaryConditionsDictList[0].values(),
                                self.DOFBoundaryPointDictList[0].values()):
            self.ebqe[('u',0)][ebNE,k]=g(x,self.vt.timeIntegration.t)
            self.isDOFBoundary[0][ebNE,k] = 1#this will get turned off if on the seepage boundary and flow is inward
        self.vt.coefficients.evaluate(self.vt.timeIntegration.t,self.ebqe)
        self.ebqe[('da',0,0,0)].flat[:]=0.0
        self.ebqe[('df',0,0)].flat[:]=0.0
        cnumericalFlux.calculateExteriorNumericalFluxRichards_sd(self.vt.coefficients.sdInfo[(0,0)][0],
                                                                 self.vt.coefficients.sdInfo[(0,0)][1],
                                                                 self.vt.coefficients.isSeepageFace,
                                                                 self.isDOFBoundary[0],
                                                                 ebqe['n'],
                                                                 self.ebqe[('u',0)],
                                                                 self.ebqe[('a',0,0)],
                                                                 ebqe[('grad(u)',0)],
                                                                 ebqe[('u',0)],
                                                                 self.ebqe[('f',0)],
                                                                 ebqe[('penalty')],
                                                                 ebqe[('diffusiveFlux',0,0)])
    def updateExteriorNumericalFluxJacobian(self,l2g,inflowFlag,q,ebqe,dphi,fluxJacobian_exterior,fluxJacobian_eb,fluxJacobian_hj):
        cnumericalFlux.calculateExteriorNumericalFluxJacobianRichards_sd(self.vt.coefficients.sdInfo[(0,0)][0],
                                                                         self.vt.coefficients.sdInfo[(0,0)][1],
                                                                         self.isDOFBoundary[0],
                                                                         ebqe['n'],
                                                                         self.ebqe[('u',0)],
                                                                         self.ebqe[('a',0,0)],
                                                                         self.ebqe[('da',0,0,0)],
                                                                         ebqe[('grad(u)',0)],
                                                                         ebqe[('grad(v)',0)],
                                                                         ebqe[('u',0)],
                                                                         self.ebqe[('df',0,0)],
                                                                         ebqe[('v',0)],
                                                                         ebqe[('penalty')],
                                                                         fluxJacobian_exterior[0][0])

class Richards_SIPG_exterior(Richards_IIPG_exterior):
    hasInterior=False
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        Richards_IIPG_exterior.__init__(self,
                                        vt,
                                        getPointwiseBoundaryConditions,
                                        getAdvectiveFluxBoundaryConditions,
                                        getDiffusiveFluxBoundaryConditions)
        self.hasInterior=False
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.penalty_constant = 100.0
        self.penalty_power = 1.0
