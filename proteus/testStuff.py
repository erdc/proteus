#! /usr/bin/env python
"""
just add a place to put things I'm testing out before they are ready for
including in proteus
"""
import Norms
import testStuffImpl
import numpy
import FemTools
#for timing
import sys,os,copy,timeit
import Profiling
from .Profiling import logEvent
TESTVPPTIMES = False

#######################################################################


class ExplicitLevelSetSolver:
    """
    Implement explicit, stabilized P^1 C^0 HJ level set algorithm from Barth and Sethian 98 (5.2)

    """
    def __init__(self,nDOF_global,RKorder=2,solnComponent=0,redistance=False,
                 enforceWeakDir=False,redEpsFact=2.0,redCFL=0.1):
        """
        RKorder              --- order of explicit RK time integration (1 or 2)
        solnComponent        --- solution component in VectorTransport to apply alg to.
        """
        self.nDOF_global= nDOF_global
        self.RKorder    = 2
        self.component  = solnComponent
        self.w          = numpy.zeros((self.nDOF_global),'d')
        self.phiP12     = numpy.zeros((self.nDOF_global),'d')
        self.phi0       = numpy.zeros((self.nDOF_global),'d')
        self.phiP1      = numpy.zeros((self.nDOF_global),'d')
        #for redistancing
        self.redPhi     = numpy.zeros((self.nDOF_global),'d')
        self.weakDirichletFlag  = numpy.zeros((self.nDOF_global),'i')
        self.redCFL     = redCFL
        self.redEpsFact = redEpsFact
        self.redistance = redistance
        self.enforceWeakDirichletBCs = enforceWeakDir
    #end init

    def advanceStage(self,phiIn,phiOut,wOut,q,ebq,l2g,elementDiameters,
                     elementQuadratureWeights,dt):
        """
        apply P^1, C^0 FE discretization with GLS style stablization and mass lumping
        to move level set function forward one RK stage
        """
        ci = self.component
        #mwf debug
        #print """advanceStage ci= %s dH shap= %s """ % (ci,q[('dH',ci,ci)].shape)
        if ('r',ci) in q.keys():
            testStuffImpl.advanceStageP1_C0_SUPG_lump(elementDiameters,
                                                      ebq['n'],
                                                      q['abs(det(J))'],
                                                      ebq['sqrt(det(g))'],
                                                      elementQuadratureWeights,
                                                      l2g,
                                                      phiIn,
                                                      q[('H',ci)],
                                                      q[('dH',ci,ci)],
                                                      q[('r',ci)],
                                                      dt,
                                                      phiOut,
                                                      wOut)
        else:
            #mwf debug
            #print """testStuff advStage SGS_lump_noSource:
#n.shape=%s q[('H',ci)].shape=%s q[('dH',ci,ci)].shape=%s phiOut.shape=%s """ % (ebq['n'].shape,
#                                                                          q[('H',ci)].shape,
#                                                                          q[('dH',ci,ci)].shape,
#                                                                          phiOut.shape)

            testStuffImpl.advanceStageP1_C0_SGS_lump_noSource(elementDiameters,
                                                              ebq['n'],
                                                              q['abs(det(J))'],
                                                              ebq['sqrt(det(g))'],
                                                              elementQuadratureWeights,
                                                              l2g,
                                                              phiIn,
                                                              q[('H',ci)],
                                                              q[('dH',ci,ci)],
                                                              dt,
                                                              phiOut,
                                                              wOut)


        return

    def computeSolution(self,vtran,dt):
        """
        compute \phi^n --> \phi^{n+1}, for \Delta t= t^{n+1}-t^n
        """
        #mwf debug
##         import math
##         for eN in range(vtran.mesh.nElements_global):
##             for ebN in range(vtran.mesh.nElementBoundaries_element):
##                 p0   = vtran.mesh.nodeArray[vtran.mesh.elementNodesArray[eN,ebN],:]
##                 pnid = [int(math.fmod(ebN+1,vtran.mesh.nElementBoundaries_element)),
##                         int(math.fmod(ebN+2,vtran.mesh.nElementBoundaries_element))]
##                 pneig= [vtran.mesh.nodeArray[vtran.mesh.elementNodesArray[eN,pnid[0]],:],
##                         vtran.mesh.nodeArray[vtran.mesh.elementNodesArray[eN,pnid[1]],:]]
##                 pm = 0.5*(pneig[0]+pneig[1])
##                 for k in range(vtran.nElementBoundaryQuadraturePoints_elementBoundary):
##                     ndot = 0.0
##                     for id in range(vtran.nSpace_global):
##                         ndot += vtran.ebq['n'][eN,ebN,k,id]*(p0[id]-pm[id])
##                     #
##                     if ndot >= 0.0:
##                         print """normCheck eN=%d ebN=%d p0=%s pneig= %s pm=%s \n n=%s ndot=%g """ % (eN,ebN,p0,pneig,pm,
##                                                                                                      vtran.ebq['n'][eN,ebN,k,:],
##                                                                                                      ndot)
##                 #k
##             #ebn
##         #en
        #mwf debug
        #print """\nexLS computeSoln dt=%g\n""" % dt

        ci = self.component
        self.advanceStage(vtran.u[ci].dof,self.phiP12,self.w,
                          vtran.q,vtran.ebq,
                          vtran.u[ci].femSpace.dofMap.l2g,
                          vtran.mesh.elementDiametersArray,
                          vtran.elementQuadratureWeights[('m',ci)],
                          dt)
        self.phiP1.flat[:] = vtran.u[ci].dof.flat[:]
        self.phiP1.flat[:]-= dt*self.phiP12.flat[:]/self.w.flat[:]
        if self.RKorder == 2:
            self.phi0.flat[:] = vtran.u[ci].dof.flat[:]
        vtran.u[ci].dof.flat[:] = self.phiP1.flat[:]
        #put in update stage?
        #vtran.updateStage()
        if self.RKorder == 2:
            #need to reevaluate coefficients here, update time level too?
            vtran.calculateCoefficients()
            self.advanceStage(self.phiP1,self.phiP12,self.w,
                              vtran.q,vtran.ebq,
                              vtran.u[ci].femSpace.dofMap.l2g,
                              vtran.mesh.elementDiametersArray,
                              vtran.elementQuadratureWeights[('m',ci)],
                              dt)
            self.phiP1.flat[:] = 0.5*(self.phiP1.flat[:]+self.phi0.flat[:])-0.5*dt*self.phiP12.flat[:]/self.w.flat[:]
            vtran.u[ci].dof.flat[:] = self.phiP1.flat[:]
            #vtran.updateStage() #use this still
        #end RK2
        #what about redistancing now
        if self.redistance:
            #evalute coefficients at new value of phi
            #assume red eqn had characteristic speed about 1
            ihmin = numpy.argmin(vtran.mesh.elementDiametersArray)
            hmin  = vtran.mesh.elementDiametersArray[ihmin]
            dtau  = self.redCFL*hmin
            nredSteps = int(4.0/self.redCFL) #approximate time to propagate 4 cells
            ihmax = numpy.argmax(vtran.mesh.elementDiametersArray)
            hmax  = vtran.mesh.elementDiametersArray[ihmax]
            eps   = self.redEpsFact*hmax
            #mwf hack
            #dtau  = min(1.0e-3,dtau)
            nredSteps = min(10,nredSteps)
            #mwf debug
            print """test redistance hmin=%g dtau=%g nredSteps=%d\n """ %(hmin,dtau,nredSteps)
            print """test redistance hmax=%g eps=%g \n """ %(hmax,eps)
            self.redPhi.flat[:] = self.phiP1.flat[:]
            if self.enforceWeakDirichletBCs:
                for i in range(nredSteps):
                    testStuffImpl.advanceStageRedistanceWeakDirP1_C0_SUPG_lump(
                        vtran.mesh.elementDiametersArray,
                        vtran.ebq['n'],
                        vtran.q['abs(det(J))'],
                        vtran.ebq['sqrt(det(g))'],
                        vtran.elementQuadratureWeights[('m',ci)],
                        vtran.u[ci].femSpace.dofMap.l2g,
                        vtran.mesh.exteriorElementBoundariesArray,
                        vtran.mesh.elementBoundaryElementsArray,
                        self.redPhi,
                        self.phiP1,
                        dtau,
                        eps,
                        self.weakDirichletFlag,
                        self.phiP12,
                        self.w)
                    self.redPhi.flat[:] -= dtau*self.phiP12.flat[:]/self.w.flat[:]
                #end for
            else:
                for i in range(nredSteps):
                    testStuffImpl.advanceStageRedistanceP1_C0_SUPG_lump(
                        vtran.mesh.elementDiametersArray,
                        vtran.ebq['n'],
                        vtran.q['abs(det(J))'],
                        vtran.ebq['sqrt(det(g))'],
                        vtran.elementQuadratureWeights[('m',ci)],
                        vtran.u[ci].femSpace.dofMap.l2g,
                        self.redPhi,
                        self.phiP1,
                        dtau,
                        eps,
                        self.phiP12,
                        self.w)
                    self.redPhi.flat[:] -= dtau*self.phiP12.flat[:]/self.w.flat[:]
                 #end for
            vtran.u[ci].dof.flat[:] = self.redPhi.flat[:]
        #end redistance
        vtran.calculateCoefficients()

        #end if
        #mwf debug
        #print """\nAfter advanceStage phiIn = %s
        #phiOut= %s
        #self.phiP12= %s\n """ % (phiIn,phiOut,self.phiP12)

        return
    #end computeSolution

import TimeIntegration

class ExplicitLevelSetIntegrator(TimeIntegration.ForwardIntegrator):
    """
    put in manual integration to play with level-set redistancing
    """
    def __init__(self,mlvtran,mlnl,dtMeth,nOptions,stepExact=True):
        TimeIntegration.ForwardIntegrator.__init__(self,mlvtran,mlnl,dtMeth,nOptions,
                                                         stepExact)
        ci = 0
        self.lsSolver = ExplicitLevelSetSolver(mlvtran.modelList[-1].u[ci].nDOF_global,RKorder=2,
                                               solnComponent=ci,redistance=True,
                                               enforceWeakDir=False,
                                               redEpsFact=20.0,redCFL=0.02)
        self.dtLast=None
        self.dtRatioMax=2.0
    #def

    def calculateSolution(self,tIn,tOut):
        """
        manually step forward with fixed steps for now
        """
        maxCFL = 1.0e-6
        ci = self.lsSolver.component
        if self.mlvTran.modelList[-1].q.has_key(('cfl',ci)):
            maxCFL = max(max(max(self.mlvTran.modelList[-1].q[('cfl',ci)])),maxCFL)
        dt = self.nOptions.runCFL/maxCFL
        if self.dtLast is None:
            self.dtLast = dt
        if dt/self.dtLast  > self.dtRatioMax:
            dt = self.dtLast*self.dtRatioMax
        self.dtLast = dt
        #mwf debug
        #print """\nexLSint tIn=%g tOut=%g dt=%g """ % (tIn,tOut,dt)
        t = tIn; failedFlag = False
        #mwf hack
        dt = 1.0e-3
        while t < tOut and failedFlag == False:
            self.mlvTran.chooseDT(dt)
            if self.stepExact and t+self.mlvTran.DT > tOut:
                self.mlvTran.resetDT()
                self.mlvTran.chooseDT(tOut-t)
            #mwf debug
            #print """\nexLSint t=%g tOut=%g DTSET=%s DT=%g """ % (t,tOut,self.DTSET,self.mlvTran.DT)
            self.writeProgress(t,self.mlvTran.DT,tOut)
            self.lsSolver.computeSolution(self.mlvTran.modelList[-1],dt)
            t += self.mlvTran.DT
            if abs(tOut-t) < 1.0e-15:
                t = tOut
            self.mlvTran.updateTimeHistory(t)

            self.tLast = t
            maxCFL=1.e-6
            if self.mlvTran.modelList[-1].q.has_key(('cfl',ci)):
                maxCFL = max(max(max(self.mlvTran.modelList[-1].q[('cfl',ci)])),maxCFL)
            dt = self.nOptions.runCFL/maxCFL
            if dt/self.dtLast  > self.dtRatioMax:
                dt = self.dtLast*self.dtRatioMax
            #if
            self.dtLast = dt
        #while

        return failedFlag,t
    #def
#class


########################################################################



########################################################################
#test projecting solutions to fine grid for computing error
########################################################################
from FemTools import *

def projectToFinestLevel(mlTransport,level,tsim=0.0):
    """
    use multilevel transport prolongation to get fine grid information
    starting at level.

    returns quadrature dictionary of projected values on fine grid
    """
    import numpy
    nLevels = len(mlTransport.uList)
    assert 0 <= level and level <= nLevels, "projectToFinestLevel range= [0,%d]" % nLevels-1

    coefficients = mlTransport.modelList[-1].coefficients
    uqprojFine = [numpy.zeros(mlTransport.modelList[-1].q[('u',ci)].shape,'d')
                  for ci in range(coefficients.nc)]
    uproj  = [[FiniteElementFunction(m.u[ci].femSpace) for m in mlTransport.modelList]
              for ci in range(coefficients.nc)]
    mFine = mlTransport.modelList[-1]
    for ci in range(coefficients.nc):
        for l in range(nLevels):
            uproj[ci][l].dof[:] = 0.0
        m = mlTransport.modelList[level]
        uproj[ci][level].dof[:] = m.u[ci].dof
        if level < nLevels-1:
            for lf in range(level,nLevels-1):
                mlTransport.meshTransfers.prolong_bcListDict[ci][lf+1].matvec(uproj[ci][lf].dof,
                                                                              uproj[ci][lf+1].dof)
                #load Dirichlet conditions in
                for dofN,g in mlTransport.modelList[lf+1].dirichletConditions[ci].DOFBoundaryConditionsDict.iteritems():
                    uproj[ci][lf+1].dof[dofN] = g(mlTransport.modelList[lf+1].dirichletConditions[ci].DOFBoundaryPointDict[dofN],tsim)
                #dirichlet conditions
            #lf up to fine
            uproj[ci][-1].getValues(mFine.q['v',ci],uqprojFine[ci])
        else: #already on fine
            uqprojFine[ci].flat[:] = mFine.q[('u',ci)].flat[:]
        #else
    #ci

    return uqprojFine

########################################################################
#try computing solution values directly on fine grid using coarse grid
#for nonconforming approximations
########################################################################
def generateParentInfo(mlMesh):
    """
    get array P[l,e] = e_c, where element e_c is the parent of element e
        P[0,:] = -1
    """
    #mwf now use new interface
    P = mlMesh.elementParentsArrayList
    import numpy
    P = {}
    nLevels = len(mlMesh.meshList)
    P[0] = numpy.ones((mlMesh.meshList[0].nElements_global,),'i')
    P[0][:]=-1
    for l in range(1,nLevels):
        P[l] = mlMesh.elementParentsArrayList[l]
    return P

def getIntegrationPointsOnCoarseGrid(xf,lf,P,lc):
    """
    given array of points N^f_e x n_q x 3 on level lf
    generate dictionary on coarse
    grid that's N^c_e x n_q^c(e_c) x 3 and holds the integration points
    assigned to the correct coarse grid element. If using uniform refinement
    should have the same number of integration points per coarse grid element
    but the number depends on the level of refinement eg  n_q^4(lf-lc).
    In general the number won't be the same for nonuniform refinement
    """
    nEf = len(P[lf])
    assert nEf == xf.shape[0], "fine grid mismatch lf=%d nEf=%d xf.shape=%s " % (lf,nEf,xf.shape)
    nq  = xf.shape[1]
    nEc = len(P[lc])
    ldiff = lf-lc
    #mwf debug
    #print """testStuff generateCoarseGridX lf=%d nEf=%d xf.shape=%s lc=%d nEc=%d ldiff=%d """ % (lf,nEf,xf.shape,lc,nEc,
    #                                                                                             ldiff)

    xc = {}
    for ec in range(nEc):
        xc[ec] = {}

    for ef in range(nEf):
        ec = ef
        for l in range(ldiff):
            ep = P[lf-l][ec]
            #mwf debug
            #print "genCoarseGridX ef=%d l=%d ec=%d ep=%d " % (ef,l,ec,ep)
            ec = ep

        for iq in range(nq):
            xc[ec][(ef,iq)] = xf[ef,iq,:]
        #iq
    #ef
    #mwf debug
    #print """getIntPointsOnCoarseGrid lf=%s lc=%s \n xf=%s""" % (lf,lc,xf)
    #for ec in range(nEc):
    #    print """xc[%s] len=%d x=%s """ % (ec,len(xc[ec]),xc[ec])


    return xc

def getCoarseGridBasisValuesOnFinePointsUniform(uc,xc,invJc,nEf,nq):
    """
    assuming have collection of points in physical space on mesh xc with same number  of points
    per element, map these back to reference element, then calculate shape function values there
    and return in a fine grid array
    """
    nEc = len(xc)
    nqc = len(xc[0]) #needs to be the same for all ec
    xArray  = numpy.zeros((nEc,nqc,3),'d')
    xiArray = numpy.zeros((nEc,nqc,3),'d')
    vcArray = numpy.zeros((nEc,nqc,uc.femSpace.max_nDOF_element),'d')
    invJ    = numpy.zeros((nEc,nqc,invJc.shape[2],invJc.shape[3]),'d')#affine

    for ec in range(nEc):
        iqc = 0
        for k,x in xc[ec].iteritems():
            #mwf debug
            #print "ec=%d iqc=%d x=%s " % (ec,iqc,x)
            xArray[ec,iqc,:] = x
            invJ[ec,iqc,:,:] = invJc[ec,0,:,:]
            iqc += 1
        #end x
    #end ec

    uc.femSpace.elementMaps.getInverseValues(invJ,
                                             xArray,
                                             xiArray)

    uc.femSpace.getBasisValuesAtArray(xiArray,vcArray)
    #mwf debug
    #print "testStuff getCoarseGridBasis xArray= %s \n xiArray=%s \n vcArray=%s " % (xArray,xiArray,vcArray)

    return vcArray,xArray,xiArray
def getCoarseGridValuesOnFineUniform(uc,vcArray,xc,nEf,nq):
    """
    given basis representation for coarse grid at fine grid points compute fem function values
    and return in fine grid array
    """
    nEc = len(xc)
    nqc = len(xc[0]) #needs to be the same for all ec
    ucvals = numpy.zeros((nEc,nqc),'d')
    uc.getValues(vcArray,ucvals)

    #now assign back to fine grid format
    ufvals = numpy.zeros((nEf,nq),'d')
    for ec in range(nEc):
        iqc = 0
        for k,x in xc[ec].iteritems():
            ef = k[0]; iqf = k[1]
            ufvals[ef,iqf] = ucvals[ec,iqc]
            iqc += 1
        #k
    #ec
    #mwf debug
    #print "getCoarseGridValues uc.dof=%s \n ucvals= %s \n ufvals= %s " % (uc.dof,ucvals,ufvals)

    return ufvals

def projectToFinestLevelNC(mlTransport,level,ci=0,tsim=0.0):
    """
    use brute force evaluation to get coarse grid quantities on fine grid
    starting at level.

    returns quadrature dictionary of projected values on fine grid
    """
    import numpy
    verbose = 0
    nLevels = len(mlTransport.modelList)

    assert 0 <= level and level < nLevels, "projectToFinestLevelNC range= [0,%d]" % nLevels-1

    P = generateParentInfo(mlTransport.mlMeshSave)
    mFine  = mlTransport.modelList[-1]
    mCoarse= mlTransport.modelList[level]
    nEf = mFine.q['x'].shape[0]; nqf = mFine.q['x'].shape[1]
    uc = mlTransport.modelList[level].u[ci]
    lf = nLevels-1 #index into list
    #if level == lf:
    if False:
        uqciprojFine = mFine.q[('u',ci)]
    else:
        xc= getIntegrationPointsOnCoarseGrid(mFine.q['x'],lf,P,level)
        vc,xArray,xiArray= getCoarseGridBasisValuesOnFinePointsUniform(uc,xc,
                                                                       mCoarse.q['inverse(J)'],nEf,nqf)
        uqciprojFine = getCoarseGridValuesOnFineUniform(uc,vc,xc,nEf,nqf)

    #mwf debug
    if verbose > 1:
        import Viewers
        if 'viewerType' in dir(Viewers) and Viewers.viewerType == 'gnuplot' and mFine.nSpace_global == 2:
            for eN in range(uqciprojFine.shape[0]):
                for k in range(uqciprojFine.shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (mFine.q['x'][eN,k,0],
                                                                       mFine.q['x'][eN,k,1],
                                                                       uqciprojFine[eN,k]))
                #
            #eN
            Viewers.datFile.write("\n \n#end uqproj ci=%d level=%d \n" % (ci,level))
            ggrid = 32
            title="uqproj ci=%d level=%d " % (ci,level)
            cmd = "set dgrid3d %d,%d,16; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (ggrid,ggrid,Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
            raw_input('press return to continue')

            for eN in range(mCoarse.q[('u',ci)].shape[0]):
                for k in range(mCoarse.q[('u',ci)].shape[1]):
                    Viewers.datFile.write("%12.5e %12.5e %12.5e \n" % (mCoarse.q['x'][eN,k,0],
                                                                       mCoarse.q['x'][eN,k,1],
                                                                       mCoarse.q[('u',ci)][eN,k]))
                #
            #eN
            Viewers.datFile.write("\n \n#end uq ci=%d level=%d \n" % (ci,level))
            ggrid = 32
            title="uc ci=%d level=%d " % (ci,level)
            cmd = "set dgrid3d %d,%d,16; set contour base; set term x11 %i; splot \'%s\' index %i with lines title \"%s\" \n" % (ggrid,ggrid,Viewers.windowNumber,
                                                                                                                                 Viewers.datFilename,
                                                                                                                                 Viewers.plotNumber,
                                                                                                                                 title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
        #viewer
    #verbose
    return uqciprojFine



def projectVelocityToFinestLevelNC(mlTransport,level,ci=0,tsim=0.0,verbose=0):
    """
    use brute force evaluation to get coarse grid quantities on fine grid
    starting at level.

    returns quadrature dictionary of projected values on fine grid
    """
    import numpy
    nLevels = len(mlTransport.modelList)

    assert 0 <= level and level < nLevels, "projectVelocityToFinestLevelNC range= [0,%d]" % nLevels-1

    mFine  = mlTransport.modelList[-1]
    mCoarse= mlTransport.modelList[level]
    if mCoarse.velocityPostProcessor is None:
        return None

    P = generateParentInfo(mlTransport.mlMeshSave)
    nEf = mFine.q['x'].shape[0]; nqf = mFine.q['x'].shape[1]
    lf = nLevels-1 #index into list
    ldiff = lf-level
    if level == lf:
        velciprojFine = mFine.q[('velocity',ci)]
        xArray        = mFine.q['x']
    else:
        #mwf debug
        import pdb
        pdb.set_trace()
        xc= getIntegrationPointsOnCoarseGrid(mFine.q['x'],lf,P,level)
        nEc = len(xc)
        #nqc = len(xc[0]) #needs to be the same for all ec
        nqc = max([len(xc[i]) for i in range(nEc)])
        xArray  = numpy.zeros((nEc,nqc,3),'d')
        for ec in range(nEc):
            iqc = 0
            for k,x in xc[ec].iteritems():
                #mwf debug
                print "ec=%d iqc=%d x=%s " % (ec,iqc,x)
                xArray[ec,iqc,:] = x
                iqc += 1
            #end x
        #end ec
        velciprojFine = numpy.zeros(mFine.q[('velocity',ci)].shape,'d')
        if mCoarse.velocityPostProcessor.postProcessingTypes[ci] == 'point-eval':
            #assume constant solution/potential gradient over coarse grid
            print "WARNING projectVelocityToFinestLevelNC type= point-eval assuming constant potential on coarse grid"

            for ef in range(nEf):
                ec = ef
                for l in range(ldiff):
                    ep = P[lf-l][ec]
                    ec = ep
                for iq in range(nqf):
                    if mFine.q.has_key(('a',ci,ci)):
                        velciprojFine[ef,iq,:] = -numpy.dot(mFine.q[('a',ci,ci)][ef,iq,:,:],
                                                              mCoarse.q[('grad(phi)',ci)][ec,0,:])
                    else:
                        velciprojFine[ef,iq,:] = 0.0
                    if mFine.q.has_key(('f',ci)):
                        velciprojFine[ef,iq,:]  += mFine.q[('f',ci)][ef,iq,:]
                #iq
            #ef
        else:
            velci0 = mCoarse.velocityPostProcessor.evaluateElementVelocityField(xArray,ci)
            for ec in range(nEc):
                iqc = 0
                for k,x in xc[ec].iteritems():
                    ef = k[0]; iqf = k[1]
                    velciprojFine[ef,iqf,:] = velci0[ec,iqc,:]
                    iqc += 1
                #x
            #ec

        #postprocessing type
    #else on level
    if verbose > 2:
        print """testStuff velocityProjNC \n xArray=%s velciprojFine= %s \n""" % (xArray,
                                                                                 velciprojFine)
    if verbose > 1:
        import Viewers
        if 'viewerType' in dir(Viewers) and Viewers.viewerType == 'gnuplot' and mFine.nSpace_global == 2:
            max_u=max(numpy.absolute(numpy.take(velciprojFine,[0],2).flat))
            max_v=max(numpy.absolute(numpy.take(velciprojFine,[1],2).flat))
            L = min((max(mFine.mesh.nodeArray[:,0]),max(mFine.mesh.nodeArray[:,1])))
            scale =10.0*max((max_u,max_v,1.0e-16))/L
            for eN in range(mFine.q['x'].shape[0]):
                for iq in range(mFine.q['x'].shape[1]):
                    x = mFine.q['x'][eN,iq,:]
                    v = velciprojFine[eN,iq,:]
                    Viewers.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],
                                                                              v[0]/scale,
                                                                              v[1]/scale))
            Viewers.datFile.write("\n \n#end velciproj ci=%d level=%d" % (ci,level))
            title = "velciproj ci=%d level=%d " % (ci,level)
            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                  Viewers.datFilename,
                                                                                                  Viewers.plotNumber,
                                                                                                  title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
            raw_input('press return to continue')

            #now try just coarse grid velocity
            max_u=max(numpy.absolute(numpy.take(mCoarse.q[('velocity',ci)],[0],2).flat))
            max_v=max(numpy.absolute(numpy.take(mCoarse.q[('velocity',ci)],[1],2).flat))
            L = min((max(mFine.mesh.nodeArray[:,0]),max(mFine.mesh.nodeArray[:,1])))
            scale =10.0*max((max_u,max_v,1.0e-16))/L
            for eN in range(mCoarse.q['x'].shape[0]):
                for iq in range(mCoarse.q['x'].shape[1]):
                    x = mCoarse.q['x'][eN,iq,:]
                    v = mCoarse.q[('velocity',ci)][eN,iq,:]
                    Viewers.datFile.write("%12.5e %12.5e %12.5e %12.5e \n" % (x[0],x[1],
                                                                              v[0]/scale,
                                                                              v[1]/scale))
            Viewers.datFile.write("\n \n#end coarse velocity ci=%d level=%d" % (ci,level))
            title = "coarse velocity ci=%d level=%d " % (ci,level)
            cmd = "set term x11 %i; plot \'%s\' index %i with vectors title \"%s\" \n" % (Viewers.windowNumber,
                                                                                          Viewers.datFilename,
                                                                                          Viewers.plotNumber,
                                                                                          title)
            Viewers.cmdFile.write(cmd)
            Viewers.viewerPipe.write(cmd)
            Viewers.newPlot()
            Viewers.newWindow()
            raw_input('press return to continue')

        #end gnuplot
    #end verbose
    return velciprojFine

########################################################################
#try tweaking default transfer operators for NC space
class MultilevelProjectionOperatorsNC:
    """
    A class that takes a hierarchical (multiLevel) mesh and generates
    the interpolation and restriction operators.

    restrictList/prolongList -- projection and restriction at only the free nodes
    restrict_bcList/prolong_bcList -- includes dirichlet boundaries as well


    By default this is set up for conforming spaces. Since the spaces
    are conforming the coarse basis functions are in the fine space so we
    need only find the coefficients of the fine space basis functions that
    yield the coarse space basis functions. This is the matrix of the
    trivial injection from coarse to fine and it is used as the projection
    operator. Restriction is taken as the matrix of the adjoint of the
    injection, which is simply the transpose of the projection matrix. These
    operators fall out if you try to solve for the error on the coarse grid:
    Starting with u_f we have a(u_f+e_f,w_f) = <f,w_f>, and we want to
    solve instead  a(e_c,w_c) = <f - a(u_f,w_f),w_c> on the coarse grid
    Using the injection i we can express this in the fine space as
    a(i e_c, i w_c) = <f - a(u_f,w_f),i w_c> writing
    this in matrix form yields p^t A_f p E = p^t R_f

    --- P1 nonconforming space ----
    Try to set up now for nonconforming P1 approximation following Chen_96b.
    Setup prolongation by evaluating coarse grid basis functions at
    fine grid interpolation condition points (face barycenters).

    Then, just need to know if fine grid interpolationCondition point falls on
    interface of coarse grid elements or not. If so, use average value of coarse
    grid quantity on fine grid. Otherwise just evaluate it

    Use this simple interpolation from coarse to fine as the projection operator.
    Restriction is taken as the matrix of the adjoint of the
    injection, which is simply the transpose of the projection matrix.

    I don't think these fall out as nicely since they're nonconforming.

    """
    ## \todo put the MultilevelProjectionOperators constructor  partially  into C
    def __init__(self,
                 multiLevelMesh,
                 femSpaceDictList,
                 offsetListList,
                 strideListList,
                 dofBoundaryConditionsDictList):
        self.nc = len(offsetListList[0])
        self.restrictList = [[]]
        self.rzvalList = [[]]
        self.restrictSumList=[[]]
        self.prolongList = [[]]
        self.pzvalList = [[]]
        self.restrict_bcListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.rbczvalListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.restrict_bcSumListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.prolong_bcListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.pbczvalListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.scaled_restrict_bcListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.scaled_rbczvalListDict = dict([(cj,[[]]) for cj in range(self.nc)])
        self.femSpaceDictList=femSpaceDictList
        self.dof_bc_DictList = dofBoundaryConditionsDictList

        #mwf get parent info
        parents = generateParentInfo(multiLevelMesh)
        usingAtLeastOneNCproj = False
        for l in range(len(multiLevelMesh.meshList)-1):
            r = {}
            p = {}
            coarse_nFreeDOF_global = 0
            fine_nFreeDOF_global=0
            for cj in range(self.nc):
                coarse_nFreeDOF_global += self.dof_bc_DictList[l][cj].nFreeDOF_global
                fine_nFreeDOF_global   += self.dof_bc_DictList[l+1][cj].nFreeDOF_global
            rSum = Vec(coarse_nFreeDOF_global)
            rColumnIndeces=[set() for row in range(coarse_nFreeDOF_global)]
            for cj in range(self.nc):
                coarseSpace = self.femSpaceDictList[l][cj]
                coarseMesh  = multiLevelMesh.meshList[l]
                coarseDOFBoundaryConditions = self.dof_bc_DictList[l][cj]
                children = multiLevelMesh.elementChildren[l]
                nChildrenMax_element =0
                for coarse_eN in range(coarseMesh.nElements_global):
                    nChildrenMax_element = max(nChildrenMax_element,len(children[coarse_eN]))
                fineSpace = self.femSpaceDictList[l+1][cj]
                fineMesh  = multiLevelMesh.meshList[l+1]
                fineDOFBoundaryConditions = self.dof_bc_DictList[l+1][cj]
                fineSpaceInterpolationFunctionals = fineSpace.referenceFiniteElement.interpolationConditions.functionalsQuadrature
                #interpolationPointsOnFineElement_reference = fineSpace.referenceFiniteElement.interpolationConditions.quadraturePointArray
                nInterpolationPoints = fineSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints
                range_nInterpolationPoints = range(nInterpolationPoints)
                referenceElement = fineSpace.referenceFiniteElement.referenceElement
                rbcSum = Vec(coarseSpace.dim)
                rbcColumnIndeces=[set() for row in range(coarseSpace.dim)]
                rbc = {}
                pbc = {}
                scaled_rbc = {}

                #mwf should tell us if interpolation condition is at coarse element interface or not
                # aka (nonconforming on coarse grid?)
                #check typecode if go to c
                isNonConformingOnCoarseGrid = numpy.zeros((fineMesh.nElements_global,nInterpolationPoints),'i')
                if referenceElement.dim  > 1 and isinstance(fineSpace,FemTools.NC_AffineLinearOnSimplexWithNodalBasis):
                    usingAtLeastOneNCproj = True
                    #need to care about interpolation conditions being conforming on coarse grid or not
                    for ebNI in range(fineMesh.nInteriorElementBoundaries_global):
                        ebN = fineMesh.interiorElementBoundariesArray[ebNI]
                        eN_left  = fineMesh.elementBoundaryElementsArray[ebN,0]
                        eN_right = fineMesh.elementBoundaryElementsArray[ebN,1]
                        ebn_left = fineMesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                        ebn_right= fineMesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                        PeN_left = parents[l+1][eN_left]
                        PeN_right= parents[l+1][eN_right]
                        if PeN_left != PeN_right:
                            #assumes a unique correspondence interpCondition <--> local element boundary
                            isNonConformingOnCoarseGrid[eN_left,ebn_left]  = 1
                            isNonConformingOnCoarseGrid[eN_right,ebn_right]= 1
                            #mwf debug
                            #print """MultilevelProjNC nc IP ebN= %d eN_left=%d PeN_left=%d eN_right=%d PeN_right=%d """ % (ebN,
                            #                                                                                               eN_left,
                            #                                                                                               PeN_left,
                            #                                                                                               eN_right,
                            #                                                                                               PeN_right)
                        #different parents
                    #interior element boundaries
                    #what about exterior ones?
                    for ebNE in range(fineMesh.nExteriorElementBoundaries_global):
                        ebN = fineMesh.exteriorElementBoundariesArray[ebNE]
                        eN_left = fineMesh.elementBoundaryElementsArray[ebN,0]
                        ebn_left = fineMesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                        isNonConformingOnCoarseGrid[eN_left,ebn_left]  = -1
                    #end exterior
                #need to care about nonconforming points
                if (not isinstance(fineSpace,FemTools.NC_AffineLinearOnSimplexWithNodalBasis) and
                    referenceElement.dim > 1 and usingAtLeastOneNCproj):
                    print """WARNING testStuff.MultilevelProjectionOperatorsNC mixing NC and C0 projection operators!!!"""
                #map reference interpolation points of fine elements to physical space
                interpolationPointsOnFineElement_physical =  fineSpace.updateInterpolationPoints()
                #copy physical space reference points on fine elements to an array for their parents
                interpolationPointsOnCoarseElement_physical = numpy.zeros((coarseMesh.nElements_global,
                                                                             nChildrenMax_element*nInterpolationPoints,
                                                                             3),
                                                                            'd')
                for coarse_eN in range(coarseMesh.nElements_global):
                    for child_N,fine_e in enumerate(children[coarse_eN]):
                        for pN in range_nInterpolationPoints:
                            interpolationPointsOnCoarseElement_physical[coarse_eN,child_N*nInterpolationPoints + pN,:] = interpolationPointsOnFineElement_physical[fine_e.N,pN]
                #map physical interpolation points on coarse elements to coarse reference coordinates
                interpolationPointsOnCoarseElement_reference = numpy.zeros((coarseMesh.nElements_global,
                                                                              nChildrenMax_element*nInterpolationPoints,
                                                                              3),
                                                                             'd')
                J = numpy.zeros((coarseMesh.nElements_global,
                                   nChildrenMax_element*nInterpolationPoints,
                                   referenceElement.dim,
                                   referenceElement.dim),
                                  'd')
                invJ = numpy.zeros((coarseMesh.nElements_global,
                                      nChildrenMax_element*nInterpolationPoints,
                                      referenceElement.dim,
                                      referenceElement.dim),
                                     'd')
                detJ = numpy.zeros((coarseMesh.nElements_global,
                                      nChildrenMax_element*nInterpolationPoints),
                                     'd')
                interpolationPointsOnCoarseElement_reference_dummy = numpy.zeros((nChildrenMax_element*nInterpolationPoints,
                                                                                    3),
                                                                                  'd')
                coarseSpace.elementMaps.getJacobianValues(interpolationPointsOnCoarseElement_reference_dummy,
                                                          J,
                                                          invJ,
                                                          detJ)
                coarseSpace.elementMaps.getInverseValues(invJ,
                                                         interpolationPointsOnCoarseElement_physical,
                                                         interpolationPointsOnCoarseElement_reference)
                #get coarse scale basis function values at these reference points
                psi = numpy.zeros((coarseMesh.nElements_global,
                                     nChildrenMax_element*nInterpolationPoints,
                                     coarseSpace.referenceFiniteElement.localFunctionSpace.dim),
                                    'd')
                coarseSpace.getBasisValuesAtArray(interpolationPointsOnCoarseElement_reference,
                                                  psi)
                for coarse_eN in range(coarseMesh.nElements_global):
                    for fine_eN,fine_e in enumerate(children[coarse_eN]):
                        for i in fineSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                            I = fineSpace.dofMap.l2g[fine_e.N,i]
                            for j in coarseSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                                J = coarseSpace.dofMap.l2g[coarse_eN,j]
                                psi_j = psi[coarse_eN,fine_eN*nInterpolationPoints:(fine_eN+1)*nInterpolationPoints,j]
                                F_ij = fineSpaceInterpolationFunctionals[i](psi_j)
                                if  abs(F_ij) > 1.0e-16:#mwf orig F_ij > 1.e-16 changed to abs(F_ij)
                                    #mwf now account for nonconforming points?
                                    if isNonConformingOnCoarseGrid[fine_e.N,i] > 0:
                                        F_ij *= 0.5
                                        #have to allow for summing up values hit multiple times
                                        #only here
                                        if rbc.has_key((J,I)):
                                            rbc[(J,I)] += F_ij
                                            pbc[(I,J)] += F_ij
                                        else:
                                            rbc[(J,I)] = F_ij
                                            pbc[(I,J)] = F_ij
                                    else:
                                        rbc[(J,I)] = F_ij
                                        pbc[(I,J)] = F_ij
                                    rbcColumnIndeces[J].add(I)
                                    if fineDOFBoundaryConditions.global2freeGlobal.has_key(I):
                                        II = fineDOFBoundaryConditions.global2freeGlobal[I]*strideListList[l+1][cj]+offsetListList[l+1][cj]
                                        if coarseDOFBoundaryConditions.global2freeGlobal.has_key(J):
                                            JJ = coarseDOFBoundaryConditions.global2freeGlobal[J]*strideListList[l][cj]+offsetListList[l][cj]
                                            rColumnIndeces[JJ].add(II)
                                            r[(JJ,II)] = rbc[(J,I)]
                                            p[(II,JJ)] = pbc[(I,J)]
                #end coarse_eN
                for I in range(coarseSpace.dim):
                    for J in rbcColumnIndeces[I]:
                        rbcSum[I] += rbc[(I,J)]
                for I in range(offsetListList[l][cj],offsetListList[l][cj]+coarseDOFBoundaryConditions.nFreeDOF_global,strideListList[l][cj]):
                    for J in rColumnIndeces[I]:
                        rSum[I] += r[(I,J)]
                for I in range(coarseSpace.dim):
                    for J in rbcColumnIndeces[I]:
                        scaled_rbc[(I,J)] = rbc[(I,J)]/rbcSum[I]
                #now make real sparse matrices
                (rbc,rbczval) = SparseMatFromDict(coarseSpace.dim,fineSpace.dim,rbc)
                (scaled_rbc,scaled_rbczval) = SparseMatFromDict(coarseSpace.dim,fineSpace.dim,scaled_rbc)
                (pbc,pbczval) = SparseMatFromDict(fineSpace.dim,coarseSpace.dim,pbc)
                self.restrict_bcSumListDict[cj].append(rbcSum)
                self.restrict_bcListDict[cj].append(rbc)
                self.scaled_restrict_bcListDict[cj].append(scaled_rbc)
                self.rbczvalListDict[cj].append(rbczval)
                self.scaled_rbczvalListDict[cj].append(scaled_rbczval)
                self.prolong_bcListDict[cj].append(pbc)
                self.pbczvalListDict[cj].append(pbczval)
            (r,rzval) = SparseMatFromDict(coarse_nFreeDOF_global,fine_nFreeDOF_global,r)
            (p,pzval) = SparseMatFromDict(fine_nFreeDOF_global,coarse_nFreeDOF_global,p)
            self.restrictSumList.append(rSum)
            self.restrictList.append(r)
            self.rzvalList.append(rzval)
            self.prolongList.append(p)
            self.pzvalList.append(pzval)


########################################################################

class AdaptiveForwardIntegrator:
    """
    class that is responsible for basic process of integrating a problem
    forward in time given a VectorTranport Problem, Nonlinear Solver,
    and a Time Integration method
    """
    def __init__(self,mlvtran,mlnl,dtMeth,nOptions,stepExact=True,stepAdapt=True,resetAfterNLfail=True):
        """
        mlvTran  --- multilevel vector transport object for system being integrated
        mlNL     --- multilevel nonlinear solver to solve discrete system
        dtMet    --- time integration method to use
        nOptions --- configuration options
        """
        self.mlvTran   = mlvtran
        self.mlNL      = mlnl
        self.dtMeth    = dtMeth
        self.tLast     = None
        self.stepExact = stepExact
        self.DTSET     = None
        self.tstring   = None
        self.nOptions  = nOptions
        self.firstStep = True
        self.stepAdapt = stepAdapt
        self.maxDTfailures = 10
        self.maxNLfailures = 10

        #reset nonlinear solver guesss after NL failure?
        self.resetAfterNLfail = resetAfterNLfail
        if self.resetAfterNLfail:
            self.uListSave = []
            for l in range(len(self.mlvTran.uList)):
                self.uListSave.append(numpy.array(self.mlvTran.uList[l]))
        else:
            self.uListSave = None
    #end init
    def initialize(self,DTSET=None,t0=0.0):
        """

        """
        self.tLast = t0
        self.DTSET = DTSET #if adapting, this is dt0 only
        if self.dtMeth == TimeIntegration.NoIntegration:
            self.DTSET = 1.0
        self.mlvTran.chooseDT(self.DTSET)
        self.mlvTran.initializeTimeIntegration(t0)
        self.firstStep = True
        if self.stepAdapt == False:
            self.maxNLfailures = 0 #allow no nl failures if can't adapt
    def calculateSolution(self,tIn,tOut):
        """
        Move forward from time tIn to time tOut
        For now doesn't worry about potential mismatch between tIn and
        last time value used by model
        """
        t = tIn
        failedFlag = False
        nNLfailures = 0
        nDTfailures = 0
        if self.dtMeth == TimeIntegration.NoIntegration:
            self.mlvTran.chooseDT(self.DTSET)
            #mwf debug
            #print """fint t=%g tOut=%g DTSET=%s DT=%g """ % (t,tOut,self.DTSET,self.mlvTran.DT)

            failedFlag=self.mlNL.solveMultilevel(uList=self.mlvTran.uList,
                                                 rList=self.mlvTran.rList)

            self.tLast = tOut
            t = tOut
            self.mlvTran.updateTimeHistory(t)
            #mwf debug
            #print """after solve failedFlag= %s """ % failedFlag

        else:
            lastStep = False
            nlSolveFailed = False
            if self.firstStep == True:
                self.mlvTran.chooseDT(self.DTSET)
                if self.stepExact and abs(t+self.mlvTran.DT) > abs(tOut - t*1.0e-8):
                    self.mlvTran.resetDT()
                    self.mlvTran.chooseDT(tOut-t)
                    lastStep = True
                    #mwf debug
                    print """adaptfint adjusting DT0 to %g """ % self.mlvTran.DT
            while t < tOut and failedFlag== False:
                #assume dt already selected
                #mwf debug
                print """\nadaptfint t=%g tOut=%g DTSET=%s DT=%g """ % (t,tOut,self.DTSET,self.mlvTran.DT)
                self.writeProgress(t,self.mlvTran.DT,tOut)
                istage = 0
                nlSolveFailed = False
                while istage < self.nOptions.nStagesTime and nlSolveFailed == False:
                    if self.resetAfterNLfail:
                        #is there a better way to reset?
                        for l in range(len(self.mlvTran.uList)):
                            self.uListSave[l].flat[:] = self.mlvTran.uList[l].flat[:]
                    #if
                    nlSolveFailed= self.mlNL.solveMultilevel(uList=self.mlvTran.uList,
                                                             rList=self.mlvTran.rList)
                    if nlSolveFailed == False:
                        self.mlvTran.updateStage()
                        istage += 1
                    #mwf debug
                    print """adaptfint t=%g istage=%d failed=%s """ % (t,istage,nlSolveFailed)
                #end stages
                if nlSolveFailed:
                    nNLfailures += 1
                stepOk,lastStep = self.chooseDT(t,tOut,nlSolveFailed=nlSolveFailed)
                if nlSolveFailed and self.resetAfterNLfail:
                    for l in range(len(self.mlvTran.uList)):
                        self.mlvTran.uList[l].flat[:] = self.uListSave[l].flat[:] #reuse last time step

                #
                if stepOk and not nlSolveFailed:
                    t += self.mlvTran.DT
                    if abs(tOut-t) < 1.0e-8*t or lastStep:
                        t=tOut
                    self.mlvTran.updateTimeHistory(t)
                    nNLfailures = 0
                    nDTfailures = 0
                #step ok
                else:
                    if nlSolveFailed == False:
                        nDTfailures += 1
                    failedFlag = (nDTfailures > self.maxDTfailures or
                                  nNLfailures > self.maxNLfailures)

            #end while t < tOut
            self.tLast = t
        #end else do time integration

        #need to put in linear solver for profiling
        #Profiling.log(self.mlNL.info(),level=5)

        return failedFlag,t

    #end calculateSolution
    def writeProgress(self,tn,dt,T):
        """
        just echo to screen what new and final time levels are
        """
        import Profiling
        if Profiling.logLevel < 2:
            eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'
            if self.tstring != None:
                sys.stdout.write(eraseTime)
        sys.stdout.write('T= %12.5e, tn = ' % T)
        self.tstring = '%12.5e' % (tn+dt,)
        if Profiling.logLevel >= 2:
            self.tstring+='\n'
        sys.stdout.write(self.tstring)
        sys.stdout.flush()
    #end writeProgress

    def chooseDT(self,t,tOut,nlSolveFailed):
        """
        try to pick dt for next time step solve. If adapting, lets
        finest level time integrator determine if error is acceptable
        and pick the output time step

        otherwise, just sets time step to be uniform one

        includes adjustment if stepping exactly to output as well

        """
        lastStep = False
        stepOk   = True
        if not self.stepAdapt:
            self.mlvTran.chooseDT(self.DTSET)
            if self.stepExact and abs(t+self.mlvTran.DT) > abs(tOut - t*1.0e-8):
                self.mlvTran.resetDT()
                self.mlvTran.chooseDT(tOut-t)
                lastStep = True
                #mwf debug
                print """adaptfint adjusting DT to %g """ % self.mlvTran.DT
            #end stepExact
            return stepOk,lastStep
        #short circuit MultilevelTransport's choose step by asking finest level
        #integrator what dt it wants. Then tell MultilevelTransport to force
        #this dt
        self.mlvTran.modelList[-1].timeIntegration.setLastSolveFailed(nlSolveFailed)
        self.mlvTran.modelList[-1].timeIntegration.chooseDT()
        #should I change integrators dt?
        stepOk     = self.mlvTran.modelList[-1].timeIntegration.lastStepErrorOk()
        self.DTSET = self.mlvTran.modelList[-1].timeIntegration.DT
        if stepOk == False or nlSolveFailed:
            self.mlvTran.resetDT()
        self.mlvTran.chooseDT(self.DTSET)
        if self.stepExact and abs(t+self.mlvTran.DT) > abs(tOut - t*1.0e-8):
            self.mlvTran.resetDT()
            self.mlvTran.chooseDT(tOut-t)
            lastStep = True
        return stepOk,lastStep

class AdaptiveBackwardEuler(TimeIntegration.BackwardEuler):
    def __init__(self,vectorTransport,runCFL=0.9,atol=1.0e-3,rtol=1.0e-3,errorType='LI_global',
                 dtRatioMax=2,dtRatioMin=0.1):
        TimeIntegration.BackwardEuler.__init__(self,vectorTransport,runCFL=runCFL)
        self.atol = atol
        self.rtol = rtol
        self.lastSolveFailed = False
        self.lastStepErrorAcc= True
        self.dtRatioMax = dtRatioMax ; self.dtRatioMin = dtRatioMin
        self.errorType = errorType #L2_global, LI_global
        self.m_pred = {}
        self.mt_last= {} ; self.mt_tmp = {}
        for ci in vectorTransport.coefficients.mass.keys():
            self.m_pred[ci] = numpy.array(vectorTransport.q[('m',ci)])
            self.mt_last[ci]= numpy.array(vectorTransport.q[('mt',ci)])
            self.mt_tmp[ci] = numpy.array(vectorTransport.q[('mt',ci)])
        #end for
        self.error_weights = {}
    def calculateElementCoefficients(self,q):
        TimeIntegration.BackwardEuler.calculateElementCoefficients(self,q)
        for ci in self.m_last.keys():
            self.mt_tmp[ci].flat[:] = q[('mt',ci)].flat[:]
            self.error_weights[ci]  = q[('dV_u',ci)]

    def chooseDT(self):
        """
        Modify self.DT
        mwf needs to be checked
        """
        if self.lastSolveFailed == True:
            self.DT *= self.dtRatioMin #0.25
            #should I reset lastSolveFailed here?
            print "AdaptBackwardEuler chooseDT NLfail DT=%s " % (self.DT)
        else:
            predFactor = 0.5

            #compute error
            error = 0.0 ;
            if self.errorType == 'L2_global':
                for ci in self.m_pred.keys():
                #compute predictor at current DT?
                    self.m_pred[ci].flat[:]  = self.mt_last[ci].flat[:]
                    self.m_pred[ci].flat[:] *= self.DT
                    self.m_pred[ci].flat[:] += self.m_last[ci].flat[:]
                    diffci = predFactor*Norms.L2errorSFEM(self.error_weights[ci],
                                                          self.m_tmp[ci],
                                                          self.m_pred[ci])
                    normci = Norms.L2normSFEM(self.error_weights[ci],
                                              self.m_tmp[ci])

                    errorci= diffci/(normci*self.rtol + self.atol)
                    error  = max(errorci,error)
                #end ci
            elif self.errorType == 'LI_global':
                for ci in self.m_pred.keys():
                    #compute predictor at current DT?
                    self.m_pred[ci].flat[:]  = self.mt_last[ci].flat[:]
                    self.m_pred[ci].flat[:] *= self.DT
                    self.m_pred[ci].flat[:] += self.m_last[ci].flat[:]
                    diffci = predFactor*Norms.LIerrorSFEM(self.error_weights[ci],
                                                          self.m_tmp[ci],
                                                          self.m_pred[ci])
                    normci = max(abs(self.m_tmp[ci].flat))

                    errorci= diffci/(normci*self.rtol + self.atol)
                    error  = max(errorci,error)
                #end ci
            else:
                assert False, "unknown error type = %s " % self.errorType
            self.lastStepErrorAcc = error < 1.0
            D = max(error,1.0e-4)
            safety = 0.9
            order  = 1.0
            tol    = 1.0#weighted type error
            dtRatN = safety*((tol/D)**(1./(order+1.)))
            dtRat  = max(self.dtRatioMin,min(dtRatN,self.dtRatioMax))
            self.DT *= dtRat

            print "AdaptBackwardEuler chooseDT error=%s dtRat=%s DT=%s " % (error,dtRat,self.DT)

    def updateTimeHistory(self,resetFromDOF=False,vectorTransport=None):
        TimeIntegration.BackwardEuler.updateTimeHistory(self,resetFromDOF,vectorTransport)
        for ci in self.m_last.keys():
            self.mt_last[ci].flat[:] = self.mt_tmp[ci].flat

    #mwf for now add separate call to notify integration method so that I don't disturb
    #original code?
    def setLastSolveFailed(self,lastSolveFailed):
        """
        tell integrator last attempted time step failed or not
        """
        self.lastSolveFailed = lastSolveFailed
    def lastStepErrorOk(self):
        """
        was the last time step acceptable or not
        """
        return self.lastStepErrorAcc

########################################################################
#begin working on limiting procedures for DG, CG
########################################################################
def stupidSort(A,hit=-12345):
    n = len(A)
    B = numpy.array(A)
    for islot in range(n):
        print "islot=%d " % (islot)
        for i in range(islot+1,n):
            print "i=%d" % (i)
            if A[i] > A[islot]:
                print "A[%d]=%s > A[%d]=%s " % (i,A[i],islot,A[islot])
                tmp = A[islot]
                A[islot]=A[i]
                A[i]= tmp
            print "after i=%d A =%s " % (i,A)
        print "after islot=%d A =%s " % (islot,A)
        #
    #islot
    return A
def checkOrder(A):
    n = len(A)
    for i in range(n-1):
        assert A[i+1] <= A[i], "out of order i=%d A[i]=%s A[i+1]=%s " % (i,A[i],A[i+1])
    #

########################################################################
#numerical fluxes
import NumericalFlux,cnumericalFlux
class RusanovNumericalFlux_Diagonal_Diffusion_LDG(NumericalFlux.Advection_DiagonalUpwind_Diffusion_LDG):
    """
    apply numerical flus :math:`f_{num}(a,b) = 1/2(f(a)+f(b)-\bar{\lambda}(b-a)` where
    :math:`\lambda >= max |f^{\prime}|` for :math:`a<= u <= b`
    this one applies flux to each component of flux separately
    """
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        NumericalFlux.Advection_DiagonalUpwind_Diffusion_LDG.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                      getAdvectiveFluxBoundaryConditions,
                                                                      getDiffusiveFluxBoundaryConditions)
        self.safetyFactor=1.1
    def calculateInteriorNumericalFlux(self,q,ebq,ebq_global):
        for ci in range(self.nc):
            #initialize  to current u so that non-dirichlet BC flux will be correct
            self.ebq[('u',ci)][:] = ebq[('u',ci)]
            for (eN,ebN,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                      self.DOFBoundaryConditionsDictList[ci].values(),
                                      self.DOFBoundaryPointDictList[ci].values()):
                self.ebq[('u',ci)][eN,ebN,k]=g(x,self.vt.T)
        self.vt.coefficients.evaluate(self.vt.T,self.ebq)
        #calculate diffusive velocity projection--requires interior and exterior loops
        for ck in self.vt.coefficients.potential.keys():
            assert (len(self.vt.coefficients.potential[ck].keys()) == 1 and
                    self.vt.coefficients.potential[ck].keys()[0] == ck), "No off-diagonal dependence in phi currently allowed in LDG"
            if self.vt.coefficients.potential[ck][ck] == 'u':
                self.ebq[('phi',ck)].flat[:]=self.ebq[('u',ck)].flat
                self.ebq[('dphi',ck,ck)].flat[:]=1.0
            self.b[ck].flat[:]=0.0
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
                                                                    self.ebq[('phi',ck)],
                                                                    ebq[('phi',ck)],
                                                                    ebq[('dphi',ck,ck)],
                                                                    self.phi_trace[ck],
                                                                    self.dphi_trace_left[ck])
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
            cfemIntegrals.calculateVelocityQuadrature_MixedForm(self.A_inv,
                                                                self.b[ck],
                                                                ebq[('v',ck)],
                                                                self.ebqV[ck],
                                                                q[('v',ck)],
                                                                self.qV[ck])
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateInteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                              self.mesh.interiorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              ebq['n'],
                                                                              ebq[('u',ci)],
                                                                              ebq[('f',ci)],
                                                                              ebq[('df',ci,ci)],
                                                                              q[('df',ci,ci)],
                                                                              ebq_global[('advectiveFlux',ci)],
                                                                              ebq_global[('dadvectiveFlux_left',ci,ci)],
                                                                              ebq_global[('dadvectiveFlux_right',ci,ci)])

            ##\todo fix numerical diffusive flux for multiple diffusion terms
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    ebq_global[('diffusiveFlux',ck,ci)].flat[:]=0.0
                    cfemIntegrals.calculateInteriorNumericalDiffusiveFlux_LDG_upwind(self.mesh.interiorElementBoundariesArray,
                                                                                     self.mesh.elementBoundaryElementsArray,
                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                     ebq['n'],
                                                                                     ebq[('u',ck)],
                                                                                     ebq[('a',ci,ck)],
                                                                                     ebq[('phi',ck)],
                                                                                     self.ebqV[ck],
                                                                                     ebq_global[('penalty')],
                                                                                     ebq_global[('diffusiveFlux',ck,ci)])
            #ck
        #ci
    #interior numerical flux
    def calculateExteriorNumericalFlux(self,inflowFlag,q,ebq,ebq_global):
        for ci in range(self.nc):
            #mwf this initialization isn't in Advection_Diagonal_LDG class
            #initialize  to current u so that non-dirichlet BC flux will be correct
            self.ebq[('u',ci)][:] = ebq[('u',ci)]
            for (eN,ebN,k),g,x in zip(self.DOFBoundaryConditionsDictList[ci].keys(),
                                      self.DOFBoundaryConditionsDictList[ci].values(),
                                      self.DOFBoundaryPointDictList[ci].values()):
                #mwf debug
                #print "Rusanov_DiagonalUpwind computing bcs eN=%d ebN=%d k=%d g=%s" % (eN,ebN,k,g(x,self.vt.T))
                self.ebq[('u',ci)][eN,ebN,k]=g(x,self.vt.T)
        self.vt.coefficients.evaluate(self.vt.T,self.ebq)
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cnumericalFlux.calculateExteriorNumericalAdvectiveFluxRusanov(self.safetyFactor,
                                                                              self.mesh.exteriorElementBoundariesArray,
                                                                              self.mesh.elementBoundaryElementsArray,
                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                              self.isDOFBoundary[ci],
                                                                              inflowFlag[ci],
                                                                              ebq['n'],
                                                                              self.ebq[('u',ci)],
                                                                              self.ebq[('f',ci)],
                                                                              self.ebq[('df',ci,ci)],
                                                                              ebq[('u',ci)],
                                                                              ebq[('f',ci)],
                                                                              ebq[('df',ci,ci)],
                                                                              q[('df',ci,ci)],
                                                                              ebq_global[('advectiveFlux',ci)],
                                                                              ebq_global[('dadvectiveFlux_left',ci,ci)])
            for ck in range(self.nc):
                if ebq.has_key(('a',ci,ck)):
                    cfemIntegrals.calculateExteriorNumericalDiffusiveFlux_LDG_upwind(self.mesh.exteriorElementBoundariesArray,
                                                                                     self.mesh.elementBoundaryElementsArray,
                                                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                     ebq['n'],
                                                                                     ebq[('u',ck)],
                                                                                     ebq[('a',ci,ck)],
                                                                                     self.ebq[('phi',ck)],
                                                                                     ebq[('phi',ck)],
                                                                                     self.ebqV[ck],
                                                                                     ebq_global[('penalty')],
                                                                                     ebq_global[('diffusiveFlux',ck,ci)])
        for ci in range(self.nc):
            if ebq.has_key(('f',ci)):
                cfemIntegrals.calculateExteriorInflowNumericalAdvectiveFlux(self.mesh.exteriorElementBoundariesArray,
                                                                            self.mesh.elementBoundaryElementsArray,
                                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                            inflowFlag[ci],
                                                                            ebq[('inflowFlux',ci)],
                                                                            ebq['n'],
                                                                            ebq[('f',ci)],
                                                                            ebq[('df',ci,ci)],
                                                                            ebq_global[('advectiveFlux',ci)],
                                                                            ebq_global[('dadvectiveFlux_left',ci,ci)])

#

########################################################################
#try special nonlinear solver to speed up rkdg
import NonlinearSolvers
import LinearSolvers
class SSPRKNewton(NonlinearSolvers.Newton):
    """
    A simple iterative solver that is Newton's method
    if you give it the right Jacobian
    """
    def __init__(self,
                 linearSolver,
                 F,J=None,du=None,par_du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True,
                 directSolver=False,
                 EWtol=True,
                 maxLSits = 100):
        import copy
        self.par_du = par_du
        if par_du != None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolvers.Newton.__init__(self,
                                         linearSolver,
                                         F,J,du,par_du,
                                         rtol_r,
                                         atol_r,
                                         rtol_du,
                                         atol_du,
                                         maxIts,
                                         norm,
                                         convergenceTest,
                                         computeRates,
                                         printInfo,
                                         fullNewton,
                                         directSolver,
                                         EWtol,
                                         maxLSits)
        self.isFactored = False

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        import Viewers
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """
        import sys
        if self.linearSolver.computeEigenvalues:
            self.u0[:]=u
        r=self.solveInitialize(u,r,b)
        if par_u != None:
            #no overlap
            #par_r.scatter_reverse_add()
            #no overlap or overlap (until we compute norms over only owned dof)
            par_r.scatter_forward_insert()
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        while (not self.converged(r) and
               not self.failed()):
            print "SSPRKNewton it", self.its, "norm(r)", self.norm_r
            if self.updateJacobian or self.fullNewton and not self.isFactored:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                print numpy.transpose(self.J)
                if self.linearSolver.computeEigenvalues:
                    self.JLast[:]=self.J
                    self.J_t_J[:]=self.J
                    self.J_t_J *= numpy.transpose(self.J)
                    self.JLsolver.prepare()
                    self.JLsolver.calculateEigenvalues()
                    self.norm_2_J_current = sqrt(max(self.JLsolver.eigenvalues_r))
                    self.norm_2_Jinv_current = 1.0/sqrt(min(self.JLsolver.eigenvalues_r))
                    self.kappa_current = self.norm_2_J_current*self.norm_2_Jinv_current
                    self.betaK_current = self.norm_2_Jinv_current
                self.linearSolver.prepare(b=r)
                self.isFactored = True
            self.du[:]=0.0
            if not self.directSolver:
                if self.EWtol:
                    self.setLinearSolverTolerance(r)
            self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
            #print self.du
            u-=self.du
            if par_u != None:
                par_u.scatter_forward_insert()
            self.computeResidual(u,r,b)
            #no overlap
            #print "local r",r
            if par_r != None:
                #no overlap
                #par_r.scatter_reverse_add()
                par_r.scatter_forward_insert()
            #print "global r",r
            if self.linearSolver.computeEigenvalues:
                #approximate Lipschitz constant of J
                self.F.getJacobian(self.dJ_t_dJ)
                self.dJ_t_dJ-=self.JLast
                self.dJ_t_dJ *= numpy.transpose(self.dJ_t_dJ)
                self.dJLsolver.prepare()
                self.dJLsolver.calculateEigenvalues()
                self.norm_2_dJ_current = sqrt(max(self.dJLsolver.eigenvalues_r))
                self.etaK_current = self.W*self.norm(self.du)
                self.gammaK_current = self.norm_2_dJ_current/self.etaK_current
                self.gammaK_max = max(self.gammaK_current,self.gammaK_max)
                self.norm_r_hist.append(self.W*self.norm(r))
                self.norm_du_hist.append(self.W*self.norm(self.du))
                if self.its  == 1:
#                     print "max(|du|) ",max(numpy.absolute(self.du))
#                     print self.du[0]
#                     print self.du[-1]
                    self.betaK_0 = self.betaK_current
                    self.etaK_0 = self.etaK_current
                if self.its  == 2:
                    self.betaK_1 = self.betaK_current
                    self.etaK_1 = self.etaK_current
                print "it = ",self.its
                print "beta(|Jinv|)  ",self.betaK_current
                print "eta(|du|)     ",self.etaK_current
                print "gamma(Lip J') ",self.gammaK_current
                print "gammaM(Lip J')",self.gammaK_max
                print "kappa(cond(J))",self.kappa_current
                if self.betaK_current*self.etaK_current*self.gammaK_current <= 0.5:
                    print "r         ",(1.0+sqrt(1.0-2.0*self.betaK_current*self.etaK_current*self.gammaK_current))/(self.betaK_current*self.gammaK_current)
                if self.betaK_current*self.etaK_current*self.gammaK_max <= 0.5:
                    print "r_max     ",(1.0+sqrt(1.0-2.0*self.betaK_current*self.etaK_current*self.gammaK_max))/(self.betaK_current*self.gammaK_max)
                print "lambda_max",max(self.linearSolver.eigenvalues_r)
                print "lambda_i_max",max(self.linearSolver.eigenvalues_i)
                print "norm_J",self.norm_2_J_current
                print "lambda_min",min(self.linearSolver.eigenvalues_r)
                print "lambda_i_min",min(self.linearSolver.eigenvalues_i)
            if self.lineSearch:
                norm_r_cur = self.norm(r)
                norm_r_last = 2.0*norm_r_cur
                ls_its = 0
                #print norm_r_cur,self.atol_r,self.rtol_r
                while ( (norm_r_cur >= 0.99 * self.norm_r + self.atol_r) and
                        (ls_its < self.maxLSits) and
                        norm_r_cur/norm_r_last < 1.0):
                    self.convergingIts = 0
                    ls_its +=1
                    self.du *= 0.5
                    u += self.du
                    if par_u != None:
                        par_u.scatter_forward_insert()
                    self.computeResidual(u,r,b)
                    #no overlap
                    if par_r != None:
                        #no overlap
                        #par_r.scatter_reverse_add()
                        par_r.scatter_forward_insert()
                    norm_r_last = norm_r_cur
                    norm_r_cur = self.norm(r)
                    print """ls #%d norm_r_cur=%s atol=%g rtol=%g""" % (ls_its,
                                                                        norm_r_cur,
                                                                        self.atol_r,
                                                                        self.rtol_r)
                if ls_its > 0:
                    logEvent("Linesearches = %i" % ls_its,level=3)
        else:
            if self.linearSolver.computeEigenvalues:
                if self.betaK_0*self.etaK_0*self.gammaK_max <= 0.5:
                    print "r_{-,0}     ",(1.0+sqrt(1.0-2.0*self.betaK_0*self.etaK_0*self.gammaK_max))/(self.betaK_0*self.gammaK_max)
                if self.betaK_1*self.etaK_1*self.gammaK_max <= 0.5 and self.its > 1:
                    print "r_{-,1}     ",(1.0+sqrt(1.0-2.0*self.betaK_1*self.etaK_1*self.gammaK_max))/(self.betaK_1*self.gammaK_max)
                print "beta0*eta0*gamma ",self.betaK_0*self.etaK_0*self.gammaK_max
                if Viewers.viewerType == 'gnuplot':
                    max_r = max(1.0,max(self.linearSolver.eigenvalues_r))
                    max_i = max(1.0,max(self.linearSolver.eigenvalues_i))
                    for lambda_r,lambda_i in zip(self.linearSolver.eigenvalues_r,self.linearSolver.eigenvalues_i):
                        Viewers.datFile.write("%12.5e %12.5e \n" % (lambda_r/max_r,lambda_i/max_i))
                    Viewers.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with points title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                      Viewers.datFilename,
                                                                                                      Viewers.plotNumber,
                                                                                                      'scaled eigenvalues')
                    Viewers.cmdFile.write(cmd)
                    Viewers.viewerPipe.write(cmd)
                    Viewers.newPlot()
                    Viewers.newWindow()
                    for it,r in zip(range(len(self.norm_r_hist)),self.norm_r_hist):
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,log(r/self.norm_r_hist[0])))
                    Viewers.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                      Viewers.datFilename,
                                                                                                      Viewers.plotNumber,
                                                                                                      'log(r)/log(r0) history')
                    Viewers.cmdFile.write(cmd)
                    Viewers.viewerPipe.write(cmd)
                    Viewers.newPlot()
                    Viewers.newWindow()
                    for it,du in zip(range(len(self.norm_du_hist)),self.norm_du_hist):
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,log(du/self.norm_du_hist[0])))
                    Viewers.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                      Viewers.datFilename,
                                                                                                      Viewers.plotNumber,
                                                                                                      'log(du) history')
                    Viewers.cmdFile.write(cmd)
                    Viewers.viewerPipe.write(cmd)
                    Viewers.newPlot()
                    Viewers.newWindow()
                #raw_input("wait")
            return self.failedFlag
        #else
    #def
    def resetFactorization(self,needToRefactor=True):
        self.isFactored = not needToRefactor
#class
