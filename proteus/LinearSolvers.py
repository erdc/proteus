r"""
A hierarchy of classes for linear algebraic system solvers.

.. inheritance-diagram:: proteus.LinearSolvers
   :parts: 1
"""
import proteus
from LinearAlgebraTools import *
import FemTools
import lapackWrappers
import superluWrappers
import TransportCoefficients
import cfemIntegrals
import Quadrature
import petsc4py
from petsc4py import PETSc as p4pyPETSc
from math import *
from .Profiling import logEvent

class LinearSolver:
    """ The base class for linear solvers. """
    def __init__(self,
                 L,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = False):
        self.solverName = "Base Class"
        self.L = L
        self.n = L.shape[0]
        self.du = Vec(self.n)
        self.rtol_r=rtol_r
        self.atol_r=atol_r
        self.rtol_du=rtol_du
        self.atol_du=atol_du
        self.maxIts=maxIts
        self.its=0
        self.solveCalls = 0
        self.recordedIts=0
        self.solveCalls_failed = 0
        self.recordedIts_failed=0
        self.rReductionFactor=0.0
        self.duReductionFactor=0.0
        self.rReductionFactor_avg=0.0
        self.duReductionFactor_avg=0.0
        self.rReductionOrder=0.0
        self.rReductionOrder_avg=0.0
        self.duReductionOrder=0.0
        self.duReductionOrder_avg=0.0
        self.ratio_r_current = 1.0
        self.ratio_r_solve = 1.0
        self.ratio_du_solve = 1.0
        self.last_log_ratio_r = 1.0
        self.last_log_ratior_du = 1.0
        self.convergenceTest = convergenceTest
        self.computeRates = computeRates
        self.computeEigenvalues=False
        self.printInfo = printInfo
        self.norm = l2Norm
        self.convergenceHistoryIsCorrupt=False
        self.norm_r0=0.0
        self.norm_r=0.0
        self.norm_du=0.0
        self.r=None
        self.leftEigenvectors=None
        self.rightEigenvectors=None
        self.eigenvalues_r=None
        self.eigenvalues_i=None
        self.work=None
        #need some information about parallel setup?
        self.par_fullOverlap = True #whether or not partitioning has full overlap
        #for petsc interface
        self.xGhosted = None
        self.b=None
    def setResTol(self,rtol,atol):
        self.rtol_r = rtol
        self.atol_r = atol
    def prepare(self,b=None):
        pass
    def solve(self,u,r=None,b=None,par_u=None,par_b=None):
        pass
    def calculateEigenvalues(self):
        pass
    def computeResidual(self,u,r,b,initialGuessIsZero=False):
        if initialGuessIsZero:
            r[:]=b
            r*=(-1.0)
        else:
            if type(self.L).__name__ == 'ndarray':
                r[:] = numpy.dot(u,self.L)
            elif type(self.L).__name__ == 'SparseMatrix':
                self.L.matvec(u,r)
            r-=b
    def solveInitialize(self,u,r,b,initialGuessIsZero=True):
        if r == None:
            if self.r is None:
                self.r = Vec(self.n)
            r=self.r
        else:
            self.r=r
        if b is None:
            if self.b is None:
                self.b = Vec(self.n)
            b=self.b
        else:
            self.b=b
        self.computeResidual(u,r,b,initialGuessIsZero)
        self.its = 0
        self.norm_r0 = self.norm(r)
        self.norm_r = self.norm_r0
        self.ratio_r_solve = 1.0
        self.ratio_du_solve = 1.0
        self.last_log_ratio_r = 1.0
        self.last_log_ratior_du = 1.0
        self.convergenceHistoryIsCorrupt=False
        return (r,b)
    def computeConvergenceRates(self):
        if self.convergenceHistoryIsCorrupt:
            return
        else:
            if self.its > 0:
                if self.norm_r < self.lastNorm_r:
                    self.ratio_r_current = self.norm_r/self.lastNorm_r
                else:
                    self.convergenceHistoryIsCorrupt=True
                    return
                if self.ratio_r_current > 1.0e-100:
                    log_ratio_r_current = log(self.ratio_r_current)
                else:
                    self.convergenceHistoryIsCorrupt
                    return
                self.ratio_r_solve *= self.ratio_r_current
                self.rReductionFactor = pow(self.ratio_r_solve,1.0/self.its)
                if self.its > 1:
                    self.rReductionOrder = log_ratio_r_current/ \
                                           self.last_log_ratio_r
                    if self.norm_du < self.lastNorm_du:
                        ratio_du_current = self.norm_du/self.lastNorm_du
                    else:
                        self.convergenceHistoryIsCorrupt=True
                        return
                    if ratio_du_current > 1.0e-100:
                        log_ratio_du_current = log(ratio_du_current)
                    else:
                        self.convergenceHistoryIsCorrupt=True
                        return
                    self.ratio_du_solve *= ratio_du_current
                    self.duReductionFactor = pow(self.ratio_du_solve,
                                                 1.0/(self.its-1))
                    if self.its > 2:
                        self.duReductionOrder = log_ratio_du_current/ \
                                                self.last_log_ratio_du
                    self.last_log_ratio_du = log_ratio_du_current
                self.last_log_ratio_r = log_ratio_r_current
                self.lastNorm_du = self.norm_du
            self.lastNorm_r = self.norm_r
    def converged(self,r):
        convergedFlag = False
        self.norm_r = self.norm(r)
        self.norm_du = self.norm(self.du)
        if self.computeRates ==  True:
            self.computeConvergenceRates()
        if self.convergenceTest == 'its':
            if self.its == self.maxIts:
                convergedFlag = True
        elif self.convergenceTest == 'r':
            if (self.its != 0 and
                self.norm_r < self.rtol_r*self.norm_r0 + self.atol_r):
                convergedFlag = True
        if convergedFlag == True and self.computeRates == True:
            self.computeAverages()
        if self.printInfo == True:
            print self.info()
        return convergedFlag
    def failed(self):
        failedFlag = False
        if self.its == self.maxIts:
            self.solveCalls_failed+=1
            self.recordedIts_failed+=self.its
            failedFlag = True
        self.its+=1
        return failedFlag
    def computeAverages(self):
        self.recordedIts+=self.its
        if self.solveCalls == 0:
            self.rReductionFactor_avg = self.rReductionFactor
            self.duReductionFactor_avg = self.duReductionFactor
            self.rReductionOrder_avg = self.rReductionOrder
            self.duReductionOrder_avg = self.duReductionOrder
            self.solveCalls+=1
        else:
            self.rReductionFactor_avg*=self.solveCalls
            self.rReductionFactor_avg+=self.rReductionFactor
            self.duReductionFactor_avg*=self.solveCalls
            self.duReductionFactor_avg+=self.duReductionFactor
            self.rReductionOrder_avg*=self.solveCalls
            self.rReductionOrder_avg+=self.rReductionOrder
            self.duReductionOrder_avg*=self.solveCalls
            self.duReductionOrder_avg+=self.duReductionOrder
            self.solveCalls +=1
            self.rReductionFactor_avg/=self.solveCalls
            self.duReductionFactor_avg/=self.solveCalls
            self.rReductionOrder_avg/=self.solveCalls
            self.duReductionOrder_avg/=self.solveCalls
    def info(self):
        self.infoString  = "************Start Linear Solver Info************\n"
        self.infoString += "its                   =  %i \n" % self.its
        self.infoString += "r reduction factor    = %12.5e\n" % self.rReductionFactor
        self.infoString += "du reduction factor   = %12.5e\n" % self.duReductionFactor
        self.infoString += "r reduction order     = %12.5e\n" % self.rReductionOrder
        self.infoString += "du reduction order    = %12.5e\n" % self.duReductionOrder
        self.infoString += "<r reduction factor>  = %12.5e\n" % self.rReductionFactor_avg
        self.infoString += "<du reduction factor> = %12.5e\n" % self.duReductionFactor_avg
        self.infoString += "<r reduction order>   = %12.5e\n" % self.rReductionOrder_avg
        self.infoString += "<du reduction order>  = %12.5e\n" % self.duReductionOrder_avg
        self.infoString += "total its             =  %i \n" % self.recordedIts
        self.infoString += "total its             =  %i \n" % self.recordedIts
        self.infoString += "solver calls          =  %i \n" % self.solveCalls
        self.infoString += "failures              =  %i \n" % self.solveCalls_failed
        self.infoString += "failed its            =  %i \n" % self.recordedIts_failed
        self.infoString += "maxIts                =  %i \n" % self.maxIts
        self.infoString += "convergenceTest       =  %s \n" % self.convergenceTest
        self.infoString += "atol_r                = %12.5e \n" % self.atol_r
        self.infoString += "rtol_r                = %12.5e \n" % self.rtol_r
        self.infoString += "norm(r0)              = %12.5e \n" % self.norm_r0
        self.infoString += "norm(r)               = %12.5e \n" % self.norm_r
        self.infoString += "norm(du)              = %12.5e \n" % self.norm_du
        if self.convergenceHistoryIsCorrupt:
            self.infoString += "HISTORY IS CORRUPT!!! \n"
        self.infoString += "************End Linear Solver Info************\n"
        return self.infoString
    def printPerformance(self):
        pass

    #petsc preconditioner interface
    def setUp(self, pc):
        self.prepare()
    def apply(self,pc,x,y):
        if self.xGhosted == None:
            self.xGhosted = self.par_b.duplicate()
            self.yGhosted = self.par_b.duplicate()
        self.xGhosted.setArray(x.getArray())
        self.xGhosted.ghostUpdateBegin(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.xGhosted.ghostUpdateEnd(p4pyPETSc.InsertMode.INSERT,p4pyPETSc.ScatterMode.FORWARD)
        self.yGhosted.zeroEntries()
        with self.yGhosted.localForm() as ylf,self.xGhosted.localForm() as xlf:
            self.solve(u=ylf.getArray(),b=xlf.getArray(),initialGuessIsZero=True)
        y.setArray(self.yGhosted.getArray())

class LU(LinearSolver):
    """
    A wrapper for pysparse's wrapper for superlu.
    """
    def __init__(self,L,computeEigenvalues=False,computeEigenvectors=None):
        import copy
        LinearSolver.__init__(self,L)
        if type(L).__name__ == 'SparseMatrix':
            self.sparseFactor = superluWrappers.SparseFactor(self.n)
        elif type(L).__name__ == 'ndarray':#mwf was array
            self.denseFactor = lapackWrappers.DenseFactor(self.n)
        self.solverName = "LU"
        self.computeEigenvalues = computeEigenvalues or (computeEigenvectors != None)
        if computeEigenvectors in ['left','both']:
            self.leftEigenvectors=numpy.zeros((self.n,self.n),'d')
            self.JOBVL='V'
        else:
            self.JOBVL='N'
        if computeEigenvectors in ['right','both']:
            self.rightEigenvectors=numpy.zeros((self.n,self.n),'d')
            self.JOBVR='V'
        else:
            self.JOBVR='N'
        if computeEigenvalues or computeEigenvectors != None:
            self.Leig=copy.deepcopy(L)
            self.work=numpy.zeros((self.n*5,),'d')
            self.eigenvalues_r = numpy.zeros((self.n,),'d')
            self.eigenvalues_i = numpy.zeros((self.n,),'d')
    def prepare(self,b=None):
        if type(self.L).__name__ == 'SparseMatrix':
            superluWrappers.sparseFactorPrepare(self.L,self.sparseFactor)
        elif type(self.L).__name__ == 'ndarray':
            if self.computeEigenvalues:
                self.Leig[:]=self.L
                self.calculateEigenvalues()
            lapackWrappers.denseFactorPrepare(self.n,
                                              self.L,
                                              self.denseFactor)
        #
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        self.du[:]=u
        self.converged(r)
        self.failed()
        u[:]=b
        if type(self.L).__name__ == 'SparseMatrix':
            superluWrappers.sparseFactorSolve(self.sparseFactor,u)
        elif type(self.L).__name__ == 'ndarray':
            lapackWrappers.denseFactorSolve(self.n,
                                            self.L,
                                            self.denseFactor,
                                            u)
        self.computeResidual(u,r,b)
        self.du -= u
        self.converged(r)
    def calculateEigenvalues(self):
        if type(self.L).__name__ == 'ndarray':
            lapackWrappers.denseCalculateEigenvalues(self.JOBVL,
                                                     self.JOBVR,
                                                     self.n,
                                                     self.Leig,
                                                     self.n,
                                                     self.eigenvalues_r,
                                                     self.eigenvalues_i,
                                                     self.leftEigenvectors,
                                                     self.n,
                                                     self.rightEigenvectors,
                                                     self.n,
                                                     self.work,
                                                     5*self.n)
            eigen_mags = numpy.sqrt(self.eigenvalues_r**2 + self.eigenvalues_i**2)
            logEvent("Minimum eigenvalue magnitude"+`eigen_mags.min()`)
            logEvent("Maximum eigenvalue magnitude"+`eigen_mags.max()`)
            logEvent("Minimum real part of eigenvalue "+`self.eigenvalues_r.min()`)
            logEvent("Maximum real part of eigenvalue "+`self.eigenvalues_r.max()`)
            logEvent("Minimum complex part of eigenvalue "+`self.eigenvalues_i.min()`)
            logEvent("Maximum complex part of eigenvalue "+`self.eigenvalues_i.max()`)
class PETSc(LinearSolver):
    def __init__(self,L,par_L,prefix=None):
        import flcbdfWrappers
        LinearSolver.__init__(self,L)
        assert type(L).__name__ == 'SparseMatrix', "PETSc can only be called with a local sparse matrix"
        self.solverName  = "PETSc"
        if prefix == None:
            self.ksp = flcbdfWrappers.KSP(par_L)
        else:
            assert isinstance(prefix,str)
            prefix.replace(' ','_')
            self.ksp = flcbdfWrappers.KSP(par_L,prefix)
        self.par_L = par_L
        self.par_fullOverlap = True
    def prepare(self,b=None):
        overlap = 1
        if self.par_fullOverlap == False:
            overlap = 0
        self.ksp.prepare(self.L,self.par_L,overlap)
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        self.ksp.solve(par_u.cparVec,par_b.cparVec)
    def useTrueResidualTest(self,par_u):
        if par_u != None:
            self.ksp.useTrueResidualConvergence(par_u.cparVec)
    def printPerformance(self):
        self.ksp.info()



class KSP_petsc4py(LinearSolver):
    """ A class that interfaces Proteus with PETSc's KSP. """
    def __init__(self,L,par_L,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 maxIts  = 100,
                 norm    = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = False,
                 prefix=None,
                 Preconditioner=None,
                 connectionList=None,
                 linearSolverLocalBlockSize=1,
                 outputResults = True,
                 bdyNullSpace=False):
        """ Initialize a petsc4py KSP object.
        
        Parameters
        -----------
        L : :class: `.superluWrappers.SparseMatrix`
        par_L :  :class: `.LinearAlgebraTools.ParMat_petsc4py`
        rtol_r : float
        atol_r : float
        maxIts : int
        norm :   norm type
        convergenceTest : 
        computeRates: bool
        printInfo : bool
        prefix : bool
        Preconditioner : :class: `.LinearSolvers.KSP_Preconditioner`
        connectionList : 
        linearSolverLocalBlockSize : int
        """
        LinearSolver.__init__(self,
                              L,
                              rtol_r=rtol_r,
                              atol_r=atol_r,
                              maxIts=maxIts,
                              norm = l2Norm,
                              convergenceTest=convergenceTest,
                              computeRates=computeRates,
                              printInfo=printInfo)
        assert type(L).__name__ == 'SparseMatrix', "petsc4py PETSc can only be called with a local sparse matrix"
        assert isinstance(par_L,ParMat_petsc4py)
        # ARB - WIP: this self.pc and self.preconditioner stuff is confusing to read
        # and not necessary.  This should be refactored before merging. 1/21/2017.  
        # initialize some class attributes
        self.pccontext = None
        self.preconditioner = None
        self.pc = None
        self.solverName  = "PETSc"
        self.par_fullOverlap = True
        self.par_firstAssembly=True
        self.par_L   = par_L
        self.petsc_L = par_L
        self.csr_rep_local = self.petsc_L.csr_rep_local
        self.csr_rep = self.petsc_L.csr_rep
        self.bdyNullSpace = bdyNullSpace

        # create petsc4py KSP object and attach operators
        self.ksp = p4pyPETSc.KSP().create()
        self._setMatOperators()
        self.ksp.setOperators(self.petsc_L,self.petsc_L)

        # set the ksp residual tolerance, options prefix and function handle for convergence message.
        self.setResTol(rtol_r,atol_r,maxIts)
        # I don't really understand this...
        convergenceTest = 'r-true'
        if convergenceTest == 'r-true':
            self.r_work = self.petsc_L.getVecLeft()
            self.rnorm0 = None
            self.ksp.setConvergenceTest(self._converged_trueRes)
        else:
            self.r_work = None
        ### DO NOT MERGE - THIS NEEDS TO BE FIXED!!!!
        if prefix != None:
           self.ksp.setOptionsPrefix(prefix)
        # set ksp preconditioner
        if Preconditioner != None:
            self._setPreconditioner(Preconditioner,par_L,prefix)
            self.ksp.setPC(self.pc)
        # set the ksp options
        self.ksp.setFromOptions()

    def setResTol(self,rtol,atol,maxIts):
        """ Set the ksp object's residual and maximum iterations. """
        self.rtol_r = rtol
        self.atol_r = atol
        self.maxIts = maxIts
        self.ksp.rtol = rtol
        self.ksp.atol = atol
        self.ksp.max_it = self.maxIts
        logEvent("KSP atol %e rtol %e" % (self.ksp.atol,self.ksp.rtol))

    def prepare(self,b=None):
        self.petsc_L.zeroEntries()
        assert self.petsc_L.getBlockSize() == 1, "petsc4py wrappers currently require 'simple' blockVec (blockSize=1) approach"
        if self.petsc_L.proteus_jacobian != None:
            self.csr_rep[2][self.petsc_L.nzval_proteus2petsc] = self.petsc_L.proteus_csr_rep[2][:]
        if self.par_fullOverlap == True:
            self.petsc_L.setValuesLocalCSR(self.csr_rep_local[0],self.csr_rep_local[1],self.csr_rep_local[2],p4pyPETSc.InsertMode.INSERT_VALUES)
        else:
            if self.par_firstAssembly:
                self.petsc_L.setOption(p4pyPETSc.Mat.Option.NEW_NONZERO_LOCATION_ERR,False)
                self.par_firstAssembly = False
            else:
                self.petsc_L.setOption(p4pyPETSc.Mat.Option.NEW_NONZERO_LOCATION_ERR,True)
            self.petsc_L.setValuesLocalCSR(self.csr_rep[0],self.csr_rep[1],self.csr_rep[2],p4pyPETSc.InsertMode.ADD_VALUES)
        self.petsc_L.assemblyBegin()
        self.petsc_L.assemblyEnd()
        self.ksp.setOperators(self.petsc_L,self.petsc_L)
        #self.ksp.setOperators(self.Lshell,self.petsc_L)
        # ARB - WIP: this self.pc and self.preconditioner stuff is confusing to read
        # and not necessary.  This should be refactored before merging. 1/21/2017.  
        if self.pc != None:
            self.pc.setOperators(self.petsc_L,self.petsc_L)
            self.pc.setUp()
            # TODO - ARB, there should be an assert here in case PCType is not set.
            if self.preconditioner.PCType=='schur':
                # this should probably live somewhere else?!
                self.preconditioner.setUp(self.ksp)
#                import pdb ; pdb.set_trace()
#                self.pc.getFieldSplitSubKSP()[1].setPCSide(1)
        self.ksp.setUp()
        self.ksp.pc.setUp()
#        self.ksp.pc.getFieldSplitSubKSP()[1].setConvergenceTest(self.preconditioner._converged_trueRes)


    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=True):
        if par_b.proteus2petsc_subdomain is not None:
            par_b.proteus_array[:] = par_b.proteus_array[par_b.petsc2proteus_subdomain]
            par_u.proteus_array[:] = par_u.proteus_array[par_u.petsc2proteus_subdomain]
        # if self.petsc_L.isSymmetric(tol=1.0e-12):
        #    self.petsc_L.setOption(p4pyPETSc.Mat.Option.SYMMETRIC, True)
        #    print "Matrix is symmetric"
        # else:
        #    print "MATRIX IS NONSYMMETRIC"
        logEvent("before ksp.rtol= %s ksp.atol= %s ksp.converged= %s ksp.its= %s ksp.norm= %s " % (self.ksp.rtol,
                                                                                                   self.ksp.atol,
                                                                                                   self.ksp.converged,
                                                                                                   self.ksp.its,
                                                                                                   self.ksp.norm))
        if self.pccontext != None:
            self.pccontext.par_b = par_b
            self.pccontext.par_u = par_u
        if self.matcontext != None:
            self.matcontext.par_b = par_b
#        import pdb ; pdb.set_trace()
        # if p4pyPETSc.COMM_WORLD.rank==0:
        #     import pdb ; pdb.set_trace()
        ### THIS IS A DEBUGGING EXPERIMENT
#        self.__createParVec(par_b)
#        petsc_view(self.tmp,'tmp')
        # petsc_view(self.vec_of_ones,'ones')
        if self.bdyNullSpace==True:
            # is there a reason par_b should not be owned by self?
#            self.__setNullSpace(par_b)
            self.__setNullSpace(par_b)
#        print '*********!!!!!!!!!!!!!!!!! par_u' + `par_u.norm()`
#        print '*********!!!!!!!!!!!!!!!!! par_u' + `self.tmp.norm()`
#        print '*********!!!!!!!!!!!!!!!!! par_b = ' + `par_b.norm()`
        self.ksp.solve(par_b,par_u)
        logEvent("after ksp.rtol= %s ksp.atol= %s ksp.converged= %s ksp.its= %s ksp.norm= %s reason = %s" % (self.ksp.rtol,
                                                                                                             self.ksp.atol,
                                                                                                             self.ksp.converged,
                                                                                                             self.ksp.its,
                                                                                                             self.ksp.norm,
                                                                                                             self.ksp.reason))
        self.its = self.ksp.its
        if self.printInfo:
            self.info()
        if par_b.proteus2petsc_subdomain is not None:
            par_b.proteus_array[:] = par_b.proteus_array[par_b.proteus2petsc_subdomain]
            par_u.proteus_array[:] = par_u.proteus_array[par_u.proteus2petsc_subdomain]
    def converged(self,r):
        return self.ksp.converged
    def failed(self):
        failedFlag = LinearSolver.failed(self)
        failedFlag = failedFlag or (not self.ksp.converged)
        return failedFlag

    def info(self):
        self.ksp.view()

    def __createParVec(self,par_b):
        """ A initial attempt to set up the null space vector
        in parallel.

        """
        #Global null space constructor
        # if p4pyPETSc.COMM_WORLD.rank==0:
        #     import pdb ; pdb.set_trace()
        #Comm information
        rank = p4pyPETSc.COMM_WORLD.rank
        size = p4pyPETSc.COMM_WORLD.size
        # information for par_vec
        par_N = par_b.getSize()
        par_n = par_b.getLocalSize()
        par_nghosts = par_b.nghosts
        subdomain2global = par_b.subdomain2global
        proteus2petsc_subdomain = par_b.proteus2petsc_subdomain
        petsc2proteus_subdomain = par_b.petsc2proteus_subdomain
        # information for pressure uknowns
        pSpace = self.par_L.pde.u[0].femSpace
        pressure_offsets = pSpace.dofMap.dof_offsets_subdomain_owned
        n_DOF_pressure = pressure_offsets[rank+1] - pressure_offsets[rank]
        N_DOF_pressure = pressure_offsets[size]
        total_pressure_unknowns = pSpace.dofMap.nDOF_subdomain
        self.tmp_array_1 = numpy.ones((par_n+par_nghosts),'d')
        self.tmp_array_2 = numpy.zeros((par_n+par_nghosts),'d')
        self.tmp_array_3 = numpy.zeros((par_n+par_nghosts),'d')
        self.vec_of_ones = ParVec_petsc4py(self.tmp_array_1,
                                           1,
                                           par_n,
                                           par_N,
                                           par_nghosts,
                                           subdomain2global,
                                           proteus2petsc_subdomain=proteus2petsc_subdomain,
                                           petsc2proteus_subdomain=petsc2proteus_subdomain)
        self.tmp = ParVec_petsc4py(self.tmp_array_2,
                                   1,
                                   par_n,
                                   par_N,
                                   par_nghosts,
                                   subdomain2global,
                                   proteus2petsc_subdomain=proteus2petsc_subdomain,
                                   petsc2proteus_subdomain=petsc2proteus_subdomain)
        # self.tmp2 = ParVec_petsc4py(tmp_array_3,
        #                            1,
        #                            par_n,
        #                            par_N,
        #                            par_nghosts,
        #                            subdomain2global,
        #                            proteus2petsc_subdomain=proteus2petsc_subdomain,
        #                            petsc2proteus_subdomain=petsc2proteus_subdomain)
        self.ksp.getOperators()[0].mult(self.vec_of_ones,self.tmp)
#        print '*************' +`self.ksp.getOperators()[0]`
        
    def __defineNullSpaceVec(self,par_b):
        """ A initial attempt to set up the null space vector
        in parallel.

        """
        #Global null space constructor
#        import pdb ; pdb.set_trace()
        #Comm information
        rank = p4pyPETSc.COMM_WORLD.rank
        size = p4pyPETSc.COMM_WORLD.size
        # information for par_vec
        par_N = par_b.getSize()
        par_n = par_b.getLocalSize()
        par_nghosts = par_b.nghosts
        subdomain2global = par_b.subdomain2global
        proteus2petsc_subdomain = par_b.proteus2petsc_subdomain
        petsc2proteus_subdomain = par_b.petsc2proteus_subdomain
        # information for pressure uknowns
        pSpace = self.par_L.pde.u[0].femSpace
        pressure_offsets = pSpace.dofMap.dof_offsets_subdomain_owned
        n_DOF_pressure = pressure_offsets[rank+1] - pressure_offsets[rank]
#        assert par_n == n_DOF_pressure
        N_DOF_pressure = pressure_offsets[size]
        total_pressure_unknowns = pSpace.dofMap.nDOF_subdomain
        t = pSpace.dofMap.nDOF_all_processes
        assert t==N_DOF_pressure
        self.tmp_array = numpy.zeros((par_n+par_nghosts),'d')
        self.tmp_array[0:n_DOF_pressure] = 1.0/(sqrt(N_DOF_pressure))
#        print 'pressure dofs global = ' +`N_DOF_pressure`
#        self.tmp_array[0:n_DOF_pressure] = 1.0/N_DOF_pressure
        # if p4pyPETSc.COMM_WORLD.rank==0:
        #     import pdb ; pdb.set_trace()
        null_space_basis = ParVec_petsc4py(self.tmp_array,
                                           1,
                                           par_n,
                                           par_N,
                                           par_nghosts,
                                           subdomain2global,
                                           proteus2petsc_subdomain=proteus2petsc_subdomain,
                                           petsc2proteus_subdomain=petsc2proteus_subdomain)
        petsc_view(null_space_basis,'pre_vec_rank_'+`p4pyPETSc.COMM_WORLD.rank`)
        null_space_basis.scatter_forward_insert()
#        print self.tmp_array
        self.tmp_array[:] = self.tmp_array[petsc2proteus_subdomain]
        from proteus import LinearAlgebraTools
#        print 'PROTEUS NORM *********' + `LinearAlgebraTools.l2Norm(self.tmp_array[0:par_n])`
#        print '*******NORM**********' + `null_space_basis.norm(1)`
        self.global_null_space = [null_space_basis]

    def __setNullSpace(self,par_b):
        """ Set up the boundary null space for the KSP solves. 

        Parameters
        ----------
        par_b : proteus.LinearAlgebraTools.ParVec_petsc4py
            The problem's RHS vector.
        """
        self.__defineNullSpaceVec(par_b)
        vecs = self.global_null_space
#        vecs = self.preconditioner.global_null_space
        # from proteus import Comm
#        print 'processor rank = ' +`p4pyPETSc.COMM_WORLD.rank` + ' has left the station'
        self.pressure_null_space = p4pyPETSc.NullSpace().create(constant = False,
                                                           vectors = vecs,
                                                           comm = p4pyPETSc.COMM_WORLD)
        self.ksp.getOperators()[0].setNullSpace(self.pressure_null_space)
#        petsc_view(par_b,'par_b_pre_remove')
#        import pdb ; pdb.set_trace()
        self.pressure_null_space.remove(par_b)
#        petsc_view(par_b,'par_b_post_remove')

    def _setMatOperators(self):
        """ Initializes python context for the ksp matrix operator """
        self.Lshell = p4pyPETSc.Mat().create()
        L_sizes = self.petsc_L.getSizes()
        L_range = self.petsc_L.getOwnershipRange()
        self.Lshell.setSizes(L_sizes)
        self.Lshell.setType('python')
        self.matcontext  = SparseMatShell(self.petsc_L.ghosted_csr_mat)
        self.Lshell.setPythonContext(self.matcontext)

    def _converged_trueRes(self,ksp,its,rnorm):
        """ Function handle to feed to ksp's setConvergenceTest  """
        ksp.buildResidual(self.r_work)
        truenorm = self.r_work.norm()
        if its == 0:
            self.rnorm0 = truenorm
            logEvent("NumericalAnalytics KSPOuterResidual: %12.5e" %(truenorm) )
            logEvent("NumericalAnalytics KSPOuterResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
                                                                                              truenorm,
                                                                                              (truenorm/ self.rnorm0),
                                                                                              ksp.atol,
                                                                                              ksp.rtol))
            return False
        else:
            logEvent("NumericalAnalytics KSPOuterResidual: %12.5e" %(truenorm) )
            logEvent("NumericalAnalytics KSPOuterResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
                                                                                              truenorm,
                                                                                              (truenorm/ self.rnorm0),
                                                                                              ksp.atol,
                                                                                              ksp.rtol))
            if truenorm < self.rnorm0*ksp.rtol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_RTOL
            if truenorm < ksp.atol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_ATOL
        return False

    def _setPreconditioner(self,Preconditioner,par_L,prefix):
        """ Sets the preconditioner type used in the KSP object """
        if Preconditioner != None:
            if Preconditioner == petsc_LU:
                logEvent("NAHeader Precondtioner LU")
                self.preconditioner = petsc_LU(par_L)
                self.pc = self.preconditioner.pc
            elif Preconditioner == petsc_ASM:
                logEvent("NAHead Preconditioner ASM")
                self.preconditioner = petsc_ASM(par_L,prefix)
                self.pc = self.preconditioner.pc
            #ARB - something else needs to be done with these commented out preconditioners,
            #I've not thought of what just yet. 1-19-17.
            # if Preconditioner == Jacobi:
            #     self.pccontext= Preconditioner(L,
            #                                    weight=1.0,
            #                                    rtol_r=rtol_r,
            #                                    atol_r=atol_r,
            #                                    maxIts=1,
            #                                    norm = l2Norm,
            #                                    convergenceTest='its',
            #                                    computeRates=False,
            #                                    printInfo=False)
            #     self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            # elif Preconditioner == GaussSeidel:
            #     self.pccontext= Preconditioner(connectionList,
            #                                    L,
            #                                    weight=1.0,
            #                                    sym=False,
            #                                    rtol_r=rtol_r,
            #                                    atol_r=atol_r,
            #                                    maxIts=1,
            #                                    norm = l2Norm,
            #                                    convergenceTest='its',
            #                                    computeRates=False,
            #                                    printInfo=False)
            #     self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            # elif Preconditioner == LU:
            #     #ARB - parallel matrices from PETSc4py don't work here...
            #     self.pccontext = Preconditioner(L)
            #     self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            # elif Preconditioner == StarILU:
            #     self.pccontext= Preconditioner(connectionList,
            #                                    L,
            #                                    weight=1.0,
            #                                    rtol_r=rtol_r,
            #                                    atol_r=atol_r,
            #                                    maxIts=1,
            #                                    norm = l2Norm,
            #                                    convergenceTest='its',
            #                                    computeRates=False,
            #                                    printInfo=False)
            #     self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            # elif Preconditioner == StarBILU:
            #     self.pccontext= Preconditioner(connectionList,
            #                                    L,
            #                                    bs=linearSolverLocalBlockSize,
            #                                    weight=1.0,
            #                                    rtol_r=rtol_r,
            #                                    atol_r=atol_r,
            #                                    maxIts=1,
            #                                    norm = l2Norm,
            #                                    convergenceTest='its',
            #                                    computeRates=False,
            #                                    printInfo=False)
            #     self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == SimpleNavierStokes3D:
                logEvent("NAHeader Preconditioner selfp" )
                self.preconditioner = SimpleNavierStokes3D(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes3D_Qp:
                logEvent("NAHeader Preconditioner Qp" )
                self.preconditioner = NavierStokes3D_Qp(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes_TwoPhaseQp:
                logEvent("NAHeader Preconditioner TwoPhaseQp" )
                self.preconditioner = NavierStokes_TwoPhaseQp(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes3D_PCD:
                logEvent("NAHeader Preconditioner PCD" )
                self.preconditioner = NavierStokes3D_PCD(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes3D_LSC:
                logEvent("NAHeader Preconditioner LSC" )
                self.preconditioner = NavierStokes3D_LSC(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes_TwoPhaseLSC:
                logEvent("NAHeader Preconditioner TwoPhaseLSC" )
                self.preconditioner = NavierStokes_TwoPhaseLSC(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokes_TwoPhasePCD:
                logEvent("NAHeader Preconditioner TwoPhasePCD" )
                self.preconditioner = NavierStokes_TwoPhasePCD(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == SimpleNavierStokes2D:
                self.preconditioner = SimpleNavierStokes2D(par_L,prefix)
                self.pc = self.preconditioner.pc
            elif Preconditioner == SimpleDarcyFC:
                self.preconditioner = SimpleDarcyFC(par_L)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokesPressureCorrection:
                self.preconditioner = NavierStokesPressureCorrection(par_L)
                self.pc = self.preconditioner.pc


class SchurOperatorConstructor:
    """ Generate matrices for use in Schur complement preconditioner operators. 

    Notes
    -----
    TODO (ARB) - This class should probably have its own test class.
    TODO (ARB) - what should be self. vs. returned?
    TODO (ARB) - Get the matrix allocation function working
    """
    def __init__(self, linear_smoother, pde_type):
        """
        Initialize a Schur Operator constructor.

        Parameters
        ----------
        linear_smoother : class
            Provides the data about the problem.
        pde_type :  str 
            Currently supports Stokes and navierStokes
        """
        if linear_smoother.PCType!='schur':
            raise Exception, 'This function only works with the' \
                'LinearSmoothers for Schur Complements.'
        self.linear_smoother=linear_smoother
        self.L = linear_smoother.L
        self.pde_type = pde_type
        if isinstance(self.L.pde,proteus.mprans.RANS2P.LevelModel):
            self.opBuilder = OperatorConstructor_rans2p(self.L.pde)
        else:
            self.opBuilder = OperatorConstructor_oneLevel(self.L.pde)
        self.Qp_built = False
        self.Ap_built = False

    def getQp(self, output_matrix=False):
        """ Return the pressure mass matrix Qp.

        Parameters
        ----------
        output_matrix : bool 
            Determines whether matrix should be exported.

        Returns
        -------
        Qp : matrix
            The pressure mass matrix.
        """
        if self.Qp_built == False:
            Qsys_petsc4py = self._massMatrix()
            self.Qp = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                                 self.linear_smoother.isp)
            if output_matrix==True:
                self._exportMatrix(self.Qp,"Qp")
            self.Qp_built = True
        return self.Qp

    def getTwoPhaseQp_rho(self,
                          output_matrix = False,
                          phase_function = None):
        """ Return the two-phase pressure mass matrix.

        Parameters
        ----------
        output_matrix : bool
            Flag for weather matrix should be output.
        phase_function : lambda
            Optional function describing the fluid phases.
        inv_scaled : bool
            Determines whether the mass matrix should be 
            scaled by the viscosity or the inverse of the
            viscosity.

        Returns
        -------
        Qp : matrix
           A two-phase pressure mass matrix.
        """
        try:
            if self.L.pde.coefficients.which_region != None:
                self._phase_function = self.L.pde.coefficients.which_region
        except AttributeError:
            pass
        Qsys_petsc4py = self._TPMassMatrix_rho()
        self.TPQp = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                               self.linear_smoother.isp)

        if output_matrix == True:
            self._exportMatrix(self.TPQp)

        return self.TPQp
    
    def getTwoPhaseInvScaledQp(self,
                               output_matrix = False,
                               phase_function = None):
        """ Return the two-phase pressure mass matrix
        scaled with the inverse of the viscosity.

        Parameters
        ----------
        output_matrix : bool
            Flag for weather matrix should be output.
        phase_function : lambda
            Optional function describing the fluid phases.
        inv_scaled : bool
            Determines whether the mass matrix should be 
            scaled by the viscosity or the inverse of the
            viscosity.

        Returns
        -------
        Qp : matrix
           A two-phase pressure mass matrix.
        """
        try:
            if self.L.pde.coefficients.which_region != None:
                self._phase_function = self.L.pde.coefficients.which_region
        except AttributeError:
            pass
        Qsys_petsc4py = self._TPInvScaledMassMatrix()
        self.TPInvScaledQp = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                                        self.linear_smoother.isp)

        if output_matrix == True:
            self._exportMatrix(self.TPInvScaledQp)

        return self.TPInvScaledQp


    def getQv(self,output_matrix=False):
        """ Return the velocity mass matrix Qv.

        Parameters
        ----------
        output_matrix : bool
            Determine whether the matrix should be exported.
        
        Returns
        -------
        Qv : matrix
            The velocity mass matrix.
        """
        Qsys_petsc4py = self._massMatrix()
        self.Qv = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                             self.linear_smoother.isv)
        if output_matrix==True:
            self._exportmatrix(self.Qv,'Qv')

        return self.Qv

    def getTwoPhaseQv_rho(self,output_matrix=False):
        """ Return the velocity mass matrix Qv.

        Parameters
        ----------
        output_matrix : bool
            Determine whether the matrix should be exported.
        
        Returns
        -------
        Qv : matrix
            The velocity mass matrix.
        """
        if self.L.pde.coefficients.which_region != None:
            self._phase_function = self.L.pde.coefficients.which_region

        Qsys_petsc4py = self._TPMassMatrix_rho()
        self.Qv = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                             self.linear_smoother.isv)
        if output_matrix==True:
            self._exportmatrix(self.Qv,'Qv')

        return self.Qv

    def getTwoPhaseQv_mu(self,output_matrix=False):
        """ Return the velocity mass matrix Qv.

        Parameters
        ----------
        output_matrix : bool
            Determine whether the matrix should be exported.
        
        Returns
        -------
        Qv : matrix
            The velocity mass matrix.
        """
        if self.L.pde.coefficients.which_region != None:
            self._phase_function = self.L.pde.coefficients.which_region

        Qsys_petsc4py = self._TPMassMatrix_mu()
        self.Qv = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                             self.linear_smoother.isv)
        if output_matrix==True:
            self._exportmatrix(self.Qv,'Qv')

        return self.Qv

    def getAp(self,output_matrix=False):
        """ Return the Laplacian pressure matrix Ap.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        Ap : matrix
            The Laplacian pressure matrix.
        """
        if self.Ap_built == False:
            Ap_sys_petsc4py = self._getLaplace()
            self.Ap = Ap_sys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                                   self.linear_smoother.isp)
            if output_matrix==True:
                self._exportmatrix(self.Ap,'Ap')
            self.Ap_built = True
        return self.Ap

    def getTwoPhaseInvScaledAp(self,
                               output_matrix = False,
                               phase_function = None):
        """ Return the two-phase pressure laplacian scaled with the
        inverse of the viscosity.

        Parameters
        ----------
        output_matrix : bool
            Flag for weather matrix should be output.
        phase_function : lambda
            Optional function describing the fluid phases.
        inv_scaled : bool
            Determines whether the mass matrix should be 
            scaled by the viscosity or the inverse of the
            viscosity.

        Returns
        -------
        Qp : matrix
           A two-phase pressure mass matrix.
        """
        try:
            if self.L.pde.coefficients.which_region != None:
                self._phase_function = self.L.pde.coefficients.which_region
        except AttributeError:
            pass
        Qsys_petsc4py = self._getTPInvScaledLaplace()
        self.TPInvScaledAp = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                                        self.linear_smoother.isp)

        if output_matrix == True:
            self._exportMatrix(self.TPInvScaledQp)

        return self.TPInvScaledAp


    def getAv(self,output_matrix=False):
        """ Return the Laplacian velocity matrix Av.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        Av : matrix
            The Laplacian velocity matrix.
        """
        Av_sys_petsc4py = self._getLaplace()
        self.Av = Av_sys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                               self.linear_smoother.isv)
        if output_matrix==True:
            self._exportMatrix(self.Av,'Av')
        return self.Av
        
    def getFp(self,output_matrix=False):
        """ Return the convection-diffusion matrix for the pressure """
        self.Fp = self.getCp()
        Ap = self.getAp()
        
        # TODO - ARB.  In the case of two-phase coefficient classes
        # with a constant viscosity, nu may not be updated to reflect
        # the viscosities nu_0 and nu_1.  I'm not sure how to fix this here.
        # In the future when we have two-phase versions of the preconditioner,
        # the two phase problems should not be allowed to use the single
        # phase versions of the PC.

        try:
            viscosity = self.opBuilder.OLT.coefficients.nu
        except ValueError:
            print "Invalid viscosity in coefficients class."
            raise
        
        self.Fp.axpy(viscosity,Ap)
        return self.Fp

    def getCp(self,output_matrix=True):
        """
        Return the convection matrix for the pressure Fp.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        Fp : matrix
             The convection-diffusion pressure matrix.

        Notes
        ----
        While I've tested this function, this may still contain a bug.
        Needs to be refactored.
        """
        #modify the diffusion term in the mass equation so the p-p block is Fp
        # First generate the advection part of Fp
        self._u = numpy.copy(self.L.pde.q[('u',1)])
        self._v = numpy.copy(self.L.pde.q[('u',2)])
        self._advectiveField = [self._u,self._v]
        Cp_sys_petsc4py = self._getAdvection()
        self.Cp = Cp_sys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                               self.linear_smoother.isp)
        if output_matrix==True:
            self._exportMatrix(self.Cp,'Cp')

        return self.Cp

    def getTwoPhaseCp_rho(self,output_matrix=True):
        """
        Return a two-phase convection matrix for the pressure.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        Cp : matrix
             The convection-diffusion pressure matrix.

        Notes
        ----
        While I've tested this function, this may still contain a bug.
        Needs to be refactored.
        """
        #modify the diffusion term in the mass equation so the p-p block is Fp
        # First generate the advection part of Fp
        self._u = numpy.copy(self.L.pde.q[('u',1)])
        self._v = numpy.copy(self.L.pde.q[('u',2)])
        self._advectiveField = [self._u,self._v]
        Cp_sys_petsc4py = self._getTwoPhaseAdvection()
        self.TPCp = Cp_sys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                                 self.linear_smoother.isp)
        if output_matrix==True:
            self._exportMatrix(self.TPCp,'Cp')

        return self.TPCp

    
    def getB(self,output_matrix=False):
        """ Create the B operator.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be stored externally.

        Returns
        -------
        B : matrix
            The operator B matrix.
        """
        Bsys_petsc4py = self._getB()
        self.B = Bsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                            self.linear_smoother.isv)
        if output_matrix==True:
            self._exportMatrix(self.B,'B')
        return self.B

    def getBt(self,output_matrix=False):
        """ Create the Bt operator.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be stored externally.

        Returns
        -------
        B : matrix
            The operator B matrix.
        """
        Bsys_petsc4py = self._getB()
        self.Bt = Bsys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                             self.linear_smoother.isp)
        if output_matrix==True:
            self._exportMatrix(self.Bt,'Bt')
        return self.Bt

    def getF(self,output_matrix=False):
        """ Return the A-block of a saddle-point system """
        rowptr, colind, nzval = self.L.pde.jacobian.getCSRrepresentation()
        self.F_rowptr = rowptr.copy()
        self.F_colind = colind.copy()
        self.F_nzval = nzval.copy()
        L_sizes = self.L.getSizes()
        nr = L_sizes[0][0]
        nc = L_sizes[1][0]
        self.F = SparseMat(nr,nc,
                           self.F_nzval.shape[0],
                           self.F_nzval,
                           self.F_colind,
                           self.F_rowptr)
        self.Fsys_petsc4py = self.L.duplicate()
        F_csr_rep_local = self.F.getSubMatCSRrepresentation(0,L_sizes[0][0])
        self.Fsys_petsc4py.setValuesLocalCSR(F_csr_rep_local[0],
                                             F_csr_rep_local[1],
                                             F_csr_rep_local[2],
                                             p4pyPETSc.InsertMode.INSERT_VALUES)
        self.Fsys_petsc4py.assemblyBegin()
        self.Fsys_petsc4py.assemblyEnd()
        self.F = self.Fsys_petsc4py.getSubMatrix(self.linear_smoother.isv,
                                                 self.linear_smoother.isv)
        if output_matrix==True:
            self._exportMatrix(self.F,'F')
        return self.F

    def getTildeA(self,output_matrix = False):
        """ Only take the symmetric parts of the A """
        self.Av = self.getAv()
        self.Av.scale(self.L.pde.coefficients.nu)
        #TODO (ARB) - needs mass matrix term for time dependent problem
        return self.Av

    def _massMatrix(self):
        """ Generates a the mass matrix.

        This function generates and returns the mass matrix for the system. This
        function is internal to the class and called by public functions which 
        take and return the relevant components (eg. the pressure or velcoity).

        Returns
        -------
        Qsys : matrix
            The system's mass matrix.
        """
        self.opBuilder.attachMassOperator()
        return superlu_2_petsc4py(self.opBuilder.MassOperator)

    def _TPMassMatrix_rho(self):
        """ Generates a two phase mass matrix.

        Returns
        -------
        Qsys : matrix
            The two-phase mass matrix.

        """
        try :
            self.opBuilder.attachTwoPhaseMassOperator_rho(phase_function = self._phase_function)
        except:
            self.opBuilder.attachTwoPhaseMassOperator_rho()
        return superlu_2_petsc4py(self.opBuilder.TPMassOperator)

    def _TPMassMatrix_mu(self):
        """ Generates a two phase mass matrix scaled by the dynamic viscosity.

        Returns
        -------
        Qsys : matrix
            The two-phase mass matrix.

        """
        self.opBuilder.attachTwoPhaseMassOperator_mu(phase_function = self._phase_function)
        return superlu_2_petsc4py(self.opBuilder.TPMassOperator)


    def _TPInvScaledMassMatrix(self):
        """ Generates a two phase mass matrix scaled by the inverse of the viscosity.

        Returns
        -------
        Qsys : matrix
            The two-phase mass matrix.

        """
        try:
            self.opBuilder.attachTwoPhaseInvScaledMassOperator(phase_function = self._phase_function)
        except AttributeError:
            self.opBuilder.attachTwoPhaseInvScaledMassOperator()
        return superlu_2_petsc4py(self.opBuilder.TPInvScaledMassOperator)
    
    def _getLaplace(self,output_matrix=False):
        """ Return the Laplacian pressure matrix Ap.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        A : matrix
            The Laplacian  matrix.

        TODO (ARB)  Write supporting tests.
        """
        self.opBuilder.attachLaplaceOperator()
        return superlu_2_petsc4py(self.opBuilder.LaplaceOperator)

    def _getTPInvScaledLaplace(self,output_matrix=False):
        """ Return the Laplacian pressure matrix Ap.

        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        A : matrix
            The Laplacian  matrix.

        """
        try:
            self.opBuilder.attachTPInvScaledLaplaceOperator(phase_function = self._phase_function)
        except AttributeError:
            self.opBuilder.attachTPInvScaledLaplaceOperator()
        return superlu_2_petsc4py(self.opBuilder.TPInvScaledLaplaceOperator)
    
    def _getAdvection(self,output_matrix=False):
        """ Return the advection matrix Fp.

        Returns
        -------
        Fp : petsc4py matrix
            A petsc4py matrix.
        """
        self.opBuilder.attachAdvectionOperator(self._advectiveField)
        return superlu_2_petsc4py(self.opBuilder.AdvectionOperator)

    def _getTwoPhaseAdvection(self,output_matrix=False):
        """ Returns the two-phase advection matrix TPCp.

        Returns
        -------
        TPCp : petsc4py matrix
        """
        try:
            self.opBuilder.attachTPAdvectionOperator(self._advectiveField, phase_function = self._phase_function)
        except AttributeError:
            self.opBuilder.attachTPAdvectionOperator()
        return superlu_2_petsc4py(self.opBuilder.TPScaledAdvectionOperator)

    def _getB(self,output_matrix=False):
        """ Return the discrete B-operator.
        
        Parameters
        ----------
        output_matrix : bool
            Determine whether matrix should be exported.

        Returns
        -------
        A : matrix
            The B operator matrix.
        """
        self.opBuilder.attachBOperator()
        return superlu_2_petsc4py(self.opBuilder.BOperator)

    def _initializeMatrix(self):
        """ Allocates memory for the matrix operators.

        Returns
        -------
        Q : sparseMat
        """
        rowptr, colind, nzval = self.L.pde.jacobian.getCSRrepresentation()
        Q_rowptr = rowptr.copy()
        Q_colind = colind.copy()
        Q_nzval  = nzval.copy()
        nr = rowptr.shape[0] - 1
        nc = nr
        Q = SparseMat(nr,nc,
                      Q_nzval.shape[0],
                      Q_nzval,
                      Q_colind,
                      Q_rowptr)        
        return Q


    def _exportMatrix(self,operator,export_name):
        """ Export the matrix operator.

        Parameters
        ----------
        operator : matrix
                   Operator to be exported.
        export_name : str
                    Export file name.
        """
        from LinearAlgebraTools import petsc_view
        petsc_view(operator, export_name) #write to export_name.m
        #running export_name.m will load into matlab  sparse matrix
        #runnig operatorf = full(Mat_...) will get the full matrix

class KSP_Preconditioner:
    """ Base class for PETSCc KSP precondtioners. """
    def __init__(self):
        pass

    def setup(self):
        pass

class petsc_ASM(KSP_Preconditioner):
    """ASM PETSc preconditioner class.

    This class provides an ASM preconditioners for PETSc4py KSP
    objects.
    """
    def __init__(self,L,prefix=None):
        """
        Initializes the ASMpreconditioner for use with PETSc.

        Parameters
        ----------
        L : the global system matrix.
        """
        self.PCType = 'asm'
        self.L = L
        self._create_preconditioner()
        self.pc.setFromOptions()

    def _create_preconditioner(self):
        """ Create the pc object. """
        self.pc = p4pyPETSc.PC().create()
        self.pc.setType('asm')

    def setUp(self,global_ksp=None):
        pass
        
class petsc_LU(KSP_Preconditioner):
    """ LU PETSc preconditioner class.
    
    This class provides an LU preconditioner for PETSc4py KSP
    objects.  Provided the LU decomposition is successful, the KSP
    iterative will converge in a single step.
    """
    def __init__(self,L,prefix=None):
        """
        Initializes the LU preconditioner for use with PETSc.

        Parameters
        ----------
        L : the global system matrix.
        """
        self.PCType = 'lu'
        self.L = L
        self._create_preconditioner()
        self.pc.setFromOptions()

    def _create_preconditioner(self):
        """ Create the ksp preconditioner object.  """
        self.pc = p4pyPETSc.PC().create()
        self.pc.setType('lu')

    def setUp(self,global_ksp=None):
        pass

    
class SchurPrecon(KSP_Preconditioner):
    """ Base class for PETSc Schur complement preconditioners.

    Notes
    -----
    TODO - needs to run for TH, Q1Q1, dim = 2 or 3 etc.
    """
    def __init__(self,L,prefix=None):
        """
        Initializes the Schur complement preconditioner for use with PETSc.

        This class creates a KSP PETSc solver object and initializes flags the
        pressure and velocity unknowns for a general saddle point problem.

        Parameters
        ----------
        L : provides the definition of the problem.
        """
        self.PCType = 'schur'
        self.L = L
        self._initializeIS(prefix)
        self.pc.setFromOptions()

    def setUp(self,global_ksp):
        pass
        
    def _initializeIS(self,prefix):
        """ Sets the index set (IP) for the pressure and velocity 
        
        Parameters
        ----------
        prefix : ???

        Notes
        -----
        TODO - this needs to set up to run for TH,Q1Q1, dim = 2 or 3 etc.
        """
        L_sizes = self.L.getSizes()
        L_range = self.L.getOwnershipRange()
        neqns = L_sizes[0][0]
        nc = self.L.pde.nc
        rank = p4pyPETSc.COMM_WORLD.rank
        size = p4pyPETSc.COMM_WORLD.size
        pSpace = self.L.pde.u[0].femSpace
        pressure_offsets = pSpace.dofMap.dof_offsets_subdomain_owned
        n_DOF_pressure = pressure_offsets[rank+1] - pressure_offsets[rank]
        N_DOF_pressure = pressure_offsets[size]
        if self.L.pde.stride[0] == 1:#assume end to end
            self.pressureDOF = numpy.arange(start=L_range[0],
                                            stop=L_range[0]+n_DOF_pressure,
                                            dtype="i")
            self.velocityDOF = numpy.arange(start=L_range[0]+n_DOF_pressure,
                                            stop=L_range[0]+neqns,
                                            step=1,
                                            dtype="i")
        else: #assume blocked
            self.pressureDOF = numpy.arange(start=L_range[0],
                                            stop=L_range[0]+neqns,
                                            step=nc,
                                            dtype="i")
            velocityDOF = []
            for start in range(1,nc):
                velocityDOF.append(numpy.arange(start=L_range[0]+start,
                                                stop=L_range[0]+neqns,
                                                step=nc,
                                                dtype="i"))
            self.velocityDOF = numpy.vstack(velocityDOF).transpose().flatten()
        self.pc = p4pyPETSc.PC().create()
        self.pc.setOptionsPrefix(prefix)
        # DO NOT MERGE - THIS NEEDS TO BE FIXED!!!!
        self.pc.setType('fieldsplit')
        self.isp = p4pyPETSc.IS()
        self.isp.createGeneral(self.pressureDOF,comm=p4pyPETSc.COMM_WORLD)
        self.isv = p4pyPETSc.IS()
        self.isv.createGeneral(self.velocityDOF,comm=p4pyPETSc.COMM_WORLD)
        self.pc.setFieldSplitIS(('velocity',self.isv),('pressure',self.isp))
        #Global null space constructor
        # if p4pyPETSc.COMM_WORLD.rank==0:
        # temp_array = numpy.zeros(shape=(neqns,1))
        # temp_array[0:n_DOF_pressure] = 1.0/(sqrt(N_DOF_pressure))
        # null_space_basis = p4pyPETSc.Vec().createWithArray(temp_array,
        #                                                    comm=p4pyPETSc.COMM_WORLD)
        # self.global_null_space = [null_space_basis]
        
class NavierStokesSchur(SchurPrecon):
    """ Schur complement preconditioners for Navier-Stokes problems.

    This class is derived from SchurPrecond and serves as the base
    class for all NavierStokes preconditioners which use the Schur complement
    method.    
    """
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        SchurPrecon.__init__(self,L,prefix)
        self.operator_constructor = SchurOperatorConstructor(self ,'navier_stokes')
        self.bdyNullSpace = bdyNullSpace
        # ARB Question - I think boundary null space should be moved into SchurPrecon

    def setUp(self,global_ksp=None):
        """
        Set up the NavierStokesSchur preconditioner.  

        Nothing needs to be done here for a generic NSE preconditioner. 
        Preconditioner arguments can be set with PETSc command line.
        """
        if self.bdyNullSpace == True:
            self._setConstantPressureNullSpace(global_ksp)
        self._setSchurlog(global_ksp)

    def _setSchurlog(self,global_ksp):
        """ Helper function that attaches a residual log to the inner solve """
        global_ksp.pc.getFieldSplitSubKSP()[1].setConvergenceTest(self._converged_trueRes)

    def _converged_trueRes(self,ksp,its,rnorm):
        """ Function handle to feed to ksp's setConvergenceTest  """
        r_work = ksp.getOperators()[1].getVecLeft()
        ksp.buildResidual(r_work)
        truenorm = r_work.norm()
        if its == 0:
            self.rnorm0 = truenorm
            logEvent("NumericalAnalytics KSPSchurResidual: %12.5e" %(truenorm) )
            logEvent("NumericalAnalytics KSPSchurResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
                                                                                              truenorm,
                                                                                              (truenorm/ self.rnorm0),
                                                                                              ksp.atol,
                                                                                              ksp.rtol))
            return False
        else:
            logEvent("NumericalAnalytics KSPSchurResidual: %12.5e" %(truenorm) )
            logEvent("NumericalAnalytics KSPSchurResidual(relative): %12.5e" %(truenorm / self.rnorm0) )
            logEvent("        KSP it %i norm(r) = %e  norm(r)/|b| = %e ; atol=%e rtol=%e " % (its,
                                                                                              truenorm,
                                                                                              (truenorm/ self.rnorm0),
                                                                                              ksp.atol,
                                                                                              ksp.rtol))
            if truenorm < self.rnorm0*ksp.rtol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_RTOL
            if truenorm < ksp.atol:
                return p4pyPETSc.KSP.ConvergedReason.CONVERGED_ATOL
        return False

    def _setConstantPressureNullSpace(self,global_ksp):
        nsp = global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[0].getNullSpace()
        global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[0].setNullSpace(nsp)

SimpleNavierStokes3D = NavierStokesSchur
NavierStokes3D = NavierStokesSchur

class NavierStokes3D_Qp(NavierStokesSchur) :
    """ A Navier-Stokes preconditioner which uses the pressure mass matrix. """
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        """
        Initializes the pressure mass matrix class.

        Parameters
        ---------
        L - petsc4py matrix
            Defines the problem's operator.
        """
        NavierStokesSchur.__init__(self,L,prefix,bdyNullSpace)

    def setUp(self,global_ksp,S_hat=False):
        """ Attaches the pressure mass matrix to PETSc KSP preconditioner.

        Parameters
        ----------
        global_ksp : PETSc KSP object
        S_hat : bool
            Flag to indicate whether the Schur operator should
            be approximated using the PETSc default or the pressure mass
            matrix. Generally this should be set to False.
        """
        # Create the pressure mass matrix and scaxle by the viscosity.
        self.Qp = self.operator_constructor.getQp()
        self.Qp.scale(1./self.L.pde.coefficients.nu)
        L_sizes = self.Qp.size
        L_range = self.Qp.owner_range

        self.S_hat = 'selfp'
        if S_hat == True:
            # Set up a PETSc shell for the Qp operator.
            self.S_hat = 'Qp'
            self.Qp_shell = p4pyPETSc.Mat().create()
            self.Qp_shell.setSizes(L_sizes)
            self.Qp_shell.setType('python')
            self.matcontext = MatrixShell(self.Qp)
            self.Qp_shell.setPythonContext(self.matcontext)
            self.Qp_shell.setUp()

        # Setup a PETSc shell for the inverse Qp operator
        self.QpInv_shell = p4pyPETSc.Mat().create()
        self.QpInv_shell.setSizes(L_sizes)
        self.QpInv_shell.setType('python')
        self.matcontext_inv = MatrixInvShell(self.Qp)
        self.QpInv_shell.setPythonContext(self.matcontext_inv)
        self.QpInv_shell.setUp()

        # Set PETSc Schur operator
        if self.S_hat == 'selfp':
            global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
            global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
            global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()
        elif self.S_hat == 'Qp':
            raise Exception, 'Currently using Qp as an approximation is not' \
                'supported.'
        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
#            import pdb ; pdb.set_trace() 
            self._setConstantPressureNullSpace(global_ksp)

class NavierStokes_TwoPhaseQp(NavierStokesSchur):
    """A two-phase pressure mass matrix preconditioner """
    def __init__(self, L, prefix=None, bdyNullSpace=False):
        """
        Initialize the two-phase pressure mass matrix preconditioner.

        Parameters
        ----------
        L - petsc4py matrix
        """
        NavierStokesSchur.__init__(self, L, prefix, bdyNullSpace)

    def setUp(self, global_ksp, S_hat = False):
        """  Attach the two-phase pressure mass matrix to preconditioner.
        
        Parameters
        ----------
        global_ksp : PETSc KSP object
        S_hat : bool
        """
        self.twoPhaseQp = self.operator_constructor.getTwoPhaseInvScaledQp()
        L_sizes = self.twoPhaseQp.size
        L_range = self.twoPhaseQp.owner_range

        self.twoPhaseQpInv = p4pyPETSc.Mat().create()
        self.twoPhaseQpInv.setSizes(L_sizes)
        self.twoPhaseQpInv.setType('python')
        self.matcontext_inv = MatrixInvShell(self.twoPhaseQp)
        self.twoPhaseQpInv.setPythonContext(self.matcontext_inv)
        self.twoPhaseQpInv.setUp()

        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()

        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            self._setConstantPressureNullSpace(global_ksp)
        
class NavierStokes3D_PCD(NavierStokesSchur) :
    def __init__(self,L,prefix=None,bdyNullSpace=False,bdyAdjust=True):
        """
        Initialize the pressure convection diffusion preconditioning class.

        Parameters
        ----------
        L :  
        prefix : 
        bdyAdjust :

        Notes
        -----
        This method runs but remains a work in progress.  Notably the
        convection diffusion operator and boundary conditions need to
        be tested and taylored to problem specific boundary conditions.
        """
        NavierStokes3D.__init__(self,L,prefix,bdyNullSpace)
        self.bdyAdjust = bdyAdjust

    def setUp(self,global_ksp):
        # Step-1: build the necessary operators
        self.Qp = self.operator_constructor.getQp()
        self.Fp = self.operator_constructor.getFp()
        self.Ap = self.operator_constructor.getAp()

        comm = Comm.get()
        if self.L.pde.forceStrongConditions == True:
            dirichlet_dof_lst = self.L.pde.dirichletConditionsForceDOF[0].DOFBoundaryConditionsDict.keys()
            if comm.size() > 1:
                self.Ap.zeroRowsColumns(ParInfo_petsc4py.subdomain2global[dirichlet_dof_lst])
                self.Fp.zeroRowsColumns(ParInfo_petsc4py.subdomain2global[dirichlet_dof_lst])
            else:
                self.Ap.zeroRowsColumns(dirichlet_dof_lst)
                self.Fp.zeroRowsColumns(dirichlet_dof_lst)
        else:
            self.Ap.zeroRowsColumns(0)
            self.Fp.zeroRowsColumns(0)

        # This can be helpful for debugging.
#        ParInfo_petsc4py.print_info()
#        petsc_view(self.tmp1,'tmp1_2')
#        petsc_view(self.Qp,'Qp_2')
#        petsc_view(self.Fp,'Fp_2')
#        petsc_view(self.Ap,'Ap_2')

        # Step-2: Set up the Shell for the  PETSc operator
        # Qp
        L_sizes = self.Qp.size
        # ??? Is L_range necessary ???
        L_range = self.Qp.owner_range
        self.PCDInv_shell = p4pyPETSc.Mat().create()
        self.PCDInv_shell.setSizes(L_sizes)
        self.PCDInv_shell.setType('python')
        
        self.matcontext_inv = PCDInv_shell(self.Qp,
                                           self.Fp,
                                           self.Ap,
                                           self.bdyAdjust)
        
        self.PCDInv_shell.setPythonContext(self.matcontext_inv)
        self.PCDInv_shell.setUp()
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()
        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            nsp = p4pyPETSc.NullSpace().create(comm=p4pyPETSc.COMM_WORLD,
                                               vectors = (),
                                               constant = True)
            global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[0].setNullSpace(nsp)
            global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[1].setNullSpace(nsp)
#            self._setConstantPressureNullSpace(global_ksp)

class NavierStokes_TwoPhasePCD(NavierStokesSchur) :
    def __init__(self, L, prefix=None, bdyNullSpace=False):
        """
        Initialize the two-phase PCD preconditioning class.

        Parameters
        ----------
        L : 
        prefix :

        """
        NavierStokes3D.__init__(self, L, prefix, bdyNullSpace)

    def setUp(self, global_ksp):
        # Step-1: build the necessary operators
        self.Np_rho = self.operator_constructor.getTwoPhaseCp_rho()
        self.Ap_invScaledRho = self.operator_constructor.getTwoPhaseInvScaledAp()
        self.Qp_rho = self.operator_constructor.getTwoPhaseQp_rho()
        self.Qp_invScaledVis = self.operator_constructor.getTwoPhaseInvScaledQp()
        L_sizes = self.Qp_rho.size
        
        L_range = self.Qp_rho.owner_range
        self.TP_PCDInv_shell = p4pyPETSc.Mat().create()
        self.TP_PCDInv_shell.setSizes(L_sizes)
        self.TP_PCDInv_shell.setType('python')
        dt = self.L.pde.timeIntegration.t - self.L.pde.timeIntegration.tLast
        self.matcontext_inv = TwoPhase_PCDInv_shell(self.Qp_invScaledVis,
                                                    self.Qp_rho,
                                                    self.Ap_invScaledRho,
                                                    self.Np_rho,
                                                    True,
                                                    dt)
        self.TP_PCDInv_shell.setPythonContext(self.matcontext_inv)
        self.TP_PCDInv_shell.setUp()
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()
        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            nsp = p4pyPETSc.NullSpace().create(comm=p4pyPETSc.COMM_WORLD,
                                               vectors = (),
                                               constant = True)
            global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[0].setNullSpace(nsp)
            global_ksp.pc.getFieldSplitSubKSP()[1].getOperators()[1].setNullSpace(nsp)
#            self._setConstantPressureNullSpace(global_ksp)



class NavierStokes3D_LSC(NavierStokesSchur) :
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        """ Initialize the least squares commutator preconditioning class.

        Parameters
        ----------
        L :  
        prefix : 
        """
        NavierStokesSchur.__init__(self,L,prefix,bdyNullSpace)

    def setUp(self,global_ksp):
        # initialize the Qv_diagonal operator
        self.Qv = self.operator_constructor.getQv()
        self.Qv_hat = p4pyPETSc.Mat().create()
        self.Qv_hat.setSizes(self.Qv.getSizes())
        self.Qv_hat.setType('aij')
        self.Qv_hat.setUp()
        self.Qv_hat.setDiagonal(self.Qv.getDiagonal())
        # initialize the B and F operators
        # ARB Note - There are two ways this can be done...first pull the
        # F and B operators from the global system.  I think this is most
        # direct and cleanest, but one could also call the operator constructors.
        self.B = global_ksp.getOperators()[0].getSubMatrix(self.isp,self.isv)
        self.F = global_ksp.getOperators()[0].getSubMatrix(self.isv,self.isv)
        self.Bt = global_ksp.getOperators()[0].getSubMatrix(self.isv,self.isp)
#        self.B = self.operator_constructor.getB()
#        self.F = self.operator_constructor.getF()
        self.matcontext_inv = LSCInv_shell(self.Qv_hat,self.B,self.Bt,self.F)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()
        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            nsp = p4pyPETSc.NullSpace().create(comm=p4pyPETSc.COMM_WORLD,
                                               vectors = (),
                                               constant = True)
            self._setConstantPressureNullSpace(global_ksp)


class NavierStokes_TwoPhaseLSC(NavierStokesSchur):
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        """ Initialize the two-phase least squares commutator preconditioning
        class.

        Parameters
        ----------
        L : 
        prefix :
        """
        NavierStokesSchur.__init__(self,L,prefix,bdyNullSpace)

    def setUp(self,global_ksp):
        self.two_phase_Qv = self.operator_constructor.getTwoPhaseQv_mu()
        self.Qv_hat = p4pyPETSc.Mat().create()
        self.Qv_hat.setSizes(self.two_phase_Qv.getSizes())
        self.Qv_hat.setType('mpiaij')
        self.Qv_hat.setUp()
        self.Qv_hat.setDiagonal(self.two_phase_Qv.getDiagonal())
        
        self.B = global_ksp.getOperators()[0].getSubMatrix(self.isp,self.isv)
        self.Bt = global_ksp.getOperators()[0].getSubMatrix(self.isv,self.isp)
        self.F = global_ksp.getOperators()[0].getSubMatrix(self.isv,self.isv)

        self.matcontext_inv = LSCInv_shell(self.Qv_hat,
                                           self.B,
                                           self.Bt,
                                           self.F)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()
        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            nsp = p4pyPETSc.NullSpace().create(comm=p4pyPETSc.COMM_WORLD,
                                               vectors = (),
                                               constant = True)
            self._setConstantPressureNullSpace(global_ksp)

class SimpleDarcyFC:
    def __init__(self,L):
        L_sizes = L.getSizes()
        L_range = L.getOwnershipRange()
        print "L_sizes",L_sizes
        neqns = L_sizes[0][0]
        print "neqns",neqns
        self.saturationDOF = numpy.arange(L_range[0],L_range[0]+neqns/2,dtype="i")
        #print "saturation",self.saturationDOF
        self.pressureDOF = numpy.arange(L_range[0]+neqns/2,L_range[0]+neqns,dtype="i")
        #print "pressure",self.pressureDOF
        self.pc = p4pyPETSc.PC().create()
        self.pc.setType('fieldsplit')
        self.isp = p4pyPETSc.IS()
        self.isp.createGeneral(self.saturationDOF,comm=p4pyPETSc.COMM_WORLD)
        self.isv = p4pyPETSc.IS()
        self.isv.createGeneral(self.pressureDOF,comm=p4pyPETSc.COMM_WORLD)
        self.pc.setFieldSplitIS(self.isp)
        self.pc.setFieldSplitIS(self.isv)
    def setUp(self):
        pass

class NavierStokes2D:
    def __init__(self,L,prefix=None):
        self.L=L
        L_sizes = L.getSizes()
        L_range = L.getOwnershipRange()
        neqns = L_sizes[0][0]
        rank = p4pyPETSc.COMM_WORLD.rank
        if self.L.pde.stride[0] == 1:#assume end to end
            pSpace = self.L.pde.u[0].femSpace
            pressure_offsets = pSpace.dofMap.dof_offsets_subdomain_owned
            nDOF_pressure = pressure_offsets[rank+1] - pressure_offsets[rank]
            self.pressureDOF = numpy.arange(start=L_range[0],
                                            stop=L_range[0]+nDOF_pressure,
                                            dtype="i")
            self.velocityDOF = numpy.arange(start=L_range[0]+nDOF_pressure,
                                            stop=L_range[0]+neqns,
                                            step=1,
                                            dtype="i")
        else: #assume blocked
            self.pressureDOF = numpy.arange(start=L_range[0],
                                            stop=L_range[0]+neqns,
                                            step=3,
                                            dtype="i")
            velocityDOF = []
            for start in range(1,3):
                velocityDOF.append(numpy.arange(start=L_range[0]+start,
                                                stop=L_range[0]+neqns,
                                                step=3,
                                                dtype="i"))
                self.velocityDOF = numpy.vstack(velocityDOF).transpose().flatten()
        self.pc = p4pyPETSc.PC().create()
        if prefix:
            self.pc.setOptionsPrefix(prefix)
        self.pc.setType('fieldsplit')
        self.isp = p4pyPETSc.IS()
        self.isp.createGeneral(self.pressureDOF,comm=p4pyPETSc.COMM_WORLD)
        self.isv = p4pyPETSc.IS()
        self.isv.createGeneral(self.velocityDOF,comm=p4pyPETSc.COMM_WORLD)
        self.pc.setFieldSplitIS(('velocity',self.isv),('pressure',self.isp))
        self.pc.setFromOptions()
    def setUp(self):
        try:
            if self.L.pde.pp_hasConstantNullSpace:
                if self.pc.getType() == 'fieldsplit':#we can't guarantee that PETSc options haven't changed the type
                    self.nsp = p4pyPETSc.NullSpace().create(constant=True,comm=p4pyPETSc.COMM_WORLD)
                    self.kspList = self.pc.getFieldSplitSubKSP()
                    self.kspList[1].setNullSpace(self.nsp)
        except:
            pass
SimpleNavierStokes2D = NavierStokes2D

class NavierStokesPressureCorrection:
    def __init__(self,L,prefix=None):
        self.L=L
        self.pc = p4pyPETSc.PC().create()
        if prefix:
            self.pc.setOptionsPrefix(prefix)
        self.pc.setFromOptions()
        self.hasNullSpace=True
        self.nsp = p4pyPETSc.NullSpace().create(constant=True,
                                                comm=p4pyPETSc.COMM_WORLD)
        self.L.setOption(p4pyPETSc.Mat.Option.SYMMETRIC, True)
        self.L.setNullSpace(self.nsp)
    def setUp(self):
        pass

class SimpleDarcyFC:
    def __init__(self,L):
        L_sizes = L.getSizes()
        L_range = L.getOwnershipRange()
        neqns = L_sizes[0][0]
        self.saturationDOF = numpy.arange(L_range[0],L_range[0]+neqns/2,dtype="i")
        self.pressureDOF = numpy.arange(L_range[0]+neqns/2,L_range[0]+neqns,dtype="i")
        self.pc = p4pyPETSc.PC().create()
        self.pc.setType('fieldsplit')
        self.isp = p4pyPETSc.IS()
        self.isp.createGeneral(self.saturationDOF,comm=p4pyPETSc.COMM_WORLD)
        self.isv = p4pyPETSc.IS()
        self.isv.createGeneral(self.pressureDOF,comm=p4pyPETSc.COMM_WORLD)
        self.pc.setFieldSplitIS(self.isp)
        self.pc.setFieldSplitIS(self.isv)
    def setUp(self):
        pass

class Jacobi(LinearSolver):
    """
    Damped Jacobi iteration.
    """
    import csmoothers
    def __init__(self,
                 L,
                 weight=1.0,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        LinearSolver.__init__(self,L,
                              rtol_r,
                              atol_r,
                              rtol_du,
                              atol_du,
                              maxIts,
                              norm,
                              convergenceTest,
                              computeRates,
                              printInfo)
        self.solverName = "Jacobi"
        self.M=Vec(self.n)
        self.w=weight
        self.node_order=numpy.arange(self.n,dtype="i")
    def prepare(self,b=None):
        if type(self.L).__name__ == 'ndarray':
            self.M = self.w/numpy.diagonal(self.L)
        elif type(self.L).__name__ == 'SparseMatrix':
            self.csmoothers.jacobi_NR_prepare(self.L,self.w,1.0e-16,self.M)
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        while (not self.converged(r) and
               not self.failed()):
            if type(self.L).__name__ == 'ndarray':
                self.du[:]=r
                self.du*=self.M
            elif type(self.L).__name__ == "SparseMatrix":
                self.csmoothers.jacobi_NR_solve(self.L,self.M,r,self.node_order,self.du)
            u -= self.du
            self.computeResidual(u,r,b)

class GaussSeidel(LinearSolver):
    """
    Damped Gauss-Seidel.
    """
    import csmoothers
    def __init__(self,
                 connectionList,
                 L,
                 weight=0.33,
                 sym=False,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        LinearSolver.__init__(self,L,
                              rtol_r,
                              atol_r,
                              rtol_du,
                              atol_du,
                              maxIts,
                              norm,
                              convergenceTest,
                              computeRates,
                              printInfo)
        self.solverName = "Gauss-Seidel"
        self.connectionList=connectionList
        self.M=Vec(self.n)
        self.node_order=numpy.arange(self.n,dtype="i")
        self.w=weight
        self.sym=sym
    def prepare(self,b=None):
        if type(self.L).__name__ == 'ndarray':
            self.M = self.w/numpy.diagonal(self.L)
        elif type(self.L).__name__ == 'SparseMatrix':
            self.csmoothers.gauss_seidel_NR_prepare(self.L,self.w,1.0e-16,self.M)
            #self.csmoothers.jacobi_NR_prepare(self.L,self.w,1.0e-16,self.M)
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        while (not self.converged(r) and
               not self.failed()):
            if type(self.L).__name__ == 'ndarray':
                self.du[:]=0.0
                for i in range(self.n):
                    rhat = r[i]
                    for j in self.connectionList[i]:
                        rhat -= self.L[j,i]*self.du[j]
                    self.du[i] = self.M[i]*rhat
                if self.sym == True:
                    u-= self.du
                    self.computeResidual(u,r,b)
                    self.du[:]=0.0
                    for i in range(self.n-1,-1,-1):
                        rhat = self.r[i]
                        for j in self.connectionList[i]:
                            rhat -= self.L[i,j]*self.du[j]
                    self.du[i] = self.M[i]*rhat
            elif type(self.L).__name__ == "SparseMatrix":
                self.csmoothers.gauss_seidel_NR_solve(self.L,self.M,r,self.node_order,self.du)
                #self.csmoothers.jacobi_NR_solve(self.L,self.M,r,self.node_order,self.du)
            u -= self.du
            self.computeResidual(u,r,b)

class StarILU(LinearSolver):
    """
    Alternating Schwarz Method on node stars.
    """
    import csmoothers
    def __init__(self,
                 connectionList,
                 L,
                 weight=1.0,
                 sym=False,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        LinearSolver.__init__(self,L,
                              rtol_r,
                              atol_r,
                              rtol_du,
                              atol_du,
                              maxIts,
                              norm,
                              convergenceTest,
                              computeRates,
                              printInfo)
        self.solverName = "StarILU"
        self.w=weight
        self.sym=sym
        if type(self.L).__name__ == 'ndarray':
            self.connectionList=connectionList
            self.subdomainIndecesList=[]
            self.subdomainSizeList=[]
            self.subdomainL=[]
            self.subdomainR=[]
            self.subdomainDU=[]
            self.subdomainSolvers=[]
            self.globalToSubdomain=[]
            for i in range(self.n):
                self.subdomainIndecesList.append([])
                connectionList[i].sort()
                self.globalToSubdomain.append(dict([(j,J+1) for J,j in
                                                    enumerate(connectionList[i])]))
                self.globalToSubdomain[i][i]=0
                nSubdomain = len(connectionList[i])+1
                self.subdomainR.append(Vec(nSubdomain))
                self.subdomainDU.append(Vec(nSubdomain))
                self.subdomainSizeList.append(len(connectionList[i]))
                self.subdomainL.append(Mat(nSubdomain,nSubdomain))
                for J,j in enumerate(connectionList[i]):
                    self.subdomainIndecesList[i].append(set(connectionList[i]) &
                                                        set(connectionList[j]))
                    self.subdomainIndecesList[i][J].update([i,j])
        elif type(L).__name__ == 'SparseMatrix':
            self.node_order=numpy.arange(self.n,dtype="i")
            self.asmFactorObject = self.csmoothers.ASMFactor(L)
    def prepare(self,b=None):
        if type(self.L).__name__ == 'ndarray':
            self.subdomainSolvers=[]
            for i in range(self.n):
                self.subdomainL[i][0,0] = self.L[i,i]
                for J,j in enumerate(self.connectionList[i]):
                    #first do row 0 (star center)
                    self.subdomainL[i][J+1,0] = self.L[j,i]
                    #now do boundary rows
                    for k in self.subdomainIndecesList[i][J]:
                        K = self.globalToSubdomain[i][k]
                        self.subdomainL[i][K,J+1]=self.L[k,j]
                self.subdomainSolvers.append(LU(self.subdomainL[i]))
                self.subdomainSolvers[i].prepare()
        elif type(self.L).__name__ == 'SparseMatrix':
            self.csmoothers.asm_NR_prepare(self.L,self.asmFactorObject)
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        while (not self.converged(r) and
               not self.failed()):
            self.du[:]=0.0
            if type(self.L).__name__ == 'ndarray':
                for i in range(self.n):
                    #load subdomain residual
                    self.subdomainR[i][0] = r[i] - self.L[i,i]*self.du[i]
                    for j in self.connectionList[i]:
                        self.subdomainR[i][0] -= self.L[j,i]*self.du[j]
                    for J,j in enumerate(self.connectionList[i]):
                        self.subdomainR[i][J+1]=r[j] - self.L[j,j]*self.du[j]
                        for k in self.connectionList[j]:
                            self.subdomainR[i][J+1] -= self.L[k,j]*self.du[k]
                    #solve
                    self.subdomainSolvers[i].solve(u=self.subdomainDU[i],
                                                   b=self.subdomainR[i])
                    #update du
                    self.subdomainDU[i]*=self.w
                    self.du[i]+=self.subdomainDU[i][0]
                    for J,j in enumerate(self.connectionList[i]):
                        self.du[j] += self.subdomainDU[i][J+1]
            elif type(self.L).__name__ == 'SparseMatrix':
                self.csmoothers.asm_NR_solve(self.L,self.w,self.asmFactorObject,self.node_order,r,self.du)
            u -= self.du
            self.computeResidual(u,r,b)

class StarBILU(LinearSolver):
    """
    Alternating Schwarz Method on 'blocks' consisting of consectutive rows in system for things like dg ...
    """
    import csmoothers
    def __init__(self,
                 connectionList,
                 L,
                 bs=1,
                 weight=1.0,
                 sym=False,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        LinearSolver.__init__(self,L,
                              rtol_r,
                              atol_r,
                              rtol_du,
                              atol_du,
                              maxIts,
                              norm,
                              convergenceTest,
                              computeRates,
                              printInfo)
        self.solverName = "StarBILU"
        self.w=weight
        self.sym=sym
        self.bs = bs
        if type(self.L).__name__ == 'ndarray':
            raise NotImplementedError
        elif type(L).__name__ == 'SparseMatrix':
            self.node_order=numpy.arange(self.n,dtype="i")
            self.basmFactorObject = self.csmoothers.BASMFactor(L,bs)
    def prepare(self,b=None):
        if type(self.L).__name__ == 'ndarray':
            raise NotImplementedError
        elif type(self.L).__name__ == 'SparseMatrix':
            self.csmoothers.basm_NR_prepare(self.L,self.basmFactorObject)
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        while (not self.converged(r) and
               not self.failed()):
            #mwf debug
            logEvent("StarBILU norm_r= %s norm_du= %s " % (self.norm_r,self.norm_du))
            self.du[:]=0.0
            if type(self.L).__name__ == 'ndarray':
                raise NotImplementedError
            elif type(self.L).__name__ == 'SparseMatrix':
                self.csmoothers.basm_NR_solve(self.L,self.w,self.basmFactorObject,self.node_order,r,self.du)
            u -= self.du
            self.computeResidual(u,r,b)
class TwoLevel(LinearSolver):
    """
    A generic two-level multiplicative Schwarz solver.
    """
    def __init__(self,
                 prolong,
                 restrict,
                 coarseL,
                 preSmoother,
                 postSmoother,
                 coarseSolver,
                 L,
                 prepareCoarse=False,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        LinearSolver.__init__(self,L,
                              rtol_r,
                              atol_r,
                              rtol_du,
                              atol_du,
                              maxIts,
                              norm,
                              convergenceTest,
                              computeRates,
                              printInfo)
        self.solverName = "TwoLevel"
        self.prolong = prolong
        self.restrict = restrict
        self.cL = coarseL
        self.preSmoother = preSmoother
        self.postSmoother = postSmoother
        self.coarseSolver = coarseSolver
        self.cb = Vec(prolong.shape[1])
        self.cr = Vec(prolong.shape[1])
        self.cdu = Vec(prolong.shape[1])
        self.prepareCoarse=prepareCoarse
    def prepare(self,b=None):
        self.preSmoother.prepare()
        self.postSmoother.prepare()
        if self.prepareCoarse is True:
            self.coarseSolver.prepare()
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        (r,b) = self.solveInitialize(u,r,b,initialGuessIsZero)
        while (not self.converged(r) and
               not self.failed()):
            self.preSmoother.solve(u,r,b,initialGuessIsZero)
            initialGuessIsZero=False
            self.restrict.matvec(r,self.cb)
            self.cdu[:]=0.0
            self.coarseSolver.solve(u=self.cdu,r=self.cr,b=self.cb,initialGuessIsZero=True)
            self.prolong.matvec(self.cdu,self.du)
            u-=self.du
            self.computeResidual(u,r,b)
            self.postSmoother.solve(u,r,b,initialGuessIsZero=False)

class MultilevelLinearSolver:
    """
    A generic multilevel solver.
    """
    def __init__(self,levelLinearSolverList,computeRates=False,printInfo=False):
        self.printInfo=printInfo
        self.solverList=levelLinearSolverList
        self.nLevels = len(self.solverList)
        self.computeEigenvalues = False
        for l in range(self.nLevels):
            levelLinearSolverList[l].computeRates=computeRates
            levelLinearSolverList[l].printInfo=self.printInfo
    def info(self):
        self.infoString="********************Start Multilevel Linear Solver Info*********************\n"
        for l in range(self.nLevels):
            self.infoString += "**************Start Level %i Info********************\n" % l
            self.infoString += self.solverList[l].info()
            self.infoString += "**************End Level %i Info********************\n" % l
        self.infoString+="********************End Multilevel Linear Solver Info*********************\n"
        return self.infoString

class MGM(MultilevelLinearSolver):
    """
    A generic multigrid W cycle.
    """
    def __init__(self,
                 prolongList,
                 restrictList,
                 LList,
                 preSmootherList,
                 postSmootherList,
                 coarseSolver,
                 mgItsList=[],
                 printInfo=False,
                 computeRates=False):
        self.printInfo=printInfo
        self.nLevels = len(LList)
        self.solverList=[coarseSolver]
        for i in range(1,len(LList)):
            if mgItsList ==[]:
                mgItsList.append(1)
            self.solverList.append(TwoLevel(prolong = prolongList[i],
                                            restrict = restrictList[i],
                                            coarseL = LList[i-1],
                                            preSmoother = preSmootherList[i],
                                            postSmoother = postSmootherList[i],
                                            coarseSolver = self.solverList[i-1],
                                            L = LList[i],
                                            maxIts = mgItsList[i],
                                            convergenceTest = 'its',
                                            computeRates=computeRates,
                                            printInfo=False))
        self.mgmSolver = self.solverList[self.nLevels-1]
        self.solverName = "TwoLevel"
    def prepare(self,b=None):
        for s in self.solverList:
            s.prepare()

    def solve(self,u,r=None,b=None,initialGuessIsZero=False):
        self.mgmSolver.solve(u,r,b,initialGuessIsZero)

class NI(MultilevelLinearSolver):
    """
    A generic nested iteration solver.
    """
    def __init__(self,
                 solverList,
                 prolongList,
                 restrictList,
                 maxIts=None,
                 tolList=None,
                 atol=None,
                 computeRates=True,
                 printInfo=False):
        self.levelSolverList=solverList
        self.solverList = [self for n in range(len(solverList))]
        MultilevelLinearSolver.__init__(self,self.solverList,computeRates=computeRates,printInfo=printInfo)
        self.prolongList = prolongList
        self.restrictList = restrictList
        self.fineMesh = self.nLevels - 1
        self.uList=[]
        self.bList=[]
        self.levelDict={}
        for l in range(self.fineMesh+1):
            n = solverList[l].n
            self.levelDict[n] = l
            self.uList.append(Vec(n))
            self.bList.append(Vec(n))
        self.levelDict[solverList[-1].n]=self.fineMesh
        self.uList.append([])
        self.bList.append([])
        self.maxIts = maxIts
        self.tolList = tolList
        self.atol_r=atol
        self.printInfo=printInfo
        self.infoString=''
    def setResTol(self,rtol,atol):
        if self.tolList != None:
            for l in range(self.nLevels):
                self.tolList[l] = rtol
            self.atol_r = atol
    def prepare(self,b=None):
        if b is not None:
            currentMesh = self.levelDict[b.shape[0]]
        else:
            currentMesh = self.fineMesh
        for s in self.levelSolverList[:currentMesh+1]:
            s.prepare()
    def solve(self,u,r=None,b=None,par_u=None,par_b=None,initialGuessIsZero=False):
        currentMesh = self.levelDict[b.shape[0]]
        if currentMesh > 0:
            self.uList[currentMesh][:] = u
            self.bList[currentMesh][:] = b
            for l in range(currentMesh,1,-1):
                if not initialGuessIsZero:
                    self.restrictList[l].matvec(self.uList[l],self.uList[l-1])
                self.restrictList[l].matvec(self.bList[l],self.bList[l-1])
            if initialGuessIsZero:
                self.uList[0][:]=0.0
        for l in range(currentMesh):
            if self.tolList != None:
                self.switchToResidualConvergence(self.levelSolverList[l],
                                                 self.tolList[l])
            self.levelSolverList[l].solve(u=self.uList[l],b=self.bList[l],initialGuessIsZero=initialGuessIsZero)
            initialGuessIsZero=False
            if self.tolList != None:
                self.revertToFixedIteration(self.levelSolverList[l])
            if l < currentMesh -1:
                self.prolongList[l+1].matvec(self.uList[l],self.uList[l+1])
            else:
                self.prolongList[l+1].matvec(self.uList[l],u)
        if self.tolList != None:
            self.switchToResidualConvergence(self.levelSolverList[currentMesh],
                                             self.tolList[currentMesh])
        self.levelSolverList[currentMesh].solve(u,r,b,initialGuessIsZero)
        self.infoString += "**************Start Level %i Info********************\n" % currentMesh
        self.infoString+=self.levelSolverList[currentMesh].info()
        self.infoString += "**************End Level %i Info********************\n" % currentMesh
        if self.tolList != None:
            self.revertToFixedIteration(self.levelSolverList[currentMesh])
    def solveMultilevel(self,bList,uList,par_bList=None,par_uList=None,initialGuessIsZero=False):
        self.infoString="*************Start Multilevel Linear Solver Info*******************\n"
        for l in range(self.fineMesh):
            if self.tolList != None:
                self.switchToResidualConvergence(self.levelSolverList[l],self.tolList[l])
            self.levelSolverList[l].solve(u=uList[l],b=bList[l],initialGuessIsZero=initialGuessIsZero)
            initialGuessIsZero=False
            if self.tolList != None:
                self.revertToFixedIteration(self.levelSolverList[l])
            self.prolongList[l+1].matvec(uList[l],uList[l+1])
            self.infoString += "**************Start Level %i Info********************\n" % l
            self.infoString+=self.levelSolverList[l].info()
            self.infoString += "**************End Level %i Info********************\n" % l
        if self.tolList != None:
            self.switchToResidualConvergence(self.levelSolverList[self.fineMesh],self.tolList[self.fineMesh])
        self.levelSolverList[self.fineMesh].solve(u=uList[self.fineMesh],b=bList[self.fineMesh],initialGuessIsZero=initialGuessIsZero)
        self.infoString += "**************Start Level %i Info********************\n" % l
        self.infoString+=self.levelSolverList[self.fineMesh].info()
        self.infoString += "**************End Level %i Info********************\n" % l
        if self.tolList != None:
            self.revertToFixedIteration(self.levelSolverList[self.fineMesh])
        self.infoString+="********************End Multilevel Linear Solver Info*********************\n"
    def info(self):
        return self.infoString
    def switchToResidualConvergence(self,solver,rtol):
        self.saved_ctest = solver.convergenceTest
        self.saved_rtol_r = solver.rtol_r
        self.saved_atol_r = solver.atol_r
        self.saved_maxIts = solver.maxIts
        self.saved_printInfo = solver.printInfo
        solver.convergenceTest = 'r'
        solver.rtol_r = rtol
        solver.atol_r = self.atol_r
        solver.maxIts = self.maxIts
        solver.printInfo = self.printInfo
    def revertToFixedIteration(self,solver):
        solver.convergenceTest = self.saved_ctest
        solver.rtol_r = self.saved_rtol_r
        solver.atol_r = self.saved_atol_r
        solver.maxIts = self.saved_maxIts
        solver.printInfo = self.saved_printInfo

    def info(self):
        return self.infoString
def multilevelLinearSolverChooser(linearOperatorList,
                                  par_linearOperatorList,
                                  multilevelLinearSolverType=NI,
                                  relativeToleranceList=None,
                                  absoluteTolerance=1.0e-8,
                                  solverConvergenceTest='r',
                                  solverMaxIts=500,
                                  printSolverInfo=False,
                                  computeSolverRates=False,
                                  levelLinearSolverType=MGM,
                                  printLevelSolverInfo=False,
                                  computeLevelSolverRates=False,
                                  smootherType=Jacobi,
                                  boundaryNullSpace=False,
                                  prolongList=None,
                                  restrictList=None,
                                  connectivityListList=None,
                                  cycles=3,
                                  preSmooths=3,
                                  postSmooths=3,
                                  printSmootherInfo=False,
                                  computeSmootherRates=False,
                                  smootherConvergenceTest='its',
                                  relaxationFactor=None,
                                  computeEigenvalues=False,
                                  parallelUsesFullOverlap = True,
                                  par_duList=None,
                                  solver_options_prefix=None,
                                  linearSolverLocalBlockSize=1):
    """ A function for setting up a multilevel linear solver."""
    logEvent("multilevelLinearSolverChooser type= %s" % multilevelLinearSolverType)
    if (multilevelLinearSolverType == PETSc or
        multilevelLinearSolverType == KSP_petsc4py or
        multilevelLinearSolverType == LU or
        multilevelLinearSolverType == Jacobi or
        multilevelLinearSolverType == GaussSeidel or
        multilevelLinearSolverType == StarILU or
        multilevelLinearSolverType == StarBILU or
        multilevelLinearSolverType == MGM):
        levelLinearSolverType = multilevelLinearSolverType
        printLevelLinearSolverInfo = printSolverInfo
        computeLevelSolverRates = computeSolverRates
    nLevels = len(linearOperatorList)
    multilevelLinearSolver = None
    levelLinearSolverList = []
    levelLinearSolver = None
    if levelLinearSolverType == MGM:
        preSmootherList=[]
        postSmootherList=[]
        mgItsList=[]
        for l in range(nLevels):
            mgItsList.append(cycles)
            if l > 0:
                if smootherType == Jacobi:
                    if relaxationFactor == None:
                        relaxationFactor = 4.0/5.0
                    preSmootherList.append(Jacobi(L=linearOperatorList[l],
                                                  weight=relaxationFactor,
                                                  maxIts=preSmooths,
                                                  convergenceTest = smootherConvergenceTest,
                                                  computeRates = computeSmootherRates,
                                                  printInfo = printSmootherInfo))
                    postSmootherList.append(Jacobi(L=linearOperatorList[l],
                                                   weight=relaxationFactor,
                                                   maxIts=postSmooths,
                                                   convergenceTest = smootherConvergenceTest,
                                                   computeRates = computeSmootherRates,
                                                   printInfo = printSmootherInfo))
                elif smootherType == GaussSeidel:
                    if relaxationFactor == None:
                        relaxationFactor = 0.33
                    preSmootherList.append(GaussSeidel(connectionList = connectivityListList[l],
                                                       L=linearOperatorList[l],
                                                       weight=relaxationFactor,
                                                       maxIts =  preSmooths,
                                                       convergenceTest = smootherConvergenceTest,
                                                       computeRates = computeSmootherRates,
                                                       printInfo = printSmootherInfo))
                    postSmootherList.append(GaussSeidel(connectionList = connectivityListList[l],
                                                        L=linearOperatorList[l],
                                                        weight=relaxationFactor,
                                                        maxIts =  postSmooths,
                                                        convergenceTest = smootherConvergenceTest,
                                                        computeRates = computeSmootherRates,
                                                        printInfo = printSmootherInfo))
                elif smootherType == StarILU:
                    if relaxationFactor == None:
                        relaxationFactor = 1.0
                    preSmootherList.append(StarILU(connectionList = connectivityListList[l],
                                                   L=linearOperatorList[l],
                                                   weight=relaxationFactor,
                                                   maxIts =  preSmooths,
                                                   convergenceTest = smootherConvergenceTest,
                                                   computeRates = computeSmootherRates,
                                                   printInfo = printSmootherInfo))
                    postSmootherList.append(StarILU(connectionList = connectivityListList[l],
                                                    L=linearOperatorList[l],
                                                    weight=relaxationFactor,
                                                    maxIts =  postSmooths,
                                                    convergenceTest = smootherConvergenceTest,
                                                    computeRates = computeSmootherRates,
                                                    printInfo = printSmootherInfo))
                elif smootherType == StarBILU:
                    if relaxationFactor == None:
                        relaxationFactor = 1.0
                    preSmootherList.append(StarBILU(connectionList = connectivityListList[l],
                                                    L=linearOperatorList[l],
                                                    bs = linearSolverLocalBlockSize,
                                                    weight=relaxationFactor,
                                                    maxIts =  preSmooths,
                                                    convergenceTest = smootherConvergenceTest,
                                                    computeRates = computeSmootherRates,
                                                    printInfo = printSmootherInfo))
                    postSmootherList.append(StarBILU(connectionList = connectivityListList[l],
                                                     L=linearOperatorList[l],
                                                     bs = linearSolverLocalBlockSize,
                                                     weight=relaxationFactor,
                                                     maxIts =  postSmooths,
                                                     convergenceTest = smootherConvergenceTest,
                                                     computeRates = computeSmootherRates,
                                                     printInfo = printSmootherInfo))
                else:
                    logEvent("smootherType unrecognized")
            else:
                preSmootherList.append([])
                postSmootherList.append([])
                coarseSolver = LU(L=linearOperatorList[l])
        levelLinearSolver = MGM(prolongList = prolongList,
                                restrictList = restrictList,
                                LList = linearOperatorList,
                                preSmootherList = preSmootherList,
                                postSmootherList = postSmootherList,
                                coarseSolver = coarseSolver,
                                mgItsList = mgItsList,
                                printInfo = printLevelSolverInfo,
                                computeRates = computeLevelSolverRates)
        levelLinearSolverList = levelLinearSolver.solverList
    elif levelLinearSolverType == LU:
        for l in range(nLevels):
            levelLinearSolverList.append(LU(linearOperatorList[l],computeEigenvalues))
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == PETSc:
        for l in range(nLevels):
            levelLinearSolverList.append(PETSc(linearOperatorList[l],par_linearOperatorList[l],
                                               prefix=solver_options_prefix))
            if solverConvergenceTest == 'r-true' and par_duList != None:
                levelLinearSolverList[-1].useTrueResidualTest(par_duList[l])
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == KSP_petsc4py:
        for l in range(nLevels):
            levelLinearSolverList.append(KSP_petsc4py(linearOperatorList[l],par_linearOperatorList[l],
                                                      maxIts = solverMaxIts,
                                                      convergenceTest = solverConvergenceTest,
                                                      rtol_r = relativeToleranceList[l],
                                                      atol_r = absoluteTolerance,
                                                      computeRates = computeLevelSolverRates,
                                                      printInfo = printLevelLinearSolverInfo,
                                                      prefix=solver_options_prefix,
                                                      Preconditioner=smootherType,
                                                      connectionList = connectivityListList[l],
                                                      linearSolverLocalBlockSize = linearSolverLocalBlockSize,
                                                      bdyNullSpace=boundaryNullSpace))
            #if solverConvergenceTest == 'r-true' and par_duList != None:
            #    levelLinearSolverList[-1].useTrueResidualTest(par_duList[l])
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == Jacobi:
        if relaxationFactor == None:
            relaxationFactor = 4.0/5.0
        for l in range(nLevels):
            levelLinearSolverList.append(Jacobi(L=linearOperatorList[l],
                                                weight=relaxationFactor,
                                                maxIts = solverMaxIts,
                                                convergenceTest = solverConvergenceTest,
                                                rtol_r = relativeToleranceList[l],
                                                atol_r = absoluteTolerance,
                                                computeRates = computeLevelSolverRates,
                                                printInfo = printLevelSolverInfo))
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == GaussSeidel:
        if relaxationFactor == None:
            relaxationFactor=0.33
        for l in range(nLevels):
            levelLinearSolverList.append(GaussSeidel(connectionList = connectivityListList[l],
                                                     L=linearOperatorList[l],
                                                     weight = relaxationFactor,
                                                     maxIts = solverMaxIts,
                                                     convergenceTest = solverConvergenceTest,
                                                     rtol_r = relativeToleranceList[l],
                                                     atol_r = absoluteTolerance,
                                                     computeRates = computeLevelSolverRates,
                                                     printInfo = printLevelSolverInfo))
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == StarILU:
        if relaxationFactor == None:
            relaxationFactor=1.0
        for l in range(nLevels):
            levelLinearSolverList.append(StarILU(connectionList = connectivityListList[l],
                                                 L=linearOperatorList[l],
                                                 weight=relaxationFactor,
                                                 maxIts = solverMaxIts,
                                                 convergenceTest = solverConvergenceTest,
                                                 rtol_r = relativeToleranceList[l],
                                                 atol_r = absoluteTolerance,
                                                 computeRates = computeLevelSolverRates,
                                                 printInfo = printLevelSolverInfo))
        levelLinearSolver = levelLinearSolverList
    elif levelLinearSolverType == StarBILU:
        if relaxationFactor == None:
            relaxationFactor=1.0
        for l in range(nLevels):
            levelLinearSolverList.append(StarBILU(connectionList = connectivityListList[l],
                                                  L=linearOperatorList[l],
                                                  bs= linearSolverLocalBlockSize,
                                                  weight=relaxationFactor,
                                                  maxIts = solverMaxIts,
                                                  convergenceTest = solverConvergenceTest,
                                                  rtol_r = relativeToleranceList[l],
                                                  atol_r = absoluteTolerance,
                                                  computeRates = computeLevelSolverRates,
                                                  printInfo = printLevelSolverInfo))
        levelLinearSolver = levelLinearSolverList
    else:
        raise RuntimeError,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Unknown level linear solver "+ levelLinearSolverType
    if multilevelLinearSolverType == NI:
        multilevelLinearSolver = NI(solverList = levelLinearSolverList,
                                    prolongList = prolongList,
                                    restrictList = restrictList,
                                    maxIts  = solverMaxIts,
                                    tolList = relativeToleranceList,
                                    atol    = absoluteTolerance,
                                    printInfo= printSolverInfo,
                                    computeRates = computeSolverRates)
    elif (multilevelLinearSolverType == PETSc or
          multilevelLinearSolverType == KSP_petsc4py or
          multilevelLinearSolverType == LU or
          multilevelLinearSolverType == Jacobi or
          multilevelLinearSolverType == GaussSeidel or
          multilevelLinearSolverType == StarILU or
          multilevelLinearSolverType == StarBILU or
          multilevelLinearSolverType == MGM):
        multilevelLinearSolver = MultilevelLinearSolver(levelLinearSolverList,
                                                        computeRates = computeSolverRates,
                                                        printInfo=printSolverInfo)
    else:
        raise RuntimeError,"Unknown linear solver %s" % multilevelLinearSolverType
    if (levelLinearSolverType == LU):
        directSolverFlag=True
    else:
        directSolverFlag=False
    for levelSolver in multilevelLinearSolver.solverList:
        levelSolver.par_fullOverlap = parallelUsesFullOverlap
    return (multilevelLinearSolver,directSolverFlag)

## @}

if __name__ == '__main__':
    from LinearAlgebra import *
    import LinearSolvers
    from LinearSolvers import *
    import Gnuplot
    from Gnuplot import *
    from math import *
    from RandomArray import *
    gf = Gnuplot.Gnuplot()
    gf("set terminal x11")
    ginit = Gnuplot.Gnuplot()
    ginit("set terminal x11")
    gsol = Gnuplot.Gnuplot()
    gsol("set terminal x11")
    gsolNI = Gnuplot.Gnuplot()
    gsolNI("set terminal x11")
    gres = Gnuplot.Gnuplot()
    gres("set terminal x11")
    levels = 7
    n=2**levels + 1
    h =1.0/(n-1.0)
    freq=10
    uFine = uniform(0,1,(n))
    uFine[0]=0.0
    uFine[n-1]=0.0
    xFine = numpy.arange(0,1.0+h,h,dtype='d')
    bFine = (freq*2*pi)**2*numpy.sin(freq*2*pi*xFine)
    gf.plot(Gnuplot.Data(xFine,bFine))
    ginit.plot(Gnuplot.Data(xFine,uFine))
    uList=[]
    bList=[]
    prolongList=[]
    restrictList=[]
    LList=[]
    LDList=[]
    hList=[]
    meshList=[]
    preSmootherList=[]
    postSmootherList=[]
    mgItsList=[]
    for l in range(levels):
        N = 2**(l+1) + 1
        L = SparseMat_old(N-2,N-2,3*(N-2),sym=True)
        LD = Mat(N-2,N-2)
        H = 1.0/(N-1.0)
        hList.append(H)
        mgItsList.append(6)
        meshList.append(numpy.arange(0,1.0+H,H,dtype='d')[1:N-1])
        u = uniform(0,1,(N))
        u[0]  = 0.0
        u[N-1] = 0.0
        b = (freq*2*pi)**2*numpy.sin(freq*2*pi*meshList[l])
        uList.append(u[1:N-1])
        bList.append(b)
        beginAssembly(L)
        for i in range(N-2):
            L[i,i] = 2.0/H**2
            LD[i,i] = 2.0/H**2
            if i > 0:
                L[i,i-1] = -1.0/H**2
                LD[i,i-1] = -1.0/H**2
            if i < N-3:
                L[i,i+1] = -1.0/H**2
                LD[i,i+1] = -1.0/H**2
            endAssembly(L)
        LList.append(L)
        LDList.append(LD)
        if l > 0:
            cN = (N - 1)/2 + 1
            restrict = SparseMat_old(cN-2,N-2,3*(N-2))
            prolong = SparseMat_old(N-2,cN-2,3*(N-2))
            for i in range(cN-2):
                restrict[i,2*i]   = 1.0/4.0
                restrict[i,2*i+1] = 2.0/4.0
                restrict[i,2*i+2] = 1.0/4.0
                prolong[2*i,i] = 1.0/2.0
                prolong[2*i+1,i]= 2.0/2.0
                prolong[2*i+2,i]= 1.0/2.0
            restrict.to_csr()
            restrictList.append(restrict)
            prolong.to_csr()
            prolongList.append(prolong)
            N = cN
            preSmootherList.append(Jacobi(L,2.0/3.0,3))
            postSmootherList.append(Jacobi(L,2.0/3.0,3))
        else:
            restrictList.append([])
            prolongList.append([])
            preSmootherList.append([])
            postSmootherList.append([])
            coarseSolver = Jacobi(L,1.0,1)
    mgm = MGM(prolongList,restrictList,LList,preSmootherList,postSmootherList,coarseSolver,mgItsList)
    mgm.prepare()
    rnorm=1.0
    mgits = 0
    while rnorm > 1.0e-8 and mgits < 20:
        mgits +=1
        mgm.solve(u=uFine[1:n-1],b=bFine[1:n-1])
        rnorm = wl2Norm(mgm.residual(),h)
    gsol.plot(Gnuplot.Data(xFine,uFine,title='numerical solution-MGM'),
              Gnuplot.Data(xFine,numpy.sin(freq*2*pi*xFine),title='exact solution'))
    #gres.plot(Gnuplot.Data(x[1:n-1],mgm.smootherList[0].res,title='final residual'))
    ni = NI(mgm.solverList,prolongList,restrictList)
    ni.prepare()
    ni.solveMultilevel(bList,uList)
    rnorm = wl2Norm(ni.residual(),h)
    gsolNI.plot(Gnuplot.Data(meshList[-1],uList[-1],
                             title='numerical solution-NI'),
                Gnuplot.Data(meshList[-1],numpy.sin(freq*2*pi*meshList[-1]),
                             title='exact solution'))
    evals=[]
    for a,b,u,h in zip(LDList,bList,uList,hList):
        lu = LU(a,computeRes=True)
        lu.prepare(b)
        lu.solve(u,b)
        dev = DenseEigenvalues(a)
        dev.computeEigenvalues()
        evals.append(dev.eigenvalues)
        ratio = (max(abs(dev.eigenvalues))/min(abs(dev.eigenvalues)))*(h**2)
        print "k*h**2 %12.5E" % ratio
    gevals = Gnuplot.Gnuplot()
    gevals("set terminal x11")
    gevals.plot(Gnuplot.Data(evals[0],title='eigenvalues'))
    for ev in evals[1:]:
        gevals.replot(Gnuplot.Data(ev,title='eigenvalues'))
    raw_input('Please press return to continue... \n')


class StorageSet(set):
    def __init__(self,initializer=[],shape=(0,),storageType='d'):
        set.__init__(self,initializer)
        self.shape = shape
        self.storageType = storageType
    def allocate(self,storageDict):
        for k in self:
            storageDict[k] = numpy.zeros(self.shape,self.storageType)

class OperatorConstructor:
    """ Base class for operator constructors. """
    def __init__(self,model):
        self.OLT = model

    def _create_inv_scaled_mass_matrix(self):
        """ TODO - generalize ? / fix names ? """
        self._mass_val = self.OLT.nzval.copy()
        self._mass_val.fill(0.)
        self.TPInvScaledMassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                                 self.OLT.nFreeVDOF_global,
                                                 self.OLT.nnz,
                                                 self._mass_val,
                                                 self.OLT.colind,
                                                 self.OLT.rowptr)
        
    def _create_mass_matrix(self):
        """ TODO - generalize ? / fix names ? """
        self._mass_val = self.OLT.nzval.copy()
        self._mass_val.fill(0.)
        self.TPMassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                        self.OLT.nFreeVDOF_global,
                                        self.OLT.nnz,
                                        self._mass_val,
                                        self.OLT.colind,
                                        self.OLT.rowptr)

    def _create_laplace_matrix(self):
        """ TODO - generalize ? / fix names ? """
        self._laplace_val = self.OLT.nzval.copy()
        self._laplace_val.fill(0.)
        self.TPInvScaledLaplaceOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                                    self.OLT.nFreeVDOF_global,
                                                    self.OLT.nnz,
                                                    self._laplace_val,
                                                    self.OLT.colind,
                                                    self.OLT.rowptr)

    def _create_advection_matrix(self):
        """ TODO - genalize? / fix names? """
        self._advection_val = self.OLT.nzval.copy()
        self._advection_val.fill(0.)
        self.TPScaledAdvectionOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                                   self.OLT.nFreeVDOF_global,
                                                   self.OLT.nnz,
                                                   self._advection_val,
                                                   self.OLT.colind,
                                                   self.OLT.rowptr)        

            
class OperatorConstructor_rans2p(OperatorConstructor):
    """ A class for building common discrete rans2p operators.

    Arguments:
    ----------
    LevelModel : :class:`proteus.mprans.RANS2P.LevelModel`
        Level transport model derived from the rans2p class.
    """
    def __init__(self,levelModel):
        OperatorConstructor.__init__(self,levelModel)

    def attachTPAdvectionOperator(self):
        """ Create a discrete two-phase advection operator matrix. """
        self._create_advection_matrix()
        self.OLT.rans2p.getTwoPhaseAdvectionOperator(self.OLT.u[0].femSpace.elementMaps.psi,
                                                     self.OLT.u[0].femSpace.elementMaps.grad_psi,
                                                     self.OLT.mesh.nodeArray,
                                                     self.OLT.mesh.elementNodesArray,
                                                     self.OLT.elementQuadratureWeights[('u',0)],
                                                     self.OLT.u[0].femSpace.psi,
                                                     self.OLT.u[0].femSpace.grad_psi,
                                                     self.OLT.u[1].femSpace.psi,
                                                     self.OLT.u[1].femSpace.grad_psi,
                                                     self.OLT.elementDiameter,
                                                     self.OLT.mesh.nodeDiametersArray,
                                                     self.OLT.mesh.nElements_global,
                                                     self.OLT.coefficients.useMetrics,
                                                     self.OLT.coefficients.epsFact_density,
                                                     self.OLT.coefficients.epsFact,
                                                     1.,
                                                     1.,
                                                     1.,
                                                     1.,
                                                     # self.OLT.coefficients.rho_0,
                                                     # self.OLT.coefficients.nu_0,
                                                     # self.OLT.coefficients.rho_1,
                                                     # self.OLT.coefficients.nu_1,
                                                     self.OLT.u[1].femSpace.dofMap.l2g,
                                                     self.OLT.u[1].dof,
                                                     self.OLT.u[2].dof,
                                                     self.OLT.coefficients.useVF,
                                                     self.OLT.coefficients.q_vf,
                                                     self.OLT.coefficients.q_phi,
                                                     self.OLT.csrRowIndeces[(0,0)],self.OLT.csrColumnOffsets[(0,0)],
                                                     self.OLT.csrRowIndeces[(1,1)],self.OLT.csrColumnOffsets[(1,1)],
                                                     self.OLT.csrRowIndeces[(2,2)],self.OLT.csrColumnOffsets[(2,2)],
                                                     self.TPScaledAdvectionOperator)

    def attachTPInvScaledLaplaceOperator(self):
        """ Create a discrete two phase laplace operator matrix. """
        self._create_laplace_matrix()
        self.OLT.rans2p.getTwoPhaseInvScaledLaplaceOperator(self.OLT.u[0].femSpace.elementMaps.psi,
                                                            self.OLT.u[0].femSpace.elementMaps.grad_psi,
                                                            self.OLT.mesh.nodeArray,
                                                            self.OLT.mesh.elementNodesArray,
                                                            self.OLT.elementQuadratureWeights[('u',0)],
                                                            self.OLT.u[0].femSpace.grad_psi,
                                                            self.OLT.u[1].femSpace.grad_psi,
                                                            self.OLT.elementDiameter,
                                                            self.OLT.mesh.nodeDiametersArray,
                                                            self.OLT.mesh.nElements_global,
                                                            self.OLT.coefficients.useMetrics,
                                                            self.OLT.coefficients.epsFact_density,
                                                            self.OLT.coefficients.epsFact,
                                                            # 1.,
                                                            # 1.,
                                                            # 1.,
                                                            # 1.,
                                                            self.OLT.coefficients.rho_0,
                                                            self.OLT.coefficients.nu_0,
                                                            self.OLT.coefficients.rho_1,
                                                            self.OLT.coefficients.nu_1,
                                                            self.OLT.u[0].femSpace.dofMap.l2g,
                                                            self.OLT.u[1].femSpace.dofMap.l2g,
                                                            self.OLT.u[0].dof,
                                                            self.OLT.u[1].dof,
                                                            self.OLT.u[2].dof,
                                                            self.OLT.coefficients.useVF,
                                                            self.OLT.coefficients.q_vf,
                                                            self.OLT.coefficients.q_phi,
                                                            self.OLT.coefficients.sdInfo[(1,1)][0],self.OLT.coefficients.sdInfo[(1,1)][1], # ARB - this should work..?
                                                            self.OLT.coefficients.sdInfo[(1,1)][0],self.OLT.coefficients.sdInfo[(1,1)][1],
                                                            self.OLT.coefficients.sdInfo[(2,2)][0],self.OLT.coefficients.sdInfo[(2,2)][1],
                                                            self.OLT.csrRowIndeces[(0,0)],self.OLT.csrColumnOffsets[(0,0)],
                                                            self.OLT.csrRowIndeces[(1,1)],self.OLT.csrColumnOffsets[(1,1)],
                                                            self.OLT.csrRowIndeces[(2,2)],self.OLT.csrColumnOffsets[(2,2)],
                                                            self.TPInvScaledLaplaceOperator)

    def attachTwoPhaseMassOperator_rho(self):
        """ Create a discrete TwoPhase Mass operator matrix. """
        self._create_mass_matrix()
        self.OLT.rans2p.getTwoPhaseScaledMassOperator(1,
                                                     self.OLT.u[0].femSpace.elementMaps.psi,
                                                     self.OLT.u[0].femSpace.elementMaps.grad_psi,
                                                     self.OLT.mesh.nodeArray,
                                                     self.OLT.mesh.elementNodesArray,
                                                     self.OLT.elementQuadratureWeights[('u',0)],
                                                     self.OLT.u[0].femSpace.psi,
                                                     self.OLT.u[0].femSpace.psi,
                                                     self.OLT.u[1].femSpace.psi,
                                                     self.OLT.u[1].femSpace.psi,
                                                     self.OLT.elementDiameter,
                                                     self.OLT.mesh.nodeDiametersArray,
                                                     self.OLT.mesh.nElements_global,
                                                     self.OLT.coefficients.useMetrics,
                                                     self.OLT.coefficients.epsFact_density,
                                                     self.OLT.coefficients.epsFact,
                                                     1.,
                                                     1.,
                                                     1.,
                                                     1.,
                                                     # self.OLT.coefficients.rho_0,
                                                     # self.OLT.coefficients.nu_0,
                                                     # self.OLT.coefficients.rho_1,
                                                     # self.OLT.coefficients.nu_1,
                                                     self.OLT.u[0].femSpace.dofMap.l2g,
                                                     self.OLT.u[1].femSpace.dofMap.l2g,
                                                     self.OLT.u[0].dof,
                                                     self.OLT.u[1].dof,
                                                     self.OLT.u[2].dof,
                                                     self.OLT.coefficients.useVF,
                                                     self.OLT.coefficients.q_vf,
                                                     self.OLT.coefficients.q_phi,
                                                     self.OLT.csrRowIndeces[(0,0)],self.OLT.csrColumnOffsets[(0,0)],
                                                     self.OLT.csrRowIndeces[(1,1)],self.OLT.csrColumnOffsets[(1,1)],
                                                     self.OLT.csrRowIndeces[(2,2)],self.OLT.csrColumnOffsets[(2,2)],
                                                     self.TPMassOperator)    
        
        

    def attachTwoPhaseInvScaledMassOperator(self):
        """Create a discrete TwoPhase Mass operator matrix. """
        self._create_inv_scaled_mass_matrix()
        self.OLT.rans2p.getTwoPhaseScaledMassOperator(0,
                                                     self.OLT.u[0].femSpace.elementMaps.psi,
                                                     self.OLT.u[0].femSpace.elementMaps.grad_psi,
                                                     self.OLT.mesh.nodeArray,
                                                     self.OLT.mesh.elementNodesArray,
                                                     self.OLT.elementQuadratureWeights[('u',0)],
                                                     self.OLT.u[0].femSpace.psi,
                                                     self.OLT.u[0].femSpace.psi,
                                                     self.OLT.u[1].femSpace.psi,
                                                     self.OLT.u[1].femSpace.psi,
                                                     self.OLT.elementDiameter,
                                                     self.OLT.mesh.nodeDiametersArray,
                                                     self.OLT.mesh.nElements_global,
                                                     self.OLT.coefficients.useMetrics,
                                                     self.OLT.coefficients.epsFact_density,
                                                     self.OLT.coefficients.epsFact,
                                                     self.OLT.coefficients.rho_0,
                                                     self.OLT.coefficients.nu_0,
                                                     self.OLT.coefficients.rho_1,
                                                     self.OLT.coefficients.nu_1,
                                                     self.OLT.u[0].femSpace.dofMap.l2g,
                                                     self.OLT.u[1].femSpace.dofMap.l2g,
                                                     self.OLT.u[0].dof,
                                                     self.OLT.u[1].dof,
                                                     self.OLT.u[2].dof,
                                                     self.OLT.coefficients.useVF,
                                                     self.OLT.coefficients.q_vf,
                                                     self.OLT.coefficients.q_phi,
                                                     self.OLT.csrRowIndeces[(0,0)],self.OLT.csrColumnOffsets[(0,0)],
                                                     self.OLT.csrRowIndeces[(1,1)],self.OLT.csrColumnOffsets[(1,1)],
                                                     self.OLT.csrRowIndeces[(2,2)],self.OLT.csrColumnOffsets[(2,2)],
                                                     self.TPInvScaledMassOperator)
        

class OperatorConstructor_oneLevel(OperatorConstructor):
    """ A class for building common discrete operators. 
    
    Arguments
    ---------
    OLT : :class:`proteus.Transport.OneLevelTransport`
        One level transport class from which operator construction 
        will be based.

    TODO - ARB: replace self.OLT with self.model.
    """
    def __init__(self,OLT):
        OperatorConstructor.__init__(self,OLT)
        self._initializeOperatorConstruction()
        self.massOperatorAttached = False
        self.Qv_constructed = False
        self.laplaceOperatorAttached = False
        self.advectionOperatorAttached = False
        self.BOperatorAttached = False

    def attachMassOperator(self, rho=1., recalculate=False):
        """Create a discrete Mass Operator matrix. """
#        import pdb ; pdb.set_trace()
        self._mass_val = self.OLT.nzval.copy()
        self._mass_val.fill(0.)
        self.MassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                      self.OLT.nFreeVDOF_global,
                                      self.OLT.nnz,
                                      self._mass_val,
                                      self.OLT.colind,
                                      self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.rho != None:
            _rho = self.OLT.coefficients.rho
        self.MassOperatorCoeff = TransportCoefficients.DiscreteMassMatrix(rho=_rho, nd=_nd)
        _t = 1.0

        Mass_q = {}
        self._allocateMassOperatorQStorageSpace(Mass_q)
        if _nd == 2:
            self.MassOperatorCoeff.evaluate(_t,Mass_q)
        self._calculateMassOperatorQ(Mass_q)
        
        Mass_Jacobian = {}
        self._allocateMatrixSpace(self.MassOperatorCoeff,
                                  Mass_Jacobian)

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak(Mass_q[('dm',ci,cj)],
                                                      Mass_q[('vXw*dV_m',cj,ci)],
                                                      Mass_Jacobian[ci][cj])
        self._createOperator(self.MassOperatorCoeff,Mass_Jacobian,self.MassOperator)
        self.massOperatorAttached = True

    def attachTwoPhaseMassOperator_mu(self,recalculate=False, phase_function = None):
        """Create a discrete Mass Operator matrix. """
        self._mass_val = self.OLT.nzval.copy()
        
        _rho_0 = self.OLT.coefficients.rho_0
        _rho_1 = self.OLT.coefficients.rho_1
        _nu_0 = self.OLT.coefficients.nu_0
        _nu_1 = self.OLT.coefficients.nu_1
        
        self._mass_val.fill(0.)
        self.TPMassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                      self.OLT.nFreeVDOF_global,
                                      self.OLT.nnz,
                                      self._mass_val,
                                      self.OLT.colind,
                                      self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.rho != None:
            _rho = self.OLT.coefficients.rho

        if phase_function == None:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseMassMatrix_mu(nd = _nd,
                                                                                         rho_0 = _rho_0,
                                                                                         nu_0 = _nu_0,
                                                                                         rho_1 = _rho_1,
                                                                                         nu_1 = _nu_1,
                                                                                         LS_model = _phase_func)
        else:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseMassMatrix_mu(nd = _nd,
                                                                                         rho_0 = _rho_0,
                                                                                         nu_0 = _nu_0,
                                                                                         rho_1 = _rho_1,
                                                                                         nu_1 = _nu_1,
                                                                                         phase_function = phase_function)
            
        _t = 1.0

        Mass_q = {}
        self._allocateTwoPhaseMassOperatorQStorageSpace(Mass_q)
        self._calculateQuadratureValues(Mass_q)
        if _nd == 2:
            self.MassOperatorCoeff.evaluate(_t,Mass_q)
        self._calculateTwoPhaseMassOperatorQ(Mass_q)
        
        Mass_Jacobian = {}
        self._allocateMatrixSpace(self.MassOperatorCoeff,
                                  Mass_Jacobian)

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak(Mass_q[('dm',ci,cj)],
                                                      Mass_q[('vXw*dV_m',cj,ci)],
                                                      Mass_Jacobian[ci][cj])

        self._createOperator(self.MassOperatorCoeff,Mass_Jacobian,self.TPMassOperator)
        self.massOperatorAttached = True


    def attachTwoPhaseMassOperator_rho(self,recalculate=False, phase_function = None):
        """Create a discrete Mass Operator matrix. """
        self._mass_val = self.OLT.nzval.copy()
        
        _rho_0 = self.OLT.coefficients.rho_0
        _rho_1 = self.OLT.coefficients.rho_1
        _nu_0 = self.OLT.coefficients.nu_0
        _nu_1 = self.OLT.coefficients.nu_1
        
        self._mass_val.fill(0.)
        self.TPMassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                      self.OLT.nFreeVDOF_global,
                                      self.OLT.nnz,
                                      self._mass_val,
                                      self.OLT.colind,
                                      self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.rho != None:
            _rho = self.OLT.coefficients.rho

        if phase_function == None:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseMassMatrix(nd = _nd,
                                                                                      rho_0 = _rho_0,
                                                                                      nu_0 = _nu_0,
                                                                                      rho_1 = _rho_1,
                                                                                      nu_1 = _nu_1,
                                                                                      LS_model = _phase_func)
        else:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseMassMatrix(nd = _nd,
                                                                                      rho_0 = _rho_0,
                                                                                      nu_0 = _nu_0,
                                                                                      rho_1 = _rho_1,
                                                                                      nu_1 = _nu_1,
                                                                                      phase_function = phase_function)
            
        _t = 1.0

        Mass_q = {}
        self._allocateTwoPhaseMassOperatorQStorageSpace(Mass_q)
        self._calculateQuadratureValues(Mass_q)
        if _nd == 2:
            self.MassOperatorCoeff.evaluate(_t,Mass_q)
        self._calculateTwoPhaseMassOperatorQ(Mass_q)
        
        Mass_Jacobian = {}
        self._allocateMatrixSpace(self.MassOperatorCoeff,
                                  Mass_Jacobian)

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak(Mass_q[('dm',ci,cj)],
                                                      Mass_q[('vXw*dV_m',cj,ci)],
                                                      Mass_Jacobian[ci][cj])

        self._createOperator(self.MassOperatorCoeff,Mass_Jacobian,self.TPMassOperator)
        self.massOperatorAttached = True

    def attachTwoPhaseInvScaledMassOperator(self, recalculate=False, phase_function = None):
        """Create a discrete Mass Operator matrix. """
        self._mass_val = self.OLT.nzval.copy()
        _rho_0 = self.OLT.coefficients.rho_0
        _rho_1 = self.OLT.coefficients.rho_1
        _nu_0 = self.OLT.coefficients.nu_0
        _nu_1 = self.OLT.coefficients.nu_1

        self._mass_val.fill(0.)
        self.TPInvScaledMassOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                                 self.OLT.nFreeVDOF_global,
                                                 self.OLT.nnz,
                                                 self._mass_val,
                                                 self.OLT.colind,
                                                 self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.rho != None:
            _rho = self.OLT.coefficients.rho

        if phase_function is None:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseInvScaledMassMatrix(nd = _nd,
                                                                                               rho_0 = _rho_0,
                                                                                               nu_0 = _nu_0,
                                                                                               rho_1 = _rho_1,
                                                                                               nu_1 = _nu_1,
                                                                                               LS_model = _phase_func)
        else:
            self.MassOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseInvScaledMassMatrix(nd = _nd,
                                                                                               rho_0 = _rho_0,
                                                                                               nu_0 = _nu_0,
                                                                                               rho_1 = _rho_1,
                                                                                               nu_1 = _nu_1,
                                                                                               phase_function = phase_function)
            
        _t = 1.0

        Mass_q = {}
        self._allocateTwoPhaseMassOperatorQStorageSpace(Mass_q)
        self._calculateQuadratureValues(Mass_q)
        if _nd == 2:
            self.MassOperatorCoeff.evaluate(_t,Mass_q)
        self._calculateTwoPhaseMassOperatorQ(Mass_q)
        
        Mass_Jacobian = {}
        self._allocateMatrixSpace(self.MassOperatorCoeff,
                                  Mass_Jacobian)

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak(Mass_q[('dm',ci,cj)],
                                                      Mass_q[('vXw*dV_m',cj,ci)],
                                                      Mass_Jacobian[ci][cj])

        self._createOperator(self.MassOperatorCoeff,Mass_Jacobian,self.TPInvScaledMassOperator)
 
    def attachLaplaceOperator(self,nu=1.0):
        """ Create a Discrete Laplace Operator matrix."""
        self._laplace_val = self.OLT.nzval.copy()
        self._laplace_val.fill(0.)
        self.LaplaceOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                         self.OLT.nFreeVDOF_global,
                                         self.OLT.nnz,
                                         self._laplace_val,
                                         self.OLT.colind,
                                         self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.nu != None:
            _nu = self.OLT.coefficients.nu
        self.LaplaceOperatorCoeff = TransportCoefficients.DiscreteLaplaceOperator(nd=_nd)
        _t = 1.0

        Laplace_phi = {}
        Laplace_dphi = {}
        self._initializeLaplacePhiFunctions(Laplace_phi,Laplace_dphi)
        self._initializeSparseDiffusionTensor(self.LaplaceOperatorCoeff)

        Laplace_q = {}
        self._allocateLaplaceOperatorQStorageSpace(Laplace_q)
        if _nd==2:
            self.LaplaceOperatorCoeff.evaluate(_t,Laplace_q)
        self._calculateLaplaceOperatorQ(Laplace_q)
        
        Laplace_Jacobian = {}
        self._allocateMatrixSpace(self.LaplaceOperatorCoeff,
                                  Laplace_Jacobian)

        for ci,ckDict in self.LaplaceOperatorCoeff.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                for cj in set(cjDict.keys()+self.LaplaceOperatorCoeff.potential[ck].keys()):
                    cfemIntegrals.updateDiffusionJacobian_weak_sd(self.LaplaceOperatorCoeff.sdInfo[(ci,ck)][0],
                                                                  self.LaplaceOperatorCoeff.sdInfo[(ci,ck)][1],
                                                                  self.OLT.phi[ck].femSpace.dofMap.l2g, #??!!??
                                                                  Laplace_q[('a',ci,ck)],
                                                                  Laplace_q[('da',ci,ck,cj)],
                                                                  Laplace_q[('grad(phi)',ck)],
                                                                  Laplace_q[('grad(w)*dV_a',ck,ci)],
                                                                  Laplace_dphi[(ck,cj)].dof,
                                                                  self._operatorQ[('v',cj)],
                                                                  self._operatorQ[('grad(v)',cj)],
                                                                  Laplace_Jacobian[ci][cj])
        self._createOperator(self.LaplaceOperatorCoeff,
                             Laplace_Jacobian,
                             self.LaplaceOperator)
        self.laplaceOperatorAttached = True

    def attachTPInvScaledLaplaceOperator(self, phase_function = None):
        """ Create a Discrete Laplace Operator matrix."""
        self._laplace_val = self.OLT.nzval.copy()
        self._laplace_val.fill(0.)

        _rho_0 = self.OLT.coefficients.rho_0
        _rho_1 = self.OLT.coefficients.rho_1
        _nu_0 = self.OLT.coefficients.nu_0
        _nu_1 = self.OLT.coefficients.nu_1

        self.TPInvScaledLaplaceOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                                    self.OLT.nFreeVDOF_global,
                                                    self.OLT.nnz,
                                                    self._laplace_val,
                                                    self.OLT.colind,
                                                    self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        if self.OLT.coefficients.nu != None:
            _nu = self.OLT.coefficients.nu

        if phase_function == None:
            self.LaplaceOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseInvScaledLaplaceOperator(nd=_nd,
                                                                                                       rho_0 = _rho_0,
                                                                                                       nu_0 = _nu_0,
                                                                                                       rho_1 = _rho_1,
                                                                                                       nu_1 = _nu_1,
                                                                                                       LS_model = _phase_func)
        else:
            self.LaplaceOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseInvScaledLaplaceOperator(nd=_nd,
                                                                                                       rho_0 = _rho_0,
                                                                                                       nu_0 = _nu_0,
                                                                                                       rho_1 = _rho_1,
                                                                                                       nu_1 = _nu_1,
                                                                                                       phase_function = phase_function)

            
        _t = 1.0

        Laplace_phi = {}
        Laplace_dphi = {}
        self._initializeLaplacePhiFunctions(Laplace_phi,Laplace_dphi)
        self._initializeSparseDiffusionTensor(self.LaplaceOperatorCoeff)

        Laplace_q = {}
        self._allocateLaplaceOperatorQStorageSpace(Laplace_q)
        self._calculateQuadratureValues(Laplace_q)
        if _nd==2:
            self.LaplaceOperatorCoeff.evaluate(_t,Laplace_q)
        self._calculateLaplaceOperatorQ(Laplace_q)
        
        Laplace_Jacobian = {}
        self._allocateMatrixSpace(self.LaplaceOperatorCoeff,
                                  Laplace_Jacobian)

        for ci,ckDict in self.LaplaceOperatorCoeff.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                for cj in set(cjDict.keys()+self.LaplaceOperatorCoeff.potential[ck].keys()):
                    cfemIntegrals.updateDiffusionJacobian_weak_sd(self.LaplaceOperatorCoeff.sdInfo[(ci,ck)][0],
                                                                  self.LaplaceOperatorCoeff.sdInfo[(ci,ck)][1],
                                                                  self.OLT.phi[ck].femSpace.dofMap.l2g, #??!!??
                                                                  Laplace_q[('a',ci,ck)],
                                                                  Laplace_q[('da',ci,ck,cj)],
                                                                  Laplace_q[('grad(phi)',ck)],
                                                                  Laplace_q[('grad(w)*dV_a',ck,ci)],
                                                                  Laplace_dphi[(ck,cj)].dof,
                                                                  self._operatorQ[('v',cj)],
                                                                  self._operatorQ[('grad(v)',cj)],
                                                                  Laplace_Jacobian[ci][cj])
        self._createOperator(self.LaplaceOperatorCoeff,
                             Laplace_Jacobian,
                             self.TPInvScaledLaplaceOperator)
        self.laplaceOperatorAttached = True



    def attachAdvectionOperator(self,advective_field):
        """Attach a Discrete Advection Operator to the Transport class.
        
        Arguments
        ---------
        advective_field : numpy array
            numpy array of the advection field.
        """
        self._advective_field = advective_field
        self._advection_val = self.OLT.nzval.copy()
        self._advection_val.fill(0.)
        self.AdvectionOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                           self.OLT.nFreeVDOF_global,
                                           self.OLT.nnz,
                                           self._advection_val,
                                           self.OLT.colind,
                                           self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        self.AdvectionOperatorCoeff = TransportCoefficients.DiscreteAdvectionOperator(self._advective_field,
                                                                                      nd=_nd)
        _t = 1.0

        Advection_q = {}
        self._allocateAdvectionOperatorQStorageSpace(Advection_q)
        if _nd==2:
            self.AdvectionOperatorCoeff.evaluate(_t,Advection_q)
        self._calculateAdvectionOperatorQ(Advection_q)

        Advection_Jacobian = {}
        self._allocateMatrixSpace(self.AdvectionOperatorCoeff,
                                  Advection_Jacobian)
        
        for ci,ckDict in self.AdvectionOperatorCoeff.advection.iteritems():
            for ck in ckDict:
                cfemIntegrals.updateAdvectionJacobian_weak_lowmem(Advection_q[('df',ci,ck)],
                                                                  self._operatorQ[('v',ck)],
                                                                  Advection_q[('grad(w)*dV_f',ci)],
                                                                  Advection_Jacobian[ci][ck])
        self._createOperator(self.AdvectionOperatorCoeff,
                             Advection_Jacobian,
                             self.AdvectionOperator)
        self.advectionOperatorAttached = True

    def attachTPAdvectionOperator(self,advective_field,phase_function=None):
        """Attach a Discrete Advection Operator to the Transport class.
        
        Arguments
        ---------
        advective_field : numpy array
            numpy array of the advection field.
        """
        self._advective_field = advective_field
        self._advection_val = self.OLT.nzval.copy()
        self._advection_val.fill(0.)

        _rho_0 = self.OLT.coefficients.rho_0
        _rho_1 = self.OLT.coefficients.rho_1
        _nu_0 = self.OLT.coefficients.nu_0
        _nu_1 = self.OLT.coefficients.nu_1
        
        self.TPAdvectionOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                             self.OLT.nFreeVDOF_global,
                                             self.OLT.nnz,
                                             self._advection_val,
                                             self.OLT.colind,
                                             self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd

        if phase_function == None:
            self.AdvectionOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseAdvectionOperator(u = self._advective_field,
                                                                                                  nd = _nd,
                                                                                                  rho_0 = _rho_0,
                                                                                                  nu_0 = _nu_0,
                                                                                                  rho_1 = _rho_1,
                                                                                                  nu_1 = _nu_1,
                                                                                                  LS_model = _phase_func)
        else:
            self.AdvectionOperatorCoeff = TransportCoefficients.DiscreteTwoPhaseAdvectionOperator(u = self._advective_field,
                                                                                                  nd = _nd,
                                                                                                  rho_0 = _rho_0,
                                                                                                  nu_0 = _nu_0,
                                                                                                  rho_1 = _rho_1,
                                                                                                  nu_1 = _nu_1,
                                                                                                  phase_function = phase_function)

        _t = 1.0

        Advection_q = {}
        self._allocateAdvectionOperatorQStorageSpace(Advection_q)
        self._calculateQuadratureValues(Advection_q)
        if _nd==2:
            self.AdvectionOperatorCoeff.evaluate(_t,Advection_q)
        self._calculateAdvectionOperatorQ(Advection_q)

        Advection_Jacobian = {}
        self._allocateMatrixSpace(self.AdvectionOperatorCoeff,
                                  Advection_Jacobian)
        
        for ci,ckDict in self.AdvectionOperatorCoeff.advection.iteritems():
            for ck in ckDict:
                cfemIntegrals.updateAdvectionJacobian_weak_lowmem(Advection_q[('df',ci,ck)],
                                                                  self._operatorQ[('v',ck)],
                                                                  Advection_q[('grad(w)*dV_f',ci)],
                                                                  Advection_Jacobian[ci][ck])
        self._createOperator(self.AdvectionOperatorCoeff,
                             Advection_Jacobian,
                             self.TPAdvectionOperator)
        self.advectionOperatorAttached = True


    def attachBOperator(self):
        """Attach a discrete B operator to the Operator Constructor """
        self._B_val = self.OLT.nzval.copy()
        self._B_val.fill(0.)
        self.BOperator = SparseMat(self.OLT.nFreeVDOF_global,
                                   self.OLT.nFreeVDOF_global,
                                   self.OLT.nnz,
                                   self._B_val,
                                   self.OLT.colind,
                                   self.OLT.rowptr)
        _nd = self.OLT.coefficients.nd
        self.BOperatorCoeff = TransportCoefficients.DiscreteBOperator(nd=_nd)
        _t = 1.0
        
        B_q = {}
        self._allocateBOperatorQStorageSpace(B_q)

        if _nd==2:
            self.BOperatorCoeff.evaluate(_t,B_q)
        self._calculateBOperatorQ(B_q)
        
        B_Jacobian = {}
        self._allocateMatrixSpace(self.BOperatorCoeff,
                                  B_Jacobian)

        for ci,cjDict in self.BOperatorCoeff.advection.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateAdvectionJacobian_weak_lowmem(B_q[('df',ci,cj)],
                                                                  self._operatorQ[('v',cj)],
                                                                  B_q[('grad(w)*dV_f',ci)],
                                                                  B_Jacobian[ci][cj])

        for ci,cjDict in self.BOperatorCoeff.hamiltonian.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(B_q[('dH',ci,cj)],
                                                                    self._operatorQ[('grad(v)',cj)],
                                                                    B_q[('w*dV_H',ci)],
                                                                    B_Jacobian[ci][cj])
        self._createOperator(self.BOperatorCoeff,
                             B_Jacobian,
                             self.BOperator)

        self.BOperatorAttached = True

    def _initializeOperatorConstruction(self):
        """ Collect basic values used by all attach operators functions. """
        self._operatorQ = {}
        self._attachJacobianInfo(self._operatorQ)
        self._attachQuadratureInfo()
        self._attachTestInfo(self._operatorQ)
        self._attachTrialInfo(self._operatorQ)

    def _attachQuadratureInfo(self):
        """Define the quadrature type used to build operators.
        
        """
        pass
#        self._elementQuadrature = self.OLT._elementQuadrature
#        self._elementBoundaryQuadrature = self.OLT._elementBoundaryQuadrature

    def _attachJacobianInfo(self,Q):
        """ This helper function attaches quadrature data related to 'J'

        Arguments
        ---------
        Q : dict
            A dictionary to store values at quadrature points.
        """
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element) ) 
        tensor_quad = StorageSet(shape={})

        tensor_quad |= set(['J',
                            'inverse(J)'])
        scalar_quad |= set(['det(J)',
                            'abs(det(J))'])

        for k in tensor_quad:
            Q[k] = numpy.zeros((self.OLT.mesh.nElements_global,
                                self.OLT.nQuadraturePoints_element,
                                self.OLT.nSpace_global,
                                self.OLT.nSpace_global),
                                'd')

        scalar_quad.allocate(Q)
        
        self.OLT.u[0].femSpace.elementMaps.getJacobianValues(self.OLT.elementQuadraturePoints,
                                                             Q['J'],
                                                             Q['inverse(J)'],
                                                             Q['det(J)'])
        Q['abs(det(J))'] = numpy.absolute(Q['det(J)'])

    def _attachTestInfo(self,Q):
        """ Attach quadrature data for test functions.
        
        Arguments
        ---------
        Q : dict
            A dictionary to store values at quadrature points.

        Notes
        -----
        TODO - This function really doesn't need to compute the whole kitchen sink.
        Find a more efficient way to handle this.
        """
        test_shape_quad = StorageSet(shape={})
        test_shapeGradient_quad = StorageSet(shape={})

        test_shape_quad |= set([('w',ci) for ci in range(self.OLT.nc)])
        test_shapeGradient_quad |= set([('grad(w)',ci) for ci in range(self.OLT.nc)])
        
        for k in test_shape_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]]),
                'd')

        for k in test_shapeGradient_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]],
                 self.OLT.nSpace_global),
                'd')

        for ci in range(self.OLT.nc):
            if Q.has_key(('w',ci)):
                self.OLT.testSpace[ci].getBasisValues(self.OLT.elementQuadraturePoints,
                                                      Q[('w',ci)])
            if Q.has_key(('grad(w)',ci)):
                self.OLT.testSpace[ci].getBasisGradientValues(self.OLT.elementQuadraturePoints,
                                                              Q[('inverse(J)')],
                                                              Q[('grad(w)',ci)])

    def _attachTrialInfo(self,Q):
        """ Attach quadrature data for trial functions.

        Arguments
        ---------
        Q : dict
            A dictionary to store values at quadrature points.
        """
        trial_shape_quad = StorageSet(shape={})
        trial_shapeGrad_quad = StorageSet(shape={})
        
        trial_shape_quad |= set([('v',ci) for ci in range(self.OLT.nc)])
        trial_shapeGrad_quad |= set([('grad(v)',ci) for ci in range(self.OLT.nc)])

        for k in trial_shape_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]]),
                'd')

        for k in trial_shapeGrad_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]],
                 self.OLT.nSpace_global),
                'd')
                                            
        for ci in range(self.OLT.nc):
            if Q.has_key(('v',ci)):
                self.OLT.testSpace[ci].getBasisValues(self.OLT.elementQuadraturePoints,
                                                      Q[('v',ci)])
            if Q.has_key(('grad(v)',ci)):
                self.OLT.testSpace[ci].getBasisGradientValues(self.OLT.elementQuadraturePoints,
                                                              Q[('inverse(J)')],
                                                              Q[('grad(v)',ci)])
        
    def _allocateMatrixSpace(self,coeff,matrixDict):
        """ Allocate space for Operator Matrix """
        for ci in range(self.OLT.nc):
            matrixDict[ci] = {}
            for cj in range(self.OLT.nc):
                if cj in coeff.stencil[ci]:
                    matrixDict[ci][cj] = numpy.zeros(
                        (self.OLT.mesh.nElements_global,
                         self.OLT.nDOF_test_element[ci],
                         self.OLT.nDOF_trial_element[cj]),
                        'd')

    def _createOperator(self,coeff,matrixDict,A):
        """ Takes the matrix dictionary and creates a CSR matrix """
        for ci in range(self.OLT.nc):
            for cj in coeff.stencil[ci]:
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.OLT.l2g[ci]['nFreeDOF'],
                                                                          self.OLT.l2g[ci]['freeLocal'],
                                                                          self.OLT.l2g[cj]['nFreeDOF'],
                                                                          self.OLT.l2g[cj]['freeLocal'],
                                                                          self.OLT.csrRowIndeces[(ci,cj)],
                                                                          self.OLT.csrColumnOffsets[(ci,cj)],
                                                                          matrixDict[ci][cj],
                                                                          A)

    def _initializeLaplacePhiFunctions(self,Laplace_phi,Laplace_dphi):
        """ Initialize the phi functions for the Laplace operator """
        for ci,space in self.OLT.testSpace.iteritems():
            Laplace_phi[ci] = FemTools.FiniteElementFunction(space)

        for ck,phi in Laplace_phi.iteritems():
            Laplace_dphi[(ck,ck)] = FemTools.FiniteElementFunction(Laplace_phi[ck].femSpace)

        for ci,dphi in Laplace_dphi.iteritems():
            dphi.dof.fill(1.0)

    def _initializeSparseDiffusionTensor(self,coeff):
        """Intialize the sparse diffusion tensor for the lapalce matrix.

        Arguments
        ---------
        coeff : `proteus:TransporCoefficients:DiscreteLaplaceOperator`
        """
        for ci,ckDict in coeff.diffusion.iteritems():
            for ck in ckDict.keys():
                if not coeff.sdInfo.has_key((ci,ck)):
                    coeff.sdInfo[(ci,ck)] = (numpy.arange(start=0,stop=self.OLT.nSpace_global**2+1,
                                                          step=self.OLT.nSpace_global,
                                                          dtype='i'),
                                             numpy.array([range(self.OLT.nSpace_global) for row in range(self.OLT.nSpace_global)],
                                                         dtype='i'))

    def _allocateMassOperatorQStorageSpace(self,Q):
        """ Allocate space for mass operator values. """
        test_shape_quad = StorageSet(shape={})
        trial_shape_quad = StorageSet(shape={})
        trial_shape_X_test_shape_quad = StorageSet(shape={})
        tensor_quad = StorageSet(shape={})
        # TODO - ARB : I don't think the 3 is necessary here...It created a
        # confusing bug in the 2-phase problem...Need to investigate.
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element,
                                        3))
  
        scalar_quad |= set([('u',ci) for ci in range(self.OLT.nc)])
        scalar_quad |= set([('m',ci) for ci in self.MassOperatorCoeff.mass.keys()])

        test_shape_quad |= set([('w*dV_m',ci) for ci in self.MassOperatorCoeff.mass.keys()])

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            trial_shape_X_test_shape_quad |= set([('vXw*dV_m',cj,ci) for cj in cjDict.keys()])

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            scalar_quad |= set([('dm',ci,cj) for cj in cjDict.keys()])

        for k in tensor_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nSpace_global,
                 self.OLT.nSpace_global),
                'd')

        for k in test_shape_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]]),
                'd')

        for k in trial_shape_X_test_shape_quad:
            Q[k] = numpy.zeros((self.OLT.mesh.nElements_global,
                                self.OLT.nQuadraturePoints_element,
                                self.OLT.nDOF_trial_element[k[1]],
                                self.OLT.nDOF_test_element[k[2]]),'d')

        scalar_quad.allocate(Q)

    def _allocateTwoPhaseMassOperatorQStorageSpace(self,Q):
        """ Allocate space for mass operator values. """
        test_shape_quad = StorageSet(shape={})
        trial_shape_quad = StorageSet(shape={})
        trial_shape_X_test_shape_quad = StorageSet(shape={})
        tensor_quad = StorageSet(shape={})
        points_quadrature = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                              self.OLT.nQuadraturePoints_element,
                                              3))
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element))

        points_quadrature |= set(['x'])
        scalar_quad |= set([('u',ci) for ci in range(self.OLT.nc)])
        scalar_quad |= set([('m',ci) for ci in self.MassOperatorCoeff.mass.keys()])
        test_shape_quad |= set([('w*dV_m',ci) for ci in self.MassOperatorCoeff.mass.keys()])
        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            trial_shape_X_test_shape_quad |= set([('vXw*dV_m',cj,ci) for cj in cjDict.keys()])
        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            scalar_quad |= set([('dm',ci,cj) for cj in cjDict.keys()])

        for k in tensor_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nSpace_global,
                 self.OLT.nSpace_global),
                'd')

        for k in test_shape_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]]),
                'd')

        for k in trial_shape_X_test_shape_quad:
            Q[k] = numpy.zeros((self.OLT.mesh.nElements_global,
                                self.OLT.nQuadraturePoints_element,
                                self.OLT.nDOF_trial_element[k[1]],
                                self.OLT.nDOF_test_element[k[2]]),'d')

        scalar_quad.allocate(Q)
        points_quadrature.allocate(Q)

    def _allocateLaplaceOperatorQStorageSpace(self,Q):
        """Initialize the storage space for the Laplace operator vals. """
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element))
        tensors_quad = StorageSet(shape={})
        vectors_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                         self.OLT.nQuadraturePoints_element,
                                         self.OLT.nSpace_global))
        gradients = StorageSet(shape={})

        points_quadrature = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                              self.OLT.nQuadraturePoints_element,
                                              3))

        points_quadrature |= set(['x'])
        scalar_quad |= set([('u',ci) for ci in range(self.OLT.nc)])
        tensors_quad |= set([('a',ci,ci) for ci in range(self.OLT.nc)])
        tensors_quad |= set([('da',ci,ci,ci) for ci in range(self.OLT.nc)])

        for ci,ckDict in self.LaplaceOperatorCoeff.diffusion.iteritems():
            gradients |= set([('grad(w)*dV_a',ck,ci) for ck in ckDict])
        
        for ci,ckDict in self.LaplaceOperatorCoeff.diffusion.iteritems():
            vectors_quad |= set([('grad(phi)',ck) for ck in ckDict.keys()])

        scalar_quad.allocate(Q)
        vectors_quad.allocate(Q)
    
        for k in tensors_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.LaplaceOperatorCoeff.sdInfo[(k[1],k[2])][0][self.OLT.nSpace_global]),
                'd')

        for k in gradients:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]],
                 self.OLT.nSpace_global),
                'd')

        points_quadrature.allocate(Q)

    def _allocateAdvectionOperatorQStorageSpace(self,Q):
        """Allocate storage space for the Advection operator values. """
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element))
        points_quadrature = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                              self.OLT.nQuadraturePoints_element,
                                              3))
        vector_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element,
                                        self.OLT.nSpace_global))
        tensor_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element,
                                        self.OLT.nSpace_global))
        gradients = StorageSet(shape={})

        points_quadrature |= set(['x'])
        scalar_quad |= set([('u',0)])
        vector_quad |= set([('f',ci) for ci in range(self.OLT.nc)])
        tensor_quad |= set([('df',0,0)])

        for i in range(self.OLT.nc):
            for j in range(1,self.OLT.nc):
                tensor_quad |= set([('df',i,j)])

        gradients |= set([('grad(w)*dV_f',ci) for ci in self.AdvectionOperatorCoeff.advection.keys()])

        scalar_quad.allocate(Q)
        vector_quad.allocate(Q)
        
        for k in tensor_quad:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nSpace_global),
                'd')

        for k in gradients:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]],
                 self.OLT.nSpace_global),
                'd')

        points_quadrature.allocate(Q)

    def _allocateBOperatorQStorageSpace(self,Q):
        """Allocate storage space for the B-operator matrix. """
        scalar_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element))
        vector_quad = StorageSet(shape=(self.OLT.mesh.nElements_global,
                                        self.OLT.nQuadraturePoints_element,
                                        self.OLT.nSpace_global))
        test_shape_quad = StorageSet(shape={})
        trial_shape_X_test_grad_quad = StorageSet(shape={})
        gradients = StorageSet(shape={})

        scalar_quad |= set([('u',ci) for ci in xrange(self.OLT.nc)])
        scalar_quad |= set([('H',ci) for ci in xrange(self.OLT.nc)])

        vector_quad |= set([('grad(u)',0)])
        vector_quad |= set([('f',0)])
        for ci,cjDict in self.BOperatorCoeff.advection.iteritems():
            vector_quad |= set([('df',ci,cj) for cj in cjDict.keys()])
        for ci,cjDict in self.BOperatorCoeff.hamiltonian.iteritems():
            vector_quad |= set([('dH',ci,cj) for cj in cjDict.keys()])

        test_shape_quad |= set([('w*dV_H',ci) for ci in self.BOperatorCoeff.hamiltonian.keys()])

        for ci,cjDict in self.BOperatorCoeff.advection.iteritems():
            trial_shape_X_test_grad_quad |= set([('v_X_grad_w_dV',cj,ci) for cj in cjDict.keys()])

        gradients |= set([('grad(w)*dV_f',ci) for ci in self.BOperatorCoeff.advection.keys()])
        
        scalar_quad.allocate(Q)
        vector_quad.allocate(Q)

        for k in test_shape_quad:
            Q[k] = numpy.zeros((self.OLT.mesh.nElements_global,
                                self.OLT.nQuadraturePoints_element,
                                self.OLT.nDOF_test_element[k[-1]]),
                               'd')
        
        for k in trial_shape_X_test_grad_quad:
            Q[k] = numpy.zeros((self.OLT.mesh.nElements_global,
                                self.OLT.nQuadraturePoints_element,
                                self.OLT.nDOF_trial_element[k[1]],
                                self.OLT.nDOF_test_element[k[2]],
                                self.OLT.coefficients.nd),'d')

        for k in gradients:
            Q[k] = numpy.zeros(
                (self.OLT.mesh.nElements_global,
                 self.OLT.nQuadraturePoints_element,
                 self.OLT.nDOF_test_element[k[-1]],
                 self.OLT.nSpace_global),
                'd')


    def _calculateMassOperatorQ(self,Q):
        """ Calculate values for mass operator. """
        elementQuadratureDict = {}

        for ci in self.MassOperatorCoeff.mass.keys():
            elementQuadratureDict[('m',ci)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        for ci in range(self.OLT.nc):
            if Q.has_key(('w*dV_m',ci)):
                cfemIntegrals.calculateWeightedShape(elementQuadratureWeights[('m',ci)],
                                                     self._operatorQ['abs(det(J))'],
                                                     self._operatorQ[('w',ci)],
                                                     Q[('w*dV_m',ci)])

        for ci in zip(range(self.OLT.nc),range(self.OLT.nc)):
                cfemIntegrals.calculateShape_X_weightedShape(self._operatorQ[('v',ci[1])],
                                                             Q[('w*dV_m',ci[0])],
                                                             Q[('vXw*dV_m',ci[1],ci[0])])

    def _calculateQuadratureValues(self,Q):
        elementQuadratureDict = {}

        elementQuadratureDict[('m',1)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)

        self.OLT.u[0].femSpace.elementMaps.getValues(elementQuadraturePoints,
                                                     Q['x'])
    
    def _calculateTwoPhaseMassOperatorQ(self,Q):
        """ Calculate values for mass operator. """
        elementQuadratureDict = {}
        for ci in self.MassOperatorCoeff.mass.keys():
            elementQuadratureDict[('m',ci)] = self._elementQuadrature
        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        for ci in range(self.OLT.nc):
            if Q.has_key(('w*dV_m',ci)):
                cfemIntegrals.calculateWeightedShape(elementQuadratureWeights[('m',ci)],
                                                     self._operatorQ['abs(det(J))'],
                                                     self._operatorQ[('w',ci)],
                                                     Q[('w*dV_m',ci)])

        for ci in zip(range(self.OLT.nc),range(self.OLT.nc)):
                cfemIntegrals.calculateShape_X_weightedShape(self._operatorQ[('v',ci[1])],
                                                             Q[('w*dV_m',ci[0])],
                                                             Q[('vXw*dV_m',ci[1],ci[0])])
        self.OLT.u[0].femSpace.elementMaps.getValues(elementQuadraturePoints,
                                                     Q['x'])

    def _calculateLaplaceOperatorQ(self,Q):
        """Calculate quadrature values for Laplace operator. """
        elementQuadratureDict = {}
        
        for ci in self.LaplaceOperatorCoeff.diffusion.keys():
            elementQuadratureDict[('a',ci)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)

        for ck in range(self.OLT.nc):
            if ['grad(phi)',ck] in Q.keys():
                Q['grad(phi)',ck].fill(1.)

        for ci,ckDict in self.LaplaceOperatorCoeff.diffusion.iteritems():
            for ck in ckDict.keys():
                cfemIntegrals.calculateWeightedShapeGradients(elementQuadratureWeights[('a',ci)],
                                                              self._operatorQ['abs(det(J))'],
                                                              self._operatorQ[('grad(w)',ci)],
                                                              Q[('grad(w)*dV_a',ck,ci)])


    def _calculateAdvectionOperatorQ(self,Q):
        """Calculate quadrature values for Advection operator. """
        elementQuadratureDict = {}

        for ci in self.AdvectionOperatorCoeff.advection.keys():
            elementQuadratureDict[('f',ci)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
            
        for ci in range(self.OLT.nc):
            cfemIntegrals.calculateWeightedShapeGradients(elementQuadratureWeights[('f',ci)],
                                                          self._operatorQ['abs(det(J))'],
                                                          self._operatorQ[('grad(w)',ci)],
                                                          Q[('grad(w)*dV_f',ci)])
        
    
    def _calculateBOperatorQ(self,Q):
        """Calculate quadrature values for B operator """
        elementQuadratureDict = {}
        
        for ci in self.BOperatorCoeff.advection.keys():
            elementQuadratureDict[('f',ci)] = self._elementQuadrature
        for ci in self.BOperatorCoeff.hamiltonian.keys():
            elementQuadratureDict[('H',ci)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        
        for ci in self.BOperatorCoeff.advection.keys():
            cfemIntegrals.calculateWeightedShapeGradients(elementQuadratureWeights[('f',ci)],
                                                          self._operatorQ['abs(det(J))'],
                                                          self._operatorQ[('grad(w)',ci)],
                                                          Q[('grad(w)*dV_f',ci)])
        
        for ci in self.BOperatorCoeff.hamiltonian.keys():
            cfemIntegrals.calculateWeightedShape(elementQuadratureWeights[('H',ci)],
                                                 self._operatorQ['abs(det(J))'],
                                                 self._operatorQ[('w',ci)],
                                                 Q[('w*dV_H',ci)])
