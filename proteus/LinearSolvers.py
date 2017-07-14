"""
A hierarchy of classes for linear algebraic system solvers.

.. inheritance-diagram:: proteus.LinearSolvers
   :parts: 1
"""
from LinearAlgebraTools import *
import FemTools
import lapackWrappers
import superluWrappers
import TransportCoefficients
import cfemIntegrals
import Quadrature
from petsc4py import PETSc as p4pyPETSc
from math import *
from .Profiling import logEvent

class LinearSolver:
    """ The base class for linear solvers.  """
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
    """ A class that interfaces Proteus with PETSc KSP. """
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

        self.setResTol(rtol_r,atol_r)

        convergenceTest = 'r-true'
        if convergenceTest == 'r-true':
            self.r_work = self.petsc_L.getVecLeft()
            self.rnorm0 = None
            self.ksp.setConvergenceTest(self._converged_trueRes)
        else:
            self.r_work = None

        if prefix != None:
            self.ksp.setOptionsPrefix(prefix)
        if Preconditioner != None:
            self._setPreconditioner(Preconditioner,par_L,prefix)
            self.ksp.setPC(self.pc)
        self.ksp.max_it = self.maxIts
        self.ksp.setFromOptions()
        
    def setResTol(self,rtol,atol):
        """ Set the ksp object's residual and maximum interations. """
        self.rtol_r = rtol
        self.atol_r = atol
        self.ksp.rtol = rtol
        self.ksp.atol = atol
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
        if self.pc != None:
            self.pc.setOperators(self.petsc_L,self.petsc_L)
            self.pc.setUp()
            if self.preconditioner:
                self.preconditioner.setUp(self.ksp)
        self.ksp.setUp()
        self.ksp.pc.setUp()

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

        if self.bdyNullSpace==True:
            self._setNullSpace(par_b)
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
        """
        decide on convention to match norms, convergence criteria
        """
        return self.ksp.converged
    def failed(self):
        failedFlag = LinearSolver.failed(self)
        failedFlag = failedFlag or (not self.ksp.converged)
        return failedFlag
    
    def info(self):
        self.ksp.view()

    def _defineNullSpaceVec(self,par_b):
        """ A initial attempt to set up the null space vector
        in parallel.
        """
        #Global null space constructor
        # Comm information
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
        null_space_basis = ParVec_petsc4py(self.tmp_array,
                                           1,
                                           par_n,
                                           par_N,
                                           par_nghosts,
                                           subdomain2global,
                                           proteus2petsc_subdomain=proteus2petsc_subdomain,
                                           petsc2proteus_subdomain=petsc2proteus_subdomain)
        null_space_basis.scatter_forward_insert()
        self.tmp_array[:] = self.tmp_array[petsc2proteus_subdomain]
        self.global_null_space = [null_space_basis]


    def _setNullSpace(self,par_b):
        """ Set up the boundary null space for the KSP solves. 

        Parameters
        ----------
        par_b : proteus.LinearAlgebraTools.ParVec_petsc4py
            The problem's RHS vector.
        """
        self._defineNullSpaceVec(par_b)
        vecs = self.global_null_space
        self.pressure_null_space = p4pyPETSc.NullSpace().create(constant=False,
                                                                vectors=vecs,
                                                                comm=p4pyPETSc.COMM_WORLD)
        self.ksp.getOperators()[0].setNullSpace(self.pressure_null_space)
        self.pressure_null_space.remove(par_b)
        
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
            if Preconditioner == Jacobi:
                self.pccontext= Preconditioner(L,
                                               weight=1.0,
                                               rtol_r=rtol_r,
                                               atol_r=atol_r,
                                               maxIts=1,
                                               norm = l2Norm,
                                               convergenceTest='its',
                                               computeRates=False,
                                               printInfo=False)
                self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == GaussSeidel:
                self.pccontext= Preconditioner(connectionList,
                                               L,
                                               weight=1.0,
                                               sym=False,
                                               rtol_r=rtol_r,
                                               atol_r=atol_r,
                                               maxIts=1,
                                               norm = l2Norm,
                                               convergenceTest='its',
                                               computeRates=False,
                                               printInfo=False)
                self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == LU:
                #ARB - parallel matrices from PETSc4py don't work here.
                self.pccontext = Preconditioner(L)
                self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == StarILU:
                self.pccontext= Preconditioner(connectionList,
                                               L,
                                               weight=1.0,
                                               rtol_r=rtol_r,
                                               atol_r=atol_r,
                                               maxIts=1,
                                               norm = l2Norm,
                                               convergenceTest='its',
                                               computeRates=False,
                                               printInfo=False)
                self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == StarBILU:
                self.pccontext= Preconditioner(connectionList,
                                               L,
                                               bs=linearSolverLocalBlockSize,
                                               weight=1.0,
                                               rtol_r=rtol_r,
                                               atol_r=atol_r,
                                               maxIts=1,
                                               norm = l2Norm,
                                               convergenceTest='its',
                                               computeRates=False,
                                               printInfo=False)
                self.pc = p4pyPETSc.PC().createPython(self.pccontext)
            elif Preconditioner == SimpleNavierStokes3D:
                logEvent("NAHeader Preconditioner selfp" )
                self.preconditioner = SimpleNavierStokes3D(par_L,prefix)
                self.pc = self.preconditioner.pc
            elif Preconditioner == SimpleNavierStokes2D:
                logEvent("NAHeader Preconditioner selfp" )
                self.preconditioner = SimpleNavierStokes2D(par_L,prefix)
                self.pc = self.preconditioner.pc
            elif Preconditioner == Schur_Qp:
                logEvent("NAHeader Preconditioner Qp" )
                self.preconditioner = Schur_Qp(par_L,prefix,self.bdyNullSpace)
                self.pc = self.preconditioner.pc
            elif Preconditioner == SimpleDarcyFC:
                self.preconditioner = SimpleDarcyFC(par_L)
                self.pc = self.preconditioner.pc
            elif Preconditioner == NavierStokesPressureCorrection:
                self.preconditioner = NavierStokesPressureCorrection(par_L)
                self.pc = self.preconditioner.pc

class SchurOperatorConstructor:
    """ 
    Generate matrices for use in Schur complement preconditioner operators. 
    """
    def __init__(self, linear_smoother, pde_type='general_saddle_point'):
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
        try:
            if isinstance(self.L.pde,proteus.mprans.RANS2P.LevelModel):
                raise Exception, 'RANS2P not supported yet.'
        except:
            self.opBuilder = OperatorConstructor_oneLevel(self.L.pde)

    def getQp(self, output_matrix=False, recalculate=False):
        """ Return the pressure mass matrix Qp.

        Parameters
        ----------
        output_matrix : bool 
            Determines whether matrix should be exported.
        recalculate : bool
            Flag indicating whther matrix should be rebuilt every time it's used

        Returns
        -------
        Qp : matrix
            The pressure mass matrix.
        """
        Qsys_petsc4py = self._massMatrix(recalculate = recalculate)
        self.Qp = Qsys_petsc4py.getSubMatrix(self.linear_smoother.isp,
                                             self.linear_smoother.isp)
        if output_matrix is True:
            self._exportMatrix(self.Qp,"Qp")
        return self.Qp

    def _massMatrix(self,recalculate=False):
        """ Generates a mass matrix.

        This function generates and returns the mass matrix for the system. This
        function is internal to the class and called by public functions which 
        take and return the relevant components (eg. the pressure or velcoity).

        Parameters
        ----------
        recalculate : bool
            Indicates whether matrix should be rebuilt everytime it's called.

        Returns
        -------
        Qsys : matrix
            The system's mass matrix.
        """
        self.opBuilder.attachMassOperator(recalculate = recalculate)
        return superlu_2_petsc4py(self.opBuilder.MassOperator)
                
class KSP_Preconditioner:
    """ Base class for PETSCc KSP precondtioners. """
    def __init__(self):
        pass

    def setup(self, global_ksp=None):
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
        prefix : str
            Prefix handle for PETSc options.
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
        prefix : str
            Prefix handle for PETSc options.
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
    """ Base class for PETSc Schur complement preconditioners. """
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        """
        Initializes the Schur complement preconditioner for use with PETSc.

        This class creates a KSP PETSc solver object and initializes flags the
        pressure and velocity unknowns for a general saddle point problem.

        Parameters
        ----------
        L : provides the definition of the problem.
        prefix : str
            Prefix identifier for command line PETSc options.
        bdyNullSpace : bool
            Flag indicating whether boundary creates nullspace.
        """
        self.PCType = 'schur'
        self.L = L
        self._initializeIS(prefix)
        self.bdyNullSpace = bdyNullSpace
        self.pc.setFromOptions()

    def setUp(self,global_ksp):
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

    def _initializeIS(self,prefix):
        """ Sets the index set (IP) for the pressure and velocity 
        
        Parameters
        ----------
        prefix : str
            Prefix identifier for command line PETSc options.
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
        self.pc.setType('fieldsplit')
        self.isp = p4pyPETSc.IS()
        self.isp.createGeneral(self.pressureDOF,comm=p4pyPETSc.COMM_WORLD)
        self.isv = p4pyPETSc.IS()
        self.isv.createGeneral(self.velocityDOF,comm=p4pyPETSc.COMM_WORLD)
        self.pc.setFieldSplitIS(('velocity',self.isv),('pressure',self.isp))

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


class NavierStokesSchur(SchurPrecon):
    """ Schur complement preconditioners for Navier-Stokes problems.

    This class is derived from SchurPrecond and serves as the base
    class for all NavierStokes preconditioners which use the Schur complement
    method.    
    """
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        SchurPrecon.__init__(self,L,prefix)
        self.operator_constructor = SchurOperatorConstructor(self,
                                                             pde_type='navier_stokes')

class Schur_Qp(SchurPrecon) :
    """ 
    A Navier-Stokes (or Stokes) preconditioner which uses the 
    viscosity scaled pressure mass matrix. 
    """
    def __init__(self,L,prefix=None,bdyNullSpace=False):
        """
        Initializes the pressure mass matrix class.

        Parameters
        ---------
        L : petsc4py matrix
            Defines the problem's operator.
        """
        SchurPrecon.__init__(self,L,prefix,bdyNullSpace)
        self.operator_constructor = SchurOperatorConstructor(self)

    def setUp(self,global_ksp):
        """ Attaches the pressure mass matrix to PETSc KSP preconditioner.

        Parameters
        ----------
        global_ksp : PETSc KSP object
        """
        # Create the pressure mass matrix and scaxle by the viscosity.
        self.Qp = self.operator_constructor.getQp()
        self.Qp.scale(1./self.L.pde.coefficients.nu)
        L_sizes = self.Qp.size
        L_range = self.Qp.owner_range

        # Setup a PETSc shell for the inverse Qp operator
        self.QpInv_shell = p4pyPETSc.Mat().create()
        self.QpInv_shell.setSizes(L_sizes)
        self.QpInv_shell.setType('python')
        self.matcontext_inv = MatrixInvShell(self.Qp)
        self.QpInv_shell.setPythonContext(self.matcontext_inv)
        self.QpInv_shell.setUp()

        # Set PETSc Schur operator
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setType('python')
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setPythonContext(self.matcontext_inv)
        global_ksp.pc.getFieldSplitSubKSP()[1].pc.setUp()

        self._setSchurlog(global_ksp)
        if self.bdyNullSpace == True:
            self._setConstantPressureNullSpace(global_ksp)
        
class NavierStokes3D:
    def __init__(self,L,prefix=None):
        self.L = L
        self.PCD=False
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
                                            step=4,
                                            dtype="i")
            velocityDOF = []
            for start in range(1,4):
                velocityDOF.append(numpy.arange(start=L_range[0]+start,
                                                stop=L_range[0]+neqns,
                                                step=4,
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

    def setUp(self, global_ksp=None):
        try:
            if self.L.pde.pp_hasConstantNullSpace:
                if self.pc.getType() == 'fieldsplit':#we can't guarantee that PETSc options haven't changed the type
                    self.nsp = p4pyPETSc.NullSpace().create(constant=True,comm=p4pyPETSc.COMM_WORLD)
                    self.kspList = self.pc.getFieldSplitSubKSP()
                    self.kspList[1].setNullSpace(self.nsp)
        except:
            pass
        if self.PCD:
            #modify the mass term in the continuity equation so the p-p block is Q
            rowptr, colind, nzval = self.L.pde.jacobian.getCSRrepresentation()
            self.Qsys_rowptr = rowptr.copy()
            self.Qsys_colind = colind.copy()
            self.Qsys_nzval = nzval.copy()
            nr = rowptr.shape[0] - 1
            nc = nr
            self.Qsys =SparseMat(nr,nc,
                                 self.Qsys_nzval.shape[0],
                                 self.Qsys_nzval,
                                 self.Qsys_colind,
                                 self.Qsys_rowptr)
            self.L.pde.q[('dm',0,0)][:] = 1.0
            self.L.pde.getMassJacobian(self.Qsys)
            self.Qsys_petsc4py = self.L.duplicate()
            L_sizes = self.L.getSizes()
            Q_csr_rep_local = self.Qsys.getSubMatCSRrepresentation(0,L_sizes[0][0])
            self.Qsys_petsc4py.setValuesLocalCSR(Q_csr_rep_local[0],
                                                 Q_csr_rep_local[1],
                                                 Q_csr_rep_local[2],
                                                 p4pyPETSc.InsertMode.INSERT_VALUES)
            self.Qsys_petsc4py.assemblyBegin()
            self.Qsys_petsc4py.assemblyEnd()
            self.Qp = self.Qsys_petsc4py.getSubMatrix(self.isp,self.isp)
            from LinearAlgebraTools import _petsc_view
            _petsc_view(self.Qp, "Qp")#write to Qp.m
            #running Qp.m will load into matlab  sparse matrix
            #runnig Qpf = full(Mat_...) will get the full matrix
            #
            #modify the diffusion term in the mass equation so the p-p block is Ap
            rowptr, colind, nzval = self.L.pde.jacobian.getCSRrepresentation()
            self.Asys_rowptr = rowptr.copy()
            self.Asys_colind = colind.copy()
            self.Asys_nzval = nzval.copy()
            L_sizes = self.L.getSizes()
            nr = L_sizes[0][0]
            nc = L_sizes[1][0]
            self.Asys =SparseMat(nr,nc,
                                 self.Asys_nzval.shape[0],
                                 self.Asys_nzval,
                                 self.Asys_colind,
                                 self.Asys_rowptr)
            self.L.pde.q[('a',0,0)][:] = self.L.pde.q[('a',1,1)][:]
            self.L.pde.q[('df',0,0)][:] = 0.0
            self.L.pde.q[('df',0,1)][:] = 0.0
            self.L.pde.q[('df',0,2)][:] = 0.0
            self.L.pde.q[('df',0,3)][:] = 0.0
            self.L.pde.getSpatialJacobian(self.Asys)#notice switched to  Spatial
            self.Asys_petsc4py = self.L.duplicate()
            A_csr_rep_local = self.Asys.getSubMatCSRrepresentation(0,L_sizes[0][0])
            self.Asys_petsc4py.setValuesLocalCSR(A_csr_rep_local[0],
                                                 A_csr_rep_local[1],
                                                 A_csr_rep_local[2],
                                                 p4pyPETSc.InsertMode.INSERT_VALUES)
            self.Asys_petsc4py.assemblyBegin()
            self.Asys_petsc4py.assemblyEnd()
            self.Ap = self.Asys_petsc4py.getSubMatrix(self.isp,self.isp)
            _petsc_view(self.Ap,"Ap")#write to A.m
            #running A.m will load into matlab  sparse matrix
            #runnig Af = full(Mat_...) will get the full matrix
            #Af(1:np,1:np) whould be Ap, the pressure diffusion matrix
            #modify the diffusion term in the mass equation so the p-p block is Ap
            #modify the diffusion term in the mass equation so the p-p block is Fp
            rowptr, colind, nzval = self.L.pde.jacobian.getCSRrepresentation()
            self.Fsys_rowptr = rowptr.copy()
            self.Fsys_colind = colind.copy()
            self.Fsys_nzval = nzval.copy()
            L_sizes = self.L.getSizes()
            nr = L_sizes[0][0]
            nc = L_sizes[1][0]
            self.Fsys =SparseMat(nr,nc,
                                 self.Fsys_nzval.shape[0],
                                 self.Fsys_nzval,
                                 self.Fsys_colind,
                                 self.Fsys_rowptr)
            self.L.pde.q[('a',0,0)][:] = self.L.pde.q[('a',1,1)][:]
            self.L.pde.q[('df',0,0)][...,0] = self.L.pde.q[('u',1)]
            self.L.pde.q[('df',0,0)][...,1] = self.L.pde.q[('u',2)]
            self.L.pde.q[('df',0,0)][...,2] = self.L.pde.q[('u',3)]
            self.L.pde.getSpatialJacobian(self.Fsys)#notice, switched  to spatial
            self.Fsys_petsc4py = self.L.duplicate()
            F_csr_rep_local = self.Fsys.getSubMatCSRrepresentation(0,L_sizes[0][0])
            self.Fsys_petsc4py.setValuesLocalCSR(F_csr_rep_local[0],
                                                 F_csr_rep_local[1],
                                                 F_csr_rep_local[2],
                                                 p4pyPETSc.InsertMode.INSERT_VALUES)
            self.Fsys_petsc4py.assemblyBegin()
            self.Fsys_petsc4py.assemblyEnd()
            self.Fp = self.Fsys_petsc4py.getSubMatrix(self.isp,self.isp)
            _petsc_view(self.Fp, "Fp")#write to F.m
            #running F.m will load into matlab  sparse matrix
            #runnig Ff = full(Mat_...) will get the full matrix
            #Ff(1:np,1:np) whould be Fp, the pressure convection-diffusion matrix
            #
            #now zero all the dummy coefficents
            #
            self.L.pde.q[('dm',0,0)][:] = 0.0
            self.L.pde.q[('df',0,0)][:] = 0.0
            self.L.pde.q[('a',0,0)][:] = 0.0
            #
            # I think the next step is to setup something that does
            #  p = ksp(Qp, Fp*ksp(Ap,b)) for some input b
            # where p represents the approximate solution of S p = b
SimpleNavierStokes3D = NavierStokes3D

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
    def setUp(self, global_ksp=None):
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
    def setUp(self, global_ksp=None):
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
    def setUp(self, global_ksp=None):
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
    def setUp(self, global_ksp=None):
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
"""
A function for setting up a multilevel linear solver.
"""
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
                                                      linearSolverLocalBlockSize = linearSolverLocalBlockSize))
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
        self.model = model

class OperatorConstructor_oneLevel(OperatorConstructor):
    """ 
    A class for building common discrete operators from one-level
    transport models.
    
    Arguments
    ---------
    OLT : :class:`proteus.Transport.OneLevelTransport`
        One level transport class from which operator construction 
        will be based.
    """
    def __init__(self,OLT):
        OperatorConstructor.__init__(self,OLT)
        self._initializeOperatorConstruction()
        self.massOperatorAttached = False

    def attachMassOperator(self, rho=1., recalculate=False):
        """Builds a discrete Mass Operator matrix.
        
        Arguments
        ---------
        rho : float
            Optional scaling factor for mass inner-products.
        recalculate : bool
            Flag indicating mass operator should be recalculated.
        """
        if self.massOperatorAttached is True and recalculate is False:
            return
        self._mass_val = self.model.nzval.copy()
        self._mass_val.fill(0.)
        self.MassOperator = SparseMat(self.model.nFreeVDOF_global,
                                      self.model.nFreeVDOF_global,
                                      self.model.nnz,
                                      self._mass_val,
                                      self.model.colind,
                                      self.model.rowptr)
        _nd = self.model.coefficients.nd
        self.MassOperatorCoeff = TransportCoefficients.DiscreteMassMatrix(rho=rho, nd=_nd)
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
        self._elementQuadrature = self.model._elementQuadrature
        self._elementBoundaryQuadrature = self.model._elementBoundaryQuadrature

    def _attachJacobianInfo(self,Q):
        """ This helper function attaches quadrature data related to 'J'

        Arguments
        ---------
        Q : dict
            A dictionary to store values at quadrature points.
        """
        scalar_quad = StorageSet(shape=(self.model.mesh.nElements_global,
                                        self.model.nQuadraturePoints_element) ) 
        tensor_quad = StorageSet(shape={})

        tensor_quad |= set(['J',
                            'inverse(J)'])
        scalar_quad |= set(['det(J)',
                            'abs(det(J))'])

        for k in tensor_quad:
            Q[k] = numpy.zeros((self.model.mesh.nElements_global,
                                self.model.nQuadraturePoints_element,
                                self.model.nSpace_global,
                                self.model.nSpace_global),
                                'd')

        scalar_quad.allocate(Q)
        
        self.model.u[0].femSpace.elementMaps.getJacobianValues(self.model.elementQuadraturePoints,
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

        test_shape_quad |= set([('w',ci) for ci in range(self.model.nc)])
        test_shapeGradient_quad |= set([('grad(w)',ci) for ci in range(self.model.nc)])
        
        for k in test_shape_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nDOF_test_element[k[-1]]),
                'd')

        for k in test_shapeGradient_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nDOF_test_element[k[-1]],
                 self.model.nSpace_global),
                'd')

        for ci in range(self.model.nc):
            if Q.has_key(('w',ci)):
                self.model.testSpace[ci].getBasisValues(self.model.elementQuadraturePoints,
                                                      Q[('w',ci)])
            if Q.has_key(('grad(w)',ci)):
                self.model.testSpace[ci].getBasisGradientValues(self.model.elementQuadraturePoints,
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
        
        trial_shape_quad |= set([('v',ci) for ci in range(self.model.nc)])
        trial_shapeGrad_quad |= set([('grad(v)',ci) for ci in range(self.model.nc)])

        for k in trial_shape_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nDOF_test_element[k[-1]]),
                'd')

        for k in trial_shapeGrad_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nDOF_test_element[k[-1]],
                 self.model.nSpace_global),
                'd')
                                            
        for ci in range(self.model.nc):
            if Q.has_key(('v',ci)):
                self.model.testSpace[ci].getBasisValues(self.model.elementQuadraturePoints,
                                                      Q[('v',ci)])
            if Q.has_key(('grad(v)',ci)):
                self.model.testSpace[ci].getBasisGradientValues(self.model.elementQuadraturePoints,
                                                              Q[('inverse(J)')],
                                                              Q[('grad(v)',ci)])

    def _allocateMassOperatorQStorageSpace(self,Q):
        """ Allocate space for mass operator values. """
        test_shape_quad = StorageSet(shape={})
        trial_shape_quad = StorageSet(shape={})
        trial_shape_X_test_shape_quad = StorageSet(shape={})
        tensor_quad = StorageSet(shape={})
        # TODO - ARB : I don't think the 3 is necessary here...It created a
        # confusing bug in the 2-phase problem...Need to investigate.
        scalar_quad = StorageSet(shape=(self.model.mesh.nElements_global,
                                        self.model.nQuadraturePoints_element,
                                        3))
  
        scalar_quad |= set([('u',ci) for ci in range(self.model.nc)])
        scalar_quad |= set([('m',ci) for ci in self.MassOperatorCoeff.mass.keys()])

        test_shape_quad |= set([('w*dV_m',ci) for ci in self.MassOperatorCoeff.mass.keys()])

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            trial_shape_X_test_shape_quad |= set([('vXw*dV_m',cj,ci) for cj in cjDict.keys()])

        for ci,cjDict in self.MassOperatorCoeff.mass.iteritems():
            scalar_quad |= set([('dm',ci,cj) for cj in cjDict.keys()])

        for k in tensor_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nSpace_global,
                 self.model.nSpace_global),
                'd')

        for k in test_shape_quad:
            Q[k] = numpy.zeros(
                (self.model.mesh.nElements_global,
                 self.model.nQuadraturePoints_element,
                 self.model.nDOF_test_element[k[-1]]),
                'd')

        for k in trial_shape_X_test_shape_quad:
            Q[k] = numpy.zeros((self.model.mesh.nElements_global,
                                self.model.nQuadraturePoints_element,
                                self.model.nDOF_trial_element[k[1]],
                                self.model.nDOF_test_element[k[2]]),'d')

        scalar_quad.allocate(Q)

    def _calculateMassOperatorQ(self,Q):
        """ Calculate values for mass operator. """
        elementQuadratureDict = {}

        for ci in self.MassOperatorCoeff.mass.keys():
            elementQuadratureDict[('m',ci)] = self._elementQuadrature

        (elementQuadraturePoints,elementQuadratureWeights,
         elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        for ci in range(self.model.nc):
            if Q.has_key(('w*dV_m',ci)):
                cfemIntegrals.calculateWeightedShape(elementQuadratureWeights[('m',ci)],
                                                     self._operatorQ['abs(det(J))'],
                                                     self._operatorQ[('w',ci)],
                                                     Q[('w*dV_m',ci)])

        for ci in zip(range(self.model.nc),range(self.model.nc)):
                cfemIntegrals.calculateShape_X_weightedShape(self._operatorQ[('v',ci[1])],
                                                             Q[('w*dV_m',ci[0])],
                                                             Q[('vXw*dV_m',ci[1],ci[0])])

    def _allocateMatrixSpace(self,coeff,matrixDict):
        """ Allocate space for Operator Matrix """
        for ci in range(self.model.nc):
            matrixDict[ci] = {}
            for cj in range(self.model.nc):
                if cj in coeff.stencil[ci]:
                    matrixDict[ci][cj] = numpy.zeros(
                        (self.model.mesh.nElements_global,
                         self.model.nDOF_test_element[ci],
                         self.model.nDOF_trial_element[cj]),
                        'd')

    def _createOperator(self,coeff,matrixDict,A):
        """ Takes the matrix dictionary and creates a CSR matrix """
        for ci in range(self.model.nc):
            for cj in coeff.stencil[ci]:
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.model.l2g[ci]['nFreeDOF'],
                                                                          self.model.l2g[ci]['freeLocal'],
                                                                          self.model.l2g[cj]['nFreeDOF'],
                                                                          self.model.l2g[cj]['freeLocal'],
                                                                          self.model.csrRowIndeces[(ci,cj)],
                                                                          self.model.csrColumnOffsets[(ci,cj)],
                                                                          matrixDict[ci][cj],
                                                                          A)
