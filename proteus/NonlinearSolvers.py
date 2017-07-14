"""
A hierarchy of classes for nonlinear algebraic system solvers.

.. inheritance-diagram:: proteus.NonlinearSolvers
   :parts: 1
"""

import numpy
import numpy as np
from math import *
import math #to disambiguate math.log and log
from LinearAlgebraTools import *
from .Profiling import *

#mwf hack for Eikonal equation solvers
import FemTools
from UnstructuredFMMandFSWsolvers import FMMEikonalSolver
from UnstructuredFMMandFSWsolvers import FSWEikonalSolver
import csmoothers

class NonlinearEquation:
    """
    The base class for nonlinear equations.
    """

    def __init__(self,dim=0,dim_proc=None):
        self.dim=dim
        if dim_proc is None:
            self.dim_proc=self.dim
        else:
            self.dim_proc = dim_proc
        #mwf decide if we can keep solver statistics here
        self.nonlinear_function_evaluations = 0
        self.nonlinear_function_jacobian_evaluations = 0

    def getResidual(u,r):
        """Evaluate the residual r = F(u)"""
        pass

    def getJacobian(jacobian,usePicard=False):
        """"""
        pass

    def resetNonlinearFunctionStatistics(self):
        self.nonlinear_function_evaluations = 0
        self.nonlinear_function_jacobian_evaluations = 0


class NonlinearSolver:
    """
    The base class for nonlinear solvers.
    """

    def __init__(self,
                 F,J=None,du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 unorm = None,
                 tol_du=0.33):
        ## @var self.F
        #NonlinearEquation
        self.F = F
        self.J = J
        if du is None:
            self.du = Vec(F.dim)
        else:
            self.du = du
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
        self.duReductionFactor=100.0
        self.rReductionFactor_avg=0.0
        self.duReductionFactor_avg=100.0
        self.rReductionOrder=0.0
        self.rReductionOrder_avg=0.0
        self.duReductionOrder=0.0
        self.duReductionOrder_avg=0.0
        self.s=100.0
        self.ratio_r_current = 1.0
        self.ratio_r_solve = 1.0
        self.ratio_du_solve = 1.0
        self.last_log_ratio_r = 1.0
        self.last_log_ratior_du = 1.0
        #mwf begin hacks for conv. rate
        self.gustafsson_alpha = -12345.0
        self.gustafsson_norm_du_last = -12345.0
        #mwf end hacks for conv. rate
        self.W = 1.0/float(self.F.dim)
        self.kappa_current = 0.0 #condition number
        self.kappa_max = 0.0
        self.norm_2_J_current = 0.0
        self.norm_2_dJ_current = 0.0
        self.betaK_current  = 0.0  #Norm of J=F'
        self.etaK_current   = 0.0  #norm of du where F' du = -F
        self.betaK_0  = 0.0  #Norm of J=F'
        self.etaK_0   = 0.0  #norm of du where F' du = -F
        self.betaK_1  = 0.0  #Norm of J=F'
        self.etaK_1   = 0.0  #norm of du where F' du = -F
        self.gammaK_current = 0.0  #Lipschitz constant of J(u) approximated as ||F'(u+du)+F'(u)||/||du||
        self.gammaK_max=0.0
        self.norm_r_hist=[]
        self.norm_du_hist=[]
        self.convergenceTest = convergenceTest
        self.computeRates = computeRates
        self.printInfo = printInfo
        self.norm_function = norm
        if unorm is not None:
            self.unorm_function = unorm
        else:
            self.unorm_function = self.norm_function
        self.tol_du = tol_du
        self.infoString=''
        self.convergenceHistoryIsCorrupt=False
        self.r=None
        self.fullResidual=True
        self.lineSearch = True
        #need some information for parallel assembly options?
        self.par_fullOverlap = True #whether or not partitioning has overlap or not
        self.linearSolverFailed = False
        self.failedFlag = False

    def norm(self,u):
        return self.norm_function(u[self.F.owned_local])

    def unorm(self,u):
        return self.unorm_function(u[self.F.owned_local])

    def fullNewtonOff(self):
        self.fullNewton=False

    def fullNewtonOn(self):
        self.fullNewton=True

    def fullResidualOff(self):
        self.fullResidual=False

    def fullResidualOn(self):
        self.fullResidual=True

    def computeResidual(self,u,r,b):
        if self.fullResidual:
            self.F.getResidual(u,r)
            if b is not None:
                r-=b
        else:
            if type(self.J).__name__ == 'ndarray':
                r[:] = numpy.dot(u,self.J)
            elif type(self.J).__name__ == 'SparseMatrix':
                self.J.matvec(u,r)
            if b is not None:
                r-=b

    def solveInitialize(self,u,r,b):
        if r is None:
            if self.r is None:
                self.r = Vec(self.F.dim)
            r=self.r
        else:
            self.r=r
        self.computeResidual(u,r,b)
        self.its = 0
        self.norm_r0 = self.norm(r)
        self.norm_r = self.norm_r0
        self.ratio_r_solve = 1.0
        self.ratio_du_solve = 1.0
        self.last_log_ratio_r = 1.0
        self.last_log_ratior_du = 1.0
        #self.convergenceHistoryIsCorrupt=False
        self.convergingIts = 0
        #mwf begin hack for conv. rate
        self.gustafsson_alpha = -12345.0
        self.gustafsson_norm_du_last = -12345.0
        #mwf end hack for conv. rate
        return r

    def computeConvergenceRates(self):
        if self.convergenceHistoryIsCorrupt:
            return
        else:
            #mwf begin hack for conv. rate
            #equation (5) in Gustafsson_Soderlind_97
            #mwf debug
            #import pdb
            #pdb.set_trace()
            if self.gustafsson_norm_du_last >= 0.0:
                tmp = self.norm_du / (self.gustafsson_norm_du_last + 1.0e-16)
                self.gustafsson_alpha = max(self.gustafsson_alpha,tmp)
            #
            if self.its > 0:
                self.gustafsson_norm_du_last = self.norm_du
            #mwf end hack for conv. rate
            if self.convergingIts > 0:
                if self.norm_r < self.lastNorm_r:
                    self.ratio_r_current = self.norm_r/self.lastNorm_r
                else:
                    logEvent("residual increase %s" % self.norm_r)
                    self.convergingIts=0
                    self.ratio_r_solve = 1.0
                    self.ratio_du_solve = 1.0
                    self.last_log_ratio_r = 1.0
                    self.last_log_ratior_du = 1.0
                    return
                if self.ratio_r_current > 1.0e-100:
                    log_ratio_r_current = math.log(self.ratio_r_current)
                else:
                    logEvent("log(ratio_r) too small ratio_r = %12.5e" % self.ratio_r_current)
                    self.convergingIts=0
                    self.ratio_r_solve = 1.0
                    self.ratio_du_solve = 1.0
                    self.last_log_ratio_r = 1.0
                    self.last_log_ratior_du = 1.0
                    return
                self.ratio_r_solve *= self.ratio_r_current
                self.rReductionFactor = pow(self.ratio_r_solve,1.0/self.convergingIts)
                if self.convergingIts > 1:
                    self.rReductionOrder = log_ratio_r_current/ \
                                           self.last_log_ratio_r
                    if self.norm_du < self.lastNorm_du:
                        ratio_du_current = self.norm_du/self.lastNorm_du
                    else:
                        logEvent("du increase norm(du_last)=%12.5e, norm(du)=%12.5e, its=%d, convergingIts=%d" % (self.lastNorm_du,self.norm_du,self.its,self.convergingIts))
                        self.convergingIts=0
                        self.ratio_r_solve = 1.0
                        self.ratio_du_solve = 1.0
                        self.last_log_ratio_r = 1.0
                        self.last_log_ratior_du = 1.0
                        return
                    if ratio_du_current > 1.0e-100:
                        log_ratio_du_current = math.log(ratio_du_current)
                    else:
                        logEvent("log(du ratio) too small to calculate ratio_du=%12.5e" % ratio_du_current)
                        self.convergingIts=0
                        self.ratio_r_solve = 1.0
                        self.ratio_du_solve = 1.0
                        self.last_log_ratio_r = 1.0
                        self.last_log_ratior_du = 1.0
                        return
                    self.ratio_du_solve *= ratio_du_current
                    self.duReductionFactor = pow(self.ratio_du_solve,
                                                 1.0/(self.convergingIts-1))
                    if self.duReductionFactor  < 1.0:
                        self.s = self.duReductionFactor/(1.0-self.duReductionFactor)
                    else:
                        self.s=100.0
                    if self.convergingIts > 2:
                        self.duReductionOrder = log_ratio_du_current/ \
                                                self.last_log_ratio_du
                    self.last_log_ratio_du = log_ratio_du_current
                self.last_log_ratio_r = log_ratio_r_current
                self.lastNorm_du = self.norm_du
            self.lastNorm_r = self.norm_r

    def converged(self,r):
        self.convergedFlag = False
        self.norm_r = self.norm(r)
        self.norm_du = self.unorm(self.du)
        if self.computeRates ==  True:
            self.computeConvergenceRates()
        if self.convergenceTest == 'its' or self.convergenceTest == 'rits':
            if self.its == self.maxIts:
                self.convergedFlag = True
        #print self.atol_r, self.rtol_r
        if self.convergenceTest == 'r' or self.convergenceTest == 'rits':
            if (self.its != 0 and
                self.norm_r < self.rtol_r*self.norm_r0 + self.atol_r):
                self.convergedFlag = True
        if self.convergenceTest == 'u':
            if (self.convergingIts != 0 and
                self.s * self.norm_du < self.tol_du):
                self.convergedFlag = True
        if self.convergedFlag == True and self.computeRates == True:
            self.computeAverages()
        if self.printInfo == True:
            print self.info()
        #print self.convergedFlag
        return self.convergedFlag

    def failed(self):
        self.failedFlag = False
        if self.linearSolverFailed == True:
            self.failedFlag = True
            return self.failedFlag
        if self.its == self.maxIts and self.convergenceTest in ['r','u']:
            self.solveCalls_failed +=1
            self.recordedIts_failed +=self.its
            self.failedFlag = True
            logEvent("   Newton it %d == maxIts FAILED convergenceTest = %s" % (self.its,self.convergenceTest))
        else:
            self.its+=1
            self.convergingIts+=1
        return self.failedFlag

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
        self.infoString =  "************Start Nonlinear Solver Info ************ \n"
        self.infoString += "its                   = %i \n" % self.its
        self.infoString += "converging its        = %i \n" % self.convergingIts
        self.infoString += "r reduction factor    = %12.5e\n" % self.rReductionFactor
        self.infoString += "du reduction factor   = %12.5e\n" % self.duReductionFactor
        self.infoString += "r reduction order     = %12.5e\n" % self.rReductionOrder
        self.infoString += "du reduction order    = %12.5e\n" % self.duReductionOrder
        self.infoString += "<r reduction factor>  = %12.5e\n" % self.rReductionFactor_avg
        self.infoString += "<du reduction factor> = %12.5e\n" % self.duReductionFactor_avg
        self.infoString += "<r reduction order>   = %12.5e\n" % self.rReductionOrder_avg
        self.infoString += "<du reduction order>  = %12.5e\n" % self.duReductionOrder_avg
        self.infoString += "total its             = %i \n" % self.recordedIts
        self.infoString += "solver calls          = %i \n" % self.solveCalls
        self.infoString += "failures              = %i \n" % self.solveCalls_failed
        self.infoString += "failed its            = %i \n" % self.recordedIts_failed
        self.infoString += "maxIts                = %i \n" % self.maxIts
        self.infoString += "convergenceTest       = %s \n" % self.convergenceTest
        self.infoString += "atol_r                = %12.5e \n" % self.atol_r
        self.infoString += "rtol_r                = %12.5e \n" % self.rtol_r
        self.infoString += "norm(r0)              = %12.5e \n" % self.norm_r0
        self.infoString += "norm(r)               = %12.5e \n" % self.norm_r
        if self.convergenceHistoryIsCorrupt:
            self.infoString += "CONVERGENCE HISTORY IS CORRUPT!!!\n"
        self.infoString += "************End Nonlinear Solver Info ************\n"
        return self.infoString


class Newton(NonlinearSolver):
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
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        self.EWtol=EWtol
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')

    def setLinearSolverTolerance(self,r):
        self.norm_r = self.norm(r)
        gamma  = 0.0001
        etaMax = 0.001
        if self.norm_r == 0.0:
            etaMin = 0.0001
        else:
            etaMin = 0.0001*(self.rtol_r*self.norm_r0 + self.atol_r)/self.norm_r
        logEvent("etaMin "+`etaMin`)
        if self.its > 1:
            etaA = gamma * self.norm_r**2/self.norm_r_last**2
            logEvent("etaA "+`etaA`)
            logEvent("gama*self.etaLast**2 "+ `gamma*self.etaLast**2`)
            if gamma*self.etaLast**2 < 0.1:
                etaC = min(etaMax,etaA)
            else:
                etaC = min(etaMax,max(etaA,gamma*self.etaLast**2))
        else:
            etaC = etaMax
        logEvent("etaC "+`etaC`)
        eta = min(etaMax,max(etaC,etaMin))
        self.etaLast = eta
        self.norm_r_last = self.norm_r
        self.linearSolver.setResTol(rtol=eta,atol=self.linearSolver.atol_r,maxIts=self.linearSolver.maxIts)
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """

        import Viewers
        memory()
        if self.linearSolver.computeEigenvalues:
            self.u0[:]=u
        r=self.solveInitialize(u,r,b)
        if par_u is not None:
            #allow linear solver to know what type of assembly to use
            self.linearSolver.par_fullOverlap = self.par_fullOverlap
            #no overlap
            if not self.par_fullOverlap:
                par_r.scatter_reverse_add()
            else:
                #no overlap or overlap (until we compute norms over only owned dof)
                par_r.scatter_forward_insert()

        self.norm_r0 = self.norm(r)
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(r) and
               not self.failed()):
            logEvent("  NumericalAnalytics NewtonIteration: %d, NewtonNorm: %12.5e"
                %(self.its-1, self.norm_r), level=1)            
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %g test=%s"
                % (self.its-1,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r)),self.convergenceTest),level=1)
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                if self.linearSolver.computeEigenvalues:
                    logEvent("Calculating eigenvalues of J^t J")
                    self.JLast[:]=self.J
                    self.J_t_J[:]=self.J
                    self.J_t_J *= numpy.transpose(self.J)
                    self.JLsolver.prepare()#eigenvalue calc happens in prepare
                    self.norm_2_J_current = sqrt(max(self.JLsolver.eigenvalues_r))
                    try:
                        self.norm_2_Jinv_current = 1.0/sqrt(min(self.JLsolver.eigenvalues_r))
                    except:
                        logEvent("Norm of J_inv_current is singular to machine prection 1/sqrt("+`min(self.JLsolver.eigenvalues_r)`+")")
                        self.norm_2_Jinv_current = np.inf
                    self.kappa_current = self.norm_2_J_current*self.norm_2_Jinv_current
                    self.betaK_current = self.norm_2_Jinv_current
                self.linearSolver.prepare(b=r)
            self.du[:]=0.0
            if not self.directSolver:
                if self.EWtol:
                    self.setLinearSolverTolerance(r)
            if not self.linearSolverFailed:
                self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                self.linearSolverFailed = self.linearSolver.failed()
            u-=self.du
            if par_u is not None:
                par_u.scatter_forward_insert()
            self.computeResidual(u,r,b)
            if par_r is not None:
                #no overlap
                if not self.par_fullOverlap:
                    par_r.scatter_reverse_add()
                else:
                    par_r.scatter_forward_insert()

            #print "global r",r
            if self.linearSolver.computeEigenvalues:
                #approximate Lipschitz constant of J
                logEvent("Calculating eigenvalues of dJ^t dJ")
                self.F.getJacobian(self.dJ_t_dJ)
                self.dJ_t_dJ-=self.JLast
                self.dJ_t_dJ *= numpy.transpose(self.dJ_t_dJ)
                self.dJLsolver.prepare()
                self.norm_2_dJ_current = sqrt(max(self.dJLsolver.eigenvalues_r))
                self.etaK_current = self.W*self.norm(self.du)
                self.gammaK_current = self.norm_2_dJ_current/self.etaK_current
                self.gammaK_max = max(self.gammaK_current,self.gammaK_max)
                self.norm_r_hist.append(self.W*self.norm(r))
                self.norm_du_hist.append(self.W*self.unorm(self.du))
                if self.its  == 1:
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
                    try:
                        print "r         ",(1.0+sqrt(1.0-2.0*self.betaK_current*self.etaK_current*self.gammaK_current))/(self.betaK_current*self.gammaK_current)
                    except:
                        pass
                if self.betaK_current*self.etaK_current*self.gammaK_max <= 0.5:
                    try:
                        print "r_max     ",(1.0+sqrt(1.0-2.0*self.betaK_current*self.etaK_current*self.gammaK_max))/(self.betaK_current*self.gammaK_max)
                    except:
                        pass
                print "lambda_max",max(self.linearSolver.eigenvalues_r)
                print "lambda_i_max",max(self.linearSolver.eigenvalues_i)
                print "norm_J",self.norm_2_J_current
                print "lambda_min",min(self.linearSolver.eigenvalues_r)
                print "lambda_i_min",min(self.linearSolver.eigenvalues_i)
            if self.lineSearch:
                norm_r_cur = self.norm(r)
                ls_its = 0
                    #print norm_r_cur,self.atol_r,self.rtol_r
    #                 while ( (norm_r_cur >= 0.99 * self.norm_r + self.atol_r) and
    #                         (ls_its < self.maxLSits) and
    #                         norm_r_cur/norm_r_last < 1.0):
                if norm_r_cur > self.rtol_r*self.norm_r0 + self.atol_r:#make sure hasn't converged already
                    while ( (norm_r_cur >= 0.9999 * self.norm_r) and
                            (ls_its < self.maxLSits)):
                        self.convergingIts = 0
                        ls_its +=1
                        self.du *= 0.5
                        u += self.du
                        if par_u is not None:
                            par_u.scatter_forward_insert()
                        self.computeResidual(u,r,b)
                        #no overlap
                        if par_r is not None:
                            #no overlap
                            if not self.par_fullOverlap:
                                par_r.scatter_reverse_add()
                            else:
                                par_r.scatter_forward_insert()
                        norm_r_cur = self.norm(r)
                        logEvent("""ls #%d norm_r_cur=%s atol=%g rtol=%g""" % (ls_its,
                                                                               norm_r_cur,
                                                                               self.atol_r,
                                                                               self.rtol_r))
                    if ls_its > 0:
                        logEvent("Linesearches = %i" % ls_its,level=3)
        else:
            if self.linearSolver.computeEigenvalues:
                try:
                    if self.betaK_0*self.etaK_0*self.gammaK_max <= 0.5:
                        print "r_{-,0}     ",(1.0+sqrt(1.0-2.0*self.betaK_0*self.etaK_0*self.gammaK_max))/(self.betaK_0*self.gammaK_max)
                    if self.betaK_1*self.etaK_1*self.gammaK_max <= 0.5 and self.its > 1:
                        print "r_{-,1}     ",(1.0+sqrt(1.0-2.0*self.betaK_1*self.etaK_1*self.gammaK_max))/(self.betaK_1*self.gammaK_max)
                except:
                    pass
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(r/self.norm_r_hist[0])))
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(du/self.norm_du_hist[0])))
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
            logEvent("  NumericalAnalytics NewtonIteration: %d, NewtonNorm: %12.5e"
                %(self.its-1, self.norm_r), level=1)                
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
                % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            logEvent(memory("Newton","Newton"),level=4)
            return self.failedFlag
        logEvent("  NumericalAnalytics NewtonIteration: %d, NewtonNorm: %12.5e"
            %(self.its-1, self.norm_r), level=1)        
        logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
            % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
        logEvent(memory("Newton","Newton"),level=4)

import deim_utils
class POD_Newton(Newton):
    """Newton's method on the reduced order system based on POD"""
    import deim_utils
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
                 maxLSits = 100,
                 use_deim=False):
        Newton.__init__(self,
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
        #setup reduced basis for solution
        self.DB = 11 #number of basis vectors for solution
        U = np.loadtxt('SVD_basis')
        self.U = U[:,0:self.DB]
        self.U_transpose = self.U.conj().T
        self.pod_J = np.zeros((self.DB,self.DB),'d')
        self.pod_linearSolver = LU(self.pod_J)
        self.J_rowptr,self.J_colind,self.J_nzval = self.J.getCSRrepresentation()
        self.pod_du = np.zeros(self.DB)
    def computeResidual(self,u,r,b):
        """
        Use DEIM algorithm to compute residual if use_deim is turned on

        Right now splits the evaluation into two 'temporal' (mass) and spatial piece

        As first step for DEIM still does full evaluation
        """
        if self.fullResidual:
            self.F.getResidual(u,r)
            if b is not None:
                r-=b
        else:
            if type(self.J).__name__ == 'ndarray':
                r[:] = numpy.dot(u,self.J)
            elif type(self.J).__name__ == 'SparseMatrix':
                self.J.matvec(u,r)
            if b is not None:
                r-=b

    def norm(self,u):
        return self.norm_function(u)
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """
        pod_u = np.dot(self.U_transpose,u)
        u[:] = np.dot(self.U,pod_u)
        r=self.solveInitialize(u,r,b)
        pod_r = np.dot(self.U_transpose,r)
        self.norm_r0 = self.norm(pod_r)
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(pod_r) and
               not self.failed()):
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %g test=%s"
                % (self.its-1,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r)),self.convergenceTest),level=1)
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.pod_J[:] = 0.0
                self.F.getJacobian(self.J)
                for i in range(self.DB):
                    for j in range(self.DB):
                        for k in range(self.F.dim):
                            for m in range(self.J_rowptr[k],self.J_rowptr[k+1]):
                                self.pod_J[i,j] += self.U_transpose[i,k]*self.J_nzval[m]*self.U[self.J_colind[m],j]
                #self.linearSolver.prepare(b=r)
                self.pod_linearSolver.prepare(b=pod_r)
            self.du[:]=0.0
            self.pod_du[:]=0.0
            if not self.linearSolverFailed:
                #self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                #self.linearSolverFailed = self.linearSolver.failed()
                self.pod_linearSolver.solve(u=self.pod_du,b=pod_r)
                self.linearSolverFailed = self.pod_linearSolver.failed()
            #pod_u-=np.dot(self.U_transpose,self.du)
            pod_u-=self.pod_du
            u[:] = np.dot(self.U,pod_u)
            #mostly for convergence norms
            self.du = np.dot(self.U,self.pod_du)
            self.computeResidual(u,r,b)
            pod_r[:] = np.dot(self.U_transpose,r)
            r[:] = np.dot(self.U,pod_r)
        else:
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
                % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            return self.failedFlag
        logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
            % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)

class POD_DEIM_Newton(Newton):
    """Newton's method on the reduced order system based on POD"""
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
                 maxLSits = 100,
                 use_deim=True):
        Newton.__init__(self,
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
        #setup reduced basis for solution
        self.DB = 43 #11 #number of basis vectors for solution
        U = np.loadtxt('SVD_basis')
        self.U = U[:,0:self.DB]
        self.U_transpose = self.U.conj().T
        #setup reduced basis for DEIM interpolants
        self.use_deim = use_deim
        self.DBf = None
        self.Uf  = None;
        self.rho_deim = None; self.Ut_Uf_PtUf_inv=None
        self.rs = None; self.rt = None
        calculate_deim_internally = True
        if self.use_deim:
            #mwf this calculates things in the code. Switch for debugging to just reading
            if calculate_deim_internally:
                Uf = np.loadtxt('Fs_SVD_basis')
                self.DBf = min(73,Uf.shape[1],self.F.dim)#debug
                self.Uf = Uf[:,0:self.DBf]
                #returns rho --> deim indices and deim 'projection' matrix
                #U(P^TU)^{-1}
                self.rho_deim,Uf_PtUf_inv = deim_utils.deim_alg(self.Uf,self.DBf)
            else:
                self.Uf = np.loadtxt('Fs_SVD_basis_truncated')
                self.DBf = self.Uf.shape[1]
                self.rho_deim = np.loadtxt('Fs_DEIM_indices_truncated',dtype='i')
                PtUf = self.Uf[self.rho_deim]
                assert PtUf.shape == (self.DBf,self.DBf)
                PtUfInv = np.linalg.inv(PtUf)
                Uf_PtUf_inv = np.dot(self.Uf,PtUfInv)
            #go ahead and left multiply projection matrix by solution basis
            #to get 'projection' from deim to coarse space
            self.Ut_Uf_PtUf_inv = np.dot(self.U_transpose,Uf_PtUf_inv)
        self.pod_J = np.zeros((self.DB,self.DB),'d')
        self.pod_Jt= np.zeros((self.DB,self.DB),'d')
        self.pod_Jtmp= np.zeros((self.DBf,self.DB),'d')
        self.pod_linearSolver = LU(self.pod_J)
        self.J_rowptr,self.J_colind,self.J_nzval = self.J.getCSRrepresentation()
        assert 'getSpatialJacobian' in dir(self.F)
        assert 'getMassJacobian' in dir(self.F)
        self.Js = self.F.initializeSpatialJacobian()
        self.Js_rowptr,self.Js_colind,self.Js_nzval = self.Js.getCSRrepresentation()
        self.Jt = self.F.initializeMassJacobian()
        self.Jt_rowptr,self.Jt_colind,self.Jt_nzval = self.Jt.getCSRrepresentation()

        self.pod_du = np.zeros(self.DB)
        self.skip_mass_jacobian_eval = True
        self.linear_reduced_mass_matrix=None
    def norm(self,u):
        return self.norm_function(u)
    #mwf add for DEIM
    def computeDEIMresiduals(self,u,rs,rt):
        """
        wrapper for computing residuals separately for DEIM
        """
        assert 'getSpatialResidual' in dir(self.F)
        assert 'getMassResidual' in dir(self.F)
        self.F.getSpatialResidual(u,rs)
        self.F.getMassResidual(u,rt)

    def solveInitialize(self,u,r,b):
        """
        if using deim modifies base initialization by
        splitting up residual evaluation into separate pieces
        interpolated by deim (right now just does 'mass' and 'space')

        NOT FINISHED
        """
        if r is None:
            if self.r is None:
                self.r = Vec(self.F.dim)
            r=self.r
        else:
            self.r=r
        self.computeResidual(u,r,b)
        if self.use_deim:
            if self.rs is None:
                self.rs = Vec(self.F.dim)
            if self.rt is None:
                self.rt = Vec(self.F.dim)
            self.computeDEIMresiduals(u,self.rs,self.rt)
        self.its = 0
        self.norm_r0 = self.norm(r)
        self.norm_r = self.norm_r0
        self.ratio_r_solve = 1.0
        self.ratio_du_solve = 1.0
        self.last_log_ratio_r = 1.0
        self.last_log_ratior_du = 1.0
        #self.convergenceHistoryIsCorrupt=False
        self.convergingIts = 0
        #mwf begin hack for conv. rate
        self.gustafsson_alpha = -12345.0
        self.gustafsson_norm_du_last = -12345.0
        #mwf end hack for conv. rate
        return r
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """
        if self.use_deim:
            return self.solveDEIM(u,r,b,par_u,par_r)
        pod_u = np.dot(self.U_transpose,u)
        u[:] = np.dot(self.U,pod_u)
        r=self.solveInitialize(u,r,b)
        pod_r = np.dot(self.U_transpose,r)
        self.norm_r0 = self.norm(pod_r)
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(pod_r) and
               not self.failed()):
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %g test=%s"
                % (self.its-1,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r)),self.convergenceTest),level=1)
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                self.pod_J[:] = 0.0
                for i in range(self.DB):
                    for j in range(self.DB):
                        for k in range(self.F.dim):
                            for m in range(self.J_rowptr[k],self.J_rowptr[k+1]):
                                self.pod_J[i,j] += self.U_transpose[i,k]*self.J_nzval[m]*self.U[self.J_colind[m],j]
                #self.linearSolver.prepare(b=r)
                self.pod_linearSolver.prepare(b=pod_r)
            self.du[:]=0.0
            self.pod_du[:]=0.0
            if not self.linearSolverFailed:
                #self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                #self.linearSolverFailed = self.linearSolver.failed()
                self.pod_linearSolver.solve(u=self.pod_du,b=pod_r)
                self.linearSolverFailed = self.pod_linearSolver.failed()
            #pod_u-=np.dot(self.U_transpose,self.du)
            pod_u-=self.pod_du
            u[:] = np.dot(self.U,pod_u)
            #mostly for norm calculations
            self.du = np.dot(self.U,self.pod_du)
            self.computeResidual(u,r,b)
            pod_r[:] = np.dot(self.U_transpose,r)
            r[:] = np.dot(self.U,pod_r)
        else:
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
                % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            return self.failedFlag
        logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
            % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
    def solveDEIM(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b

        using DEIM
        Start with brute force just testing things
        """
        assert self.use_deim
        pod_u = np.dot(self.U_transpose,u)
        u[:] = np.dot(self.U,pod_u)
        #evaluate fine grid residuals directly
        r=self.solveInitialize(u,r,b)
        #mwf debug
        tmp = r-self.rt-self.rs
        assert np.absolute(tmp).all() < 1.0e-12

        #r_deim = self.rt[self.rho_deim].copy()
        r_deim = self.rs[self.rho_deim]
        pod_rt = np.dot(self.U_transpose,self.rt)
        pod_r  = np.dot(self.Ut_Uf_PtUf_inv,r_deim)
        pod_r += pod_rt
        #mwf debug
        #import pdb
        #pdb.set_trace()
        assert not numpy.isnan(pod_r).any()
        #
        self.norm_r0 = self.norm(pod_r)
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(pod_r) and
               not self.failed()):
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %g test=%s"
                % (self.its-1,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r)),self.convergenceTest),level=1)
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                #go ahead and evaluate spatial grid on fine grid for now
                self.F.getSpatialJacobian(self.Js)
                assert not numpy.isnan(self.Js_nzval).any()
                self.F.getMassJacobian(self.Jt)
                assert not numpy.isnan(self.Jt_nzval).any()
                #mwf hack, speed up mass matrix calculation
                if self.skip_mass_jacobian_eval and self.linear_reduced_mass_matrix is None:
                    #now this holds U^T Jt U
                    self.pod_Jt[:] = 0.0
                    for i in range(self.DB):
                        for j in range(self.DB):
                            for k in range(self.F.dim):
                                for m in range(self.Jt_rowptr[k],self.Jt_rowptr[k+1]):
                                    self.pod_Jt[i,j] += self.U_transpose[i,k]*self.Jt_nzval[m]*self.U[self.Jt_colind[m],j]
                    #combined DEIM, coarse grid projection
                    #self.pod_Jt = np.dot(self.Ut_Uf_PtUf_inv,self.pod_Jtmp)
                    assert not numpy.isnan(self.pod_Jt).any()
                    self.linear_reduced_mass_matrix  = self.pod_Jt.copy()
                #have to scale by dt in general shouldn't affect constant mass matrix
                self.Jt_nzval /= self.F.timeIntegration.dt
                #mwf debug
                self.F.getJacobian(self.J)
                tmp = self.Jt_nzval+self.Js_nzval-self.J_nzval
                assert numpy.absolute(tmp).all() < 1.0e-12
                #now this holds P^T J_s U
                self.pod_Jtmp[:] = 0.0
                for i in range(self.DBf):
                    deim_i = self.rho_deim[i]
                    for j in range(self.DB):
                        for m in range(self.Js_rowptr[deim_i],self.Js_rowptr[deim_i+1]):
                            self.pod_Jtmp[i,j] += self.Js_nzval[m]*self.U[self.Js_colind[m],j]
                #combined DEIM, coarse grid projection
                tmp = np.dot(self.Ut_Uf_PtUf_inv,self.pod_Jtmp)
                self.pod_J.flat[:] = tmp.flat[:]
                assert not numpy.isnan(self.pod_J).any()
                if not self.skip_mass_jacobian_eval:
                    #now this holds U^T Jt U
                    self.pod_Jt[:] = 0.0
                    for i in range(self.DB):
                        for j in range(self.DB):
                            for k in range(self.F.dim):
                                for m in range(self.Jt_rowptr[k],self.Jt_rowptr[k+1]):
                                    self.pod_Jt[i,j] += self.U_transpose[i,k]*self.Jt_nzval[m]*self.U[self.Jt_colind[m],j]
                    #combined DEIM, coarse grid projection
                    #self.pod_Jt = np.dot(self.Ut_Uf_PtUf_inv,self.pod_Jtmp)
                    assert not numpy.isnan(self.pod_Jt).any()
                else:
                    #mwf hack just copy over from precomputed mass matrix
                    self.pod_Jt.flat[:] = self.linear_reduced_mass_matrix.flat[:]
                    self.pod_Jt /= self.F.timeIntegration.dt
                self.pod_J += self.pod_Jt
                #self.linearSolver.prepare(b=r)
                self.pod_linearSolver.prepare(b=pod_r)
            self.du[:]=0.0
            self.pod_du[:]=0.0
            if not self.linearSolverFailed:
                #self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                #self.linearSolverFailed = self.linearSolver.failed()
                self.pod_linearSolver.solve(u=self.pod_du,b=pod_r)
                self.linearSolverFailed = self.pod_linearSolver.failed()
            assert not self.linearSolverFailed
            #pod_u-=np.dot(self.U_transpose,self.du)
            assert not numpy.isnan(self.pod_du).any()
            pod_u-=self.pod_du
            u[:] = np.dot(self.U,pod_u)
            #mostly for norm calculations
            self.du = np.dot(self.U,self.pod_du)
            #self.computeResidual(u,r,b)
            self.computeDEIMresiduals(u,self.rs,self.rt)
            #mwf debug
            self.F.getResidual(u,r)
            tmp = r-self.rt-self.rs
            assert np.absolute(tmp).all() < 1.0e-12
            #mwf debug
            #r_deim = self.rt[self.rho_deim].copy()
            r_deim = self.rs[self.rho_deim]
            pod_rt = np.dot(self.U_transpose,self.rt)
            pod_r = np.dot(self.Ut_Uf_PtUf_inv,r_deim)
            pod_r += pod_rt
            assert not numpy.isnan(pod_r).any()
            r[:] = np.dot(self.U,pod_r)
        else:
            logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
                % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            return self.failedFlag
        logEvent("   Newton it %d norm(r) = %12.5e  \t\t norm(r)/(rtol*norm(r0)+atol) = %12.5e"
            % (self.its,self.norm_r,(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)


class NewtonNS(NonlinearSolver):
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
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        #mwf turned back on self.lineSearch = False
        self.EWtol=EWtol
        #mwf added
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')

    def setLinearSolverTolerance(self,r):
        self.norm_r = self.norm(r)
        gamma  = 0.01
        etaMax = 0.01
        if self.norm_r == 0.0:
            etaMin = 0.01*self.atol_r
        else:
            etaMin = 0.01*(self.rtol_r*self.norm_r0 + self.atol_r)/self.norm_r
        if self.its > 1:
            etaA = gamma * self.norm_r**2/self.norm_r_last**2
            if self.its > 2:
                if gamma*self.etaLast**2 < 0.01:
                    etaC = min(etaMax,etaA)
                else:
                    etaC = min(etaMax,max(etaA,gamma*self.etaLast**2))
            else:
                etaC = min(etaMax,etaA)
        else:
            etaC = etaMax
        eta = min(etaMax,max(etaC,etaMin))
        self.etaLast = eta
        self.norm_r_last = self.norm_r
        self.linearSolver.setResTol(rtol=eta,atol=0.0)


    def converged(self,r):
        self.convergedFlag = False
        self.norm_r = self.norm(r)
        self.norm_cont_r = self.norm_function(r[:self.F.dim_proc/4])
        self.norm_mom_r  = self.norm_function(r[self.F.dim_proc/4:self.F.dim_proc])

        #self.norm_cont_r = self.norm(r[:r.shape[0]/4])
        #self.norm_mom_r  = self.norm(r[r.shape[0]/4:])
        self.norm_du= 1.0/float(self.its+2)
        if self.computeRates ==  True:
            self.computeConvergenceRates()
        if self.convergenceTest == 'its' or self.convergenceTest == 'rits':
            if self.its == self.maxIts:
                self.convergedFlag = True
        #print self.atol_r, self.rtol_r
        if self.convergenceTest == 'r' or self.convergenceTest == 'rits':
            if (self.its != 0 and
                self.norm_cont_r < self.atol_r and #cek enforce mass conservation using atol_r only
                self.norm_mom_r  < self.rtol_r*self.norm_mom_r0  + self.atol_r):
                self.convergedFlag = True

        if self.convergedFlag == True and self.computeRates == True:
            self.computeAverages()
        if self.printInfo == True:
            print self.info()
        #print self.convergedFlag
        return self.convergedFlag


    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """

        import Viewers
        memory()
        if self.linearSolver.computeEigenvalues:
            self.u0[:]=u
        r=self.solveInitialize(u,r,b)

        if par_u is not None:
            #allow linear solver to know what type of assembly to use
            self.linearSolver.par_fullOverlap = self.par_fullOverlap
            #no overlap
            if not self.par_fullOverlap:
                par_r.scatter_reverse_add()
            else:
                #no overlap or overlap (until we compute norms over only owned dof)
                par_r.scatter_forward_insert()


        self.norm_cont_r0 = self.norm_function(r[:self.F.dim_proc/4])
        self.norm_mom_r0  = self.norm_function(r[self.F.dim_proc/4:self.F.dim_proc])

        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(r) and
               not self.failed()):
            logEvent("   Newton it %d Mom.  norm(r) = %12.5e   tol = %12.5e" % (self.its-1,self.norm_mom_r,self.atol_r),level=1)
            logEvent("   Newton it %d Cont. norm(r) = %12.5e   tol = %12.5e" % (self.its-1,self.norm_cont_r,self.rtol_r*self.norm_mom_r0  + self.atol_r),level=1)

            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                logEvent("Start assembling jacobian",level=4)
                self.F.getJacobian(self.J)
                logEvent("Done  assembling jacobian",level=4)

                if self.linearSolver.computeEigenvalues:

                    logEvent("Performing eigen analyses",level=4)
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
            self.du[:]=0.0
            if not self.directSolver:
                if self.EWtol:
                    self.setLinearSolverTolerance(r)

            logEvent("Start linear solve",level=4)
            if not self.linearSolverFailed:
                self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                self.linearSolverFailed = self.linearSolver.failed()
            self.linearSolver.printPerformance()
            #print self.du
            #if par_du is not None:
            #    par_du.scatter_forward_insert()
            u-=self.du
            if par_u is not None:
                par_u.scatter_forward_insert()
            self.computeResidual(u,r,b)
            #no overlap
            #print "local r",r
            if par_r is not None:
                #no overlap
                if not self.par_fullOverlap:
                    par_r.scatter_reverse_add()
                else:
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
                self.norm_du_hist.append(self.W*self.unorm(self.du))
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
#                 while ( (norm_r_cur >= 0.99 * self.norm_r + self.atol_r) and
#                         (ls_its < self.maxLSits) and
#                         norm_r_cur/norm_r_last < 1.0):
                while ( (norm_r_cur >= 0.9999 * self.norm_r) and
                        (ls_its < self.maxLSits)):
                    self.convergingIts = 0
                    ls_its +=1
                    self.du *= 0.5
                    u += self.du
                    if par_u is not None:
                        par_u.scatter_forward_insert()
                    self.computeResidual(u,r,b)
                    #no overlap
                    if par_r is not None:
                        #no overlap
                        if not self.par_fullOverlap:
                            par_r.scatter_reverse_add()
                        else:
                            par_r.scatter_forward_insert()
                    norm_r_last = norm_r_cur
                    norm_r_cur = self.norm(r)
                    logEvent("""ls #%d norm_r_cur=%s atol=%g rtol=%g""" % (ls_its,
                                                                           norm_r_cur,
                                                                           self.atol_r,
                                                                           self.rtol_r))
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(r/self.norm_r_hist[0])))
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(du/self.norm_du_hist[0])))
                    Viewers.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                      Viewers.datFilename,
                                                                                                      Viewers.plotNumber,
                                                                                                      'log(du) history')
                    Viewers.cmdFile.write(cmd)
                    Viewers.viewerPipe.write(cmd)
                    Viewers.newPlot()
                    Viewers.newWindow()

        logEvent("   Final       Mom.  norm(r) = %12.5e   %12.5e" % (self.norm_mom_r,self.rtol_r*self.norm_mom_r0  + self.atol_r),level=1)
        logEvent("   Final       Cont. norm(r) = %12.5e   %12.5e" % (self.norm_cont_r,self.rtol_r*self.norm_mom_r0  + self.atol_r),level=1)

        logEvent(memory("NSNewton","NSNewton"),level=4)

class SSPRKNewton(Newton):
    """
    Version of Newton for SSPRK so doesn't refactor unnecessarily
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
        self.par_du = par_du
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        Newton.__init__(self,
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
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """

        import Viewers

        if self.linearSolver.computeEigenvalues:
            self.u0[:]=u
        r=self.solveInitialize(u,r,b)
        if par_u is not None:
            #no overlap
            #par_r.scatter_reverse_add()
            #no overlap or overlap (until we compute norms over only owned dof)
            par_r.scatter_forward_insert()
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(r) and
               not self.failed()):
            logEvent("SSPRKNewton it "+`self.its`+" norm(r) " + `self.norm_r`,level=3)
            if self.updateJacobian or self.fullNewton and not self.isFactored:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                #print numpy.transpose(self.J)
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
            if not self.linearSolverFailed:
                self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                self.linearSolverFailed = self.linearSolver.failed()
            #print self.du
            u-=self.du
            if par_u is not None:
                par_u.scatter_forward_insert()
            self.computeResidual(u,r,b)
            #no overlap
            #print "local r",r
            if par_r is not None:
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
                    if par_u is not None:
                        par_u.scatter_forward_insert()
                    self.computeResidual(u,r,b)
                    #no overlap
                    if par_r is not None:
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
            return self.failedFlag

    def resetFactorization(self,needToRefactor=True):
        self.isFactored = not needToRefactor

class PicardNewton(Newton):
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
        Newton.__init__(self,
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
        self.picardIts = 1
        self.picardTol = 1000.0
        self.usePicard = True

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        """
        Solve F(u) = b

        b -- right hand side
        u -- solution
        r -- F(u) - b
        """

        import Viewers

        if self.linearSolver.computeEigenvalues:
            self.u0[:]=u
        r=self.solveInitialize(u,r,b)
        if par_u is not None:
            #allow linear solver to know what type of assembly to use
            self.linearSolver.par_fullOverlap = self.par_fullOverlap
            #no overlap
            if not self.par_fullOverlap:
                par_r.scatter_reverse_add()
            else:
                #no overlap or overlap (until we compute norms over only owned dof)
                par_r.scatter_forward_insert()

        self.norm_r0 = self.norm(r)
        self.norm_r_hist = []
        self.norm_du_hist = []
        self.gammaK_max=0.0
        self.linearSolverFailed = False
        while (not self.converged(r) and
               not self.failed()):
            if self.maxIts>1:
                logEvent("   Newton it %d norm(r) = %12.5e  %12.5g \t\t norm(r)/(rtol*norm(r0)+atol) = %g"
                            % (self.its-1,self.norm_r,100*(self.norm_r/self.norm_r0),(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                if self.usePicard and (self.its < self.picardIts or self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r) > self.picardTol):
                    print "Picard iteration"
                    self.F.getJacobian(self.J,self.usePicard)
                else:
                    self.F.getJacobian(self.J)
                #mwf commented out print numpy.transpose(self.J)
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
            self.du[:]=0.0
            if not self.directSolver:
                if self.EWtol:
                    self.setLinearSolverTolerance(r)
            if not self.linearSolverFailed:
                self.linearSolver.solve(u=self.du,b=r,par_u=self.par_du,par_b=par_r)
                self.linearSolverFailed = self.linearSolver.failed()
            #print self.du
            u-=self.du
            if par_u is not None:
                par_u.scatter_forward_insert()
            self.computeResidual(u,r,b)
            #no overlap
            #print "local r",r
            if par_r is not None:
                #no overlap
                if not self.par_fullOverlap:
                    par_r.scatter_reverse_add()
                else:
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
                self.norm_du_hist.append(self.W*self.unorm(self.du))
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
                ls_its = 0
                #print norm_r_cur,self.atol_r,self.rtol_r
#                 while ( (norm_r_cur >= 0.99 * self.norm_r + self.atol_r) and
#                         (ls_its < self.maxLSits) and
#                         norm_r_cur/norm_r_last < 1.0):
                while ( (norm_r_cur >= 0.9999 * self.norm_r) and
                        (ls_its < self.maxLSits)):
                    self.convergingIts = 0
                    ls_its +=1
                    self.du *= 0.5
                    u += self.du
                    if par_u is not None:
                        par_u.scatter_forward_insert()
                    self.computeResidual(u,r,b)
                    #no overlap
                    if par_r is not None:
                        #no overlap
                        if not self.par_fullOverlap:
                            par_r.scatter_reverse_add()
                        else:
                            par_r.scatter_forward_insert()
                    norm_r_cur = self.norm(r)
                    logEvent("""ls #%d norm_r_cur=%s atol=%g rtol=%g""" % (ls_its,
                                                                           norm_r_cur,
                                                                           self.atol_r,
                                                                           self.rtol_r))
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(r/self.norm_r_hist[0])))
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
                        Viewers.datFile.write("%12.5e %12.5e \n" % (it,math.log(du/self.norm_du_hist[0])))
                    Viewers.datFile.write("\n \n")
                    cmd = "set term x11 %i; plot \'%s\' index %i with linespoints title \"%s\" \n" % (Viewers.windowNumber,
                                                                                                      Viewers.datFilename,
                                                                                                      Viewers.plotNumber,
                                                                                                      'log(du) history')
                    Viewers.cmdFile.write(cmd)
                    Viewers.viewerPipe.write(cmd)
                    Viewers.newPlot()
                    Viewers.newWindow()
            if self.maxIts>1:
                logEvent("   Newton it %d norm(r) = %12.5e  %12.5g \t\t norm(r)/(rtol*norm(r0)+atol) = %g"
                             % (self.its-1,self.norm_r,100*(self.norm_r/self.norm_r0),(self.norm_r/(self.rtol_r*self.norm_r0+self.atol_r))),level=1)
            return self.failedFlag

class NLJacobi(NonlinearSolver):
    """
    Nonlinear Jacobi iteration.
    """

    def __init__(self,
                 F,J,du,
                 weight=4.0/5.0,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True):
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.linearSolver = LinearSolver(J)#dummy
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.M=Vec(self.F.dim)
        self.w=weight
        self.node_order=numpy.arange(self.F.dim,dtype='i')
        self.node_order=numpy.arange(self.F.dim-1,-1,-1,dtype='i')

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        r=self.solveInitialize(u,r,b)
        while (not self.converged(r) and
               not self.failed()):
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                if type(self.J).__name__ == 'ndarray':
                    self.M = self.w/numpy.diagonal(self.J)
                elif type(self.J).__name__ == 'SparseMatrix':
                    csmoothers.jacobi_NR_prepare(self.J,self.w,1.0e-16,self.M)
            if type(self.J).__name__ == 'ndarray':
                self.du[:]=r
                self.du*=self.M
            elif type(self.J).__name__ == "SparseMatrix":
                csmoothers.jacobi_NR_solve(self.J,self.M,r,self.node_order,self.du)
            u -= self.du
            self.computeResidual(u,r,b)
        else:
            return self.failedFlag

class NLGaussSeidel(NonlinearSolver):
    """
    Nonlinear Gauss-Seidel.
    """

    def __init__(self,
                 connectionList,
                 F,J,du,
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
                 printInfo = True,
                 fullNewton=True):
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.linearSolver = LinearSolver(J)#dummy
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.w=weight
        self.connectionList=connectionList
        self.M=Vec(self.F.dim)
        self.node_order=numpy.arange(self.F.dim,dtype='i')
        self.sym = sym
        self.node_order=numpy.arange(self.F.dim,dtype='i')
        self.node_order=numpy.arange(self.F.dim-1,-1,-1,dtype='i')

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        r=self.solveInitialize(u,r,b)
        while (not self.converged(r) and
               not self.failed()):
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                if type(self.J).__name__ == 'ndarray':
                    self.M = self.w/numpy.diagonal(self.J)
                elif type(self.J).__name__ == 'SparseMatrix':
                    dtol = min(numpy.absolute(r))*1.0e-8
                    csmoothers.gauss_seidel_NR_prepare(self.J,self.w,dtol,self.M)
            if type(self.J).__name__ == 'ndarray':
                self.du[:]=0.0
                for i in range(self.F.dim):
                    rhat = r[i]
                    for j in self.connectionList[i]:
                        rhat -= self.J[j,i]*self.du[j]
                    self.du[i] = self.M[i]*rhat
                if self.sym == True:
                    u-= self.du
                    self.computeResidual(u,r,b)
                    self.du[:]=0.0
                    for i in range(self.n-1,-1,-1):
                        rhat = self.r[i]
                        for j in self.connectionList[i]:
                            rhat -= self.J[i,j]*self.du[j]
                    self.du[i] = self.M[i]*rhat
            elif type(self.J).__name__ == "SparseMatrix":
                csmoothers.gauss_seidel_NR_solve(self.J,self.M,r,self.node_order,self.du)
            u -= self.du
            self.computeResidual(u,r,b)
        else:
            return self.failedFlag

class NLStarILU(NonlinearSolver):
    """
    Nonlinear alternating Schwarz on node stars.
    """

    def __init__(self,
                 connectionList,
                 F,J,du,
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
                 printInfo = True,
                 fullNewton=True):
        NonlinearSolver.__init__(self,
                                 F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.linearSolver = LinearSolver(J)#dummy
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.w=weight
        self.sym=sym
        if type(self.J).__name__ == 'ndarray':
            self.connectionList=connectionList
            self.subdomainIndecesList=[]
            self.subdomainSizeList=[]
            self.subdomainJ=[]
            self.subdomainR=[]
            self.subdomainDU=[]
            self.subdomainSolvers=[]
            self.globalToSubdomain=[]
            for i in range(self.F.dim):
                self.subdomainIndecesList.append([])
                connectionList[i].sort()
                self.globalToSubdomain.append(dict([(j,J+1) for J,j in
                                                    enumerate(connectionList[i])]))
                self.globalToSubdomain[i][i]=0
                nSubdomain = len(connectionList[i])+1
                self.subdomainR.append(Vec(nSubdomain))
                self.subdomainDU.append(Vec(nSubdomain))
                self.subdomainSizeList.append(len(connectionList[i]))
                self.subdomainJ.append(Mat(nSubdomain,nSubdomain))
                for J,j in enumerate(connectionList[i]):
                    self.subdomainIndecesList[i].append(set(connectionList[i]) &
                                                        set(connectionList[j]))
                    self.subdomainIndecesList[i][J].update([i,j])
        elif type(self.J).__name__ == 'SparseMatrix':
            self.node_order=numpy.arange(self.F.dim-1,-1,-1,dtype='i')
            self.asmFactorObject = self.csmoothers.ASMFactor(self.J)

    def prepareSubdomains(self):
        if type(self.J).__name__ == 'ndarray':
            self.subdomainSolvers=[]
            for i in range(self.F.dim):
                self.subdomainJ[i][0,0] = self.J[i,i]
                for J,j in enumerate(self.connectionList[i]):
                    #first do row 0 (star center)
                    self.subdomainJ[i][J+1,0] = self.J[j,i]
                    #now do boundary rows
                    for k in self.subdomainIndecesList[i][J]:
                        K = self.globalToSubdomain[i][k]
                        self.subdomainJ[i][K,J+1]=self.J[k,j]
                self.subdomainSolvers.append(LU(self.subdomainJ[i]))
                self.subdomainSolvers[i].prepare()
        elif type(self.J).__name__ == 'SparseMatrix':
            self.csmoothers.asm_NR_prepare(self.J,self.asmFactorObject)

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        r=self.solveInitialize(u,r,b)
        while (not self.converged(r) and
               not self.failed()):
            if self.updateJacobian or self.fullNewton:
                self.updateJacobian = False
                self.F.getJacobian(self.J)
                self.prepareSubdomains()
            self.du[:]=0.0
            if type(self.J).__name__ == 'ndarray':
                for i in range(self.F.dim):
                    #load subdomain residual
                    self.subdomainR[i][0] = r[i] - self.J[i,i]*self.du[i]
                    for j in self.connectionList[i]:
                        self.subdomainR[i][0] -= self.J[j,i]*self.du[j]
                    for J,j in enumerate(self.connectionList[i]):
                        self.subdomainR[i][J+1]=r[j] - self.J[j,j]*self.du[j]
                        for k in self.connectionList[j]:
                            self.subdomainR[i][J+1] -= self.J[k,j]*self.du[k]
                    #solve
                    self.subdomainSolvers[i].solve(u=self.subdomainDU[i],
                                                   b=self.subdomainR[i])
                    #update du
                    self.subdomainDU[i]*=self.w
                    self.du[i]+=self.subdomainDU[i][0]
                    for J,j in enumerate(self.connectionList[i]):
                        self.du[j] += self.subdomainDU[i][J+1]
            elif type(self.J).__name__ == 'SparseMatrix':
                self.csmoothers.asm_NR_solve(self.J,self.w,self.asmFactorObject,self.node_order,r,self.du)
            u -= self.du
            self.computeResidual(u,r,b)
        else:
            return self.failedFlag

class FasTwoLevel(NonlinearSolver):
    """
    A generic nonlinear two-level solver based on the full approximation scheme. (FAS).
    """
    def __init__(self,
                 prolong,
                 restrict,
                 restrictSum,
                 coarseF,
                 preSmoother,
                 postSmoother,
                 coarseSolver,
                 F,J=None,du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True):
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)

        class dummyL:

            def __init__(self):
                self.shape = [F.dim,F.dim]

        self.linearSolver = LinearSolver(dummyL())#dummy
        self.prolong = prolong
        self.restrict = restrict
        self.restrictSum = restrictSum
        self.cF = coarseF
        self.preSmoother = preSmoother
        self.postSmoother = postSmoother
        self.coarseSolver = coarseSolver
        self.cb = Vec(prolong.shape[1])
        self.cr = Vec(prolong.shape[1])
        self.crFAS = Vec(prolong.shape[1])
        self.cdu = Vec(prolong.shape[1])
        self.cu = Vec(prolong.shape[1])
        self.cError = 'FAS-chord'#'linear'#'FAS'#'FAS-chord'#'linear'

    def fullNewtonOff(self):
        self.preSmoother.fullNewtonOff()
        self.postSmoother.fullNewtonOff()
        self.coarseSolver.fullNewtonOff()

    def fullNewtonOn(self):
        self.preSmoother.fullNewtonOn()
        self.postSmoother.fullNewtonOn()
        self.coarseSolver.fullNewtonOn()

    def fullResidualOff(self):
        self.preSmoother.fullResidualOff()
        self.postSmoother.fullResidualOff()
        self.coarseSolver.fullResidualOff()

    def fullResidualOn(self):
        self.preSmoother.fullResidualOn()
        self.postSmoother.fullResidualOn()
        self.coarseSolver.fullResidualOn()

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        r=self.solveInitialize(u,r,b)
        while (not self.converged(r) and
               not self.failed()):
            self.preSmoother.solve(u,r,b)
            self.restrict.matvec(r,self.cb)
            if self.cError == 'linear':
                self.coarseSolver.fullResidualOff()
                self.coarseSolver.fullNewtonOff()
                lsSave = self.coarseSolver.lineSearch
                self.coarseSolver.lineSearch = False
                self.coarseSolver.solve(u=self.cdu,r=self.crFAS,b=self.cb)
                self.coarseSolver.lineSearch = lsSave
                self.coarseSolver.fullNewtonOn()
                self.coarseSolver.fullResidualOn()
            elif self.cError == 'FAS-chord':
                self.crFAS[:]= -self.cb
                self.restrict.matvec(u,self.cu)
                #change restriction to
                #be weighted average for u
                for i in range(self.cu.shape[0]):
                    self.cu[i]/=self.restrictSum[i]
                self.cF.getResidual(self.cu,self.cr)
                self.cb+=self.cr
                self.cdu[:]=self.cu
                self.coarseSolver.fullNewtonOff()
                self.coarseSolver.solve(u=self.cdu,r=self.crFAS,b=self.cb)
                self.coarseSolver.fullNewtonOn()
                self.cdu-=self.cu
            else:
                self.crFAS[:]= -self.cb
                self.restrict.matvec(u,self.cu)
                #change restriction to
                #be weighted average for u
                for i in range(self.cu.shape[0]):
                    self.cu[i]/=self.restrictSum[i]
                self.cF.getResidual(self.cu,self.cr)
                self.cb+=self.cr
                self.cdu[:]=self.cu
                self.coarseSolver.solve(u=self.cdu,r=self.crFAS,b=self.cb)
                self.cdu-=self.cu
            self.prolong.matvec(self.cdu,self.du)
            u-=self.du
            self.postSmoother.solve(u,r,b)
        else:
            return self.failedFlag

class FAS:
    """
    A generic nonlinear multigrid W cycle using the Full Approximation Scheme (FAS).
    """

    def __init__(self,
                 prolongList,
                 restrictList,
                 restrictSumList,
                 FList,
                 preSmootherList,
                 postSmootherList,
                 coarseSolver,
                 mgItsList=[],
                 printInfo=True):
        self.solverList=[coarseSolver]
        for i in range(1,len(FList)):
            if mgItsList ==[]:
                mgItsList.append(1)
            self.solverList.append(
                FasTwoLevel(prolong = prolongList[i],
                            restrict = restrictList[i],
                            restrictSum = restrictSumList[i],
                            coarseF = FList[i-1],
                            preSmoother = preSmootherList[i],
                            postSmoother = postSmootherList[i],
                            coarseSolver = self.solverList[i-1],
                            F = FList[i],
                            maxIts = mgItsList[i],
                            convergenceTest = 'its',
                            printInfo=False))
        self.fasSolver = self.solverList[len(FList)-1]

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        self.fasSolver.solve(u,r,b)
        return self.fasSolver.failedFlag


class MultilevelNonlinearSolver:
    """
    A generic multilevel solver.
    """

    def __init__(self,
                 fList,
                 levelNonlinearSolverList,
                 computeRates=False,
                 printInfo=False):
        self.fList = fList
        self.solverList=levelNonlinearSolverList
        self.nLevels = len(self.solverList)
        for l in range(self.nLevels):
            self.solverList[l].computeRates = computeRates
            self.solverList[l].printInfo = printInfo

    def solveMultilevel(self,uList,rList,bList=None,par_uList=None,par_rList=None):
        if bList is None:
            bList = [None for r in rList]
        for l in range(self.nLevels):
            if par_uList is not None and len(par_uList) > 0:
                par_u=par_uList[l]
                par_r=par_rList[l]
            else:
                par_u=None
                par_r=None
            logEvent("  NumericalAnalytics Newton iteration for level " + `l`, level = 0)                
            self.solverList[l].solve(u = uList[l],
                                     r = rList[l],
                                     b = bList[l],
                                     par_u = par_u,
                                     par_r = par_r)
        return self.solverList[-1].failedFlag

    def updateJacobian(self):
        for l in range(self.nLevels):
            self.solverList[l].updateJacobian = True

    def info(self):
        self.infoString="********************Start Multilevel Nonlinear Solver Info*********************\n"
        for l in range(self.nLevels):
            self.infoString += "**************Start Level %i Info********************\n" % l
            self.infoString += self.solverList[l].info()
            self.infoString += "**************End Level %i Info********************\n" % l
        self.infoString+="********************End Multilevel Nonlinear Solver Info*********************\n"
        return self.infoString


class EikonalSolver:
    """
    Simple wrapper for special purpose Eikonal equation solvers on a single level.
    Current types allowed::

       FMMEikonalSolver
       FSWEikonalSolver

    """
#    TODO Feb 20
#      Debug change in use truncation approach for positive and negative solutions
#      Debug local Reconstruction
    def __init__(self,
                 levelSolverType,
                 F,
                 relativeTolerance = 0.0,
                 absoluteTolerance = 1.0e-8,
                 maxSolverIts      = 100,
                 frontTolerance    = 1.0e-4,
                 frontInitType     = 'magnitudeOnly',#'magnitudeOnly','frontIntersection'
                 bandTolerance     = -1.0e-2,#copy over everything from Eikonal solution whose mag. >= tol
                 eikonalVariable   = 0, #which variable in F is to be used for ic
                 localReconstruction = None,#'localPWL'
                 printInfo=False):
        self.linearSolver = None # unneeded
        self.solverType = levelSolverType
        assert self.solverType == FMMEikonalSolver or self.solverType == FSWEikonalSolver, "unrecognized solver type"
        self.F = F
        self.eikonalVariable = eikonalVariable
        #allow limited set of fem spaces for eikonal variable
        self.allowedEikonalSpaces = [FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                                     FemTools.DG_AffineLinearOnSimplexWithNodalBasis,
                                     FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis]
        self.eikonalVariableFemSpace = None
        for space in self.allowedEikonalSpaces:
            if isinstance(self.F.u[self.eikonalVariable].femSpace,space):
                self.eikonalVariableFemSpace = space
        assert self.eikonalVariableFemSpace is not None, "allowed spaces= %s" % self.allowedEikonalSpaces

        #for determining if a point is on front or not
        self.frontTolerance = frontTolerance

        #allow some local reconstruction of front rather than just freezing values
        #PWL option can move zero level set location, None should result in just accepting "known"
        #values. This should freeze 0 level set if criteria for tagging known cells is broad enough
        #
        self.localReconstruction = localReconstruction
        self.utmp = None
        if (self.localReconstruction == 'localPWL' and
            self.eikonalVariableFemSpace == FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
            self.utmp = numpy.zeros(self.F.u[self.eikonalVariable].dof.shape,'d')

        if (self.localReconstruction == 'localPWL' and
            self.eikonalVariableFemSpace != FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
            logEvent("""WARNING EikonalSolver localReconstruction = %s only allowed with P1C0 for now""" % self.localReconstruction)
            self.localReconstruction = None
        #if
        #do not overwrite values from input if magnitude <  bandTolerance
        self.bandTolerance = bandTolerance
        #iteration params not used for FMM
        self.relativeTolerance = relativeTolerance
        self.absoluteTolerance = absoluteTolerance
        self.maxSolverIts      = maxSolverIts
        #
        self.printInfo         = printInfo
        self.baseSolver = None
        #set eikonal solver to use node number global dof
        if self.solverType == FSWEikonalSolver:
            #need to add relative tol
            self.baseSolver = FSWEikonalSolver(F.mesh,F.mesh.elementNodesArray,F.nSpace_global,
                                               iterAtol=self.absoluteTolerance,
                                               iterRtol=self.relativeTolerance,
                                               maxIts =self.maxSolverIts,
                                               frontInitType=frontInitType)
        else:
            self.baseSolver = FMMEikonalSolver(F.mesh,F.mesh.elementNodesArray,F.nSpace_global,
                                               frontInitType=frontInitType)

        #for holding initial condition with plus and minus direction results
        self.phi0p  = numpy.zeros((self.F.mesh.nNodes_global,),'d')
        self.phi0m  = numpy.zeros((self.F.mesh.nNodes_global,),'d')
        #go ahead and use a C0 FemSpace for solution variable to facilitate switching between arbitrary input spaces
        #and expected C0 P1 rep
        self.phidof = numpy.zeros((self.F.mesh.nNodes_global,),'d')
        if self.eikonalVariableFemSpace == FemTools.C0_AffineLinearOnSimplexWithNodalBasis:
            self.C0P1space = self.F.u[self.eikonalVariable].femSpace
        else:
            self.C0P1space = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(self.F.mesh,self.F.nSpace_global)
            #compute shape function values at interpolation points for actual space to save
            #time in transfering between C0 representation
            refPoints = self.F.u[self.eikonalVariable].femSpace.referenceFiniteElement.interpolationConditions.quadraturePointArray
            assert refPoints.shape[0] == self.F.u[self.eikonalVariable].femSpace.referenceFiniteElement.localFunctionSpace.dim
            self.C0P1shapeValues = numpy.zeros((self.F.mesh.nElements_global,
                                                  self.F.u[self.eikonalVariable].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                                  self.C0P1space.referenceFiniteElement.localFunctionSpace.dim),
                                                 'd')
            self.eikonalVarInterpCondVals = numpy.zeros((self.F.mesh.nElements_global,
                                                           self.F.u[self.eikonalVariable].femSpace.referenceFiniteElement.localFunctionSpace.dim),
                                                          'd')

            self.C0P1space.getBasisValues(refPoints,self.C0P1shapeValues)

        #for compatibility with generic multilevel solver, these are ignored
        self.updateJacobian = False
        self.computeRates = False

    def solve(self,u,r,b=None):
        if (self.localReconstruction == 'localPWL' and
            self.eikonalVariableFemSpace == FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
            self.baseSolver.localPWLreconstruction(self.F.nSpace_global,self.F.mesh.cmesh,u,self.utmp,self.frontTolerance,0)
            u.flat[:] = self.utmp.flat[:]
            assert False, "Needs to be debgged more Feb 20"
        #allowed reconstruction only for C0P1 for now
        self.convertToC0P1Rep(u,self.phidof)
        self.phi0p.flat[:] = numpy.maximum(self.phidof, 0.0)
        #mwf Feb 24 insert abs, need to be consistent in FSW and FMM on assumption for input values
        self.phi0m.flat[:] = numpy.abs(numpy.minimum(self.phidof, 0.0))
        #mwf debug
        #print "entering EikonalSolve u=%s \n phi0p=%s \n phi0m=%s " % (u,self.phi0p,self.phi0m)
        #negative direction first
        failedFlagM = self.baseSolver.solve(self.phi0m,self.phidof,zeroTol=self.frontTolerance,verbose=0)
        #save for later
        self.phi0m.flat[:] = self.phidof.flat

        #positive direction
        failedFlagP = self.baseSolver.solve(self.phi0p,self.phidof,zeroTol=self.frontTolerance,verbose=0)

        #T = T^+ - T^-
        self.phidof.flat[:] -= self.phi0m.flat
        self.convertFromC0P1Rep(self.phidof,u)
        #compute residual to update models copy of u
        self.F.getResidual(u,r)
        #mwf add updateTimeHistory? or call from somewhere else?
        self.F.updateTimeHistory(self.F.timeIntegration.t,resetFromDOF=False)
        failedFlag = failedFlagP or failedFlagM
        return failedFlag

    def info(self):
        myinfo = """EikonalSolver type=%s
        atol=%s rtol=%s maxits=%s eikonalVariable=%s """ % (self.solverType,self.absoluteTolerance,self.relativeTolerance,
                                                            self.maxIts,self.eikonalVariable)
        return myinfo

    def convertToC0P1Rep(self,dofin,dofout,interpType='min'):
        """
        allow input FEM space to be something besides C0P1 and then convert to expected
        representation here
        """
        if self.eikonalVariableFemSpace == FemTools.C0_AffineLinearOnSimplexWithNodalBasis:
            dofout.flat[:] =  dofin.flat #decide if needs to be deep copy or not
        elif (self.eikonalVariableFemSpace == FemTools.DG_AffineLinearOnSimplexWithNodalBasis or
              self.eikonalVariableFemSpace == FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis):
            #take first nd+1 dofs which are associated with nodes
            if interpType == 'min':
                dofout.flat[:] = 1.2345e28
                for eN in range(self.F.mesh.nElements_global):
                    for iv in range(self.F.nSpace_global+1):
                        IG = self.F.u[self.eikonalVariable].femSpace.dofMap.l2g[eN,iv]
                        nN = self.F.mesh.elementNodesArray[eN,iv]
                        dofout[nN] = min(dofout[nN],dofin[IG])
                    #iv
                #eN
            elif interpType == 'max':
                dofout.flat[:] = -1.2345e28
                for eN in range(self.F.mesh.nElements_global):
                    for iv in range(self.F.nSpace_global+1):
                        IG = self.F.u[self.eikonalVariable].femSpace.dofMap.l2g[eN,iv]
                        nN = self.F.mesh.elementNodesArray[eN,iv]
                        dofout[nN] = max(dofout[nN],dofin[IG])
                    #iv
                #eN
            else: #average by default
                dofout.flat[:] = 0.0
                for eN in range(self.F.mesh.nElements_global):
                    for iv in range(self.F.nSpace_global+1):
                        IG = self.F.u[self.eikonalVariable].femSpace.dofMap.l2g[eN,iv]
                        nN = self.F.mesh.elementNodesArray[eN,iv]
                        dofout[nN] += dofin[IG]
                    #iv
                #eN
                for nN in range(self.F.mesh.nNodes_global):
                    dofout[nN] /= self.F.mesh.nElements_node[nN]
            #avg
        #DG P[1,2] --> C0 P1
        else:
            assert False, "space= %s not supported" % self.eikonalVariableFemSpace

    def convertFromC0P1Rep(self,dofin,dofout):
        """
        allow output FEM space to be something besides C0P1 and convert from C0P1 to expected
        representation here
        """

        #mwf debug
        #print """into convert from c0p1 rep dofin.shape=%s dofout.shape=%s """ %(dofin.shape,dofout.shape)
        if self.eikonalVariableFemSpace == FemTools.C0_AffineLinearOnSimplexWithNodalBasis:
            #could check to make sure l2g <==> elementNodesArray
            #dofout.flat[:] = dofin.flat #decide if need deep copy or not
            FMMEikonalSolver.cfmmfsw.copyOverOutsideBand(self.bandTolerance,dofin,dofout)
        elif self.eikonalVariableFemSpace == FemTools.DG_AffineLinearOnSimplexWithNodalBasis:
            for eN in range(self.F.mesh.nElements_global):
                for iv in range(self.F.nSpace_global+1):
                    IG = self.F.u[self.eikonalVariable].femSpace.dofMap.l2g[eN,iv]
                    nN = self.F.mesh.elementNodesArray[eN,iv]
                    if abs(dofout[nN]) > self.bandTolerance:
                        dofout[IG] = dofin[nN]
                #iv
            #eN
        #DG P1 --> C0 P1
        else:
            #assuming nodal interpolation conditions for now
            #compute C0P1 values at output space interpolation points and then project
            #output dofs from interpolation points
            print "WARNING, default Eikonal convertFromC0P1 rep doesn't respect band width"
            inFEM = FemTools.FiniteElementFunction(self.C0P1space,dof=dofin)
            inFEM.getValues(self.C0P1shapeValues,self.eikonalVarInterpCondVals)
            outFEM= FemTools.FiniteElementFunction(self.F.u[self.eikonalVariable].femSpace,dof=dofout)
            outFEM.projectFromInterpolationConditions(self.eikonalVarInterpCondVals)

class MultilevelEikonalSolver:
    """
    Attempt a wrapper for multilevel Eikonal equation solves with restricted interface
    """
    def __init__(self,
                 levelSolverList,
                 printInfo=False):
        self.solverList=levelSolverList
        self.nLevels = len(self.solverList)
        self.printInfo = printInfo

    def solveMultilevel(self,uList,rList,bList=None,par_uList=None,par_rList=None):
        if bList is None:
            bList = [None for r in rList]
        failedFlag = False
        for l in range(self.nLevels):
            #mwf debug
            #print "entering MLEik l=%s uList[l]= %s " % (l,uList[l])
            failedFlag = self.solverList[l].solve(uList[l],rList[l])
            #mwf debug
            #print "failedFlag = %s leaving MLEik l=%s uList[l]= %s " % (failedFlag,l,uList[l])

        return failedFlag
    def updateJacobian(self):
        pass
    def info(self):
        self.infoString="********************Start Multilevel Eikonal Solver Info*********************\n"
        for l in range(self.nLevels):
            self.infoString += "**************Start Level %i Info********************\n" % l
            #self.infoString += self.solverList[l].info()
            self.infoString += "**************End Level %i Info********************\n" % l
        self.infoString+="********************End Multilevel Eikonal Solver Info*********************\n"
        return self.infoString


class NLNI(MultilevelNonlinearSolver):
    """
    Nonlinear nested iteration.
    """

    def __init__(self,
                 fList=[],
                 solverList=[],
                 prolongList=[],
                 restrictList=[],
                 restrictSumList=[],
                 maxIts = None,
                 tolList=None,
                 atol=None,
                 computeRates=True,
                 printInfo=True):
        MultilevelNonlinearSolver.__init__(self,fList,solverList,computeRates,printInfo)
        self.prolongList = prolongList
        self.restrictList = restrictList
        self.restrictSumList = restrictSumList
        self.fineLevel = self.nLevels - 1
        self.levelDict={}
        self.uList=[]
        self.rList=[]
        self.bList=[]
        for l in range(self.fineLevel+1):
            n = solverList[l].F.dim
            self.levelDict[n] = l
            self.uList.append(Vec(n))
            self.rList.append(Vec(n))
            self.bList.append(Vec(n))
        self.levelDict[solverList[-1].F.dim]=self.fineLevel
        self.maxIts = maxIts
        self.tolList = tolList
        self.atol = atol
        self.printInfo = printInfo
        self.infoString=''

    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        #\todo this is broken because prolong isn't right
        currentMesh = self.levelDict[u.shape[0]]
        if currentMesh > 0:
            self.uList[currentMesh][:] = u
            self.rList[currentMesh][:] = r
            self.bList[currentMesh][:] = b
            for l in range(currentMesh,1,-1):
                self.restrictList[l].matvec(self.uList[l],self.uList[l-1])
                if b is not None:
                    self.restrictList[l].matvec(self.bList[l],self.bList[l-1])
                for i in range(self.uList[l-1].shape[0]):
                    self.uList[l-1][i]/=self.restrictSumList[l][i]
            for l in range(currentMesh):
                if self.tolList is not None:
                    self.switchToResidualConvergence(self.solverList[l],
                                                     self.tolList[l])
                self.solverList[l].solve(u=self.uList[l],r=self.rList[l],b=self.bList[l],par_u=self.par_uList[l],par_r=self.par_rList[l])
                if self.tolList is not None:
                    self.revertToFixedIteration(self.solverList[l])
            if l < currentMesh -1:
                self.prolongList[l+1].matvec(self.uList[l],self.uList[l+1])
            else:
                self.prolongList[l+1].matvec(self.uList[l],u)
        if self.tolList is not None:
            self.switchToResidualConvergence(self.solverList[currentMesh],
                                             self.tolList[currentMesh])
        self.solverList[currentMesh].solve(u,r,b)
        if self.tolList is not None:
            self.revertToFixedIteration(self.solverList[currentMesh])
        return self.solverList[currentMesh].failedFlag

    def solveMultilevel(self,uList,rList,bList=None,par_uList=None,par_rList=None):
        if bList is None:
            bList = [None for r in rList]
        self.infoString="********************Start Multilevel Nonlinear Solver Info*********************\n"
        for l in range(self.fineLevel):
            if self.tolList is not None:
                self.switchToResidualConvergence(self.solverList[l],self.tolList[l])
            self.solverList[l].solve(u=uList[l],r=rList[l],b=bList[l],par_u=par_uList[l],par_r=par_rList[l])
            self.infoString+="****************Start Level %i Info******************\n" %l
            self.infoString+=self.solverList[l].info()
            self.infoString+="****************End Level %i Info******************\n" %l
            if self.tolList is not None:
                self.revertToFixedIteration(self.solverList[l])
            #\todo see if there's a better way to do this
            #copy user u,r into internal
            #because fas will over write u and r with ones
            #that aren't solutions and we need to save
            #the true solutions on this level
            self.uList[l][:]=uList[l]
            self.rList[l][:]=rList[l]
            for ci,p in self.prolongList.iteritems():
                p[l+1].matvec(self.solverList[l].F.u[ci].dof,self.solverList[l+1].F.u[ci].dof)
            self.solverList[l+1].F.setFreeDOF(uList[l+1])
        if self.tolList is not None:
            self.switchToResidualConvergence(self.solverList[self.fineLevel],self.tolList[self.fineLevel])
        self.solverList[self.fineLevel].solve(u=uList[self.fineLevel],
                                              r=rList[self.fineLevel],
                                              b=bList[self.fineLevel],
                                              par_u=par_uList[self.fineLevel],
                                              par_r=par_rList[self.fineLevel])
        self.infoString+="****************Start Level %i Info******************\n" %self.fineLevel
        self.infoString+=self.solverList[self.fineLevel].info()
        self.infoString+="****************End Level %i Info******************\n" %self.fineLevel
        if self.tolList is not None:
            self.revertToFixedIteration(self.solverList[self.fineLevel])
        #reset u and r on other levels:
        for l in range(self.fineLevel):
            self.fList[l].setUnknowns(self.uList[l])
            rList[l][:]=self.rList[l]
        self.infoString+="********************End Multilevel Nonlinear Solver Info*********************\n"
        return self.solverList[-1].failedFlag

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
        solver.atol_r = self.atol
        solver.maxIts = self.maxIts
        solver.printInfo = self.printInfo

    def revertToFixedIteration(self,solver):
        solver.convergenceTest = self.saved_ctest
        solver.rtol_r = self.saved_rtol_r
        solver.atol_r = self.saved_atol_r
        solver.maxIts = self.saved_maxIts
        solver.printInfo = self.saved_printInfo


def multilevelNonlinearSolverChooser(nonlinearOperatorList,
                                     jacobianList,
                                     par_jacobianList,
                                     duList=None,
                                     par_duList=None,
                                     relativeToleranceList=None,
                                     absoluteTolerance=1.0e-8,
                                     multilevelNonlinearSolverType = NLNI,#Newton,NLJacobi,NLGaussSeidel,NLStarILU,
                                     computeSolverRates=False,
                                     printSolverInfo=False,
                                     linearSolverList=None,
                                     linearDirectSolverFlag=False,
                                     solverFullNewtonFlag=True,
                                     maxSolverIts=500,
                                     solverConvergenceTest='r',
                                     levelNonlinearSolverType=Newton,#NLJacobi,NLGaussSeidel,NLStarILU
                                     levelSolverFullNewtonFlag=True,
                                     levelSolverConvergenceTest='r',
                                     computeLevelSolverRates=False,
                                     printLevelSolverInfo=False,
                                     relaxationFactor=None,
                                     connectionListList=None,
                                     smootherType = 'Jacobi',
                                     prolong_bcList = None,
                                     restrict_bcList = None,
                                     restrict_bcSumList = None,
                                     prolongList = None,
                                     restrictList = None,
                                     restrictionRowSumList = None,
                                     preSmooths=3,
                                     postSmooths=3,
                                     cycles=3,
                                     smootherConvergenceTest='its',
                                     computeSmootherRates=False,
                                     printSmootherInfo=False,
                                     smootherFullNewtonFlag=True,
                                     computeCoarseSolverRates=False,
                                     printCoarseSolverInfo=False,
                                     EWtol=True,
                                     maxLSits=100,
                                     parallelUsesFullOverlap = True,
                                     nonlinearSolverNorm = l2Norm):
    if (multilevelNonlinearSolverType == Newton or
        multilevelNonlinearSolverType == NLJacobi or
        multilevelNonlinearSolverType == NLGaussSeidel or
        multilevelNonlinearSolverType == NLStarILU):
        levelNonlinearSolverType = multilevelNonlinearSolverType
    nLevels = len(nonlinearOperatorList)
    multilevelNonlinearSolver=None
    levelNonlinearSolverList=[]
    if levelNonlinearSolverType == FAS:
        preSmootherList=[]
        postSmootherList=[]
        mgItsList=[]
        for l in range(nLevels):
            mgItsList.append(cycles)
            if l > 0:
                if smootherType == NLJacobi:
                    if relaxationFactor is None:
                        relaxationFactor = 2.0/5.0#4.0/5.0
                    preSmootherList.append(NLJacobi(F=nonlinearOperatorList[l],
                                                    J=jacobianList[l],
                                                    du=duList[l],
                                                    weight=relaxationFactor,
                                                    maxIts=preSmooths,
                                                    convergenceTest=smootherConvergenceTest,
                                                    computeRates = computeSmootherRates,
                                                    printInfo=printSmootherInfo,
                                                    fullNewton=smootherFullNewtonFlag))
                    postSmootherList.append(NLJacobi(F=nonlinearOperatorList[l],
                                                     J=jacobianList[l],
                                                     du=duList[l],
                                                     weight=relaxationFactor,
                                                     maxIts=postSmooths,
                                                     convergenceTest=smootherFullNewtonFlag,
                                                     computeRates = computeSmootherRates,
                                                     printInfo=printSmootherInfo,
                                                     fullNewton=smootherFullNewtonFlag))
                elif smootherType == NLGaussSeidel:
                    if relaxationFactor is None:
                        relaxationFactor = 3.0/5.0
                    preSmootherList.append(NLGaussSeidel(connectionList = connectionListList[l],
                                                         F=nonlinearOperatorList[l],
                                                         J=jacobianList[l],
                                                         du=duList[l],
                                                         weight=relaxationFactor,
                                                         maxIts=preSmooths,
                                                         convergenceTest=smootherConvergenceTest,
                                                         computeRates = computeSmootherRates,
                                                         printInfo=printSmootherInfo,
                                                         fullNewton=smootherFullNewtonFlag))
                    postSmootherList.append(NLGaussSeidel(connectionList = connectionListList[l],
                                                          F=nonlinearOperatorList[l],
                                                          J=jacobianList[l],
                                                          du=duList[l],
                                                          weight=relaxationFactor,
                                                          maxIts=postSmooths,
                                                          convergenceTest=smootherConvergenceTest,
                                                          computeRates = computeSmootherRates,
                                                          printInfo=printSmootherInfo,
                                                          fullNewton=smootherFullNewtonFlag))
                elif smootherType == NLStarILU:
                    if relaxationFactor is None:
                        relaxationFactor = 2.0/5.0
                    preSmootherList.append(NLStarILU(connectionList = connectionListList[l],
                                                     F = nonlinearOperatorList[l],
                                                     J = jacobianList[l],
                                                     du=duList[l],
                                                     weight = relaxationFactor,
                                                     maxIts = preSmooths,
                                                     convergenceTest = smootherConvergenceTest,
                                                     computeRates = computeSmootherRates,
                                                     printInfo = printSmootherInfo,
                                                     fullNewton = smootherFullNewtonFlag))
                    postSmootherList.append(NLStarILU(connectionList = connectionListList[l],
                                                      F=nonlinearOperatorList[l],
                                                      J=jacobianList[l],
                                                      du=duList[l],
                                                      weight=relaxationFactor,
                                                      maxIts=postSmooths,
                                                      convergenceTest=smootherConvergenceTest,
                                                      computeRates = computeSmootherRates,
                                                      printInfo=printSmootherInfo,
                                                      fullNewton=smootherFullNewtonFlag))
                else:
                    raise RuntimeError("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!smootherType unrecognized")
            else:
                if smootherType == NLJacobi:
                    if relaxationFactor is None:
                        relaxationFactor = 4.0/5.0
                    coarseSolver = NLJacobi(F=nonlinearOperatorList[l],
                                            J=jacobianList[l],
                                            du=duList[l],
                                            weight=relaxationFactor,
                                            maxIts=postSmooths,
                                            convergenceTest=smootherFullNewtonFlag,
                                            computeRates = computeSmootherRates,
                                            printInfo=printSmootherInfo,
                                            fullNewton=smootherFullNewtonFlag,
                                            norm = nonlinearSolverNorm)
                elif smootherType == NLGaussSeidel:
                    if relaxationFactor is None:
                        relaxationFactor = 3.0/5.0
                    coarseSolver = NLGaussSeidel(connectionList = connectionListList[l],
                                                 F=nonlinearOperatorList[l],
                                                 J=jacobianList[l],
                                                 du=duList[l],
                                                 weight=relaxationFactor,
                                                 maxIts=postSmooths,
                                                 convergenceTest=smootherConvergenceTest,
                                                 computeRates = computeSmootherRates,
                                                 printInfo=printSmootherInfo,
                                                 fullNewton=smootherFullNewtonFlag,
                                                 norm = nonlinearSolverNorm)
                elif smootherType == NLStarILU:
                    if relaxationFactor is None:
                        relaxationFactor = 2.0/5.0
                    coarseSolver = NLStarILU(connectionList = connectionListList[l],
                                             F = nonlinearOperatorList[l],
                                             J = jacobianList[l],
                                             du=duList[l],
                                             weight = relaxationFactor,
                                             maxIts = preSmooths,
                                             convergenceTest = smootherConvergenceTest,
                                             computeRates = computeSmootherRates,
                                             printInfo = printSmootherInfo,
                                             fullNewton = smootherFullNewtonFlag,
                                             norm = nonlinearSolverNorm)
                else:
                    raise RuntimeError("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!smootherType unrecognized")
                preSmootherList.append([])
                postSmootherList.append([])

        levelNonlinearSolver = FAS(prolongList = prolongList,
                                   restrictList = restrictList,
                                   restrictSumList = restrictionRowSumList,
                                   FList = nonlinearOperatorList,
                                   preSmootherList = preSmootherList,
                                   postSmootherList = postSmootherList,
                                   coarseSolver = coarseSolver,
                                   mgItsList = mgItsList,
                                   printInfo=printLevelSolverInfo)
        levelNonlinearSolverList = levelNonlinearSolver.solverList
    elif levelNonlinearSolverType == Newton:
        for l in range(nLevels):
            if par_duList is not None and len(par_duList) > 0:
                par_du=par_duList[l]
            else:
                par_du=None
            levelNonlinearSolverList.append(Newton(linearSolver=linearSolverList[l],
                                                   F=nonlinearOperatorList[l],
                                                   J=jacobianList[l],
                                                   du=duList[l],
                                                   par_du=par_du,
                                                   rtol_r=relativeToleranceList[l],
                                                   atol_r=absoluteTolerance,
                                                   maxIts=maxSolverIts,
                                                   norm = nonlinearSolverNorm,
                                                   convergenceTest = levelSolverConvergenceTest,
                                                   computeRates = computeLevelSolverRates,
                                                   printInfo=printLevelSolverInfo,
                                                   fullNewton=levelSolverFullNewtonFlag,
                                                   directSolver=linearDirectSolverFlag,
                                                   EWtol=EWtol,
                                                   maxLSits=maxLSits ))
    elif levelNonlinearSolverType in [POD_Newton,POD_DEIM_Newton]:
        for l in range(nLevels):
            if par_duList is not None and len(par_duList) > 0:
                par_du=par_duList[l]
            else:
                par_du=None
            levelNonlinearSolverList.append(levelNonlinearSolverType(linearSolver=linearSolverList[l],
                                                                     F=nonlinearOperatorList[l],
                                                                     J=jacobianList[l],
                                                                     du=duList[l],
                                                                     par_du=par_du,
                                                                     rtol_r=relativeToleranceList[l],
                                                                     atol_r=absoluteTolerance,
                                                                     maxIts=maxSolverIts,
                                                                     norm = nonlinearSolverNorm,
                                                                     convergenceTest = levelSolverConvergenceTest,
                                                                     computeRates = computeLevelSolverRates,
                                                                     printInfo=printLevelSolverInfo,
                                                                     fullNewton=levelSolverFullNewtonFlag,
                                                                     directSolver=linearDirectSolverFlag,
                                                                     EWtol=EWtol,
                                                                     maxLSits=maxLSits ))
    elif levelNonlinearSolverType == NewtonNS:
        for l in range(nLevels):
            if par_duList is not None and len(par_duList) > 0:
                par_du=par_duList[l]
            else:
                par_du=None
            levelNonlinearSolverList.append(NewtonNS(linearSolver=linearSolverList[l],
                                                   F=nonlinearOperatorList[l],
                                                   J=jacobianList[l],
                                                   du=duList[l],
                                                   par_du=par_du,
                                                   rtol_r=relativeToleranceList[l],
                                                   atol_r=absoluteTolerance,
                                                   maxIts=maxSolverIts,
                                                   norm = nonlinearSolverNorm,
                                                   convergenceTest = levelSolverConvergenceTest,
                                                   computeRates = computeLevelSolverRates,
                                                   printInfo=printLevelSolverInfo,
                                                   fullNewton=levelSolverFullNewtonFlag,
                                                   directSolver=linearDirectSolverFlag,
                                                   EWtol=EWtol,
                                                   maxLSits=maxLSits ))
    elif levelNonlinearSolverType == NLJacobi:
        if relaxationFactor is None:
            relaxationFactor = 4.0/5.0
        for l in range(nLevels):
            levelNonlinearSolverList.append(NLJacobi(F=nonlinearOperatorList[l],
                                                     J=jacobianList[l],
                                                     du=duList[l],
                                                     rtol_r=relativeToleranceList[l],
                                                     atol_r=absoluteTolerance,
                                                     norm = nonlinearSolverNorm,
                                                     maxIts=maxSolverIts,
                                                     convergenceTest = levelSolverConvergenceTest,
                                                     weight=relaxationFactor,
                                                     computeRates = computeLevelSolverRates,
                                                     printInfo=printLevelSolverInfo,
                                                     fullNewton=levelSolverFullNewtonFlag))
    elif levelNonlinearSolverType == NLGaussSeidel:
        if relaxationFactor is None:
            relaxationFactor = 4.0/5.0
        for l in range(nLevels):
            levelNonlinearSolverList.append(NLGaussSeidel(F=nonlinearOperatorList[l],
                                                          J=jacobianList[l],
                                                          du=duList[l],
                                                          connectionList = connectionListList[l],
                                                          rtol_r=relativeToleranceList[l],
                                                          atol_r=absoluteTolerance,
                                                          maxIts=maxSolverIts,
                                                          convergenceTest = levelSolverConvergenceTest,
                                                          weight=relaxationFactor,
                                                          computeRates = computeLevelSolverRates,
                                                          printInfo=printLevelSolverInfo,
                                                          fullNewton=levelSolverFullNewtonFlag))
    elif levelNonlinearSolverType == NLStarILU:
        if relaxationFactor is None:
            relaxationFactor = 3.0/5.0
        for l in range(nLevels):
            levelNonlinearSolverList.append(NLStarILU(F=nonlinearOperatorList[l],
                                                      J=jacobianList[l],
                                                      du=duList[l],
                                                      connectionList = connectionListList[l],
                                                      rtol_r=relativeToleranceList[l],
                                                      atol_r=absoluteTolerance,
                                                      maxIts=maxSolverIts,
                                                      norm = nonlinearSolverNorm,
                                                      convergenceTest = levelSolverConvergenceTest,
                                                      weight=relaxationFactor,
                                                      computeRates = computeLevelSolverRates,
                                                      printInfo=printLevelSolverInfo,
                                                      fullNewton=levelSolverFullNewtonFlag))
    #mwf what about wrapping Eikonal equation solvers as specific types of nonlinear solvers?
    elif levelNonlinearSolverType == FMMEikonalSolver or levelNonlinearSolverType == FSWEikonalSolver:
        assert multilevelNonlinearSolverType == MultilevelEikonalSolver,  "level solver requires MultilevelEikonalSolver "
        for l in range(nLevels):
            levelNonlinearSolverList.append(EikonalSolver(levelNonlinearSolverType,
                                                          nonlinearOperatorList[l],
                                                          relativeTolerance = relativeToleranceList[l],
                                                          absoluteTolerance = absoluteTolerance,
                                                          maxSolverIts      = maxSolverIts))
                                                          #can set front tolerance and eikonal variable index too


    elif levelNonlinearSolverType == SSPRKNewton:
        for l in range(nLevels):
            if par_duList is not None and len(par_duList) > 0:
                par_du=par_duList[l]
            else:
                par_du=None
            levelNonlinearSolverList.append(SSPRKNewton(linearSolver=linearSolverList[l],
                                                        F=nonlinearOperatorList[l],
                                                        J=jacobianList[l],
                                                        du=duList[l],
                                                        par_du=par_du,
                                                        rtol_r=relativeToleranceList[l],
                                                        atol_r=absoluteTolerance,
                                                        maxIts=maxSolverIts,
                                                        norm = nonlinearSolverNorm,
                                                        convergenceTest = levelSolverConvergenceTest,
                                                        computeRates = computeLevelSolverRates,
                                                        printInfo=printLevelSolverInfo,
                                                        fullNewton=levelSolverFullNewtonFlag,
                                                        directSolver=linearDirectSolverFlag,
                                                        EWtol=EWtol,
                                                        maxLSits=maxLSits ))

    else:
        try:
            for l in range(nLevels):
                if par_duList is not None and len(par_duList) > 0:
                    par_du=par_duList[l]
                else:
                    par_du=None
                levelNonlinearSolverList.append(levelNonlinearSolverType(linearSolver=linearSolverList[l],
                                                       F=nonlinearOperatorList[l],
                                                       J=jacobianList[l],
                                                       du=duList[l],
                                                       par_du=par_du,
                                                       rtol_r=relativeToleranceList[l],
                                                       atol_r=absoluteTolerance,
                                                       maxIts=maxSolverIts,
                                                       norm = nonlinearSolverNorm,
                                                       convergenceTest = levelSolverConvergenceTest,
                                                       computeRates = computeLevelSolverRates,
                                                       printInfo=printLevelSolverInfo,
                                                       fullNewton=levelSolverFullNewtonFlag,
                                                       directSolver=linearDirectSolverFlag,
                                                       EWtol=EWtol,
                                                       maxLSits=maxLSits ))
        except:
            raise RuntimeError,"Unknown level nonlinear solver "+ levelNonlinearSolverType
    if multilevelNonlinearSolverType == NLNI:
        multilevelNonlinearSolver = NLNI(fList = nonlinearOperatorList,
                                         solverList = levelNonlinearSolverList,
                                         prolongList = prolong_bcList,
                                         restrictList = restrict_bcList,
                                         restrictSumList = restrict_bcSumList,
                                         maxIts = maxSolverIts,
                                         tolList = relativeToleranceList,
                                         atol=absoluteTolerance,
                                         computeRates = computeSolverRates,
                                         printInfo=printSolverInfo)
    elif (multilevelNonlinearSolverType == Newton or
          multilevelNonlinearSolverType == POD_Newton or
          multilevelNonlinearSolverType == POD_DEIM_Newton or
          multilevelNonlinearSolverType == NewtonNS or
          multilevelNonlinearSolverType == NLJacobi or
          multilevelNonlinearSolverType == NLGaussSeidel or
          multilevelNonlinearSolverType == NLStarILU):
        multilevelNonlinearSolver = MultilevelNonlinearSolver(nonlinearOperatorList,
                                                              levelNonlinearSolverList,
                                                              computeRates = computeSolverRates,
                                                              printInfo = printSolverInfo)
    elif multilevelNonlinearSolverType == MultilevelEikonalSolver:
        #should I take care of assignment here
        logEvent("Warning Using Multilevel Eikonal Equation Solver, hope equation is correct!")
        multilevelNonlinearSolver = MultilevelEikonalSolver(levelNonlinearSolverList,
                                                            printInfo = printSolverInfo)

    else:
        raise RuntimeError,"!!!!!!!!!!!!!!!!!!!!!!!!!Unknown multilevelNonlinearSolverType " + multilevelNonlinearSolverType

    #add minimal configuration for parallel?
    for levelSolver in multilevelNonlinearSolver.solverList:
        levelSolver.par_fullOverlap = parallelUsesFullOverlap
    return multilevelNonlinearSolver
