from LinearAlgebra import *
from TimeIntegrationTools import *
from NonlinearSolvers import *
from ScalarTransport import *
import numpy

class PhiDummy(object):
    """
    dummy class for passing off FEM potential function to TimeIntegration
    classes

    """
    def __init__(self,dim):
        self.dim=dim
        self.dof=dim
    #end init
#end PhiDummy

class dUdtEqMinusLambda(object):
    """
    test problem for TimeIntegration classes. Integrates the
    canonical test problem

       du/dt = -lambda . u

    from t=0. Here, lambda is a positive real scalar. We'll call the
    right hand side -r(u) to be consistent with the ScalarTransport classes

    Rather than take a method of lines approach, a Rothe paradigm
    is followed, so the problem class is 'in control' of the integration.
    It queries the TimeIntegrator for the step size, and provides the
    interface for the nonlinear equation solvers to solve the resulting
    discrete problem for the solution at the next time level or stage.

    The NonlinearEquation interface is

    self.dim
    getResidual(u,r)
    getJacobian(jacobian)

    The TimeIntegration classes expect that this class has a
    dictionary of quantities q
    q usually correspond to numerical integration
    points on a spatial mesh.

    """
    def __init__(self,u0,TimeIntegrationClass=BackwardEuler,
                 cfl=Numeric.ones(1,Numeric.Float),
                 lam=Numeric.ones(1,Numeric.Float),
                 nstages=1,order=1):
        self.dim = Numeric.size(lam,0)
        self.lam= lam
        #element quadrature terms ala ScalarTransport
        self.q = {}
        #numerical flux quadrature terms ala ScalarTransport
        self.ebq_global = {}
        # #  set initial conditions # #
        for term in ['u','m','mt','dm','dmt','f','div(f)',
                     'a','grad(u)','r','dr','cfl']:
            self.q[term] = Numeric.zeros(self.dim,Numeric.Float)
        #end for
        self.q['cfl'].flat[:] = cfl.flat[:]
        #potential is m=u, make a dummy
        self.phi = PhiDummy(self.dim)
        self.nFreeDOF_global=self.dim
        #solution mass, is m=u
        #initial condition is u=1
        self.q['u'].flat[:] = u0.flat[:]
        self.q['m'].flat[:] = self.q['u'].flat[:]
        #reaction term is actually where right hand side is implemented
        self.q['r']= lam*self.q['m']
        # #  end array initialization
        if TimeIntegrationClass == SSPRKintegration:
            self.timeIntegrator = TimeIntegrationClass(self,order)
        else:
            self.timeIntegrator = TimeIntegrationClass(self)
        #end if
    #end init
    def setUnknowns(self,u):
        """
        copy over solution values u into dictionary holding unknowns
        """
        self.q['u'].flat[:] = u.flat[:]
    #end setUnknowns
    def updateCoefficients(self):
        """
        just evalute reaction term and mass term
        """
        self.q['m'].flat[:] = self.q['u'].flat[:]
        self.q['dm'].flat[:]= 1.0
        self.q['r'].flat[:] = self.lam*self.q['u'].flat[:]
        self.q['dr'].flat[:]= self.lam

    #end updateCoefficients
    def initializeTimeIntegration(self):
        """
        get ready for next time step
        """
        self.updateCoefficients()
        self.timeIntegrator.setInitialStageValues()
        self.timeIntegrator.updateTimeHistory()
    #end initializeTimeIntegration

    def getResidual(self,u,r):
        """
        evaluate residual r = du/dt + r(u)
                            = du/dt + lam*u
        """
        r.flat[:] = 0.0
        self.timeIntegrator.calculateU(u)
        self.setUnknowns(u)
        self.updateCoefficients()
        self.timeIntegrator.updateMass(self.q['m'],self.q['mt'],
                                        self.q['dm'],self.q['dmt'])
        self.timeIntegrator.updateReaction(self.q['r'],self.q['dr'])
        r.flat[:] = self.q['mt'] + self.q['r']

    #end get residual

    def getJacobian(self,jacobian):
        """
        build jacobian

          d/du(dm/dt) + dr/du

        where

          dr/du = -lambda
        """
        beginAssembly(jacobian)
        for i in range(self.dim):
            for j in range(self.dim):
                jacobian[i,j] = 0.0
            #end j
        #end i
        #get mass derivative from time integration or assume
        #updateMass has already been called
        for i in range(self.dim):
            #mwf debug
            #print """dUdt...getJac(%d,%d) dmt=%g dr=%g""" % (i,i,
            #                                                 self.q['dmt'][i],
            #                                                 self.q['dr'][i])
            jacobian[i,i] += self.q['dmt'][i] + self.q['dr'][i]
        #end i
        endAssembly(jacobian)
#end dUdt


if __name__ == '__main__':
    import sys
    import numpy
    from TimeIntegrationTools import *
    from LinearSolvers import *

    from optparse import OptionParser
    parser = OptionParser()
    #options controlling simulation behavior

    parser.add_option('-C','--cfl',
                      default=0.1,
                      help="""target cfl number [0.1]""")

    parser.add_option('-L','--lamScale',
                      default=1.0,
                      help="""basic \lambda in \dot{\vec u}=-\mat{\Lamda}\vec u
                      where \mat{\Lamda}_ii = -\lambda*i
                      [1]""")

    parser.add_option('-m','--method',
                      default=0,
                      help="""which integration method to use [0]
                              0 -- Backward Euler
                              1 -- Forward Euler
                              2 -- SSPRK order 1--5
                       """)

    parser.add_option('-n','--neq',
                      default=1,
                      help="""number of equations in system
                              \dot{\vec u}=-\mat{\Lamda}\vec u [1]""")

    parser.add_option('-O','--rkOrder',
                      default=1,
                      help="""order for SSPRK integrator [1]""")

    parser.add_option('--ntMax',
                      default=100,
                      help="""max number of steps allowed. If adaptDT is false
                      then takes ntMax steps [100] """)

    parser.add_option('--adaptDT',
                      default=False,
                      action='store_true',
                      help="""use time step adaption and error control? [False]
                      """)

    parser.add_option('-T','--tend',
                      default=1.0,
                      help="""stopping time in simulation [1.0]""")



    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")
    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args
    #
    verbose = int(options.verbose)

    #problem size
    dim = int(options.neq)
    #need positive real scalar value
    lamScale = float(options.lamScale)
    lam = (1.0+Numeric.arange(dim,dtype=Numeric.Float))*lamScale
    lamMax = max(lam)
    T   = float(options.tend)
    ntMax = int(options.ntMax)
    adaptDT = bool(options.adaptDT)
    methodFlag = int(options.method)
    targetCFL  = float(options.cfl)
    #pick a time integrator
    nstages = 1
    order   = 1
    TimeIntegrationClass = BackwardEuler
    if methodFlag == 1:
        TimeIntegrationClass = ForwardEuler
    elif methodFlag == 2:
        TimeIntegrationClass = SSPRKintegration
        order   = int(options.rkOrder)
        nstages = order
    #end
    #estimate a stable maximum time step?

    #initial time
    t0= 0.0

    #initial condition and solution
    u0 = Numeric.zeros(dim,Numeric.Float)
    u  = Numeric.zeros(dim,Numeric.Float)
    u0.flat[:] = 1.0
    u.flat[:]  = 1.0

    #use lambda to determine max allowable step
    cfl= Numeric.ones(dim,Numeric.Float)
    cfl.flat[:] = lam.flat[:]

    ode = dUdtEqMinusLambda(u0,TimeIntegrationClass,cfl,lam,nstages,order)

    ode.updateCoefficients()
    ode.timeIntegrator.chooseDT()

    r0 = Numeric.ones(dim,Numeric.Float)
    r  = Numeric.ones(dim,Numeric.Float)

    ode.getResidual(u0,r0)
    ode.getResidual(u0,r)

    if verbose > 1:
        print('u0= ',u0)
        print('initial residual = ',r0)
        print('initial dt= ',ode.timeIntegrator.DT)
    #end verbose

    if verbose > 5:
        print('ode.q[u] = ',ode.q['u'])
        print('ode.q[m] = ',ode.q['m'])
        print('ode.q[r] = ',ode.q['r'])
        print('ode.q[mt]= ',ode.q['mt'])
        print('ode.q[dm]= ',ode.q['dm'])
        print('ode.q[dr]= ',ode.q['dr'])
    #end if


    # # create nonlinear system for solving problem
    if verbose > 0:
        print('creating nonlinear system')
    #end verbose
    MatType = Mat
    matType = 'dense'
    linearSolverType=  'DenseLU'

    jacobian = MatType(ode.dim,ode.dim,ode.dim)

    # # create linear system for solving problem
    if verbose > 5:
        print('creating linear system')
        print('initial jacobian mat =\n',jacobian)
    #end verbose

    linearSolver = DenseLU(jacobian)
    printNLinfo = False
    if verbose > 10:
        printNLinfo = True
    nlSolver = Newton(linearSolver,ode,jacobian,printInfo=printNLinfo)

    if verbose > 5:
        print('ode u0= ',u0)
        ode.getJacobian(jacobian)
        print('ode residual= ',r0)
        print('system jacobian =\n',jacobian)
    # # solving problem

    #now try a time loop?
    ode.timeIntegrator.runCFL = targetCFL
    t = t0
    nsteps = 0
    dtFix = T/ntMax
    L2err = 0.0
    LIerr = 0.0
    errT  = 0.0
    while t < T and nsteps < ntMax:

        ode.timeIntegrator.chooseDT()
        dtMin = min(T-t,ode.timeIntegrator.DT)
        if not adaptDT:
            dtMin = min(dtMin,dtFix)
        #end not adapting

        ode.timeIntegrator.DT=dtMin
        if nsteps == 0:
            ode.initializeTimeIntegration()
        else:
            ode.timeIntegrator.updateTimeHistory()
        #end if
        t += ode.timeIntegrator.DT #do this before or after solve?

        for i in range(nstages):
            ode.getResidual(u,r)
            #solve system for new time level
            nlSolver.solve(u,r)
            #have to have this here for multistage values?
            #ode.updateCoefficients()
            ode.timeIntegrator.updateStage()
        #end stage loop
        nsteps += 1

        #exact solution?
        uex = u0*Numeric.exp(-lam*t)
        err = u-uex
        errT= l2Norm(err)
        L2err += errT**2*ode.timeIntegrator.DT
        LIerr = max(LIerr,max(abs(err)))
        if verbose > 2:
            print("""t= %g \n\tu  = %s \n\tuex= %s
        err= %g LIe= %g""" % (t,u,uex,errT,LIerr))
        #end verbose

        #make sure solution is kept by ode
        ode.setUnknowns(u)

    #end time loop
    #go ahead and get coefficients at final time solution
    ode.updateCoefficients()

    L2err = sqrt(L2err)

    print("""reached t= %g \n\tu  = %s \n\tuex= %s
    err= %g L2err= %g LIe= %g""" % (t,u,uex,errT,L2err,LIerr))
