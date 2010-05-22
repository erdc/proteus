from math import *
import cPickle
import os
from MeshTools import *
from FemTools import *
from QuadTools import *
from LinearSolvers import *
from NonlinearSolvers import *
#import Gnuplot
from NormTools import *
from AnalyticalSolutions import *
from ScalarTransport import *
from LatexReport import *
## \defgroup ScalarTransportTests ScalarTransportTests
#
# @{

"""
A module for testing models of nonlinear advection-diffusion-reaction equations.
"""

def runTests():
    #First define the mesh independent problem definition
    #consisting of the following information
    #see examples below for definitions
    nd={}  #number of spatial dimensions
    getDirichletConditions={}
    fluxBoundaryConditions={}
    coefficients={}
    analyticalSolution={}
    timeIntegration = {}
    getInitialConditions ={}
    T = {}
    fullNewtonFlag = {}
    try:
        os.mkdir('tmp')
    except:
        pass
    report = openLatexReport('tmp/scalarTransportReport.tex','Scalar Transport Report')
    #Define test problems
    testProblems = []
    #1D
    #
    # (bu - a u_x)_x = 0; u(0) = 1; u(1) = 0
    #
    test='Linear-AD-SS-1D'
    testProblems.append(test)
    timeIntegration[test] = NoIntegration
    nd[test]=1
    def getDBC(x):
        if x == 0.0:
            return lambda x,t: 1.0
        if x == 1.0:
            return lambda x,t: 0.0
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    #a0=1.0e-4
    #a0=5.0e-2
    a0=0.1
    #a0=0.0
    A0=Numeric.array([[a0]])
    #b0=1.0
    b0=1.0
    B0=Numeric.array([b0])
    C0=1.0
    M0=0.0        
    coefficients[test] = LinearADR_ConstantCoefficients(M=0.0,A=A0,B=B0,C=0.0)
    if a0 > 0.0 and (b0/a0 < 1.0/2.5e-2):
        analyticalSolution[test] = LinearAD_SteadyState(b=b0,a=a0)
    else:
        analyticalSolution[test] = None
    getInitialConditions[test]=None
    fullNewtonFlag[test] = False
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    #define nonlinear coefficients for
    # M*u^p_t + \div( B*u^q - Au^t \grad u^r) + C*u^s = 0
    #
    class NonlinearADR_pqrst(ScalarTransportCoefficients):
        from transportCoefficients import nonlinearADR_pqrstEvaluate
        def __init__(self,M,A,B,C,
                     p=1.0,q=1.0,r=1.0,s=1.0,t=0.0):
            self.M = M
            self.A = A
            self.B = B
            self.C = C
            self.p = p
            self.q = q
            self.r = r
            self.s = s
            self.t = t
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.nonlinearADR_pqrstEvaluate(u.shape[0],
                                                f.shape[1],
                                                self.M,
                                                self.A,
                                                self.B,
                                                self.C,
                                                self.p,
                                                self.q,
                                                self.r,
                                                self.s,
                                                self.t,
                                                t,
                                                x,
                                                u,
                                                m,dm,
                                                f,df,
                                                a,da,
                                                phi,dphi,
                                                r,dr)
            else:
                maxu0 = Numeric.maximum(u,0.0)

                if self.p > 1.0:
                    Numeric.multiply(Numeric.power(maxu0,self.p,m),self.M,m)
                    Numeric.multiply(Numeric.power(maxu0,self.p-1,dm),self.M*self.p,dm)
                else:
                    m[:] = u
                    m *= self.M
                    dm[:] = self.M
                if self.q > 1.0:
                    Numeric.power(maxu0,self.q,f[:,0])
                    for ix in range(1,a.shape[2]-1):
                        f[:,ix]=f[:,0]
                    Numeric.multiply(f,self.B,f)
                    Numeric.power(maxu0,self.q-1,df[:,0])
                    for ix in range(1,a.shape[2]-1):
                        df[:,ix]=df[:,0]
                    Numeric.multiply(df,self.B*self.q,df)
                else:
                    Numeric.multiply(self.B,u[:,Numeric.newaxis],f)
                    df[:] = self.B
                if self.t > 0.0:
                    Numeric.power(maxu0,self.t,a[:,0,0])
                    Numeric.power(maxu0,self.t-1,da[:,0,0])
                    for ix in range(1,a.shape[2]-1):
                        a[:,ix,ix]=a[:,0,0]
                        da[:,ix,ix]=da[:,0,0]
                    Numeric.multiply(a,self.A,a)
                    Numeric.multiply(da,self.A*self.t,da)
                else:
                    a[:] = self.A
                    da[:,:,:] = 0.0
                if self.r > 1.0:
                    Numeric.power(maxu0,self.r,phi)
                    Numeric.multiply(Numeric.power(maxu0,self.r-1,dphi),self.r,dphi)
                else:
                    phi[:]=u
                    dphi[:]=1.0
                if self.s > 1.0:
                    Numeric.multiply(Numeric.power(maxu0,self.s,r),self.C,r)
                    Numeric.multiply(Numeric.power(maxu0,self.s-1,dr),self.C*self.s,dr)
                else:
                    r[:] = u
                    r*=self.C
                    dr[:] = self.C

    #
    # (v max(u,0)^2 - a u_x)_x = 0; u(0) = 1; u(1) = 0
    #
    test='NonlinearAD-q2r1-SS-1D'
    testProblems.append(test)
    timeIntegration[test] = NoIntegration
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    coefficients[test] = NonlinearADR_pqrst(M=0.0,A=A0,B=B0,C=0.0,q=2,t=1)
    if a0 != 0.0 and b0/a0 < 1.0/2.5e-2:
        analyticalSolution[test] = NonlinearAD_SteadyState(b=abs(b0),
                                                           q=2,
                                                           a=a0,
                                                           r=1)
    else:
        analyticalSolution[test] = None
    getInitialConditions[test]=None
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].reaction = None
    #
    # [vu - a (max(u,0)^2)_x]_x = 0; u(0) = 1; u(1) = 0
    #
    #just need to redefine f and phi from last problem
    test='NonlinearAD-q1r2-SS-1D'
    nd[test]=1
    testProblems.append(test)
    timeIntegration[test] = NoIntegration
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    coefficients[test] = NonlinearADR_pqrst(M=0.0,A=A0,B=B0,C=0.0,q=1,r=2)
    if a0 != 0.0 and b0/a0 < 1.0/2.5e-2:
        analyticalSolution[test] = NonlinearAD_SteadyState(b=abs(b0),q=1,a=a0,r=2)
    else:
        analyticalSolution[test] = None
    getInitialConditions[test]=None
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    # [vu - a (max(u,0)^2) u_x]_x = 0; u(0) = 1; u(1) = 0
    #
    test='NonlinearAD-q1r1t2-SS-1D'
    nd[test]=1
    testProblems.append(test)
    timeIntegration[test] = NoIntegration
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    diff_t = 2
    coefficients[test] = NonlinearADR_pqrst(M=0.0,A=A0,B=B0,C=0.0,q=1,r=1,t=diff_t)
    analyticalSolution[test] = None
    getInitialConditions[test]=None
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    #define nonlinear coefficients for
    # M*u^p_t + \div( B*u^q - A \grad (1 - max(1-u)^r) + C*u^s = 0
    #
    class NonlinearADR_pqrstDual(NonlinearADR_pqrst):
        from transportCoefficients import nonlinearADR_pqrstDualEvaluate
        def __init__(self,M,A,B,C,
                     p1=1.0,q1=1.0,r1=1.0,s1=1.0,t1=0.0,
                     p2=1.0,q2=1.0,r2=1.0,s2=1.0,t2=0.0):
            NonlinearADR_pqrst.__init__(self,M=M,A=A,B=B,C=C,
                                        p=p1,q=q1,r=r1,s=s1,t=t1)
            self.p2 = p2
            self.q2 = q2
            self.r2 = r2
            self.s2 = s2
            self.t2 = t2
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.nonlinearADR_pqrstDualEvaluate(u.shape[0],
                                                    f.shape[1],
                                                    self.M,
                                                    self.A,
                                                    self.B,
                                                    self.C,
                                                    self.p,
                                                    self.q,
                                                    self.r,
                                                    self.s,
                                                    self.t,
                                                    self.p2,
                                                    self.q2,
                                                    self.r2,
                                                    self.s2,
                                                    self.t2,
                                                    t,
                                                    x,
                                                    u,
                                                    m,dm,
                                                    f,df,
                                                    a,da,
                                                    phi,dphi,
                                                    r,dr)
            else:
                NonlinearADR_pqrst.evaluate(self,t,x,u,m,dm,f,df,a,da,phi,dphi,r,dr)
                max_1mu_0 = Numeric.maximum(1.0 - u,0.0)
                if self.p2 > 1.0:
                    m[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.p2),m)
                    dm[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.p2-1),self.p2*dm)
                if self.q2 > 1.0:
                    f[:] = Numeric.multiply(Numeric.power(max_1mu_0,self.q2)[:,Numeric.newaxis],f)
                    df[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.q2-1)[:,Numeric.newaxis],self.q2*df)
                if self.t2 > 0.0:
                    atmp = Numeric.power(max_1mu_0,self.t2)
                    Numeric.multiply(atmp,a[:,0,0],a[:,0,0])
                    datmp = Numeric.power(max_1mu_0,self.t2-1)
                    Numeric.multiply(datmp,self.t2*da[:,0,0],da[:,0,0])
                    for ix in range(1,a.shape[2]-1):
                        Numeric.multiply(atmp,a[:,ix,ix],a[:,ix,ix])
                        Numeric.multiply(datmp,self.t2*da[:,ix,ix],da[:,ix,ix])
                if self.r2 > 1.0:
                    phi[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.r2),phi)
                    dphi[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.r2-1),self.r2*dphi)
                if self.s2 > 1.0:
                    r[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.s2*r),r)
                    dr[:]=Numeric.multiply(Numeric.power(max_1mu_0,self.s2-1),self.s2*dr)
    #
    #
    # [v max(u,0)^2 - a (1.0 - max(1.0-u,0)^2)_x] = 0; u(0) = 1; u(1) = 0
    #
    #just need to redefine f and phi from last problem
    test='NonlinearAD-simpRE-SS-1D'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    coefficients[test] = NonlinearADR_pqrstDual(M=0.0,A=A0,B=B0,C=0.0,q1=1,q2=2,r1=2,r2=1)
    analyticalSolution[test] = None
    getInitialConditions[test]=None
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    #
    # head based Richards' with Mualem-van Genuchten
    #
    #just need to redefine f and phi from last problem
    test='NonlinearAD-RE-MvG'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    #
    #porous media properties
    #
    lengthScale   = 0.1     #m
    permeability  = 1.0e-7  #m^2
    viscosity     = 8.9e-4  #kg/(m*s)
    density       = 997.0   #kg/m^3
    gravity       = 9.8     #m/s^2
    thetaS        = 0.301   #-
    thetaR        = 0.093   #-
    mvg_alpha     = 5.47    #1/m
    mvg_n         = 4.264
    mvg_m         = 1.0 - 1.0/mvg_n
    #make non-dimensional
    dimensionless_conductivity  = density*gravity*permeability/(viscosity*sqrt(gravity*lengthScale))
    dimensionless_density  = 1.0
    dimensionless_gravity  = Numeric.array([1.0,
                                            0.0,
                                            0.0])
    dimensionless_alpha    = mvg_alpha*lengthScale
    def getDBC_Richards(x):
        if x == 0.0:
            return lambda x,t: 0.1/lengthScale   #0.1 m of  pressure head
        if x == 1.0:
            return lambda x,t: -10.0/lengthScale # -10.0 m of pressure head 
    getDirichletConditions[test]=getDBC_Richards
    class ConservativeHeadRichardsMualemVanGenuchten:
        from transportCoefficients import  conservativeHeadRichardsMualemVanGenuchtenHomEvaluate
        def __init__(self,
                     hydraulicConductivity,
                     gravity,
                     density,
                     thetaS,
                     thetaR,
                     alpha,
                     n,
                     m):
            self.Ks = hydraulicConductivity
            self.gravity=gravity
            self.rho = density
            self.thetaS = thetaS
            self.thetaR = thetaR
            self.thetaSR = thetaS - thetaR
            self.alpha = alpha
            self.n = n
            self.m = m
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr): 
            phi.flat[:]=u.flat
            dphi.flat[:]=1.0
            r.flat[:]=0.0
            dr.flat[:]=0.0
            self.conservativeHeadRichardsMualemVanGenuchtenHomEvaluate(u.shape[0],
                                                                       f.shape[1],
                                                                       self.rho,
                                                                       self.gravity,
                                                                       self.alpha,
                                                                       self.n,
                                                                       self.m,
                                                                       self.thetaR,
                                                                       self.thetaSR,
                                                                       self.Ks,
                                                                       u,
                                                                       m,
                                                                       dm,
                                                                       f,
                                                                       df,
                                                                       a,
                                                                       da)
##             a[:,0,0]=1.0e-3
##             da[:,0,0]=0.0
##             print 'm',m
##             print 'dm',dm
##             print 'a',a
##             print 'da',da
##             print 'f',f
##             print 'df',df
##             print 'phi',phi
##             print 'dphi',dphi
    coefficients[test] = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                    gravity=dimensionless_gravity,
                                                                    density=dimensionless_density,
                                                                    thetaS=thetaS,
                                                                    thetaR=thetaR,
                                                                    alpha= dimensionless_alpha,
                                                                    n = mvg_n,
                                                                    m = mvg_m)
#     xTest = Numeric.array([i/1000.0 for i in range(1001)])
#     uTest = Numeric.array([(i/500.0 - 1.0) for i in range(1001)])
#     m = Numeric.zeros((1001,),Numeric.Float)
#     dm = Numeric.zeros((1001,),Numeric.Float)
#     f = Numeric.zeros((1001,1),Numeric.Float)
#     df = Numeric.zeros((1001,1),Numeric.Float)
#     a = Numeric.zeros((1001,1,1),Numeric.Float)
#     da = Numeric.zeros((1001,1,1),Numeric.Float)
#     phi = Numeric.zeros((1001,),Numeric.Float)
#     dphi = Numeric.zeros((1001,),Numeric.Float)
#     r = Numeric.zeros((1001,),Numeric.Float)
#     dr = Numeric.zeros((1001,),Numeric.Float)
#     coefficients[test].evaluate(0.0,
#                                 xTest,
#                                 uTest,
#                                 m,dm,
#                                 f,df,
#                                 a,da,
#                                 phi,dphi,
#                                 r,dr)
#     mPlot = Gnuplot.Gnuplot()
#     dmPlot = Gnuplot.Gnuplot()
#     fPlot = Gnuplot.Gnuplot()
#     dfPlot = Gnuplot.Gnuplot()
#     aPlot = Gnuplot.Gnuplot()
#     daPlot = Gnuplot.Gnuplot()
#     charPlot = Gnuplot.Gnuplot()
#     mPlot.plot(Gnuplot.Data(uTest,
#                             m,
#                             with='lines',
#                             title='m'))
#     dmPlot.plot(Gnuplot.Data(uTest,
#                             dm,
#                             with='lines',
#                             title='dm'))
#     fPlot.plot(Gnuplot.Data(uTest,
#                             f[:,0],
#                             with='lines',
#                             title='f'))
#     dfPlot.plot(Gnuplot.Data(uTest,
#                             df[:,0],
#                             with='lines',
#                             title='df'))
#     aPlot.plot(Gnuplot.Data(uTest,
#                             a[:,0,0],
#                             with='lines',
#                             title='a'))
#     daPlot.plot(Gnuplot.Data(uTest,
#                             da[:,0,0],
#                             with='lines',
#                             title='da'))
#     charPlot.plot(Gnuplot.Data(uTest,
#                                Numeric.array([dfi/(dmi+1.0e-16) for  dfi,dmi in zip(df[:,0],dm)]),
#                                with='lines',
#                                title='df/dm'))
    analyticalSolution[test] = None
    class LinearIC:
        def __init__(self,dbc):
            self.uLeft = dbc(0.0)(0.0,0.0)
            self.uRight = dbc(1.0)(0.0,0.0)
        def uOfXT(self,x,t):
            return x[0]*self.uRight + (1.0-x[0])*self.uLeft
    getInitialConditions[test]=LinearIC(getDBC_Richards)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # RE-unsat
    test='NonlinearAD-RE-MvG-unsat'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    def getDBC_Richards_unsat(x):
        if x == 0.0:
            return lambda x,t: -0.5/lengthScale  # -0.5 m of  pressure head
        if x == 1.0:
            return lambda x,t: -10.0/lengthScale # -10.0 m of pressure head 
    getDirichletConditions[test]=getDBC_Richards_unsat
    coefficients[test] = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                    gravity=dimensionless_gravity,
                                                                    density=dimensionless_density,
                                                                    thetaS=thetaS,
                                                                    thetaR=thetaR,
                                                                    alpha= dimensionless_alpha,
                                                                    n = mvg_n,
                                                                    m = mvg_m)
    analyticalSolution[test] = None
    getInitialConditions[test]=LinearIC(getDBC_Richards_unsat)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # RE-nearly-sat
    test='NonlinearAD-RE-MvG-nearly-sat'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    def getDBC_Richards_nearly_sat(x):
        if x == 0.0:
            return lambda x,t:  0.1/lengthScale # 0.1 m of  pressure head
        if x == 1.0:
            return lambda x,t: -0.5/lengthScale # -0.5 m of pressure head 
    getDirichletConditions[test]=getDBC_Richards_nearly_sat
    coefficients[test] = ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                    gravity=dimensionless_gravity,
                                                                    density=dimensionless_density,
                                                                    thetaS=thetaS,
                                                                    thetaR=thetaR,
                                                                    alpha= dimensionless_alpha,
                                                                    n = mvg_n,
                                                                    m = mvg_m)
    analyticalSolution[test] = None
    getInitialConditions[test]=LinearIC(getDBC_Richards_nearly_sat)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    #
    # saturation based Richards' with Mualem-van Genuchten
    #
    #just need to redefine f and phi from last problem
    re_sat_eps=1.0e-2
    re_unsat_eps=1.0e-4
    test='NonlinearAD-satRE-MvG'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    def getDBC_satRichards(x):
        if x == 0.0:
            return lambda x,t: 1.0-re_sat_eps
        if x == 1.0:
            return lambda x,t: re_unsat_eps
    getDirichletConditions[test]=getDBC_satRichards
    class ConservativeSatRichardsMualemVanGenuchten:
        from transportCoefficients import  conservativeSatRichardsMualemVanGenuchtenHomEvaluate
        def __init__(self,
                     hydraulicConductivity,
                     gravity,
                     density,
                     thetaS,
                     thetaR,
                     alpha,
                     n,
                     m):
            self.Ks = hydraulicConductivity
            self.gravity=gravity
            self.rho = density
            self.thetaS = thetaS
            self.thetaR = thetaR
            self.thetaSR = thetaS - thetaR
            self.alpha = alpha
            self.n = n
            self.m = m
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr): 
            r.flat[:]=0.0
            dr.flat[:]=0.0
            self.conservativeSatRichardsMualemVanGenuchtenHomEvaluate(u.shape[0],
                                                                       f.shape[1],
                                                                       self.rho,
                                                                       self.gravity,
                                                                       self.alpha,
                                                                       self.n,
                                                                       self.m,
                                                                       self.thetaR,
                                                                       self.thetaSR,
                                                                       self.Ks,
                                                                       u,
                                                                       m,
                                                                       dm,
                                                                       f,
                                                                       df,
                                                                       a,
                                                                       da,
                                                                       phi,
                                                                       dphi)
    coefficients[test] = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                   gravity=dimensionless_gravity,
                                                                   density=dimensionless_density,
                                                                   thetaS=thetaS,
                                                                   thetaR=thetaR,
                                                                   alpha= dimensionless_alpha,
                                                                   n = mvg_n,
                                                                   m = mvg_m)
 ##    xTest = Numeric.array([i/1000.0 for i in range(1001)])
##     uTest = Numeric.array([i/1000.0 for i in range(1001)])
##     m = Numeric.zeros((1001,),Numeric.Float)
##     dm = Numeric.zeros((1001,),Numeric.Float)
##     f = Numeric.zeros((1001,1),Numeric.Float)
##     df = Numeric.zeros((1001,1),Numeric.Float)
##     a = Numeric.zeros((1001,1,1),Numeric.Float)
##     da = Numeric.zeros((1001,1,1),Numeric.Float)
##     phi = Numeric.zeros((1001,),Numeric.Float)
##     dphi = Numeric.zeros((1001,),Numeric.Float)
##     r = Numeric.zeros((1001,),Numeric.Float)
##     dr = Numeric.zeros((1001,),Numeric.Float)
##     coefficients[test].evaluate(0.0,
##                                 xTest,
##                                 uTest,
##                                 m,dm,
##                                 f,df,
##                                 a,da,
##                                 phi,dphi,
##                                 r,dr)
##     mPlot = Gnuplot.Gnuplot()
##     dmPlot = Gnuplot.Gnuplot()
##     phiPlot = Gnuplot.Gnuplot()
##     dphiPlot = Gnuplot.Gnuplot()
##     fPlot = Gnuplot.Gnuplot()
##     dfPlot = Gnuplot.Gnuplot()
##     aPlot = Gnuplot.Gnuplot()
##     daPlot = Gnuplot.Gnuplot()
##     charPlot = Gnuplot.Gnuplot()
##     mPlot.plot(Gnuplot.Data(uTest,
##                             m,
##                             with='lines',
##                             title='m'))
##     dmPlot.plot(Gnuplot.Data(uTest,
##                             dm,
##                             with='lines',
##                             title='dm'))
##     phiPlot.plot(Gnuplot.Data(uTest,
##                               phi,
##                             with='lines',
##                             title='phi'))
##     dphiPlot.plot(Gnuplot.Data(uTest,
##                                dphi,
##                             with='lines',
##                             title='dphi'))
##     fPlot.plot(Gnuplot.Data(uTest,
##                             f[:,0],
##                             with='lines',
##                             title='f'))
##     dfPlot.plot(Gnuplot.Data(uTest,
##                             df[:,0],
##                             with='lines',
##                             title='df'))
##     aPlot.plot(Gnuplot.Data(uTest,
##                             a[:,0,0],
##                             with='lines',
##                             title='a'))
##     daPlot.plot(Gnuplot.Data(uTest,
##                             da[:,0,0],
##                             with='lines',
##                             title='da'))
##     charPlot.plot(Gnuplot.Data(uTest,
##                                Numeric.array([dfi/(dmi+1.0e-16) for  dfi,dmi in zip(df[:,0],dm)]),
##                                with='lines',
##                                title='df/dm'))
    analyticalSolution[test] = None
    getInitialConditions[test]=LinearIC(getDBC_satRichards)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    # satRE-unsat
    test='NonlinearAD-satRE-MvG-unsat'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    def getDBC_satRichards_unsat(x):
        if x == 0.0:
            return lambda x,t: 0.5
        if x == 1.0:
            return lambda x,t: re_unsat_eps
    getDirichletConditions[test]=getDBC_satRichards_unsat
    coefficients[test] = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                   gravity=dimensionless_gravity,
                                                                   density=dimensionless_density,
                                                                   thetaS=thetaS,
                                                                   thetaR=thetaR,
                                                                   alpha= dimensionless_alpha,
                                                                   n = mvg_n,
                                                                   m = mvg_m)
    analyticalSolution[test] = None
    getInitialConditions[test]=LinearIC(getDBC_satRichards_unsat)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    # sat RE-nearly-sat
    test='NonlinearAD-satRE-MvG-nearly-sat'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    nd[test]=1
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    def getDBC_satRichards_nearly_sat(x):
        if x == 0.0:
            return lambda x,t:  1.0-re_sat_eps
        if x == 1.0:
            return lambda x,t:  0.5
    getDirichletConditions[test]=getDBC_satRichards_nearly_sat
    coefficients[test] = ConservativeSatRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                   gravity=dimensionless_gravity,
                                                                   density=dimensionless_density,
                                                                   thetaS=thetaS,
                                                                   thetaR=thetaR,
                                                                   alpha= dimensionless_alpha,
                                                                   n = mvg_n,
                                                                   m = mvg_m)
    analyticalSolution[test] = None
    getInitialConditions[test]=LinearIC(getDBC_satRichards_nearly_sat)
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    # (vu - au_x) + cu + d= 0; u(0) = 1; u(1) = 0
    #
    # u = sin(n pi (x + 1/(2n)))
    #need to redefine phi and r
    test='LinearADR-Sine'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='noFlow'
    nd[test]=1
    N=3.0
    analyticalSolution[test] = LinearADR_Sine(
        b=B0,
        a=A0,
        c=C0,
        omega=Numeric.array([(2.0*N+1.0)*pi/2.0]),
        omega0=pi/2.0)
    class LinearAD_UserR_Coefficients(ScalarTransportCoefficients):
        from transportCoefficients import linearADR_ConstantCoefficientsEvaluate
        def __init__(self,
                     M,
                     A,
                     B,
                     rFunc=LinearADR_Sine()):
            self.M = M
            self.A = A
            self.B = B
            self.rFunc=rFunc
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if  self.useC:
                self.linearADR_ConstantCoefficientsEvaluate(x.shape[0],
                                                            f.shape[1],
                                                            self.M,
                                                            self.A,
                                                            self.B,
                                                            0.0, #C=0
                                                            t,
                                                            x,
                                                            u,
                                                            m,dm,
                                                            f,df,
                                                            a,da,
                                                            phi,dphi,
                                                            r,dr)
            else:
                m[:] = u
                m *= self.M
                dm[:] = self.M

                f[:]=u[:,Numeric.newaxis]
                f*=self.B
                df[:] = self.B

                a[:,:] = self.A
                da[:,:] = 0.0

                phi[:]=u
                dphi[:]=1.0
            for i in range(u.shape[0]):
                r[i] = self.rFunc.rOfUX(u[i],x[i])
                dr[i] = self.rFunc.drOfUX(u[i],x[i])
    coefficients[test]=LinearAD_UserR_Coefficients(
        M=0.0,
        A=A0,
        B=B0,
        rFunc=analyticalSolution[test])
    getInitialConditions[test]=None
    fullNewtonFlag[test] = False
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = 'linear'
    #
    # u_t + (bu - a u_x)_x = 0; u(0) = 0
    #
    test='LinearAD-DiracIC'
    testProblems.append(test)
    timeIntegration[test] = ForwardEuler
    #timeIntegration[test] = OuterTheta
    #timeIntegration[test] = IMEX
    #timeIntegration[test] = BackwardEuler
    nd[test]=1
    def getDBC_hom(x):
        if x == 0.0:
            return lambda x,t: 0.0
        #if x == 1.0:
        #    return lambda x,t: 0.0
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = LinearADR_ConstantCoefficients(M=1.0,A=A0,B=B0,C=0.0)
    analyticalSolution[test] = LinearAD_DiracIC(b=B0,a=a0,tStart=0.25)
    T[test]=1.0
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 0
    #
    test='BurgersDiracIC'
    testProblems.append(test)
    timeIntegration[test]=BackwardEuler
    nd[test]=1
    def getDBC_hom(x):
        if x == 0.0:
            return lambda x,t: 0.0
#          if x == 1.0:
#              return lambda x,t: 0.0
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=0.0,q=2)
    analyticalSolution[test] = None
    T[test]=0.3
    getInitialConditions[test] = LinearAD_DiracIC(b=B0,a=a0,tStart=0.25)
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 0
    #
    test='Degenerate_BurgersDiracIC'
    testProblems.append(test)
    timeIntegration[test]=BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=0.0,q=2,r=6)
    analyticalSolution[test] = None
    T[test]=0.3
    getInitialConditions[test] = LinearAD_DiracIC(b=B0,a=a0,tStart=0.25)
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 0
    #
    test='Degenerate_Burgers2DiracIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    diff_t = 5.0
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=(diff_t+1.0)*A0,B=B0,C=0.0,q=2,r=1,t=diff_t)
    analyticalSolution[test] = None
    T[test]=0.3
    getInitialConditions[test] = LinearAD_DiracIC(b=B0,a=a0,tStart=0.25)
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu - a u_x)_x = 0; u(0) = 1
    #
    test='LinearAD_ShockIC'
    testProblems.append(test)
    timeIntegration[test] = ForwardEuler
    nd[test]=1
    class ShockIC:
        def uOfXT(self,x,t):
            if x <= 0.0:
                return 1.0
            else:
                return 0.0
    def getDBC_out(x):
        if x == 0.0:
            return lambda x,t: 1.0
    getDirichletConditions[test]=getDBC_out
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = LinearADR_ConstantCoefficients(M=1.0,A=A0,B=B0,C=0.0)
    analyticalSolution[test] = None
    T[test]=0.5
    getInitialConditions[test] = ShockIC()
    fullNewtonFlag[test] = False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 1
    #
    test='Burgers-ShockIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=0.0,q=2)
    analyticalSolution[test] = None
    T[test]=0.5
    getInitialConditions[test] = ShockIC()
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 1
    #
    test='Degenerate-Burgers_ShockIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=0.0,q=2.0,r=6.0)
    analyticalSolution[test] = None
    T[test]=0.5
    getInitialConditions[test] = ShockIC()
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'nonlinear'
    coefficients[test].reaction = None
    #
    # u_t + (bu^2 - a u_x)_x = 0; u(0) = 1
    #
    test='Degenerate_Burgers2_ShockIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='outFlow'
    diff_t = 5.0
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=(diff_t+1.0)*A0,B=B0,C=0.0,q=2,r=1,t=diff_t)
    analyticalSolution[test] = None
    T[test]=0.5
    getInitialConditions[test] = ShockIC()
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u^2_t + (bu - a u_x)_x = 0; u(0) = 1
    #
    test='Degenerate_Burgers3_ShockIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=0.0,p=2)
    analyticalSolution[test] = None
    T[test]=0.3
    getInitialConditions[test] = ShockIC()
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'nonlinear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # Richards
    #
    test='RE-MvG-ShockIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    pondingPressure=0.01
    def getDBC_Richards_Shock(x):
        if x == 0.0:
            return lambda x,t: pondingPressure
        if x == 1.0:
            return lambda x,t: -1.0+1.0*dimensionless_gravity[0]*dimensionless_density
    getDirichletConditions[test]=getDBC_Richards_Shock
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] =  ConservativeHeadRichardsMualemVanGenuchten(hydraulicConductivity=dimensionless_conductivity,
                                                                    gravity=Numeric.array([dimensionless_gravity]),
                                                                    density=dimensionless_density,
                                                                    thetaS=thetaS,
                                                                    thetaR=thetaR,
                                                                    alpha= mvg_alpha*lengthScale,
                                                                    n = mvg_n,
                                                                    m = mvg_m)
    analyticalSolution[test] = None
    T[test]=1.0/(dimensionless_conductivity*( (1.0+pondingPressure)/0.1 + dimensionless_gravity[0]*dimensionless_density))
    class ShockIC_Richards:
        def uOfXT(self,x,t):
            return -1.0+x[0]*dimensionless_gravity[0]*dimensionless_density
    getInitialConditions[test] = ShockIC_Richards()
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'nonlinear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'nonlinear'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + (bu - a u_x)_x + cu = 0; u(0) = 1
    #
    test='LinearADR_DecayDiracIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = LinearADR_ConstantCoefficients(M=1.0,A=A0,B=B0,C=C0)
    analyticalSolution[test] = LinearADR_Decay_DiracIC(b=B0,a=a0,c=C0,tStart=0.25)
    T[test]=0.3
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = 'linear'
    #
    # u_t + (bu - a u_x)_x + cu^s= 0; u(0) = 1;
    #
    test='NonlinearADR_DecayDiracIC'
    testProblems.append(test)
    timeIntegration[test] = BackwardEuler
    nd[test]=1
    getDirichletConditions[test]=getDBC_hom
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test] = NonlinearADR_pqrst(M=1.0,A=A0,B=B0,C=C0,s=4)
    analyticalSolution[test] = NonlinearADR_Decay_DiracIC(b=B0,a=a0,c=C0,d=4,tStart=0.25)
    T[test]=0.3
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = 'nonlinear'
    #
    # u_t + \deld bu = 0  on (0,1)x(0,1) 
    #
    # u(x,y,0) = 1/4(1+cos pi X)(1+cox pi Y) X^2+Y^2 <= 1
    #
    # X = (x-1/2), Y = (y-1/2)
    # b = (Y,-X)
    test = 'RotatingCone2D'
    testProblems.append(test)
    timeIntegration[test]=BackwardEuler#SSPRKintegration#BackwardEuler
    def getNoDBC(x):
        pass
    getDirichletConditions[test]=getNoDBC
    fluxBoundaryConditions[test]='outFlow'
    nd[test]=2
    N=3.0
    class RotatingCone2D:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX = 0.25*sin(2*pi*t) + 0.5
            centerY = 0.25*cos(2*pi*t) + 0.5
            coneX = x[0] - centerX
            coneY = x[1] - centerY
            if math.sqrt(coneX**2 + coneY**2) < self.radius:
                return 0.25*(1.0 + math.cos(math.pi*coneX/self.radius))*(1.0+math.cos(math.pi*coneY/self.radius))
            else:
                return 0.0
    analyticalSolution[test] = RotatingCone2D(1.0/8.0)
    class UnitSquareRotation(ScalarTransportCoefficients):
	from transportCoefficients import unitSquareRotationEvaluate
        def __init__(self):
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.unitSquareRotationEvaluate(u.shape[0],
                                                f.shape[1],
                                                x,
                                                u,
                                                m,dm,
                                                f,df)
            else:
                m[:] = u
                dm[:] = 1.0
                for i in range(u.shape[0]):
                    vx = 2*math.pi*(x[i][1] - 0.5)
                    vy = 2*math.pi*(0.5     - x[i][0]) 
                    f[i][0] = vx*u[i]
                    f[i][1] = vy*u[i]
                    df[i][0] = vx
                    df[i][1] = vy
    coefficients[test]=UnitSquareRotation()
    T[test]=1.0
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = None
    coefficients[test].potential = None
    coefficients[test].reaction = None
    #
    # u_t + \deld (bu - a \grad u)= 0  on (0,1)
    #
    # initial and boundary conditions given by u(x,y,t) below
    test = 'RotatingPulse2D-BE'
    testProblems.append(test)
    timeIntegration[test]=BackwardEuler
    nd[test]=2
    N=3.0
    class RotatingPulse2D:
        def __init__(self,sigma,xc,yc,a):
            self.two_sigma_squared = 2.0*sigma**2
            self.XC = xc
            self.YC = yc
            self.a = a
        def uOfXT(self,x,t):
            centerX = x[0] - 0.5
            centerY = x[1] - 0.5
            barX = centerX*cos(4.0*t) + centerY*sin(4.0*t)
            barY = -centerX*sin(4.0*t) + centerY*cos(4.0*t)
            den = self.two_sigma_squared + 4.0*self.a*t
            return self.two_sigma_squared *exp(- ((barX-self.XC)**2 + (barY-self.YC)**2)/den)/den
    pulseSol=analyticalSolution[test] = RotatingPulse2D(sigma=0.0447,xc=-0.25,yc=0.0,a=a0)
    def getRotatingPulseDBC2D(x):
	pass
#        if (x[X] == 0.0
#            or x[X] == 1.0
#            or x[Y] == 0.0
#            or x[Y] == 1.0):
#            return pulseSol.uOfXT
    getDirichletConditions[test]=getRotatingPulseDBC2D
    fluxBoundaryConditions[test]='outFlow'
    class RotatingPulseVel(ScalarTransportCoefficients):
        from transportCoefficients import rotatingPulseVelEvaluate
        def __init__(self,a):
            self.a=a
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.rotatingPulseVelEvaluate(u.shape[0],
                                              f.shape[1],
                                              self.a,
                                              x,
                                              u,
                                              m,dm,
                                              f,df,
                                              a,da,
                                              phi,dphi)
            else:
                m[:] = u
                dm[:] = 1.0
            
            for ix in range(a.shape[2]):
                a[:,ix,ix] = self.a
            da[:,:,:] = 0.0
            
	    phi[:]=u
	    dphi[:]=1.0
	    
            for i in range(u.shape[0]):
                centerX = x[i][0] - 0.5
                centerY = x[i][1] - 0.5
                vx = -4.0*centerY
                vy = 4.0*centerX
                f[i][0] = vx*u[i]
                f[i][1] = vy*u[i]
                df[i][0] = vx
                df[i][1] = vy
    coefficients[test]=RotatingPulseVel(a=a0)
    T[test]=1.0
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True#False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    pulseProblem = 'RotatingPulse2D-BE'
    #same problem with crank-nicholson
    test = 'RotatingPulse2D-CN'
    testProblems.append(test)
    timeIntegration[test]=OuterTheta
    nd[test]=nd[pulseProblem]
    analyticalSolution[test] = analyticalSolution[pulseProblem]
    getDirichletConditions[test]=getDirichletConditions[pulseProblem] 
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients[pulseProblem]
    T[test]=T[pulseProblem]
    getInitialConditions[test] = getInitialConditions[pulseProblem]
    fullNewtonFlag[test] = fullNewtonFlag[pulseProblem]
    coefficients[test].mass = coefficients[pulseProblem].mass
    coefficients[test].advection = coefficients[pulseProblem].advection
    coefficients[test].diffusion = coefficients[pulseProblem].diffusion
    coefficients[test].potential = coefficients[pulseProblem].potential 
    coefficients[test].reaction  = coefficients[pulseProblem].reaction
    #same problem with forward euler
    test = 'RotatingPulse2D-FE'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    nd[test]=nd[pulseProblem]
    analyticalSolution[test] = analyticalSolution[pulseProblem]
    getDirichletConditions[test]=getDirichletConditions[pulseProblem] 
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients[pulseProblem]
    T[test]=T[pulseProblem]
    getInitialConditions[test] = getInitialConditions[pulseProblem]
    fullNewtonFlag[test] = fullNewtonFlag[pulseProblem]
    coefficients[test].mass = coefficients[pulseProblem].mass
    coefficients[test].advection = coefficients[pulseProblem].advection
    coefficients[test].diffusion = coefficients[pulseProblem].diffusion
    coefficients[test].potential = coefficients[pulseProblem].potential 
    coefficients[test].reaction  = coefficients[pulseProblem].reaction
    #same problem with forward euler
    test = 'RotatingPulse2D-IMEX'
    testProblems.append(test)
    timeIntegration[test]=IMEX
    nd[test]=nd[pulseProblem]
    analyticalSolution[test] = analyticalSolution[pulseProblem]
    getDirichletConditions[test]=getDirichletConditions[pulseProblem] 
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients[pulseProblem]
    T[test]=T[pulseProblem]
    getInitialConditions[test] = getInitialConditions[pulseProblem]
    fullNewtonFlag[test] = fullNewtonFlag[pulseProblem]
    coefficients[test].mass = coefficients[pulseProblem].mass
    coefficients[test].advection = coefficients[pulseProblem].advection
    coefficients[test].diffusion = coefficients[pulseProblem].diffusion
    coefficients[test].potential = coefficients[pulseProblem].potential 
    coefficients[test].reaction  = coefficients[pulseProblem].reaction
    #
    # u_t + \deld (bu - a \grad u)= 0  on (0,1)
    #
    # initial and boundary conditions given by u(x,y,t) below
    test = 'DiscontinuousRotatingPulse2D'
    testProblems.append(test)
    #timeIntegration[test]=BackwardEuler
    timeIntegration[test]=ForwardEuler
    #timeIntegration[test]=OuterTheta
    nd[test]=2
    N=3.0
    getDirichletConditions[test]=getRotatingPulseDBC2D
    fluxBoundaryConditions[test]='outFlow'
    analyticalSolution[test] = pulseSol
    class DiscontinuousRotatingPulseVel(ScalarTransportCoefficients):
        from transportCoefficients import disRotatingPulseVelEvaluate
        def __init__(self,a):
            self.a=a
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.disRotatingPulseVelEvaluate(u.shape[0],
                                                 f.shape[1],
                                                 self.a,
                                                 x,
                                                 u,
                                                 m,dm,
                                                 f,df,
                                                 a,da,
                                                 phi,dphi)
            else:
                m[:] = u
                dm[:] = 1.0
            
            for ix in range(a.shape[2]):
                a[:,ix,ix] = self.a
            da[:,:,:] = 0.0
            phi[:]=u
	    dphi[:]=1.0
            for i in range(u.shape[0]):
                centerX = x[i][0] - 0.5
                centerY = x[i][1] - 0.5
                vx = -4.0*centerY
                vy = 4.0*centerX
                f[i][0] = vx*u[i]
                f[i][1] = vy*u[i]
                df[i][0] = vx
                df[i][1] = vy
                if sqrt(centerX**2 + centerY**2) < 0.25:
                    f[i][0] *= 0.001
                    f[i][1] *=  0.001
                    df[i][0] *= 0.001
                    df[i][1] *= 0.001
    coefficients[test]=DiscontinuousRotatingPulseVel(a=a0)
    T[test]=0.3
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + \deld (bu - a \grad u)= 0  on (0,1)
    #
    # initial and boundary conditions given by u(x,y,t) below
    test = 'DiscontinuousSlug2D'
    testProblems.append(test)
    #timeIntegration[test]=BackwardEuler
    timeIntegration[test]=ForwardEuler
    #timeIntegration[test]=OuterTheta
    nd[test]=2
    N=3.0
    def Slug2D(x):
        if x[0] == 0.0:
            return lambda x,t: 1.0
    getDirichletConditions[test]=Slug2D
    fluxBoundaryConditions[test]='outFlow'
    analyticalSolution[test] = None
    class DiscontinuousVel(ScalarTransportCoefficients):
        from transportCoefficients import disVelEvaluate
        def __init__(self,a):
            self.a=a
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.disVelEvaluate(u.shape[0],
                                    f.shape[1],
                                    self.a,
                                    x,
                                    u,
                                    m,dm,
                                    f,df,
                                    a,da,
                                    phi,dphi)
            else:
                m[:] = u
                dm[:] = 1.0
            
            for ix in range(a.shape[2]):
                a[:,ix,ix] = self.a
            da[:,:,:] = 0.0
            phi[:]=0.0
	    dphi[:]=1.0
            for i in range(u.shape[0]):
                f[i,0] = u[i]
                f[i,1] = 0.0
                df[i,0] = 1.0
                df[i,1] = 0.0
                if x[i,1] > 0.5:
                    f[i,0] *= 0.25
                    df[i,0] *= 0.25
    coefficients[test]=DiscontinuousVel(a=a0)
    T[test]=0.5
    class Slug2DIC:
        def uOfXT(self,x,t):
            if x[0] <= 0.0:
                return 1.0
            else:
                return 0.0
    getInitialConditions[test] = Slug2DIC()
    fullNewtonFlag[test] = False
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # u_t + \deld (bu^2 - a \grad u) = 0  on (0,1)x(0,1)
    #
    # intial and boundary conditions given by u(x,y,t) below
    test = 'BurgersDiagonal2D-BE'
    testProblems.append(test)
    timeIntegration[test]=BackwardEuler
    nd[test]=2
    N=3.0
    class BurgersDiagonal2D:
        def __init__(self,a):
            self.a = a
        def uOfXT(self,x,t):
            e_arg = (x[0] + x[1] - (t+0.25))/self.a
            if abs(e_arg) < 5.0e2:
                return 1.0/(1.0 + exp(e_arg))
            elif e_arg > 0.0:
                return 0.0
            else:
                return 1.0
                
    burgersSol=analyticalSolution[test] = BurgersDiagonal2D(a=a0)
    def getBurgersDiagonalDBC2D(x):
        if (x[X] == 0.0
            or x[X] == 1.0
            or x[Y] == 0.0
            or x[Y] == 1.0):
            return burgersSol.uOfXT
    getDirichletConditions[test]=getBurgersDiagonalDBC2D
    fluxBoundaryConditions[test]='outFlow'
    class BurgersDiagonalVel(ScalarTransportCoefficients):
        from transportCoefficients import burgersDiagonalVelEvaluate
        def __init__(self,a):
            self.a=a
            self.useC=True
        def evaluate(self,
                     t,
                     x,
                     u,
                     m,dm,
                     f,df,
                     a,da,
                     phi,dphi,
                     r,dr):
            if self.useC:
                self.burgersDiagonalVelEvaluate(u.shape[0],
                                                f.shape[1],
                                                self.a,
                                                u,
                                                m,dm,
                                                f,df,
                                                a,da,
                                                phi,dphi)
            else:
                m[:] = u
                dm[:] = 1.0

            for ix in range(a.shape[2]):
                a[:,ix,ix] = self.a
            da[:,:,:] = 0.0
            phi[:]=u
	    dphi[:]=1.0
            maxu0 = Numeric.maximum(u,0.0)
            Numeric.power(maxu0,2.0,f[:,0])
            for ix in range(1,f.shape[1]):
                f[:,ix]=f[:,0]
            Numeric.multiply(f,0.5,f)
            for ix in range(f.shape[1]):
                df[:,ix]=maxu0
    #coefficients[test]=BurgersDiagonalVel(a=a0)
    coefficients[test]=NonlinearADR_pqrst(M=1.0,A=A0,B=Numeric.array([1.0,1.0,1.0]),C=0.0,q=2)
    T[test]=0.5
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # do FE on same problem
    #
    test = 'BurgersDiagonal2D-FE'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    nd[test]=nd['BurgersDiagonal2D-BE']
    analyticalSolution[test] = analyticalSolution['BurgersDiagonal2D-BE']
    getDirichletConditions[test]=getBurgersDiagonalDBC2D
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients['BurgersDiagonal2D-BE']
    T[test]=T['BurgersDiagonal2D-BE']
    getInitialConditions[test] = getInitialConditions['BurgersDiagonal2D-BE']
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # do CN on same problem
    #
    test = 'BurgersDiagonal2D-CN'
    testProblems.append(test)
    timeIntegration[test]=OuterTheta
    nd[test]=nd['BurgersDiagonal2D-BE']
    analyticalSolution[test] = analyticalSolution['BurgersDiagonal2D-BE']
    getDirichletConditions[test]=getBurgersDiagonalDBC2D
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients['BurgersDiagonal2D-BE']
    T[test]=T['BurgersDiagonal2D-BE']
    getInitialConditions[test] = getInitialConditions['BurgersDiagonal2D-BE']
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # do IMEX on same problem
    #
    test = 'BurgersDiagonal2D-IMEX'
    testProblems.append(test)
    timeIntegration[test]=IMEX
    nd[test]=nd['BurgersDiagonal2D-BE']
    analyticalSolution[test] = analyticalSolution['BurgersDiagonal2D-BE']
    getDirichletConditions[test]=getBurgersDiagonalDBC2D
    fluxBoundaryConditions[test]='outFlow'
    coefficients[test]=coefficients['BurgersDiagonal2D-BE']
    T[test]=T['BurgersDiagonal2D-BE']
    getInitialConditions[test] = getInitialConditions['BurgersDiagonal2D-BE']
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].advection = 'nonlinear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = None
    #
    # Owen's problem in 2D
    #
    #need to redefine phi and r
    test='PoissonsEquation2D'
    testProblems.append(test)
    timeIntegration[test]=NoIntegration
    def getHomogeneousDBC2D(x):
        if x[X] == 0.0:
            if x[Y] <= 0.5:
                return lambda x,t: 1.0
            else:
                return lambda x,t: 0.0
        #elif x[X] == 1.0:
        #    if x[Y] >= 0.5:
        #        return lambda x,t: 0.0
        elif x[Y] == 0.0:
            if x[X] <= 0.5:
                return lambda x,t: 1.0
            else:
                return lambda x,t: 0.0
        #elif x[Y] == 1.0:
        #    if x[X] >= 0.5:
        #        return lambda x,t: 0.0
    getDirichletConditions[test]=getHomogeneousDBC2D
    #fluxBoundaryConditions[test]='noFlow'
    fluxBoundaryConditions[test]='outFlow'
    nd[test]=2
    N=3.0
    analyticalSolution[test] = PoissonsEquation(10,nd=nd[test])
    coefficients[test]=LinearAD_UserR_Coefficients(
        M=0.0,
        A=Numeric.array([[0.0,0.0],[0.0,0.0]]),
        B=Numeric.array([100000.0,100000.0]),
#        B=Numeric.zeros((2,),Numeric.Float),
        rFunc=analyticalSolution['PoissonsEquation2D'])
    getInitialConditions[test]=None
    fullNewtonFlag[test] = True
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = 'linear'
    #
    # Owen's problem
    #
    #need to redefine phi and r
    test='PoissonsEquation3D'
    testProblems.append(test)
    timeIntegration[test] = NoIntegration
    def getHomogeneousDBC3D(x):
        if x[X] == 0.0:
            return lambda x,t: 0.0
        elif x[X] == 1.0:
            return lambda x,t: 0.0
        elif x[Y] == 0.0:
            return lambda x,t: 0.0
        elif x[Y] == 1.0:
            return lambda x,t: 0.0
        elif x[Z] == 0.0:
            return lambda x,t: 0.0
        elif x[Z] == 1.0:
            return lambda x,t: 0.0
    getDirichletConditions[test]=getHomogeneousDBC3D
    fluxBoundaryConditions[test]='noFlow'
    nd[test]=3
    N=3.0
    analyticalSolution[test] = PoissonsEquation(10)
    coefficients[test]=LinearAD_UserR_Coefficients(
        M=0.0,
        A=Numeric.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]),
        B=Numeric.zeros((3,),Numeric.Float),
        rFunc=analyticalSolution['PoissonsEquation3D'])
    getInitialConditions[test]=None
    fullNewtonFlag[test] = False
    coefficients[test].mass = None
    coefficients[test].advection = 'linear'
    coefficients[test].diffusion = 'constant'
    coefficients[test].potential = 'linear'
    coefficients[test].reaction = 'linear'
#          #
#          # two-phase flow (saturation equation)
#          #
#          test='TwoPhase1D'
#          testProblems.append(test) isTimeDependent[test] = False
#          #use boundary conditions from other 1D
#          getDirichletConditions[test]=getDBC
#          nd[test]=1
#          N=3.0
#          analyticalSolution[test] = None
#          a0 =1.0e-4
#          A = ETen(EVec(a0,0.0,0.0),EVec(0.0,a0,0.0),EVec(0.0,0.0,a0))
#          DA = ZEROVEC
#          b0=1.0
#          B0=EVec(b0,0.0,0.0)
#          def m(u,x,t):
#              omega = 1.0
#              rhow = 1.0
#              s = u
#              m = omega*rhow*s
#              dm = omega*rhow
#              return (m,dm)
#          def a(u,x,t):
#              return (A0,dA0)
#          def phi(u,x,t):
#              return (u,1.0)
#          def f(u,x,t):
#              s=u
#              pc = u
#              qt = B0
#              gradz = ZEROVEC
#  #              qt = ZEROVEC
#  #              gradz = B0
#              K = 1.0
#              rhon = 0.5
#              rhow = 1.0
#              muw = 1.0
#              mun = 1.0
#              b = rhon/rhow
#              krw = s**2
#              dkrw = 2.0*s
#              krn = 0.5*(1.0-s)**2
#              dkrn = (s - 1.0)
#              lambdaw = rhow*krw/muw
#              dlambdaw = rhow*dkrw/muw
#              lambdan = rhon*krn/mun
#              dlambdan = rhon*dkrn/mun
#              lambdat = lambdaw + lambdan
#              dlambdat = dlambdaw + dlambdan
#              fn = lambdan/lambdat
#              dfn = (lambdat*dlambdan - lambdan*dlambdat)/(lambdat**2)
#              fw = lambdaw/lambdat
#              dfw = (lambdat*dlambdaw - lambdaw*dlambdat)/(lambdat**2)
#              f = fw*qt+K*lambdaw*fn*(b*rhon - rhow)*gradz
#              df = dfw*qt + (dlambdaw*fn + lambdaw*dfn)*K*(b*rhon-rhow)*gradz
#              return (f,df)
#          def r(u,x,t):
#              R = 0.0
#              DR = 0.0
#              return (R,DR)
#          mFunc[test]=m
#          aFunc[test]=a
#          phiFunc[test]=phi
#          fFunc[test]=f
#          rFunc[test]=r
#          #loop over test problems and solve
    for nTest,test in enumerate(testProblems):
        print nTest,test
    #
    #uncomment  next line to test performance 
    #
    #for test in testProblems[11:12]:
    #
    #uncomment next line to verify discretization
    #
    #for test in [testProblems[0],testProblems[1],testProblems[5],testProblems[19],testProblems[26],testProblems[28]]:
    #
    #
    #for test in [testProblems[5],testProblems[28]]:#,testProblems[29]]:
    #
    #all the  tests
    #
    #for test in testProblems[25:26]:
    for test in testProblems[1:2]:
        report.write('\section{'+test+'}\n')
        tolFac = 0.01
        atol = 1.0e-8
        linTolFac = 0.001
        runCFL = 0.01
        DG = False
	polyOrder = 1
        tOrder = 1
        if timeIntegration[test] == SSPRKintegration:
            tOrder=3
            nStages = tOrder
        else:
            tOrder = 1
            nStages = 1
        if DG:
	    if polyOrder == 1:
		FemSpaceType = DG_AffineLinearOnSimplexWithNodalBasis
	    elif polyOrder == 2:
		FemSpaceType = DG_AffineQuadraticOnSimplexWithNodalBasis
	    else:
		print "polyOrder should be 1 or 2"
		return
            conservativeFlux = None
            numericalFlux = True
            stabilization=None
            shockCapturing=None
            shockCapturingDiffusion=0.0
        else:
	    if polyOrder == 1:
		FemSpaceType = C0_AffineLinearOnSimplexWithNodalBasis
	    elif polyOrder == 2:
		FemSpaceType = C0_AffineQuadraticOnSimplexWithNodalBasis
	    else:
		print "polyOrder should be 1 or 2"
		return
            conservativeFlux = None#'pwl'
            numericalFlux = None
            stabilization=None#'2'
            shockCapturing= None#'2'# #Remember that shock capturing makes  the problem nonlinear or at least makes diffusion time-dependent
            shockCapturingDiffusion = 0.9
        quadratureOrder=3
        massLumping=False#True
        preSmooths = 2
        postSmooths = 2
        cycles = 2
        nLevels=5
        computeSpaceTimeError=False
        nn=3 #number of nodes on the coarsest mesh
        print "Starting Test "+`test`
        print "Setting up quadrature"
        quadrature={}
        gq = SimplexGaussQuadrature(nd[test])
        #gq = SimplexLobattoQuadrature(nd[test])
        gq.setOrder(quadratureOrder)
        for integral in OneLevelScalarTransport.integralKeys:
            quadrature[integral] = gq 
        if stabilization != None:
            quadrature['stab'] = gq 
        if shockCapturing != None:
            quadrature['numDiff'] = gq
        if massLumping == True:
            quadrature['m'] = SimplexLobattoQuadrature(nd[test])
        elementBoundaryQuadrature={}
        ebgq = SimplexGaussQuadrature(nd[test]-1)
        ebgq.setOrder(quadratureOrder)
        for elementBoundaryIntegral in OneLevelScalarTransport.elementBoundaryIntegralKeys:
            elementBoundaryQuadrature[elementBoundaryIntegral] = ebgq 
        #
        #define the mesh hierarchy
        #
        print "Setting up MultilevelMesh"
        mlMesh = None
        mlMeshFileName = "tmp/mlMesh%dD.%d" % (nd[test],nLevels)
        try:
            mlMeshFile = open(mlMeshFileName,'rb')
            print "reading mesh"
            mlMesh = cPickle.load(mlMeshFile)
            print "done reading mesh"
        except:
            print "generating mesh"
            mlMeshFile = open(mlMeshFileName,'wb')
            if nd[test]==1:
                mlMesh = MultilevelEdgeMesh(nn,1,1,refinementLevels=nLevels)
            elif nd[test]==2:
                mlMesh = MultilevelTriangularMesh(nn,nn,1,
                                                  refinementLevels=nLevels)
            elif nd[test]==3:
                mlMesh = MultilevelTetrahedralMesh(nn,nn,nn,
                                                   refinementLevels=nLevels)
            cPickle.dump(mlMesh,mlMeshFile,protocol=cPickle.HIGHEST_PROTOCOL)
            print "done generating mesh"
        print "Setting up MultilevelScalarTransport"
        tolList=[]
        linTolList=[]
        for l in range(nLevels):
            mlMesh.meshList[l].computeGeometricInfo()
            tolList.append(tolFac*(mlMesh.meshList[l].h**2))
            linTolList.append(linTolFac*(mlMesh.meshList[l].h**2))
	#matType = 'dense'
	matType = 'csr'
        mlScalarTransport = MultilevelScalarTransport(
            nd[test],
            mlMesh,
            FemSpaceType,
            FemSpaceType,
            matType,
            getDirichletConditions[test],
            coefficients[test],
            quadrature,
            elementBoundaryQuadrature,
            fluxBoundaryConditions[test],
            stabilization,
            shockCapturing,
            shockCapturingDiffusion,
            conservativeFlux,
            numericalFlux,
            timeIntegration[test],
            tOrder)
	#setting up some storage for spatial error calculations
        utmp = []
        uqtmp = Numeric.zeros(mlScalarTransport.modelList[-1].q['u'].shape,Numeric.Float)
        asolq = Numeric.zeros(mlScalarTransport.modelList[-1].q['u'].shape,Numeric.Float)
        for l in range(nLevels):
            utmp.append(FiniteElementFunction(mlScalarTransport.modelList[l].trialSpace))
        print "Setting up LinearSolver"
        computeEigenvalues=False
        if computeEigenvalues:
            gevals = Gnuplot.Gnuplot()
            kplot = Gnuplot.Gnuplot()
        if matType == 'dense':
            multilevelLinearSolverType = levelLinearSolverType = 'DenseLU'
        else:
            multilevelLinearSolverType  = levelLinearSolverType = 'SparseLU'
            #multilevelLinearSolverType = 'NI'
            #levelLinearSolverType = 'MGM'
	(multilevelLinearSolver,directSolverFlag) = multilevelLinearSolverChooser(
	    linearOperatorList = mlScalarTransport.jacobianList,
	    multilevelLinearSolverType = multilevelLinearSolverType,
            computeSolverRates=True,
            printSolverInfo=False,
	    levelLinearSolverType = levelLinearSolverType,
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
	    smootherType = 'GaussSeidel',#'StarILU'
            computeSmootherRates=True,
            printSmootherInfo=False,
	    prolongList = mlScalarTransport.meshTransfers.prolongList,
	    restrictList = mlScalarTransport.meshTransfers.restrictList,
	    connectivityListList = [mlScalarTransport.modelList[l].freeNodeStarList for l in range(nLevels)],
	    relativeToleranceList = linTolList,
	    absoluteTolerance = min(linTolList),
	    solverMaxIts = 500,
	    cycles=3,
	    preSmooths=3,
	    postSmooths=3)
        print "Setting up NonlinearSolver"
	multilevelNonlinearSolver = multilevelNonlinearSolverChooser(
	    nonlinearOperatorList = mlScalarTransport.modelList,
	    jacobianList = mlScalarTransport.jacobianList,
            multilevelNonlinearSolverType = 'NLNI',#'Newton','FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            #multilevelNonlinearSolverType = 'Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeSolverRates=True,
            printSolverInfo=False,
	    relativeToleranceList = tolList,
	    absoluteTolerance = atol,
	    levelNonlinearSolverType='Newton',#'FAS','NLJacobi','NLGaussSeidel','NLJacobi'
            computeLevelSolverRates=True,
	    printLevelSolverInfo=False,
	    smootherType = None,#'NLJacobi','NLGaussSeidel','NLStarILU'
            computeSmootherRates=True,
	    printSmootherInfo=False,
	    preSmooths=3,
	    postSmooths=3,
	    cycles=3,
	    maxSolverIts=100,
	    prolong_bcList = mlScalarTransport.meshTransfers.prolong_bcList,
	    restrict_bcList = mlScalarTransport.meshTransfers.restrict_bcList,
	    restrict_bcSumList = mlScalarTransport.meshTransfers.restrict_bcSumList,
	    prolongList = mlScalarTransport.meshTransfers.prolongList,
	    restrictList = mlScalarTransport.meshTransfers.restrictList,
	    restrictionRowSumList = mlScalarTransport.meshTransfers.restrictSumList,
	    connectionListList=[mlScalarTransport.modelList[l].freeNodeStarList for l in range(nLevels)],
	    linearSolverList=multilevelLinearSolver.solverList,
	    linearDirectSolverFlag=directSolverFlag,
	    solverFullNewtonFlag=fullNewtonFlag[test],
	    levelSolverFullNewtonFlag=fullNewtonFlag[test],
	    smootherFullNewtonFlag=fullNewtonFlag[test])
        print "Running Solver"
	try:
	    os.mkdir('tmp/'+test)
	except:
	    pass
	os.chdir('tmp/'+test)
        if timeIntegration[test] == NoIntegration:
            T[test] = 1.0
        tn = 0.0
        nSteps = 0
        if getInitialConditions[test] != None:
            mlScalarTransport.setInitialConditions(getInitialConditions[test],tn)
#          nx=[]
#          ny=[]
#          x=[]
#          y=[]
#          aSol=[]
#          solPlot = Gnuplot.Gnuplot()
#          solPlot("set terminal x11")
#          aSolPlot = Gnuplot.Gnuplot()
#          aSolPlot("set terminal x11")
#          for l in range(nLevels):
#              if DG != True:
#                  if nd[test]==1:
#                      solPlot.title(test)
#                      aSolPlot.title(test)
#                      nap=101
#                      dxap=Numeric.array([1.0/(nap - 1.0),0.0,0.0])
#                      P = [(i*dxap) for i in range(nap)]
#                      Px = [x[0] for x in P]
#                      solPlot.plot(Gnuplot.Data(mlMesh.meshList[l].nodeArray[:,0],
#                                                mlScalarTransport.modelList[l].u.dof,
#                                                with='linespoints',
#                                                title='numerical solution (initial condition/guess)'))
#                      if analyticalSolution[test] != None:
#                          aSolPlot.plot(Gnuplot.Data(Px,
#                                                     [analyticalSolution[test].uOfXT(x,tn) for x in P],
#                                                     with='lines',
#                                                     title='analytical solution'))
#                      else:
#                          aSolPlot=solPlot
#                  elif nd[test]==2:
#                      nx.append((nn-1)*(2**l)+1)
#                      ny.append(nx[l])
#                      x.append(Numeric.arange(nx[l])/float(nx[l]-1))
#                      y.append(Numeric.arange(nx[l])/float(nx[l]-1))
#                      nSol = Numeric.reshape(mlScalarTransport.modelList[l].u.dof,(nx[l],ny[l]))
#                      solPlot('set parametric')
#                      solPlot('set data style lines')
#                      solPlot('set hidden')
#                      solPlot('set contour base')
#                      solPlot('set cntrparam levels incremental 0.1,0.1,1.0')
#                      solPlot.xlabel('x')
#                      solPlot.ylabel('y')
#                      solPlot.splot(Gnuplot.GridData(nSol,
#                                                     x[l],
#                                                     y[l],
#                                                     binary=0,
#                                                     inline=0
#                                                     ))
#                      if analyticalSolution[test] != None:
#                          aSol.append(Numeric.zeros((nx[l],ny[l]),Numeric.Float))
#                          for i in range(nx[l]):
#                              for j in range(ny[l]):
#                                  aSol[l][i,j] = analyticalSolution[test].uOfXT(Numeric.array([x[l][i],y[l][j],0.0]),tn)
#                          aSolPlot('set parametric')
#                          aSolPlot('set data style lines')
#                          aSolPlot('set hidden')
#                          aSolPlot('set contour base')
#                          aSolPlot('set cntrparam levels incremental 0.1,0.1,1.0')
#                          aSolPlot.xlabel('x')
#                          aSolPlot.ylabel('y')
#                          aSolPlot.splot(Gnuplot.GridData(aSol[l],
#                                                         x[l],
#                                                         y[l],
#                                                         binary=0,
#                                                         inline=0
#                                         ))
        mlScalarTransport.modelList[-1].timeIntegration.runCFL = runCFL
        if int(T[test]/mlMesh.meshList[-1].h) != 0.0:
            DTSET = 0.01*T[test]/int(T[test]/mlMesh.meshList[-1].h)
        else:
            DTSET = T[test]
        #DTSET=None
        eSpace = {}
        eSpaceTime = {}
        import sys
        tstring=None
        eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b'
#        mlScalarTransport.modelList[-1].trialSpace.writeFunctionGnuplot(mlScalarTransport.modelList[-1].u,'initial  conditions')
	if polyOrder == 1:
	    mlMesh.meshList[-1].writeMeshEnsight(test,test)
	elif polyOrder == 2:
	    mlScalarTransport.modelList[-1].trialSpace.writeMeshEnsight(test,test)
        mlScalarTransport.modelList[-1].u.name='u'
        #mlScalarTransport.modelList[-1].trialSpace.writeFunctionEnsight(mlScalarTransport.modelList[-1].u,test,append=False)
	timeValues = [tn]
        failedFlag=False
        while (tn < T[test] and not failedFlag==True):
            if timeIntegration[test] != NoIntegration:
                mlScalarTransport.chooseDT(DTSET)
                if nSteps == 0:
                    mlScalarTransport.initializeTimeIntegration()
                tn += mlScalarTransport.DT
                if tstring != None:
                    sys.stdout.write(eraseTime)
                else:
                    sys.stdout.write('T = %12.5e, tn = ' % T[test])
                tstring='%12.5e' % (tn,)
                sys.stdout.write(tstring)
                sys.stdout.flush()
                nSteps += 1
                testOut = test + ('%4.4i' % nSteps)
            else:
                mlScalarTransport.DT = 1.0
                tn=1.0
                nSteps +=1
                testOut = test
            for i in range(nStages):
                failedFlag=multilevelNonlinearSolver.solveMultilevel(uList   = 
                                                                     mlScalarTransport.uList,
                                                                     rList   = 
                                                                     mlScalarTransport.rList)
                mlScalarTransport.updateStage()
            #raw_input('Please press return to continue... \n')
            timeValues.append(tn)
#             for l in range(nLevels):
#                 mlScalarTransport.modelList[l].trialSpace.writeFunctionEnsight(mlScalarTransport.modelList[l].u,test,append=True)
#                 mlScalarTransport.modelList[l].trialSpace.writeFunctionGnuplot(mlScalarTransport.modelList[l].u,'solution')
            #mlScalarTransport.modelList[-1].trialSpace.writeFunctionEnsight(mlScalarTransport.modelList[-1].u,test,append=True)
            #mlScalarTransport.modelList[-1].trialSpace.writeFunctionGnuplot(mlScalarTransport.modelList[-1].u,'solution')
	    mlScalarTransport.updateTimeHistory()
            if conservativeFlux == 'pwc':
                mlScalarTransport.modelList[-1].getConservationFluxPWC()
            elif conservativeFlux == 'pwl':
                mlScalarTransport.modelList[-1].getConservationFluxPWL()
            #print mlScalarTransport.modelList[-1].e['conservationResidual']
#              elif numericalFlux != None:
#                  mlScalarTransport.modelList[-1].e['conservationResidual'].flat[:]=0.0
#                  for eN in range(mlScalarTransport.modelList[-1].mesh.nElements_global):
#                      for i in range(mlScalarTransport.modelList[-1].nDOF_element):
#                          mlScalarTransport.modelList[-1].e['conservationResidual'][eN]+=mlScalarTransport.modelList[-1].elementResidual[eN,i]
	    #if conservativeFlux == 'pwc' or conservativeFlux == 'pwl' or numericalFlux != None:
		#print "Max mass conservation error "+`max(abs(mlScalarTransport.modelList[-1].e['conservationResidual']))`
            #rnorm = wl2Norm(ni.resid ual(),h)
            #print "rnorm"+`rnorm`
#              for l in range(nLevels):
#                  if DG != True:
    #                     nSolPlot.hardcopy(testOut+'_sol.eps', eps=1,enhanced=1,color=1)
    #                      nResPlot = Gnuplot.Gnuplot()
    #                      nResPlot("set terminal x11")
    #                      nResPlot.title(testOut)
    #                      nResPlot.plot(
    #                          Gnuplot.Data(mlScalarTransport.rList[-1],
    #                                       with='linespoints',
    #                                       title='numerical residual'))
        #  #                  nResPlot.hardcopy(testOut+'_res.eps', eps=1,enhanced=1,color=1)
 #                     if nd[test]==1:
#                          solPlot.title(testOut)
#                          solPlot.plot(Gnuplot.Data(mlMesh.meshList[l].nodeArray[:,0],
#                                                    mlScalarTransport.modelList[l].u.dof,
#                                                    with='linespoints',
#                                                    title='numerical solution'))
#                          if analyticalSolution[test] != None:
#                              aSolPlot.title(testOut)
#                              aSolPlot.plot(Gnuplot.Data(Px,
#                                                         [analyticalSolution[test].uOfXT(x,tn) for x in P],
#                                                         with='lines',
#                                                         title='analytical solution'))
#                      elif nd[test]==2:
#                          nSol = Numeric.reshape(mlScalarTransport.modelList[l].u.dof,(nx[l],ny[l]))
#                          solPlot('set parametric')
#                          solPlot('set data style lines')
#                          solPlot('set hidden')
#                          solPlot('set contour base')
#                          solPlot('set cntrparam levels incremental 0.1,0.1,1.0')
#                          solPlot.xlabel('x')
#                          solPlot.ylabel('y')
#                          solPlot.splot(Gnuplot.GridData(nSol,
#                                                         x[l],
#                                                         y[l],
#                                                         binary=0,
#                                                         inline=0,
#                                                         ))
#                          if analyticalSolution[test] != None:
#                              for i in range(nx[l]):
#                                  for j in range(ny[l]):
#                                      aSol[l][i,j] = analyticalSolution[test].uOfXT(Numeric.array([x[l][i],y[l][j],0.0]),tn)
#                              aSolPlot('set parametric')
#                              aSolPlot('set data style lines')
#                              aSolPlot('set hidden')
#                              aSolPlot('set contour base')
#                              aSolPlot('set cntrparam levels incremental 0.1,0.1,1.0')
#                              aSolPlot.xlabel('x')
#                              aSolPlot.ylabel('y')
#                              aSolPlot.splot(Gnuplot.GridData(aSol[l],
#                                                             x[l],
#                                                             y[l],
#                                                             binary=0,
#                                                             inline=0
#                                                              ))
            if computeEigenvalues:
                gevals("set terminal x11")
                gevals.title(testOut+' eigenvalues');
                gevals.xlabel(r'real(\lambda)')
                gevals.ylabel(r'imag(\lambda)')
                gevals.plot()
                kplot("set terminal x11")            
                kplot.plot()
                kplot.title(testOut+r' ln(\kappa) vs. ln(1/h)')
                kplot.xlabel(r'ln(\kappa)')
                kplot.xlabel('ln(1/h)')
                kList=[]
            if failedFlag == True:
                tn = T[test]
            if computeSpaceTimeError or tn >= T[test]:
                eCoarse=1.0
                eFine=1.0
                hCoarse=1.0
                hFine=1.0
                mFine = mlScalarTransport.modelList[-1]
                for l in range(nLevels):
                    utmp[l].dof[:]=0.0
                if analyticalSolution[test] != None:
                    for eN in range(mFine.q['x'].shape[0]):
                        for k in range(mFine.q['x'].shape[1]):
                            asolq[eN,k] = analyticalSolution[test].uOfXT(mFine.q['x'][eN,k],tn)
                for m,jac,mesh,l in zip(mlScalarTransport.modelList,
                                        mlScalarTransport.jacobianList,
                                        mlMesh.meshList,
                                        range(nLevels)):
                    utmp[l].dof[:] = m.u.dof
                    if l < nLevels-1:
                        for lf in range(l,nLevels-1):
                            mlScalarTransport.meshTransfers.prolong_bcList[lf+1].matvec(utmp[lf].dof,utmp[lf+1].dof)
                            #Load the Dirichlet conditions into the projection 
                            for dofN,g in mlScalarTransport.modelList[lf+1].dirichletConditions.DOFBoundaryConditionDict.iteritems():
                                utmp[lf+1].dof[dofN] = g(mlScalarTransport.modelList[lf+1].dirichletConditions.DOFBoundaryPointDict[dofN],tn)
                        #solPlot.replot(Gnuplot.Data(mlMesh.meshList[-1].nodeArray[:,0],
                        #                          utmp[-1].dof,
                        #                          with='linespoints',
                        #                          title='numerical solution'))
                        utmp[-1].getValues(mFine.q['v'],uqtmp)
                    else:
                        uqtmp[:]=mFine.q['u']
                    if analyticalSolution[test] == None:
                        asolq[:]=mFine.q['u']
                    eCoarse=eFine
                    hCoarse=hFine
                    hFine = mesh.h
                    eFine = L2errorSFEM(mFine.q['dx_a'],asolq,uqtmp)
                    if eSpaceTime.has_key(hFine):
                        eSpaceTime[hFine] += mlScalarTransport.DT*eFine**2
                    else:
                        eSpaceTime[hFine] = mlScalarTransport.DT*eFine**2
                    eSpace[hFine] = eFine
                    if computeEigenvalues:
##                     dev = DenseEigenvalues(jac)
##                     dev.computeEigenvalues()
##                     try:
##                         gevals.replot(
##                             Gnuplot.Data([l for l in dev.eigenvalues],
##                                          [0.0 for l in dev.eigenvalues],
##                                          title='h = %12.5e' % hFine,with='points'))
##                     except TypeError:
##                         gevals.replot(
##                             Gnuplot.Data([l.real for l in dev.eigenvalues],
##                                          [l.imag for l in dev.eigenvalues],
##                                          title='h = %12.5e' % hFine,
##                                          with='points'))
##                     k = max(abs(dev.eigenvalues))/min(abs(dev.eigenvalues))
##                     kList.append(k)
##                     ratio = k*(hFine**2)
##                     print "k*h**2 %12.5E" % ratio
                        multilevelLinearSolver.solverList[l].calculateEigenvalues()
                        gevals.replot(
                            Gnuplot.Data(multilevelLinearSolver.solverList[l].eigenvalues_r,
                                         multilevelLinearSolver.solverList[l].eigenvalues_i,
                                         title='h = %12.5e' % hFine,
                                         with='points'))
                        k = max(abs(multilevelLinearSolver.solverList[l].eigenvalues_r))/min(abs(multilevelLinearSolver.solverList[l].eigenvalues_r))
                        kList.append(k)
                        ratio = k*(hFine**2)
                        print "k*h**2 %12.5E" % ratio
        hFine = 0
        errors='$\|error \|_{L_2(\Omega)}$'
        errorsSpaceTime='$\| \|error \|_{L_2 (\Omega) } \|_{L_2(T)}$ '
        orders='spatial order'
        columns='c|'
        hs="$\\Delta t = %4.2e$, $h=$" % mlScalarTransport.DT
        for mesh in mlMesh.meshList:
            columns+='c'
            hCoarse=hFine
            hFine = mesh.h
            hs += "& %4.2e" % hFine
            if hCoarse != 0:
                if eSpace[hFine] != 0.0 and eSpace[hCoarse] != 0.0:
                    p = (log(eSpace[hFine]) - log(eSpace[hCoarse]))/(log(hFine) - log(hCoarse))
                else:
                    p=0
            else:
                p = 0
            errors+="& %4.2e" % eSpace[hFine]
            orders+="& %4.2e" % p
            errorsSpaceTime+="& %4.2e" % sqrt(eSpaceTime[hFine])
        hs += '\\\ \n \hline'
        errors += '\\\ \n'
        orders += '\\\ \n'
        errorsSpaceTime += '\\\ \n'
        report.write('\\begin{center} \n \\begin{tabular}{'+columns+'}\n')
        report.write(hs)
        report.write(errors)
        report.write(orders)
        report.write(errorsSpaceTime)
        report.write('\\end{tabular}\n \\end{center}\n')
        report.flush()
        print errors
        print orders
        if computeEigenvalues:
            kplot.replot(
                Gnuplot.Data([log(1.0/mesh.h) for mesh in mlMesh.meshList],
                             [log(k) for k in kList],
                             with='linespoints'))
            gevals.hardcopy(testOut+'_eig.eps',eps=1,enhanced=1,color=1)
            kplot.hardcopy(testOut+'_cond.eps',eps=1,enhanced=1,color=1)
	print multilevelNonlinearSolver.info()
        #print multilevelLinearSolver.info()
        print "nsteps %i" % nSteps
        if DG != True:
            nsolFile = 'tmp'+testOut+'_sol.eps'
            asolFile = 'tmp'+testOut+'_asol.eps'
            #solPlot.hardcopy(nsolFile, eps=1,enhanced=1,color=1)
            #aSolPlot.hardcopy(asolFile, eps=1,enhanced=1,color=1)
            #report.write("\\begin{figure}{%s}\n\epsfig{file=%s,scale=0.65}\epsfig{file=%s,scale=0.65}\n \\end{figure}\n" % (test,nsolFile,asolFile))
            report.flush()
        #mlScalarTransport.modelList[-1].writeBoundaryTermsEnsight(test)
        raw_input('Please press return to continue... \n')
	mlScalarTransport.modelList[-1].trialSpace.endTimeSeriesEnsight(timeValues,test,test)
	os.chdir('../../')
    closeLatexReport(report)
    #os.spawnlp(os.P_WAIT,'latex','latex','scalarTransportReport.tex')
    #os.spawnlp(os.P_WAIT,'latex','latex','scalarTransportReport.tex')
    #os.spawnlp(os.P_WAIT,'xdvi','xdvi','scalarTransportReport.dvi')

## @}

if __name__ == '__main__':
    runTests()
