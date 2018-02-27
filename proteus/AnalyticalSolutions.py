"""
Classes representing analytical solutions of differential and partial differential equations.

.. inheritance-diagram:: proteus.AnalyticalSolutions
   :parts: 1
"""
from math import *
from EGeometry import *
from LinearAlgebraTools import *
from Profiling import logEvent

class AS_base:
    """
    The base class for general analytical solutions, u(x,y,z,t)

    For basic error calculations only the :func:`uOfXT` member need be
    overridden. The vectorized member functions such as
    :func:`getUOfXT` are provided for future optimizations but are
    implemented in :class:`AS_base` in a simple, inefficient way.
    """
    def uOfXT(self,X,T):
        """
        Return the solution at a 3D point, X, and a time, T
        """
        pass
    def duOfXT(self,X,T):
        """
        Return the gradient of the solution at a 3D point, X, and a time, T

        .. todo:: Change du to Du
        """
        pass
    def rOfUX(self,u,x):
        """
        Return the reaction term for the (manufactured) solution at a 3D point, X, and a time, T
        """
        pass
    def drOfUX(self,u,x):
        """
        Return the derivative of the reaction term for the (manufactured) solution at a 3D point, X, and a time, T
        """
        pass
    def advectiveFluxOfX(self,x):
        """
        Return the advective flux at a 3D point, X

        .. todo:: Change advective flux to depend on T
        """
        pass
    def diffusiveFluxOfX(self,x):
        """
        Return the diffusive flux at a 3D point, X

        .. todo:: Change advective flux to depend on T
        """
        pass
    def totalFluxOfX(self,x):
        """
        Return the total flux at a 3D point, X

        .. todo:: Decide if we need total flux
        """
        pass
    def getUOfXT(self,X,T,U):
        """
        Set the solution, U, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.uOfXT,X,T,U)
    def getDUOfXT(self,X,T,DU):
        """
        Set the gradient of the solution, DU, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.duOfXT,X,T,DU)
    def getROfUX(self,X,T,R):
        """
        Set the reaction term, R, at an array of 3D points, X, and a time, T
        """
        self.coeff_p_to_v(self.rOfXT,U,X,T,R)
    def getDROfUX(self,X,T,DR):
        """
        Set the derivative of the reaction term, DR, at an array of 3D points, X, and a time, T
        """
        self.coeff_p_to_v(self.drOfXT,U,X,T,DR)
    def getAdvectiveFluxOfXT(self,X,T,F):
        """
        Set the advective flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.advectiveFluxOfXT,X,T,F)
    def getDiffusiveFluxOfX(self,X,T,F):
        """
        Set the diffusive flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.diffusiveFluxOfXT,X,T,F)
    def getTotalFluxOfX(self,X,T,F):
        """
        Set the total flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.totalFluxOfXT,X,T,F)
    def p_to_v(self,func,X,T,funcValVec):
        """
        Calculate a point-wise function over and array of points.
        """
        nPoints=1
        for d in X.shape[:-1]:
            nPoints*=d
        for i in range(nPoints):
            funcValVec[i]=funcValVec(X.flat[i*3:i*3+3],T)

class DAE_base(AS_base):
    """
    The base class for differential-algebraic equation solutions

    To use this class override :func:`uOfT`
    """
    def uOfXT(self,X,T):
        return self.uOfT(self,T)
    def uOfT(self,T):
        pass

class NonlinearDAE(DAE_base):
    r"""
    The exact solution of the simple nonlinear DAE

    .. math: u_t  = - a \max(u,0)^p \quad u(0) = 1
    """
    def __init__(self,a,p):
        """
        Set the coefficients a and p
        """
        self.a_ = a
        self.p_ = p
        if self.p_ == 1:
            self.func = lambda t : exp(-self.a_*t)
        else:
            q = 1.0/(1.0 - self.p_)
            self.func = lambda t : max(1.0 - (1.0 - self.p_)*self.a_*t,0.0)**q
    def uOfT(self,t):
        return self.func(t)
    def fOfUT(self,u,t):
        return -self.a_*max(u,0.0)**self.p_

class SteadyState(AS_base):
    """
    Based class for steady-state solutions u(x)

    Override :func:`uOfX` to define steady state solutions.
    """
    def uOfX(self,X):
        pass
    def uOfXT(self,X,T=None):
        return self.uOfX(X)

class LinearAD_SteadyState(SteadyState):
    r"""
    The solution of the steady linear advection-diffusion equation

    The boundary value problem is

    .. math: (bu - au_x)_x = 0 \quad u(0) = 1 \quad u(1) = 0
    """
    def __init__(self,b=1.0,a=5.0e-1):
        self.b_ = b
        self.a_ = a
        if b!=0.0:
            self.D_ = (1.0/(exp(b/a)-1.0))
        else:
            self.D_ = 0.0
        self.C_ = -self.D_*exp(b/a)
    def uOfX(self,X):
        x=X[0]
        if self.D_ !=0.0:
            return -self.D_*exp(self.b_*x/self.a_) - self.C_
        else:
            return 1.0-x

class NonlinearAD_SteadyState(LinearAD_SteadyState):
    r"""
    The solution of a steady nonlinear advection-diffusion equation

    The boundary value problem is

    .. math: (bu^q - a(u^r)_x)_x = 0 \quad u(0) = 1 \quad u(1) = 0

    Currently :math:`q=1,r=2`, and :math:`r=2,q=1` are implemented
    """
    def __init__(self,b=1.0,q=2,a=5.0e-1,r=1):
        LinearAD_SteadyState.__init__(self,b,a)
        self.r_ = r
        self.q_ = q
        if (q==2 and r==1):
            if b!= 0.0:
                def f(rtmC):
                    return rtmC*tanh(b*rtmC/a) - 1.0
                def df(rtmC):
                    return rtmC*(1.0/cosh(b*rtmC/a)**2)*(b/a) + tanh(b*rtmC/a)
                logEvent("Solving for sqrt(-C) for q=2,r=1")
                rtmC = sqrt(1.5)
                while abs(f(rtmC)) > 1.0e-8:
                    rtmC -= f(rtmC)/df(rtmC)
                logEvent("sqrt(-C)="+`rtmC`)
                self.rtmC_ = rtmC
                self.sol_ = lambda x : self.rtmC_* \
                            tanh((-self.b_*self.rtmC_/self.a_)*(x - 1.0))
            else:
                self.sol_ = lambda x: 1.0 - x
        elif (q==1 and r==2):
            logEvent("Solving for C in q=1,r=2")
            def f(C):
                return 2.0*C*(log(C-1.0) - log(C)) + 2.0 + self.b_/self.a_
            def df(C):
                return 2.0*(log(C-1.0) - log(C)) + 2.0*C*(1.0/(C-1.0) - 1.0/C)
            C = 1.0 + 1.0e-10
            f0 = f(C)
            print f0
            while abs(f(C)) > (1.0e-7*abs(f0) + 1.0e-7):
                dC = -f(C)/df(C)
                logEvent("dc")
                print dC
                Ctmp = C + dC
                while (abs(f(Ctmp)) > 0.99*abs(f(C))
                       or Ctmp <= 1.0):
                    print f(Ctmp)
                    print f(C)
                    logEvent("ls")
                    dC*=0.9
                    Ctmp = C + dC
                logEvent("out")
                print Ctmp
                print f(Ctmp)
                print df(Ctmp)
                C=Ctmp
            logEvent("C="+`C`)
            self.nlC_ = C
            self.nlD_ = 0.5*(2.0*C*log(C*(C-1)) - \
                             4.0*C + 2.0 - self.b_/self.a_)
            logEvent("D="+`self.nlD_`)
#              def f(C):
#                  return (2.0*self.a_/self.b_)*(1.0 +
#                                                C*log((C-1.0)/C))
#              def df(C):
#                  return (2.0*self.a_/self.b_)*(log((C-1.0)/C) +
#                                                1.0/(C-1.0))
#              C = 1.8*self.b_
#              f0 = f(C)
#              print f0
#              while abs(f(C)) > 1.0e-8*f0 + 1.0e-8:
#                  C -= f(C)/df(C)
#                  print C
#                  print f(C)
#              self.nlC_ = C
        else:
            logEvent("q,r not implemented")
    def uOfX(self,X):
        x=X[0]
        if (self.q_==2 and self.r_==1):
            return self.sol_(x)
        elif (self.q_==1 and self.r_==2):
            def f(u):
                return (2.0*(self.nlC_*log(self.nlC_-u) - \
                        (self.nlC_-u)) -
                        self.nlD_) - self.b_*x/self.a_
            def df(u):
                return (2.0*self.nlC_/(u-self.nlC_)+2.0)
            u=LinearAD_SteadyState.uOfX(self,X)
            f0 = f(u)
            while abs(f(u)) > 1.0e-6*abs(f0) + 1.0e-6:
                u-=f(u)/df(u)
            return u
#              def f(u):
#                  return ((2.0*self.a_/self.b_**2)*
#                          (self.nlC*log(self.nlC-self.b_*u) -
#                           (self.nlC - self.b_*u) -
#                           self.nlD) - x)
#              def df(u):
#                  return ((2.0*self.a_/self.b_**2)*
#                          (-self.nlC*self.b_/(self.nlC-self.b_*u) +
#                           self.b_))
#              u=LinearAD_SteadyState.uOfX(self,X)
#              f0 = abs(f(u))
#              while abs(f(u)) > 1.0e-10*f0 +1.0e-10:
#                  u -= f(u)/df(u)
#              return u
        else:
            logEvent("q,r not implemented")

class LinearADR_Sine(SteadyState):
    r"""
    An exact solution and source term for

    .. math:: \nabla \cdot (\mathbf{b} u - \mathbf{a} \nabla u) + c u + d = 0

    The solution and source are

    .. math::
       :nowrap:

       \begin{eqnarray}
       u(x)   &=& \sin[\mathbf{\omega} \cdot \mathbf{x} + \omega_0) ] = \sin(Ax - b) \\
       r(u,x) &=& - ((\mathbf{a} \mathbf{\omega}) \cdot \mathbf{\omega}) u \\
              & & - (\mathbf{b} \cdot \omega) \cos(Ax - b) \\
              &=& c u + D \cos(Ax - b) \\
              &=& cu + d
       \end{eqnarray}

    also returns the advective, diffusive, and total flux at a point.
    """

    def __init__(self,
                 b=numpy.array((1.0,0.0,0.0)),
                 a=numpy.array(((1.0e-2,0.0,0.0),
                                  (0.0,1.0e-2,0.0),
                                  (0.0,0.0,1.0e-2))),
                 c=1.0,
                 omega=numpy.array((2*pi,0.0,0.0)),
                 omega0=0.0):
        """
        Set the coefficients a,b,c,omega, and omega0
        """
        self.b_ = b
        self.a_ = a
        self.c_ = c
        self.omega_=omega
        self.omega0_=omega0
        self.D_ = - numpy.dot(b,omega)
        self.E_ = - numpy.dot(numpy.dot(a,omega),omega) - c

    def uOfX(self,x):
        return sin(numpy.dot(self.omega_,x[:self.omega_.shape[0]]) + self.omega0_)

    def duOfX(self,x):
        return self.omega_ * cos(numpy.dot(self.omega_,x[:self.omega_.shape[0]])+self.omega0_)

    def rOfUX(self,u,x):
        return self.c_*u + self.D_*cos(numpy.dot(self.omega_,x[:self.omega_.shape[0]])+self.omega0_) + self.E_*sin(numpy.dot(self.omega_,x[:self.omega_.shape[0]]) + self.omega0_)

    def drOfUX(self,u,x):
        return self.c_

    def advectiveFluxOfX(self,x):
        return self.b_*self.uOfX(x)

    def diffusiveFluxOfX(self,x):
        return -numpy.dot(self.a_,self.duOfX(x))

    def totalFluxOfX(self,x):
        return self.advectiveFluxOfX(x) + self.diffusiveFluxOfX(x)

class PoissonsEquation(SteadyState):
    r"""
    Manufactured solution of Poisson's equation

    .. math:: - \Delta u - f = 0

    where u and f are given by

    .. math::
       :nowrap:

       \begin{eqnarray}
       u &=& K x(1-x)y(1-y)z(1-z)e^{x^2 + y^2 + z^2} \\
       f &=& -K\{[y(1-y)z(1-z)][4x^3 - 4x^2 + 6x - 2]+ \\
           &&[x(1-x)z(1-z)][4y^3 - 4y^2 + 6y - 2]+ \\
           &&[x(1-x)y(1-y)][4z^3 - 4z^2 + 6z - 2]\}e^{x^2 + y^2 + z^2}
       \end{eqnarray}
    """
    def __init__(self,K=10.0,nd=3):
        self.K_=K
        self.nd_=nd
    def uOfX(self,X):
        if self.nd_==1:
            x = X[0]
            return self.K_*x*(1.0-x)*exp(x**2)
        elif self.nd_==2:
            x = X[0]
            y = X[1]
            return self.K_*x*(1.0-x)*y*(1.0-y)*exp(x**2 + y**2)
        else:
            x = X[0]
            y = X[1]
            z = X[2]
            return self.K_*x*(1.0-x)*y*(1.0-y)*z*(1.0-z)*exp(x**2 + y**2 + z**2)
    def rOfUX(self,u,X):
        if self.nd_==1:
            x = X[0]
            return  self.K_*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*(1.0-x)*x - 2.0)*exp(x**2)
        elif self.nd_==2:
            x = X[0]
            y = X[1]
            return  self.K_*(y*(1.0-y)*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*x*(1.0-x) - 2.0)+
                             x*(1.0-x)*(4.0*(1.0-y)*y**3 - 4.0*y**2 + 6.0*y*(1.0-y) - 2.0))*exp(x**2 + y**2)
        else:
            x = X[0]
            y = X[1]
            z = X[2]
            return  self.K_*(y*(1.0-y)*z*(1.0-z)*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*x*(1.0-x) - 2.0)+
                             x*(1.0-x)*z*(1.0-z)*(4.0*(1.0-y)*y**3 - 4.0*y**2 + 6.0*y*(1.0-y) - 2.0)+
                             x*(1.0-x)*y*(1.0-y)*(4.0*(1.0-z)*z**3 - 4.0*z**2 + 6.0*z*(1.0-z) - 2.0))*exp(x**2 + y**2 + z**2)
    def drOfUX(self,u,x):
        return 0.0

class LinearAD_DiracIC(AS_base):
    r"""
    The exact solution of

    .. math:: u_t + \nabla \cdot (b u - a \nabla u) = 0

    on an infinite domain with dirac initial data

    .. math:: u0 = \int u0 \delta(x - x0)

    also returns advective, diffusive, and total flux
    """

    def __init__(self,
                 n=1.0,
                 b=numpy.array((1.0,0.0,0.0)),
                 a=1.0e-2,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0,0.0,0.0))):
        self.n_  = n
        self.u0_ = u0
        self.x0_ = x0
        self.b_  = b
        self.a_  = a
        self.tStart = tStart
    def uOfXT(self,x,T):
        t = T + self.tStart
        y = (x[:self.b_.shape[0]] - self.x0_[:self.b_.shape[0]] - self.b_*t)
        exp_arg=numpy.dot(y,y) / (4.0 * self.a_ * t)
        if exp_arg > 100:
            return 0.0
        else:
            return self.u0_*exp(-exp_arg) / (4.0*self.a_*pi*t)**(self.n_/2.0)

    def duOfXT(self,x,T):
        t = T + self.tStart
        y = (x[:self.b_.shape[0]] - self.x0_[:self.b_.shape[0]] - self.b_*t)
        return self.uOfXT(x,T)*2.0*y/(4.0*self.a_*t)

    def advectiveFluxOfXT(self,x,T):
        return self.b_*self.uOfXT(x,T)

    def diffusiveFluxOfXT(self,x,T):
        return -self.a_*self.duOfXT(x,T)

    def totalFluxOfXT(self,x,T):
        return self.advectiveFluxOfXT(x,T) + self.diffusiveFluxOfXT(x,T)

class LinearADR_Decay_DiracIC(LinearAD_DiracIC):
    r"""
    The exact solution of

    .. math:: u_t + \nabla \cdot (bu - a \nabla u) + cu= 0

    on an infinite domain with Dirac initial data.  Also
    returns the fluxes (by inheritance).
    """
    def __init__(self,
                 n=1.0,
                 b=numpy.array((1.0,0.0,0.0)),
                 a=1.0e-2,
                 c=1.0,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0,0.0,0.0))):
        LinearAD_DiracIC.__init__(self,n,b,a,tStart,u0,x0)
        self.c_ = c

    def uOfXT(self,x,T):
        t = T + self.tStart
        return LinearAD_DiracIC.uOfXT(self,x,T)*exp(- self.c_*t)

    def rOfUXT(self,u,x,T):
        return self.c_*u

    def drOfUXT(self,u,x,T):
        return self.c_

class NonlinearADR_Decay_DiracIC(LinearADR_Decay_DiracIC):
    r"""
    The approximate analytical solution of

    .. math:: u_t + \nabla \cdot (bu - a \nabla u) + cu^d= 0

    on an infinite domain with Dirac initial data.
    Also returns the fluxes (by inheritance).
    """

    def __init__(self,
                 n=1.0,
                 b=numpy.array((1.0,0.0,0.0)),
                 a=1.0e-2,
                 c=1.0,
                 d=2.0,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0,0.0,0.0))):
        LinearADR_Decay_DiracIC.__init__(self,n,b,a,c,tStart,u0,x0)
        self.d_=d

    def uOfXT(self,x,T):
        t = T + self.tStart
        u1=LinearAD_DiracIC.uOfXT(self,x,T)
        if u1 > 0.0:
            return u1*exp(-(2.0*self.c_*t*pow(u1,self.d_-1.0))/(self.d_+1.0))
        else:
            return u1

    def rOfUXT(self,u,x,T):
        return self.c_*(u**self.d_)

    def drOfUXT(self,u,x,T):
        return self.d_*self.c_*(u**(self.d_-1))

class PlaneCouetteFlow_u(SteadyState):
    """
    The exact solution for the u component of  velocity in plane Couette Flow
    """
    from canalyticalSolutions import PlaneCouetteFlow_u
    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0,0.0,0.0]):
        self.iwork = numpy.zeros((1,),'i')
        self.rwork = numpy.array([upperPlateVelocity,
                                  plateSeperation,
                                  origin[0],
                                  origin[1],
                                  origin[2]])
    def uOfX(self,x):
        uList=numpy.array([[0.0]])
        xList=numpy.array([x])
        t=0.0
        nPoints=1
        self.PlaneCouetteFlow_u(self.iwork,self.rwork,t,xList,uList)
        return uList.flat[0]

class PlaneCouetteFlow_v(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Couette Flow
    """
    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0,0.0,0.0]):
        pass
    def uOfX(self,x):
        return 0.0

class PlaneCouetteFlow_p(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Couette Flow
    """
    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0,0.0,0.0]):
        pass
    def uOfX(self,x):
        return 0.0

class PlanePoiseuilleFlow_u(SteadyState):
    from canalyticalSolutions import PlanePoiseuilleFlow_u
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q = 0.0,
                 origin=[0.0,0.0,0.0]):
        #allocated iwork and rwork arrays and load in values
        self.iwork = numpy.zeros((1,),'i') #( (shape_0,shape_1,...), typecode)
        self.rwork = numpy.array([plateSeperation,
                                    mu,
                                    grad_p,
                                    q,
                                    origin[0],
                                    origin[1],
                                    origin[2]]) # numpy.array([component_0,component_2,...]) builds a double* initialized to values in list [...]
    #wrapped "vectorized" analytical solutions that Moira wrote with a point-wise evaluation
    def uOfX(self,x):
        uList=numpy.array([0.0]) #array of doubles for storing solution (just the solution at a single point)
        xList=numpy.array([x])   #array of double*'s for storing spatial locations (just a single point in 3D in this case)
        t=0.0     #time at which solution is evaluated (this is a steady state problem)
        nPoints=1 #just a single point
        self.PlanePoiseuilleFlow_u(self.iwork,self.rwork,t,xList,uList) #make call to wrapped C function
        return uList[0] #return value  of solution at this point

class PlanePoiseuilleFlow_v(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """
    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q = 0.0,
                 origin=[0.0,0.0,0.0]):
        pass
    def uOfX(self,x):
        return 0.0

class PlanePoiseuilleFlow_p(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """
    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q = 0.0,
                 origin=[0.0,0.0,0.0]):
        self.grad_p = grad_p
    def uOfX(self,x):
        return self.grad_p*x[0]

#encapsulate analytical solver for BuckleyLeverett Riemann problems
from proteus.ObjectiveFunctions import OsherFuncCoef
from proteus.Optimizers import fminbound
class Buckley_Leverett_RiemannSoln(AS_base):
    def __init__(self,coefficients,uLeft=1.0,uRight=0.0,t0=0.0,x0=0.0,T=0.5,ftol=1.0e-8,
                 useShallowCopyCoef=True):
        #ought to start with range and spline for eval later I guess
        self.coefficients = coefficients
        self.uL = uLeft; self.uR = uRight;
        self.t0 = t0; self.x0 = x0;
        self.T  = T;  self.ftol = ftol
        self.riemF = OsherFuncCoef(self.uL,self.uR,self.coefficients,self.T-self.t0,self.x0,
                                   useShallowCopy=useShallowCopyCoef)
        self.solver= fminbound(self.riemF,self.ftol)
        #mwf if want to look at coefficients
        #self.generateAplot(x0,x0+1.0,101,T)
    def uOfXT(self,x,t):
        if abs(t-self.t0) < 1.0e-7:
            if x[0]-self.x0 <= 0.0:
                return self.uL
            return self.uR
        self.riemF.xi = (x[0]-self.x0)/(t-self.t0)
        u,f = self.solver.solve(Guess_x = 0.5*(self.uL+self.uR))
        return u
    def uOfX(self,x):
        return self.uOfXT(x,self.T)
    def generateAplot(self,xLeft,xRight,nnx,t,filename='Buckley_Leverett_ex.dat'):
        """
        save exact solution to a file
        """
        dx = (xRight-xLeft)/(nnx-1.)
        fout = open(filename,'w')
        for i in range(nnx):
            x = (xLeft+dx*i,0.0,0.0)
            u = self.uOfXT(x,t)
            fout.write('%12.5e  %12.5e \n' % (x[0],u))
        #
        fout.close()
class PlaneBase(SteadyState):
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0,1.0,1.0]):
        self.plane_n = numpy.array([cos(plane_theta)*sin(plane_phi),
                                    sin(plane_theta)*sin(plane_phi),
                                    cos(plane_phi)])
        if (plane_phi > -pi/2.0 and plane_phi < pi/2.0):
            if plane_phi == 0.0:
                logEvent("plate is in x-y plane")
                v_phi = pi/2.0
            else:
                v_phi = atan(-1.0/(tan(plane_phi)*cos(v_theta-plane_theta)))
        else:
            logEvent("plate is is parallel to z-axis, using angle of velocity with z-axis instead of angle with x-axis")
            v_theta = plane_theta - pi/2.0
        self.v_n = numpy.array([cos(v_theta)*sin(v_phi),
                                sin(v_theta)*sin(v_phi),
                                cos(v_phi)])
        minZ = (0.0,numpy.array([0.0,0.0,0.0]))
        maxZ = (0.0,numpy.array([0.0,0.0,0.0]))
        for i in [0,1]:
            for j in [0,1]:
                for k in [0,1]:
                    x = numpy.array([i*L[0],j*L[1],k*L[2]])
                    Z = numpy.dot(x,self.plane_n)
                    if Z < minZ[0]:
                        minZ = (Z,x)
                    if Z > maxZ[0]:
                        maxZ = (Z,x)
        self.Zshift = minZ[0]
        self.Zsep = maxZ[0] - minZ[0]
        print self.Zshift
        #allocated iwork and rwork arrays and load in values
#         self.iwork = numpy.zeros((1,),'i') #( (shape_0,shape_1,...), typecode)
#         self.rwork = numpy.array([mu,
#                                   grad_p,
#                                   v_norm,
#                                   self.plane_n[0],
#                                   self.plane_n[1],
#                                   self.plane_n[2],
#                                   self.v_n[0],
#                                   self.v_n[1],
#                                   self.v_n[2],
#                                   minZ[1][0],
#                                   minZ[1][1],
#                                   minZ[1][2]])
        self.mu = mu
        self.grad_p = grad_p
        self.v_norm = v_norm
    def Z(self,x):
        return numpy.dot(x,self.plane_n) - self.Zshift
    def X(self,x):
        return numpy.dot(x,self.v_n)
    def U(self,x):
        Z = self.Z(x)
        return (Z*self.v_norm/self.Zsep) - self.grad_p*Z*(self.Zsep-Z)/(2.0*self.mu)

class PlanePoiseuilleFlow_u2(PlaneBase):
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0,1.0,1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)
    def uOfX(self,x):
        return self.U(x)*self.v_n[0]

class PlanePoiseuilleFlow_v2(PlaneBase):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """
    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0,1.0,1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)
    def uOfX(self,x):
        return self.U(x)*self.v_n[1]

class PlanePoiseuilleFlow_w2(PlaneBase):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """
    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0,1.0,1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)
    def uOfX(self,x):
        return self.U(x)*self.v_n[2]

class PlanePoiseuilleFlow_p2(PlaneBase):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """
    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0,1.0,1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)
    def uOfX(self,x):
        return self.grad_p*self.X(x)

class VortexDecay_u(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n=n
        self.Re=Re
    def uOfXT(self,x,t):
        #return -cos(x[0])*sin(x[1])*exp(-2.0*t)
        return -cos(self.n*pi*x[0])*sin(self.n*pi*x[1])*exp(-2.0 * (self.n*pi)**2 * t / self.Re)

class VortexDecay_v(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n=n
        self.Re=Re
    def uOfXT(self,x,t):
        #return sin(x[0])*cos(x[1])*exp(-2.0*t)
        return sin(self.n*pi*x[0])*cos(self.n*pi*x[1])*exp(-2.0 * (self.n*pi)**2 * t / self.Re)

class VortexDecay_p(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    #constructor or initializer function, reads in parameters needed in the iwork and rwork arrays
    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n=n
        self.Re=Re
    def uOfXT(self,x,t):
        #return -self.Re*0.25*(cos(2.0*x[0])+cos(2.0*x[1]))*exp(-4.0*t)
        return -0.25*(cos(2.0*self.n*pi*x[0]) + cos(2.0*self.n*pi*x[1]))*exp(-4.0 * (self.n*pi)**2 * t / self.Re)

## @}

if __name__ == '__main__':
    for i in range(5):
        for j in range(5):
            for k in range(5):
                for l in range(5):
                    p = PlaneBase(plane_theta=i*math.pi/4.0,
                                  plane_phi=j*math.pi/4.0-pi/2.0,
                                  v_theta=k*math.pi/4.0,
                                  v_phi=l*math.pi/4.0-pi/2.0,
                                  v_norm=1.0,
                                  mu=1.0,
                                  grad_p=1.0)
#     import Gnuplot
#     import MeshTools
#     import AnalyticalSolutions
#     import canalyticalSolutions
#     g = Gnuplot.Gnuplot(debug=1)
#     grid = MeshTools.RectangularGrid(500,1,1,Lx=1.0)
#     X=numpy.array([n.p for n in grid.nodeList])
#     Xx=numpy.array([n.p[0] for n in grid.nodeList])
#     logEvent("Testing Solutions in 1D")
#     logEvent("NonlinearDAE")
#     sol = NonlinearDAE(5.0,1.0)
#     y1=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,2.0)
#     y2=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,3.0)
#     y3=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,0.75)
#     y075=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,0.25)
#     y025=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,-0.25)
#     yn025=[sol.uOfT(x) for x in Xx]
#     sol = NonlinearDAE(5.0,-1.25)
#     yn125=[sol.uOfT(x) for x in Xx]
#     g.plot(Gnuplot.Data(Xx,y1,title='Solution,q=1'),
#            Gnuplot.Data(Xx,y2,title='Solution,q=2'),
#            Gnuplot.Data(Xx,y3,title='Solution,q=3'),
#            Gnuplot.Data(Xx,y075,title='Solution,q=0.75'),
#            Gnuplot.Data(Xx,y025,title='Solution,q=0.25'),
#            Gnuplot.Data(Xx,yn025,title='Solution,q=-0.25'),
#            Gnuplot.Data(Xx,yn125,title='Solution,q=-1.25'))
#     #raw_input('Please press return to continue... \n')
#     logEvent("AD_SteadyState")
#     sol=LinearAD_SteadyState()
#     linearSS_data = Gnuplot.Data(Xx,[sol.uOfX(x) for x in X],
#                         title='Solution,q=1,r=1')
#     logEvent("NonlinearAD_SteadyState")
#     sol=NonlinearAD_SteadyState(q=2,r=1)
#     nonlinearSS_datap1q2 = Gnuplot.Data(Xx,
#                                         [sol.uOfX(x) for x in X],
#                                         title='Solution,q=2,r=1')
#     sol=NonlinearAD_SteadyState(q=1,r=2)
#     nonlinearSS_datap2q1 = Gnuplot.Data(Xx,
#                                         [sol.uOfX(x) for x in X],
#                                         title='Solution,q=1,r=2')
#     g.plot(linearSS_data,nonlinearSS_datap1q2,nonlinearSS_datap2q1)
#     #raw_input('Please press return to continue... \n')
#     logEvent("LinearADR_Sine")
#     sol=LinearADR_Sine()
#     g.plot(Gnuplot.Data(Xx,[sol.uOfX(x) for x in X],
#                         title='Solution'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,[sol.duOfX(x)[0] for x in X],
#                         title='Gradient'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfX(x)[0] for x in X],
#                         title='Advective Flux'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfX(x)[0] for x in X],
#                         title='Diffusive Flux'))
#     #raw_input('Please press return to continue... \n')

#     g.plot(Gnuplot.Data(Xx,
#                         [sol.totalFluxOfX(x)[0] for x in X],
#                         title='Total Flux'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.rOfUX(sol.uOfX(x),x) for x in X],
#                         title='reaction'))
#     #raw_input('Please press return to continue... \n')
#     logEvent("LinearAD_DiracIC")
#     sol=LinearAD_DiracIC()
#     g.plot(Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.25) for x in X],
#                         title='Solution,t=0.25'),
#            Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.75) for x in X],
#                         title='Solution,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.25)[0] for x in X],
#                         title='Gradient,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.75)[0] for x in X],
#                         title='Gradient,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Advective Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Advective Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Diffusive Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Diffusive Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Total Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Total Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     logEvent("LinearADR_Decay_DiracIC")
#     sol=LinearADR_Decay_DiracIC()
#     g.plot(Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.25) for x in X],
#                         title='Solution,t=0.25'),
#            Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.75) for x in X],
#                         title='Solution,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.25)[0] for x in X],
#                         title='Gradient,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.75)[0] for x in X],
#                         title='Gradient,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Advective Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Advective Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Diffusive Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Diffusive Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Total Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Total Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.rOfUXT(sol.uOfXT(x,T=0.25),x,T=0.25) for x in X],
#                         title='Reaction,T=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.rOfUXT(sol.uOfXT(x,T=0.75),x,T=0.75) for x in X],
#                         title='Reaction,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     logEvent("NonlinearADR_Decay_DiracIC")
#     sol=NonlinearADR_Decay_DiracIC()
#     g.plot(Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.25) for x in X],
#                         title='Solution,t=0.25'),
#            Gnuplot.Data(Xx,[sol.uOfXT(x,T=0.75) for x in X],
#                         title='Solution,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.25)[0] for x in X],
#                         title='Gradient,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.duOfXT(x,T=0.75)[0] for x in X],
#                         title='Gradient,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Advective Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.advectiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Advective Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Diffusive Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.diffusiveFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Diffusive Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.25)[0] for x in X],
#                         title='Total Flux,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.totalFluxOfXT(x,T=0.75)[0] for x in X],
#                         title='Total Flux,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     g.plot(Gnuplot.Data(Xx,
#                         [sol.rOfUXT(sol.uOfXT(x,T=0.25),x,T=0.25) for x in X],
#                         title='Reaction,t=0.25'),
#            Gnuplot.Data(Xx,
#                         [sol.rOfUXT(sol.uOfXT(x,T=0.75),x,T=0.75) for x in X],
#                         title='Reaction,t=0.75'))
#     #raw_input('Please press return to continue... \n')
#     logEvent("-------------------------------------------")
#     logEvent("-------------------------------------------")
#     logEvent("Starting test of canalyticalSolutionsModule")
#     logEvent("-------------------------------------------")
#     logEvent("-------------------------------------------")
#     iwork = numpy.zeros((2,),'i')
#     iwork[0]=3
#     iwork[1]=4
#     rwork = numpy.zeros((1,),'d')
#     nPoints_x = 33
#     nPoints_y = 33
#     nPoints_z = 34
#     t=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#                     for k in range(nPoints_z):
#                             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#                             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#                             x[i,j,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros(x.shape[:-2],'d')
#     #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((2,),'d')
#     nPoints = 101
#     t=1.0
#     x = numpy.zeros((nPoints,3),'d')
#     for i in range(nPoints):
#         x[i,0] = i*(1.0/(nPoints-1.0))
#     gnuplot = Gnuplot.Gnuplot()
#     gnuplot('set multiplot')
#     gnuplot('set yrange [0:1]')
#     gnuplot('set key on graph .93,.97')
#     rwork[0]=5.0
#     rwork[1]=1.0
#     u = numpy.zeros(x.shape[0],'d')
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=1'))
#     rwork[1]=2.0
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.93')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=2'))
#     rwork[1]=3.0
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.89')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=3'))
#     rwork[1]=0.75
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.85')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=0.75'))
#     rwork[1]=0.25
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.81')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=0.25'))
#     rwork[1]=-0.25
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.77')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=-0.25'))
#     rwork[1]=-1.25
#     canalyticalSolutions.NonlinearDAE(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.73')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearDAE, q=-1.25'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------python code----------------------------------------
#     # upy = numpy.zeros(x.shape[0],'d')
#     # pySol = AnalyticalSolutions.LinearADR_Sine(b=numpy.array((1.0,0.0,0.0)),
#     #                                            a=numpy.array(((1.0e-2,0.0,0.0),
#     #                                                             (0.0,1.0e-2,0.0),
#     #                                                             (0.0,0.0,1.0e-2))),
#     #                                            c=1.0,
#     #                                            omega=numpy.array((iwork[0]*2*math.pi,0.0,0.0)),
#     #                                            omega0=0.0)
#     # for i in range(nPoints):
#     #     upy[i]=pySol.uOfX(x[i,:])
#     # gnuplot.plot(Gnuplot.Data(x[:,0],
#     #                           upy,
#     #
#     #                           title='NonlinearDAE'))
#     # #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((2,),'d')
#     rwork[0]=5.0
#     rwork[1]=2.0
#     t=1.0
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                 x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                 x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                 x[i,j,k,2] = k*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     # x is the initial u.
#     u = x.copy()
#     canalyticalSolutions.NonlinearDAE_f(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     slice = nPoints_z/2
#     gnuplot.splot(Gnuplot.GridData(u[:,:,slice,0],
#                                    x[:,0,slice,0],
#                                    x[0,:,slice,1],

#                                    title='NonlinearDAE_f'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((2,),'d')
#     rwork[0]=1.0
#     rwork[1]=0.5
#     nPoints = 101
#     t=0.0
#     x = numpy.zeros((nPoints,3),'d')
#     for i in range(nPoints):
#         x[i,0] = i*(1.0/(nPoints-1.0))
#     u = numpy.zeros(x.shape[0],'d')
#     canalyticalSolutions.LinearAD_SteadyState(iwork,rwork,t,x,u)
#     # gnuplot.clear()
#     gnuplot('set multiplot')
#     gnuplot('set yrange [0:1]')
#     gnuplot('set key on graph .93,.97')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='LinearAD_SteadyState: q=1, r=1'))
#     # #raw_input('Please press return to continue... \n')
#     #---------------------------------------------------------------
#     # rwork = numpy.zeros((2,),'d')
#     # rwork[0]=1.0
#     # rwork[1]=0.5
#     # nPoints = 101
#     # t=0.0
#     # x = numpy.zeros((nPoints,3),'d')
#     # for i in range(nPoints):
#     #     x[i,0] = i*(1.0/(nPoints-1.0))
#     # u = numpy.zeros(x.shape[0],'d')
#     # gnuplot.clear()
#     # gnuplot('set multiplot')
#     # gnuplot('set yrange [0:1]')
#     iwork = numpy.zeros((2,),'i')
#     iwork[0]=2
#     iwork[1]=1
#     canalyticalSolutions.NonlinearAD_SteadyState(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.93')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearAD_SteadyState: q=2, r=1'))
#     iwork[0]=1
#     iwork[1]=2
#     canalyticalSolutions.NonlinearAD_SteadyState(iwork,rwork,t,x,u)
#     gnuplot('set key on graph .93,.89')
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,
#                               title='NonlinearAD_SteadyState: q=1, r=2'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((4,),'d')
#     rwork[0]=2.0*3.14159
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=0.0
#     t=0.0
#     nPoints_x = 51
#     nPoints_y = 51
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.LinearADR_Sine(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u,
#                                    x[:,0,0],
#                                    x[0,:,1],

#                                    title='LinearADR_Sine'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((4,),'d')
#     rwork[0]=2.0*3.14159
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=0.0
#     t=0.0
#     nPoints_x=34
#     nPoints_y=34
#     nPoints_z=34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#              for k in range(nPoints_z):
#                   x[i,j,k,0] = i/(nPoints_x-1.0)
#                   x[i,j,k,1] = j/(nPoints_y-1.0)
#                   x[i,j,k,2] = k/(nPoints_z-1.0)
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     canalyticalSolutions.LinearADR_Sine_du(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                    x[:,0,0,0],
#                                    x[0,:,0,1],

#                                    title='LinearADR_Sine_du: Gradient'))
#     #raw_input('Please press return to continue... \n')
#     # #-----------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((17,),'d')
#     rwork[0]=2.0*3.14159
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=0.0
#     rwork[4]=1.0
#     rwork[5]=0.0
#     rwork[6]=0.0
#     rwork[7]=0.01
#     rwork[8]=0.0
#     rwork[9]=0.0
#     rwork[10]=0.0
#     rwork[11]=0.01
#     rwork[12]=0.0
#     rwork[13]=0.0
#     rwork[14]=0.0
#     rwork[15]=0.01
#     rwork[16]=1.0
#     t=0.0
#     nPoints_x=34
#     nPoints_y=34
#     nPoints_z=34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#              for k in range(nPoints_z):
#                   x[i,j,k,0] = i/(nPoints_x-1.0)
#                   x[i,j,k,1] = j/(nPoints_y-1.0)
#                   x[i,j,k,2] = k/(nPoints_z-1.0)
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     canalyticalSolutions.LinearADR_Sine_advectiveVelocity(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                    x[:,0,0,0],
#                                    x[0,:,0,1],

#                                    title='LinearADR_Sine_AvectiveVelocity'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((17,),'d')
#     rwork[0]=2.0*3.14159
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=0.0
#     rwork[4]=1.0
#     rwork[5]=0.0
#     rwork[6]=0.0
#     rwork[7]=0.01
#     rwork[8]=0.0
#     rwork[9]=0.0
#     rwork[10]=0.0
#     rwork[11]=0.01
#     rwork[12]=0.0
#     rwork[13]=0.0
#     rwork[14]=0.0
#     rwork[15]=0.01
#     rwork[16]=1.0
#     t=0.0
#     nPoints_x=34
#     nPoints_y=34
#     nPoints_z=34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#              for k in range(nPoints_z):
#                   x[i,j,k,0] = i/(nPoints_x-1.0)
#                   x[i,j,k,1] = j/(nPoints_y-1.0)
#                   x[i,j,k,2] = k/(nPoints_z-1.0)
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     canalyticalSolutions.LinearADR_Sine_diffusiveVelocity(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                    x[:,0,0,0],
#                                    x[0,:,0,1],

#                                    title='LinearADR_Sine_DiffusiveVelocity'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((17,),'d')
#     rwork[0]=2.0*3.14159
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=0.0
#     rwork[4]=1.0
#     rwork[5]=0.0
#     rwork[6]=0.0
#     rwork[7]=0.01
#     rwork[8]=0.0
#     rwork[9]=0.0
#     rwork[10]=0.0
#     rwork[11]=0.01
#     rwork[12]=0.0
#     rwork[13]=0.0
#     rwork[14]=0.0
#     rwork[15]=0.01
#     rwork[16]=1.0
#     t=0.0
#     nPoints_x=34
#     nPoints_y=34
#     nPoints_z=34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#              for k in range(nPoints_z):
#                   x[i,j,k,0] = i/(nPoints_x-1.0)
#                   x[i,j,k,1] = j/(nPoints_y-1.0)
#                   x[i,j,k,2] = k/(nPoints_z-1.0)
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     canalyticalSolutions.LinearADR_Sine_totalVelocity(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                    x[:,0,0,0],
#                                    x[0,:,0,1],

#                               title='LinearADR_Sine_TotalVelocity'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((17,),'d')
# #     rwork[0]=2.0*3.14159
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=0.0
# #     rwork[4]=1.0
# #     rwork[5]=0.0
# #     rwork[6]=0.0
# #     rwork[7]=0.01
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=0.0
# #     rwork[11]=0.01
# #     rwork[12]=0.0
# #     rwork[13]=0.0
# #     rwork[14]=0.0
# #     rwork[15]=0.01
# #     rwork[16]=1.0
# #     t=0.0
# #     nPoints_x=51
# #     nPoints_y=51
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #                 x[i,j,0] = i/(nPoints_x-1.0)
# #                 x[i,j,1] = j/(nPoints_y-1.0)
# #     u = numpy.zeros((nPoints_x,nPoints_y,),'d')
# #     canalyticalSolutions.LinearADR_Sine_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     gnuplot('set xlabel "x"')
# #     gnuplot('set ylabel "y"')
# #     gnuplot('set zlabel "z"')
# #     gnuplot.splot(Gnuplot.GridData(u[:,:],
# #                                    x[:,0,0],
# #                                    x[0,:,1],
# #
# #                                    title='LinearADR_Sine_r: reaction'))
# #     #raw_input('Please press return to continue... \n')
#     # #---------------LinearADR_SINE_r-3D-----------------------------
#     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((17,),'d')
# #     rwork[0]=2.0*3.14159
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=0.0
# #     rwork[4]=1.0
# #     rwork[5]=0.0
# #     rwork[6]=0.0
# #     rwork[7]=0.01
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=0.0
# #     rwork[11]=0.01
# #     rwork[12]=0.0
# #     rwork[13]=0.0
# #     rwork[14]=0.0
# #     rwork[15]=0.01
# #     rwork[16]=1.0
# #     nPoints=101
# #     t=0.0
# #     x = numpy.zeros((nPoints,3),'d')
# #     for i in range(nPoints):
# #             x[i,0] = i*(1.0/(nPoints-1.0))
# #     u = numpy.zeros((nPoints,),'d')
# #     canalyticalSolutions.LinearADR_Sine_dr(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot.plot(Gnuplot.Data(x[:,0],
# #                               u,
# #
# #                               title='LinearADR_Sine_dr'))
# #     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((1,),'d')
#     rwork[0]=5
#     nPoints = 101
#     t=0.0
#     x = numpy.zeros((nPoints,3),'d')
#     for i in range(nPoints):
#         x[i,0] = i*(1.0/(nPoints-1.0))
#     u = numpy.zeros(x.shape[0],'d')
#     canalyticalSolutions.poissonsEquationExp1D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,

#                               title='poissonsEquationExp1D'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((1,),'d')
#     rwork[0]=5
#     nPoints_x = 51
#     nPoints_y = 51
#     t=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.poissonsEquationExp2D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u,
#                                    x[:,0,0],
#                                    x[0,:,1],
#                                    title='poissonsEquationExp2D'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((1,),'d')
#     rwork[0]=5
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.poissonsEquationExp3D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     slice=nPoints_x/2
#     gnuplot.splot(Gnuplot.GridData(u[slice,:,:],
#                                    x[slice,:,0,1],
#                                    x[slice,0,:,2],
#                                    title='poissonsEquationExp3D-x/2'))
#     gnuplot1 = Gnuplot.Gnuplot()
#     gnuplot1('set parametric')
#     gnuplot1('set data style lines')
#     gnuplot1('set hidden')
#     gnuplot1('set contour base')
#     slice=nPoints_y/2
#     gnuplot1.splot(Gnuplot.GridData(u[:,slice,:],
#                                    x[:,slice,0,0],
#                                    x[0,slice,:,2],
#                                    title='poissonsEquationExp3D-y/2'))
#     gnuplot2 = Gnuplot.Gnuplot()
#     gnuplot2('set parametric')
#     gnuplot2('set data style lines')
#     gnuplot2('set hidden')
#     gnuplot2('set contour base')
#     slice=nPoints_z/2
#     gnuplot2.splot(Gnuplot.GridData(u[:,:,slice],
#                                    x[:,0,slice,0],
#                                    x[0,:,slice,1],
#                                    title='poissonsEquationExp3D-z/2'))
#     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((1,),'d')
# #     rwork[0]=5
# #     nPoints = 101
# #     t=0.0
# #     x = numpy.zeros((nPoints,3),'d')
# #     for i in range(nPoints):
# #         x[i,0] = i*(1.0/(nPoints-1.0))
# #     u = numpy.zeros(x.shape[0],'d')
# #     canalyticalSolutions.poissonsEquationExp1D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot.plot(Gnuplot.Data(x[:,0],
# #                               u,
# #
# #                               title='poissonsEquationExp1D_r'))
# #     #raw_input('Please press return to continue... \n')
# #     # # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((1,),'d')
# #     rwork[0]=5
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     t=0.0
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     canalyticalSolutions.poissonsEquationExp2D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     gnuplot.splot(Gnuplot.GridData(u,
# #                                    x[:,0,0],
# #                                    x[0,:,1],
# #                                    title='poissonsEquationExp2D_r'))
# #     #raw_input('Please press return to continue... \n')
# #     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((1,),'d')
# #     rwork[0]=5
# #     nPoints_x = 34
# #     nPoints_y = 34
# #     nPoints_z = 34
# #     t=0.0
# #     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #                     for k in range(nPoints_z):
# #                             x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
# #                             x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
# #                             x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     canalyticalSolutions.poissonsEquationExp3D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     slice=nPoints_x/2
# #     gnuplot.splot(Gnuplot.GridData(u[slice,:,:],
# #                                    x[slice,:,0,1],
# #                                    x[slice,0,:,2],
# #                                    title='poissonsEquationExp3D_r-x/2'))
# #     gnuplot1 = Gnuplot.Gnuplot()
# #     gnuplot1('set parametric')
# #     gnuplot1('set data style lines')
# #     gnuplot1('set hidden')
# #     gnuplot1('set contour base')
# #     slice=nPoints_y/2
# #     gnuplot1.splot(Gnuplot.GridData(u[:,slice,:],
# #                                    x[:,slice,0,0],
# #                                    x[0,slice,:,2],
# #                                    title='poissonsEquationExp3D_r-y/2'))
# #     gnuplot2 = Gnuplot.Gnuplot()
# #     gnuplot2('set parametric')
# #     gnuplot2('set data style lines')
# #     gnuplot2('set hidden')
# #     gnuplot2('set contour base')
# #     slice=nPoints_z/2
# #     gnuplot2.splot(Gnuplot.GridData(u[:,:,slice],
# #                                    x[:,0,slice,0],
# #                                    x[0,:,slice,1],
# #                                    title='poissonsEquationExp3D_r-z/2'))
# #     #raw_input('Please press return to continue... \n')
# #     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((1,),'d')
# #     rwork[0]=5
# #     nPoints_x = 101
# #     t=0.0
# #     x = numpy.zeros((nPoints_x,3),'d')
# #     for i in range(nPoints_x):
# #         x[i,0] = i*(1.0/(nPoints_x-1.0))
# #     u = numpy.zeros(x.shape[0],'d')
# #     canalyticalSolutions.poissonsEquationExp3D_dr(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot.plot(Gnuplot.Data(x[:,0],
# #                               u,
# #                               title='poissonsEquationExp3D_dr'))
# #     #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     logEvent("LinearAD_DiracIC - solution")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((10),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     nPoints_x = 51
#     nPoints_y = 51
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#                       x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     tt =[0.25, 0.75]
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearAD_DiracIC(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u,
#                                        x[:,0,0],
#                                        x[0,:,1],

#                                        title='LinearAD_DiracIC'))
#         #raw_input('Please press return to continue... \n')
#     #---------------------------------------------------------------
#     logEvent("Gradient")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((10),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     tt =[0.25, 0.75]
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearAD_DiracIC_du(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                        x[:,0,0,0],
#                                        x[0,:,0,1],

#                                        title='LinearAD_DiracIC_du: Gradient'))
#         #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     logEvent("Advective Velocity")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((10),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     tt =[0.25, 0.75]
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearAD_DiracIC_advectiveVelocity(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                        x[:,0,0,0],
#                                        x[0,:,0,1],

#                                        title='LinearAD_DiracIC_AdvectiveVelocity'))
#         #raw_input('Please press return to continue... \n')
#     # # #---------------------------------------------------------------
#     logEvent("Diffusive Velocity")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((10),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     tt =[0.25, 0.75]
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearAD_DiracIC_diffusiveVelocity(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                        x[:,0,0,0],
#                                        x[0,:,0,1],

#                                        title='LinearAD_DiracIC_DiffusiveVelocity'))
#         #raw_input('Please press return to continue... \n')
#     # # #---------------------------------------------------------------
#     logEvent("Total Velocity")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((10),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     u = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     tt =[0.25, 0.75]
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearAD_DiracIC_totalVelocity(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u[:,:,0,0],
#                                        x[:,0,0,0],
#                                        x[0,:,0,1],

#                                        title='LinearAD_DiracIC_TotalVelocity'))
#         #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
#     # #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((11),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     rwork[10]=1.0
#     nPoints_x = 51
#     nPoints_y = 51
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#                       x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     tt =[0.25, 0.75]
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.LinearADR_Decay_DiracIC(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u,
#                                        x[:,0,0],
#                                        x[0,:,1],

#                                        title='LinearADR_Decay_DiracIC'))
#         #raw_input('Please press return to continue... \n')
#     # #---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((11),'d')
# #     rwork[0]=1.0
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=1.0
# #     rwork[4]=0.01
# #     rwork[5]=0.0
# #     rwork[6]=0.1
# #     rwork[7]=0.0
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=1.0
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #                             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #                             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     tt =[0.25, 0.75]
# #     for t in tt:
# #         logEvent("time t = "+`t`)
# #         canalyticalSolutions.LinearADR_Decay_DiracIC_r(iwork,rwork,t,x,u)
# #         gnuplot.splot(Gnuplot.GridData(u,
# #                                        x[:,0,0],
# #                                        x[0,:,1],
# #
# #                                        title='LinearADR_Decay_DiracIC_r'))
# #         #raw_input('Please press return to continue... \n')
# #     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((11),'d')
# #     rwork[0]=1.0
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=1.0
# #     rwork[4]=0.01
# #     rwork[5]=0.0
# #     rwork[6]=0.1
# #     rwork[7]=0.0
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=1.0
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     t=0.25
# #     canalyticalSolutions.LinearADR_Decay_DiracIC_dr(iwork,rwork,t,x,u)
# #     gnuplot.splot(Gnuplot.GridData(u,
# #                                    x[:,0,0],
# #                                    x[0,:,1],
# #
# #                                    title='LinearADR_Decay_DiracIC_dr'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     # ---------------------------------------------------------------
#     logEvent("NonlinearADR_Decay_DiracIC")
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((12),'d')
#     rwork[0]=1.0
#     rwork[1]=0.0
#     rwork[2]=0.0
#     rwork[3]=1.0
#     rwork[4]=0.01
#     rwork[5]=0.0
#     rwork[6]=0.1
#     rwork[7]=0.0
#     rwork[8]=0.0
#     rwork[9]=0.0
#     rwork[10]=1.0
#     rwork[11]=2.0
#     nPoints_x = 51
#     nPoints_y = 51
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#                       x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     tt =[0.25, 0.75]
#     for t in tt:
#         logEvent("time t = "+`t`)
#         canalyticalSolutions.NonlinearADR_Decay_DiracIC(iwork,rwork,t,x,u)
#         gnuplot.splot(Gnuplot.GridData(u,
#                                        x[:,0,0],
#                                        x[0,:,1],

#                                        title='NonlinearADR_Decay_DiracIC'))
#         #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((12),'d')
# #     rwork[0]=1.0
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=1.0
# #     rwork[4]=0.01
# #     rwork[5]=0.0
# #     rwork[6]=0.1
# #     rwork[7]=0.0
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=1.0
# #     rwork[11]=2.0
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     tt =[0.25, 0.75]
# #     for t in tt:
# #         logEvent("time t = "+`t`)
# #         canalyticalSolutions.NonlinearADR_Decay_DiracIC_r(iwork,rwork,t,x,u)
# #         gnuplot.splot(Gnuplot.GridData(u[:,:],
# #                                        x[:,0,0],
# #                                        x[0,:,1],
# #
# #                                        title='NonlinearADR_Decay_DiracIC_r'))
#         #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     iwork = numpy.zeros((1,),'i')
# #     rwork = numpy.zeros((12),'d')
# #     rwork[0]=1.0
# #     rwork[1]=0.0
# #     rwork[2]=0.0
# #     rwork[3]=1.0
# #     rwork[4]=0.01
# #     rwork[5]=0.0
# #     rwork[6]=0.1
# #     rwork[7]=0.0
# #     rwork[8]=0.0
# #     rwork[9]=0.0
# #     rwork[10]=1.0
# #     rwork[11]=2.0
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')

# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     tt =[0.25, 0.75]
# #     for t in tt:
# #         logEvent("time t = "+`t`)
# #         canalyticalSolutions.NonlinearADR_Decay_DiracIC_dr(iwork,rwork,t,x,u)
# #         gnuplot.splot(Gnuplot.GridData(u[:,:],
# #                                        x[:,0,0],
# #                                        x[0,:,1],
# #
# #                                        title='NonlinearADR_Decay_DiracIC_dr'))
#         #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     iwork[0]=1
#     rwork = numpy.zeros((1,),'d')
#     nPoints = 101
#     t=0.0
#     x = numpy.zeros((nPoints,3),'d')
#     for i in range(nPoints):
#         x[i,0] = i*(1.0/(nPoints-1.0))
#     u = numpy.zeros(x.shape[0],'d')
#     canalyticalSolutions.diffusionSin1D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot.plot(Gnuplot.Data(x[:,0],
#                               u,

#                               title='diffusionSin1D'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((2,),'i')
#     iwork[0]=1
#     iwork[1]=1
#     rwork = numpy.zeros((1,),'d')
#     nPoints_x = 51
#     nPoints_y = 51
#     t=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
#             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
#     u = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.diffusionSin2D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot.splot(Gnuplot.GridData(u,
#                                    x[:,0,0],
#                                    x[0,:,1],
#                                    title='diffusionSin2D'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((3,),'i')
#     iwork[0]=1
#     iwork[1]=1
#     iwork[2]=1
#     rwork = numpy.zeros((1,),'d')
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     import time
#     a = time.clock()
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#               for k in range(nPoints_z):
#                       x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
#                       x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
#                       x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
#     b=time.clock()
#     print b- a
#     a = time.clock()
#     u = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.diffusionSin3D(iwork,rwork,t,x,u)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     slice=nPoints_x/2
#     gnuplot.splot(Gnuplot.GridData(u[slice,:,:],
#                                    x[slice,:,0,1],
#                                    x[slice,0,:,2],
#                                    title='diffusionSin3D-x/2'))
#     gnuplot1 = Gnuplot.Gnuplot()
#     gnuplot1('set parametric')
#     gnuplot1('set data style lines')
#     gnuplot1('set hidden')
#     gnuplot1('set contour base')
#     slice=nPoints_y/2
#     gnuplot1.splot(Gnuplot.GridData(u[:,slice,:],
#                                    x[:,slice,0,0],
#                                    x[0,slice,:,2],
#                                    title='diffusionSin3D-y/2'))
#     gnuplot2 = Gnuplot.Gnuplot()
#     gnuplot2('set parametric')
#     gnuplot2('set data style lines')
#     gnuplot2('set hidden')
#     gnuplot2('set contour base')
#     slice=nPoints_z/2
#     gnuplot2.splot(Gnuplot.GridData(u[:,:,slice],
#                                    x[:,0,slice,0],
#                                    x[0,:,slice,1],
#                                    title='diffusionSin3D-z/2'))
#     b=time.clock()
#     print b- a
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((1,),'i')
# #     iwork[0]=1
# #     rwork = numpy.zeros((1,),'d')
# #     nPoints = 101
# #     t=0.0
# #     x = numpy.zeros((nPoints,3),'d')
# #     for i in range(nPoints):
# #         x[i,0] = i*(1.0/(nPoints-1.0))
# #     u = numpy.zeros(x.shape[0],'d')
# #     canalyticalSolutions.diffusionSin1D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot.plot(Gnuplot.Data(x[:,0],
# #                               u,
# #
# #                               title='diffusionSin1D_r'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((2,),'i')
# #     iwork[0]=1
# #     iwork[1]=1
# #     rwork = numpy.zeros((1,),'d')
# #     nPoints_x = 51
# #     nPoints_y = 51
# #     t=0.0
# #     x = numpy.zeros((nPoints_x,nPoints_y,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #             x[i,j,0] = i*(1.0/(nPoints_x-1.0))
# #             x[i,j,1] = j*(1.0/(nPoints_y-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     canalyticalSolutions.diffusionSin2D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     # print u.shape
# #     # print x.shape
# #     # print x[:,0,0]
# #     # print x[0,:,1]
# #     # print u
# #     gnuplot.splot(Gnuplot.GridData(u,
# #                                    x[:,0,0],
# #                                    x[0,:,1],
# #                                    title='diffusionSin2D_r'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
# #     iwork = numpy.zeros((3,),'i')
# #     iwork[0]=1
# #     iwork[1]=1
# #     iwork[2]=1
# #     rwork = numpy.zeros((1,),'d')
# #     nPoints_x = 34
# #     nPoints_y = 34
# #     nPoints_z = 34
# #     t=0.0
# #     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
# #     for i in range(nPoints_x):
# #         for j in range(nPoints_y):
# #                     for k in range(nPoints_z):
# #                             x[i,j,k,0] = i*(1.0/(nPoints_x-1.0))
# #                             x[i,j,k,1] = j*(1.0/(nPoints_y-1.0))
# #                             x[i,j,k,2] = k*(1.0/(nPoints_z-1.0))
# #     u = numpy.zeros(x.shape[:-1],'d')
# #     canalyticalSolutions.diffusionSin3D_r(iwork,rwork,t,x,u)
# #     gnuplot.clear()
# #     gnuplot('set parametric')
# #     gnuplot('set data style lines')
# #     gnuplot('set hidden')
# #     gnuplot('set contour base')
# #     slice=nPoints_x/2
# #     gnuplot.splot(Gnuplot.GridData(u[slice,:,:],
# #                                    x[slice,:,0,1],
# #                                    x[slice,0,:,2],
# #                                    title='diffusionSin3D_r-x/2'))
# #     gnuplot1 = Gnuplot.Gnuplot()
# #     gnuplot1('set parametric')
# #     gnuplot1('set data style lines')
# #     gnuplot1('set hidden')
# #     gnuplot1('set contour base')
# #     slice=nPoints_y/2
# #     gnuplot1.splot(Gnuplot.GridData(u[:,slice,:],
# #                                    x[:,slice,0,0],
# #                                    x[0,slice,:,2],
# #                                    title='diffusionSin3D_r-y/2'))
# #     gnuplot2 = Gnuplot.Gnuplot()
# #     gnuplot2('set parametric')
# #     gnuplot2('set data style lines')
# #     gnuplot2('set hidden')
# #     gnuplot2('set contour base')
# #     slice=nPoints_z/2
# #     gnuplot2.splot(Gnuplot.GridData(u[:,:,slice],
# #                                    x[:,0,slice,0],
# #                                    x[0,:,slice,1],
# #                                    title='diffusionSin3D_r-z/2'))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((8,),'d')
#     rwork[0]=0.1
#     rwork[1]=0.0
#     rwork[2]=0.0
#     nPoints_x = 31
#     nPoints_y = 31
#     nPoints_z = 31
#     t=0.0
#     dx = 1.0/(nPoints_x-1.0)
#     dy = 1.0/(nPoints_y-1.0)
#     dz = 1.0/(nPoints_z-1.0)
#     rwork[3]=0.1
#     rwork[4]=0.5
#     rwork[5]=0.5
#     rwork[6]=0.5
#     rwork[7]=1.003e-3
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                           x[i,j,k,0] = i*dx
#                           x[i,j,k,1] = j*dy
#                           x[i,j,k,2] = k*dz
#     p = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.STflowSphere_P(iwork,rwork,t,x,p)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     slice=nPoints_x/2
#     gnuplot.splot(Gnuplot.GridData(p[slice,:,:],
#                                    x[slice,:,0,1],
#                                    x[slice,0,:,2],
#                                    title='STflowSphere_P: yz plane'))
#     #gnuplot.hardcopy('STflowSphere_P_x0.5.eps', eps=1,enhanced=1,color=1)
#     gnuplot1 = Gnuplot.Gnuplot()
#     gnuplot1('set parametric')
#     gnuplot1('set data style lines')
#     gnuplot1('set hidden')
#     gnuplot1('set contour base')
#     slice=nPoints_y/2
#     gnuplot1.splot(Gnuplot.GridData(p[:,slice,:],
#                                     x[:,slice,0,0],
#                                     x[0,slice,:,2],
#                                     title='STflowSphere_P: xz plane'))
#     #gnuplot1.hardcopy('STflowSphere_P_y0.5.eps', eps=1,enhanced=1,color=1)
#     gnuplot2 = Gnuplot.Gnuplot()
#     gnuplot2('set parametric')
#     gnuplot2('set data style lines')
#     gnuplot2('set hidden')
#     gnuplot2('set contour base')
#     slice=nPoints_z/2
#     gnuplot2.splot(Gnuplot.GridData(p[:,:,slice],
#                                     x[:,0,slice,0],
#                                     x[0,:,slice,1],
#                                     title='STflowSphere_P: xy plane'))
#     #gnuplot2.hardcopy('STflowSphere_P_z0.5.eps', eps=1,enhanced=1,color=1)
#     #raw_input('Please press return to continue... \n')
#     #
#     ux = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.STflowSphere_Vx(iwork,rwork,t,x,ux)
#     uy = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.STflowSphere_Vy(iwork,rwork,t,x,uy)
#     uz = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.STflowSphere_Vz(iwork,rwork,t,x,uz)
#     gnuplot2 = Gnuplot.Gnuplot()
#     gnuplot2('set parametric')
#     gnuplot2('set mxtics 5')
#     gnuplot2('set mytics 5')
#     gnuplot2('set grid xtics ytics mxtics mytics')
#     gnuplot2('set xlabel "x"')
#     gnuplot2('set ylabel "y"')
#     slice=nPoints_z/2
#     gnuplot2.plot(Gnuplot.Data(numpy.reshape(x[:,:,slice,0],(nPoints_x*nPoints_y,)),
#                                numpy.reshape(x[:,:,slice,1],(nPoints_x*nPoints_y,)),
#                                numpy.reshape(ux[:,:,slice],(nPoints_x*nPoints_y,)),
#                                numpy.reshape(uy[:,:,slice],(nPoints_x*nPoints_y,)),
#                                title='STflowSphere : xy plane',
#                                ))
#     #gnuplot2.hardcopy('STflowSphere_V_z0.5.eps', eps=1,enhanced=1,color=1)
#     #raw_input('Please press return to continue... \n')
#     gnuplot2('set xlabel "x"')
#     gnuplot2('set ylabel "z"')
#     slice=nPoints_y/2
#     gnuplot2.plot(Gnuplot.Data(numpy.reshape(x[:,slice,:,0],(nPoints_x*nPoints_z,)),
#                                numpy.reshape(x[:,slice,:,2],(nPoints_x*nPoints_z,)),
#                                numpy.reshape(ux[:,slice,:],(nPoints_x*nPoints_z,)),
#                                numpy.reshape(uz[:,slice,:],(nPoints_x*nPoints_z,)),
#                                title='STflowSphere : xz plane',
#                                ))
#     #gnuplot2.hardcopy('STflowSphere_V_y0.5.eps', eps=1,enhanced=1,color=1)
#     #raw_input('Please press return to continue... \n')
#     gnuplot2('set xlabel "y"')
#     gnuplot2('set ylabel "z"')
#     slice=nPoints_x/2
#     gnuplot2.plot(Gnuplot.Data(numpy.reshape(x[slice,:,:,1],(nPoints_z*nPoints_y,)),
#                                numpy.reshape(x[slice,:,:,2],(nPoints_z*nPoints_y,)),
#                                numpy.reshape(uy[slice,:,:],(nPoints_z*nPoints_y,)),
#                                numpy.reshape(uz[slice,:,:],(nPoints_z*nPoints_y,)),
#                                title='STflowSphere : yz plane',
#                                ))
#     #gnuplot2.hardcopy('STflowSphere_V_x0.5.eps', eps=1,enhanced=1,color=1)
#     #raw_input('Please press return to continue... \n')
#     #---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((5,),'d')
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     dx = 1.0/(nPoints_x-1.0)
#     dy = 1.0/(nPoints_y-1.0)
#     dz = 1.0/(nPoints_z-1.0)
#     rwork[0]=0.1
#     rwork[1]=dy*(nPoints_y-1.0)
#     rwork[2]=0.0
#     rwork[3]=0.0
#     rwork[4]=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                 x[i,j,k,0] = i*dx
#                 x[i,j,k,1] = j*dy
#                 x[i,j,k,2] = k*dz
#     ux = numpy.zeros(x.shape[:-1],'d')
#     uy = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.CouetteFlow(iwork,rwork,t,x,ux)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 3')
#     gnuplot('set ytics nomirror')
#     gnuplot('set xtics nomirror')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "h"')
#     gnuplot('set zlabel "z"')
#     slice=nPoints_z/2
#     gnuplot.splot(Gnuplot.GridData(ux[:,:,slice],
#                                     x[:,0,slice,0],
#                                     x[0,:,slice,1],
#                                     title='Couette Flow: xy plane'))
#     #raw_input('Please press return to continue... \n')
#     gnuplot.clear()
#     gnuplot('set title "Couette Flow"')
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 3')
#     gnuplot('set ytics nomirror')
#     gnuplot('set xtics nomirror')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "h"')
#     gnuplot('set zlabel "z"')
#     gnuplot('set arrow 1 from graph 0,1.02 to 1.15,1.02 back head size graph 0.03,30.0 lw 3 lt -1')
#     gnuplot.plot(Gnuplot.Data(numpy.reshape(x[:,:,slice,0],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(x[:,:,slice,1],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(ux[:,:,slice], (nPoints_x*nPoints_y,)),
#                               numpy.reshape(uy[:,:,slice],(nPoints_x*nPoints_y,)),
#                               ))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((7,),'d')
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     dx = 1.0/(nPoints_x-1.0)
#     dy = 1.0/(nPoints_y-1.0)
#     dz = 1.0/(nPoints_z-1.0)
#     rwork[0]=dy*(nPoints_y-1.0)
#     rwork[1]=0.001
#     rwork[2]=-1.0
#     # rwork[3]=83.3
#     rwork[4]=0.0
#     rwork[5]=0.0
#     rwork[6]=0.0
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                            x[i,j,k,0] = i*dx
#                            x[i,j,k,1] = j*dy
#                            x[i,j,k,2] = k*dz
#     ux = numpy.zeros(x.shape[:-1],'d')
#     uy = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.PoiseuilleFlow(iwork,rwork,t,x,ux)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 3')
#     gnuplot('set ytics nomirror')
#     gnuplot('set xtics nomirror')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "h"')
#     gnuplot('set zlabel "z"')
#     slice=nPoints_z/2
#     gnuplot.splot(Gnuplot.GridData(ux[:,:,slice],
#                                     x[:,0,slice,0],
#                                     x[0,:,slice,1],
#                                     title='Poiseuille Flow'))
#     #raw_input('Please press return to continue... \n')
#     gnuplot.clear()
#     gnuplot('set title "   Poiseuille Flow"')
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 3')
#     gnuplot('set ytics nomirror')
#     gnuplot('set xtics nomirror')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "h"')
#     gnuplot('set zlabel "z"')
#     gnuplot.plot(Gnuplot.Data(numpy.reshape(x[:,:,slice,0],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(x[:,:,slice,1],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(ux[:,:,slice], (nPoints_x*nPoints_y,)),
#                               numpy.reshape(uy[:,:,slice],(nPoints_x*nPoints_y,)),
#                               ))
#     #raw_input('Please press return to continue... \n')
#     # ---------------------------------------------------------------
#     # ---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((7,),'d')
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     dx = 1.0/(nPoints_x-1.0)
#     dy = 1.0/(nPoints_y-1.0)
#     dz = 1.0/(nPoints_z-1.0)
#     rwork[0]=dy*(nPoints_y-1.0)/2.0
#     rwork[1]=0.001
#     rwork[2]=-1.0
#     # rwork[3]=24.5
#     rwork[4]=0.5
#     rwork[5]=0.5
#     rwork[6]=0.5
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                            x[i,j,k,0] = i*dx
#                            x[i,j,k,1] = j*dy
#                            x[i,j,k,2] = k*dz
#     ux = numpy.zeros(x.shape[:-1],'d')
#     uy = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.PoiseuillePipeFlow(iwork,rwork,t,x,ux)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 7')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "Diameter"')
#     gnuplot('set zlabel "z"')
#     slice=nPoints_z/2
#     gnuplot.splot(Gnuplot.GridData(ux[:,:,slice],
#                                     x[:,0,slice,0],
#                                     x[0,:,slice,1],
#                                     title='Poiseuille Pipe Flow'))
#     #raw_input('Please press return to continue... \n')
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 7')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "Diameter"')
#     gnuplot('set zlabel "z"')
#     gnuplot('set ytics nomirror')
#     gnuplot('set xtics nomirror')
#     gnuplot('set title "Poiseuille Pipe Flow"')
#     gnuplot.plot(Gnuplot.Data(numpy.reshape(x[:,:,slice,0],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(x[:,:,slice,1],(nPoints_x*nPoints_y,)),
#                               numpy.reshape(ux[:,:,slice], (nPoints_x*nPoints_y,)),
#                               numpy.reshape(uy[:,:,slice],(nPoints_x*nPoints_y,)),
#                               ))
#     #raw_input('Please press return to continue... \n')
#     # # ---------------------------------------------------------------
#     # # ---------------------------------------------------------------
#     iwork = numpy.zeros((1,),'i')
#     rwork = numpy.zeros((8,),'d')
#     nPoints_x = 34
#     nPoints_y = 34
#     nPoints_z = 34
#     t=0.0
#     dx = 1.0/(nPoints_x-1.0)
#     dy = 1.0/(nPoints_y-1.0)
#     dz = 1.0/(nPoints_z-1.0)
#     rwork[0]=dy*(nPoints_y-1.0)/2.0
#     rwork[1]=0.001
#     rwork[2]=-1.0
#     # rwork[3]=24.5
#     rwork[4]=1.0
#     rwork[5]=0.5
#     rwork[6]=0.5
#     rwork[7]=0.5
#     x = numpy.zeros((nPoints_x,nPoints_y,nPoints_z,3),'d')
#     for i in range(nPoints_x):
#         for j in range(nPoints_y):
#             for k in range(nPoints_z):
#                            x[i,j,k,0] = i*dx
#                            x[i,j,k,1] = j*dy
#                            x[i,j,k,2] = k*dz
#     ux = numpy.zeros(x.shape[:-1],'d')
#     uy = numpy.zeros(x.shape[:-1],'d')
#     canalyticalSolutions.PoiseuillePipeFlow_P(iwork,rwork,t,x,ux)
#     gnuplot.clear()
#     gnuplot('set parametric')
#     gnuplot('set multiplot')
#     gnuplot('set data style lines')
#     gnuplot('set hidden')
#     gnuplot('set contour base')
#     gnuplot('set border 7')
#     gnuplot('set xlabel "x"')
#     gnuplot('set ylabel "Diameter"')
#     gnuplot('set zlabel "z"')
#     slice=nPoints_z/2
#     gnuplot.splot(Gnuplot.GridData(ux[:,:,slice],
#                                     x[:,0,slice,0],
#                                     x[0,:,slice,1],
#                                     title='Poiseuille Pipe Flow: Pressure'))
#     #raw_input('Please press return to continue... \n')
#     # # ---------------------------------------------------------------
#     # # ---------------------------------------------------------------
