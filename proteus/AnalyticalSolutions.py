"""
Classes representing analytical solutions of differential and partial differential equations.

.. inheritance-diagram:: proteus.AnalyticalSolutions
   :parts: 1
"""
from math import *
from .EGeometry import *
from .LinearAlgebraTools import *
from .Profiling import logEvent
from . import canalyticalSolutions

class AS_base(object):
    """
    The base class for general analytical solutions, u(x,y,z,t)

    For basic error calculations only the :func:`uOfXT` member need be
    overridden. The vectorized member functions such as
    :func:`getUOfXT` are provided for future optimizations but are
    implemented in :class:`AS_base` in a simple, inefficient way.
    """

    def uOfXT(self, X, T):
        """
        Return the solution at a 3D point, X, and a time, T
        """
        pass

    def duOfXT(self, X, T):
        """
        Return the gradient of the solution at a 3D point, X, and a time, T

        .. todo:: Change du to Du
        """
        pass

    def rOfUX(self, u, x):
        """
        Return the reaction term for the (manufactured) solution at a 3D point, X, and a time, T
        """
        pass

    def drOfUX(self, u, x):
        """
        Return the derivative of the reaction term for the (manufactured) solution at a 3D point, X, and a time, T
        """
        pass

    def advectiveFluxOfX(self, x):
        """
        Return the advective flux at a 3D point, X

        .. todo:: Change advective flux to depend on T
        """
        pass

    def diffusiveFluxOfX(self, x):
        """
        Return the diffusive flux at a 3D point, X

        .. todo:: Change advective flux to depend on T
        """
        pass

    def totalFluxOfX(self, x):
        """
        Return the total flux at a 3D point, X

        .. todo:: Decide if we need total flux
        """
        pass

    def getUOfXT(self, X, T, U):
        """
        Set the solution, U, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.uOfXT, X, T, U)

    def getDUOfXT(self, X, T, DU):
        """
        Set the gradient of the solution, DU, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.duOfXT, X, T, DU)

    def getROfUX(self, X, T, R):
        """
        Set the reaction term, R, at an array of 3D points, X, and a time, T
        """
        self.coeff_p_to_v(self.rOfXT, U, X, T, R)

    def getDROfUX(self, X, T, DR):
        """
        Set the derivative of the reaction term, DR, at an array of 3D points, X, and a time, T
        """
        self.coeff_p_to_v(self.drOfXT, U, X, T, DR)

    def getAdvectiveFluxOfXT(self, X, T, F):
        """
        Set the advective flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.advectiveFluxOfXT, X, T, F)

    def getDiffusiveFluxOfX(self, X, T, F):
        """
        Set the diffusive flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.diffusiveFluxOfXT, X, T, F)

    def getTotalFluxOfX(self, X, T, F):
        """
        Set the total flux, F, at an array of 3D points, X, and a time, T
        """
        self.p_to_v(self.totalFluxOfXT, X, T, F)

    def p_to_v(self, func, X, T, funcValVec):
        """
        Calculate a point-wise function over and array of points.
        """
        nPoints = 1
        for d in X.shape[:-1]:
            nPoints *= d
        for i in range(nPoints):
            funcValVec[i] = funcValVec(X.flat[i*3:i*3+3], T)


class DAE_base(AS_base):
    """
    The base class for differential-algebraic equation solutions

    To use this class override :func:`uOfT`
    """

    def uOfXT(self, X, T):
        return self.uOfT(self, T)

    def uOfT(self, T):
        pass


class NonlinearDAE(DAE_base):
    r"""
    The exact solution of the simple nonlinear DAE

    .. math: u_t  = - a \max(u,0)^p \quad u(0) = 1
    """

    def __init__(self, a, p):
        """
        Set the coefficients a and p
        """
        self.a_ = a
        self.p_ = p
        if self.p_ == 1:
            self.func = lambda t: exp(-self.a_*t)
        else:
            q = 1.0/(1.0-self.p_)
            self.func = lambda t: max(1.0 - (1.0 - self.p_)*self.a_*t, 0.0)**q

    def uOfT(self, t):
        return self.func(t)

    def fOfUT(self, u, t):
        return -self.a_*max(u, 0.0)**self.p_


class SteadyState(AS_base):
    """
    Based class for steady-state solutions u(x)

    Override :func:`uOfX` to define steady state solutions.
    """

    def uOfX(self, X):
        pass

    def uOfXT(self, X, T=None):
        return self.uOfX(X)


class LinearAD_SteadyState(SteadyState):
    r"""
    The solution of the steady linear advection-diffusion equation

    The boundary value problem is

    .. math: (bu - au_x)_x = 0 \quad u(0) = 1 \quad u(1) = 0
    """

    def __init__(self, b=1.0, a=5.0e-1):
        self.b_ = b
        self.a_ = a
        if b != 0.0:
            self.D_ = (1.0/(exp(b/a)-1.0))
        else:
            self.D_ = 0.0
        self.C_ = -self.D_*exp(b/a)

    def uOfX(self, X):
        x = X[0]
        if self.D_ != 0.0:
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

    def __init__(self, b=1.0, q=2, a=5.0e-1, r=1):
        LinearAD_SteadyState.__init__(self, b, a)
        self.r_ = r
        self.q_ = q
        if (q == 2 and r == 1):
            if b != 0.0:
                def f(rtmC):
                    return rtmC*tanh(b*rtmC/a) - 1.0

                def df(rtmC):
                    return rtmC*(1.0/cosh(b*rtmC/a)**2)*(b/a) + tanh(b*rtmC/a)
                logEvent("Solving for sqrt(-C) for q=2,r=1")
                rtmC = sqrt(1.5)
                while abs(f(rtmC)) > 1.0e-8:
                    rtmC -= f(rtmC)/df(rtmC)
                logEvent("sqrt(-C)="+repr(rtmC))
                self.rtmC_ = rtmC
                self.sol_ = lambda x: self.rtmC_ * \
                    tanh((-self.b_*self.rtmC_/self.a_)*(x - 1.0))
            else:
                self.sol_ = lambda x: 1.0 - x
        elif (q == 1 and r == 2):
            logEvent("Solving for C in q=1,r=2")

            def f(C):
                return 2.0*C*(log(C-1.0) - log(C)) + 2.0 + self.b_/self.a_

            def df(C):
                return 2.0*(log(C-1.0) - log(C)) + 2.0*C*(1.0/(C-1.0) - 1.0/C)
            C = 1.0 + 1.0e-10
            f0 = f(C)
            print(f0)
            while abs(f(C)) > (1.0e-7*abs(f0) + 1.0e-7):
                dC = -f(C)/df(C)
                logEvent("dc")
                print(dC)
                Ctmp = C + dC
                while (abs(f(Ctmp)) > 0.99*abs(f(C))
                       or Ctmp <= 1.0):
                    print(f(Ctmp))
                    print(f(C))
                    logEvent("ls")
                    dC *= 0.9
                    Ctmp = C + dC
                logEvent("out")
                print(Ctmp)
                print(f(Ctmp))
                print(df(Ctmp))
                C = Ctmp
            logEvent("C="+repr(C))
            self.nlC_ = C
            self.nlD_ = 0.5*(2.0*C*log(C*(C-1)) -
                             4.0*C + 2.0 - self.b_/self.a_)
            logEvent("D="+repr(self.nlD_))
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

    def uOfX(self, X):
        x = X[0]
        if (self.q_ == 2 and self.r_ == 1):
            return self.sol_(x)
        elif (self.q_ == 1 and self.r_ == 2):
            def f(u):
                return (2.0*(self.nlC_*log(self.nlC_-u) -
                             (self.nlC_-u)) -
                        self.nlD_) - self.b_*x/self.a_

            def df(u):
                return (2.0*self.nlC_/(u-self.nlC_)+2.0)
            u = LinearAD_SteadyState.uOfX(self, X)
            f0 = f(u)
            while abs(f(u)) > 1.0e-6*abs(f0) + 1.0e-6:
                u -= f(u)/df(u)
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
                 b=numpy.array((1.0, 0.0, 0.0)),
                 a=numpy.array(((1.0e-2, 0.0, 0.0),
                                (0.0, 1.0e-2, 0.0),
                                (0.0, 0.0, 1.0e-2))),
                 c=1.0,
                 omega=numpy.array((2*pi, 0.0, 0.0)),
                 omega0=0.0):
        """
        Set the coefficients a,b,c,omega, and omega0
        """
        self.b_ = b
        self.a_ = a
        self.c_ = c
        self.omega_ = omega
        self.omega0_ = omega0
        self.D_ = - numpy.dot(b, omega)
        self.E_ = - numpy.dot(numpy.dot(a, omega), omega) - c

    def uOfX(self, x):
        return sin(numpy.dot(self.omega_, x[:self.omega_.shape[0]]) + self.omega0_)

    def duOfX(self, x):
        return self.omega_ * cos(numpy.dot(self.omega_, x[:self.omega_.shape[0]])+self.omega0_)

    def rOfUX(self, u, x):
        return self.c_*u + self.D_*cos(numpy.dot(self.omega_, x[:self.omega_.shape[0]])+self.omega0_) + self.E_*sin(numpy.dot(self.omega_, x[:self.omega_.shape[0]]) + self.omega0_)

    def drOfUX(self, u, x):
        return self.c_

    def advectiveFluxOfX(self, x):
        return self.b_*self.uOfX(x)

    def diffusiveFluxOfX(self, x):
        return -numpy.dot(self.a_, self.duOfX(x))

    def totalFluxOfX(self, x):
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

    def __init__(self, K=10.0, nd=3):
        self.K_ = K
        self.nd_ = nd

    def uOfX(self, X):
        if self.nd_ == 1:
            x = X[0]
            return self.K_*x*(1.0-x)*exp(x**2)
        elif self.nd_ == 2:
            x = X[0]
            y = X[1]
            return self.K_*x*(1.0-x)*y*(1.0-y)*exp(x**2 + y**2)
        else:
            x = X[0]
            y = X[1]
            z = X[2]
            return self.K_*x*(1.0-x)*y*(1.0-y)*z*(1.0-z)*exp(x**2 + y**2 + z**2)

    def rOfUX(self, u, X):
        if self.nd_ == 1:
            x = X[0]
            return self.K_*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*(1.0-x)*x - 2.0)*exp(x**2)
        elif self.nd_ == 2:
            x = X[0]
            y = X[1]
            return self.K_*(y*(1.0-y)*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*x*(1.0-x) - 2.0) +
                            x*(1.0-x)*(4.0*(1.0-y)*y**3 - 4.0*y**2 + 6.0*y*(1.0-y) - 2.0))*exp(x**2 + y**2)
        else:
            x = X[0]
            y = X[1]
            z = X[2]
            return self.K_*(y*(1.0-y)*z*(1.0-z)*(4.0*(1.0-x)*x**3 - 4.0*x**2 + 6.0*x*(1.0-x) - 2.0) +
                            x*(1.0-x)*z*(1.0-z)*(4.0*(1.0-y)*y**3 - 4.0*y**2 + 6.0*y*(1.0-y) - 2.0) +
                            x*(1.0-x)*y*(1.0-y)*(4.0*(1.0-z)*z**3 - 4.0*z**2 + 6.0*z*(1.0-z) - 2.0))*exp(x**2 + y**2 + z**2)

    def drOfUX(self, u, x):
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
                 b=numpy.array((1.0, 0.0, 0.0)),
                 a=1.0e-2,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0, 0.0, 0.0))):
        self.n_ = n
        self.u0_ = u0
        self.x0_ = x0
        self.b_ = b
        self.a_ = a
        self.tStart = tStart

    def uOfXT(self, x, T):
        t = T + self.tStart
        y = (x[:self.b_.shape[0]] - self.x0_[:self.b_.shape[0]] - self.b_*t)
        exp_arg = numpy.dot(y,y)/(4.0*self.a_*t)
        if exp_arg > 100:
            return 0.0
        else:
            return self.u0_*exp(-exp_arg) / (4.0*self.a_*pi*t)**(self.n_/2.0)

    def duOfXT(self, x, T):
        t = T + self.tStart
        y = (x[:self.b_.shape[0]] - self.x0_[:self.b_.shape[0]] - self.b_*t)
        return self.uOfXT(x, T)*2.0*y/(4.0*self.a_*t)

    def advectiveFluxOfXT(self, x, T):
        return self.b_*self.uOfXT(x, T)

    def diffusiveFluxOfXT(self, x, T):
        return -self.a_*self.duOfXT(x, T)

    def totalFluxOfXT(self, x, T):
        return self.advectiveFluxOfXT(x, T) + self.diffusiveFluxOfXT(x, T)


class LinearADR_Decay_DiracIC(LinearAD_DiracIC):
    r"""
    The exact solution of

    .. math:: u_t + \nabla \cdot (bu - a \nabla u) + cu= 0

    on an infinite domain with Dirac initial data.  Also
    returns the fluxes (by inheritance).
    """

    def __init__(self,
                 n=1.0,
                 b=numpy.array((1.0, 0.0, 0.0)),
                 a=1.0e-2,
                 c=1.0,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0, 0.0, 0.0))):
        LinearAD_DiracIC.__init__(self, n, b, a, tStart, u0, x0)
        self.c_ = c

    def uOfXT(self, x, T):
        t = T + self.tStart
        return LinearAD_DiracIC.uOfXT(self, x, T)*exp(- self.c_*t)

    def rOfUXT(self, u, x, T):
        return self.c_*u

    def drOfUXT(self, u, x, T):
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
                 b=numpy.array((1.0, 0.0, 0.0)),
                 a=1.0e-2,
                 c=1.0,
                 d=2.0,
                 tStart=0.0,
                 u0=0.1,
                 x0=numpy.array((0.0, 0.0, 0.0))):
        LinearADR_Decay_DiracIC.__init__(self, n, b, a, c, tStart, u0, x0)
        self.d_ = d

    def uOfXT(self, x, T):
        t = T + self.tStart
        u1 = LinearAD_DiracIC.uOfXT(self, x, T)
        if u1 > 0.0:
            return u1*exp(-(2.0*self.c_*t*pow(u1,self.d_-1.0))/(self.d_+1.0))
        else:
            return u1

    def rOfUXT(self, u, x, T):
        return self.c_*(u**self.d_)

    def drOfUXT(self, u, x, T):
        return self.d_*self.c_*(u**(self.d_-1))


class PlaneCouetteFlow_u(SteadyState):
    """
    The exact solution for the u component of  velocity in plane Couette Flow
    """
    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0, 0.0, 0.0]):
        self.iwork = numpy.zeros((1,), 'i')
        self.rwork = numpy.array([upperPlateVelocity,
                                  plateSeperation,
                                  origin[0],
                                  origin[1],
                                  origin[2]])

    def uOfX(self, x):
        uList = numpy.array([[0.0]])
        xList = numpy.array([x])
        t = 0.0
        nPoints = 1
        canalyticalSolutions.PlaneCouetteFlow_u(self.iwork, self.rwork, t, xList, uList)
        return uList.flat[0]


class PlaneCouetteFlow_v(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Couette Flow
    """

    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0, 0.0, 0.0]):
        pass

    def uOfX(self, x):
        return 0.0


class PlaneCouetteFlow_p(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Couette Flow
    """

    def __init__(self,
                 plateSeperation=1.0,
                 upperPlateVelocity=0.01,
                 origin=[0.0, 0.0, 0.0]):
        pass

    def uOfX(self, x):
        return 0.0


class PlanePoiseuilleFlow_u(SteadyState):
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q=0.0,
                 origin=[0.0, 0.0, 0.0]):
        # allocated iwork and rwork arrays and load in values
        # ( (shape_0,shape_1,...), typecode)
        self.iwork = numpy.zeros((1,), 'i')
        self.rwork = numpy.array([plateSeperation,
                                  mu,
                                  grad_p,
                                  q,
                                  origin[0],
                                  origin[1],
                                  origin[2]])  # numpy.array([component_0,component_2,...]) builds a double* initialized to values in list [...]
    # wrapped "vectorized" analytical solutions that Moira wrote with a point-wise evaluation

    def uOfX(self, x):
        # array of doubles for storing solution (just the solution at a single point)
        uList = numpy.array([0.0])
        # array of double*'s for storing spatial locations (just a single point in 3D in this case)
        xList = numpy.array([x])
        # time at which solution is evaluated (this is a steady state problem)
        t = 0.0
        nPoints = 1  # just a single point
        # make call to wrapped C function
        canalyticalSolutions.PlanePoiseuilleFlow_u(self.iwork, self.rwork, t, xList, uList)
        return uList[0]  # return value  of solution at this point


class PlanePoiseuilleFlow_v(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """

    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q=0.0,
                 origin=[0.0, 0.0, 0.0]):
        pass

    def uOfX(self, x):
        return 0.0


class PlanePoiseuilleFlow_p(SteadyState):
    """
    The exact solution for the v component of  velocity in plane Poiseuille Flow
    """

    def __init__(self,
                 plateSeperation=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 q=0.0,
                 origin=[0.0, 0.0, 0.0]):
        self.grad_p = grad_p

    def uOfX(self, x):
        return self.grad_p*x[0]


# encapsulate analytical solver for BuckleyLeverett Riemann problems
from proteus.ObjectiveFunctions import OsherFuncCoef
from proteus.Optimizers import fminbound


class Buckley_Leverett_RiemannSoln(AS_base):
    def __init__(self, coefficients, uLeft=1.0, uRight=0.0, t0=0.0, x0=0.0, T=0.5, ftol=1.0e-8,
                 useShallowCopyCoef=True):
        # ought to start with range and spline for eval later I guess
        self.coefficients = coefficients
        self.uL = uLeft
        self.uR = uRight
        self.t0 = t0
        self.x0 = x0
        self.T = T
        self.ftol = ftol
        self.riemF = OsherFuncCoef(self.uL, self.uR, self.coefficients, self.T-self.t0, self.x0,
                                   useShallowCopy=useShallowCopyCoef)
        self.solver = fminbound(self.riemF, self.ftol)
        # mwf if want to look at coefficients
        # self.generateAplot(x0,x0+1.0,101,T)

    def uOfXT(self, x, t):
        if abs(t-self.t0) < 1.0e-7:
            if x[0]-self.x0 <= 0.0:
                return self.uL
            return self.uR
        self.riemF.xi = (x[0]-self.x0)/(t-self.t0)
        u, f = self.solver.solve(Guess_x=0.5*(self.uL+self.uR))
        return u

    def uOfX(self, x):
        return self.uOfXT(x, self.T)

    def generateAplot(self, xLeft, xRight, nnx, t, filename='Buckley_Leverett_ex.dat'):
        """
        save exact solution to a file
        """
        dx = (xRight-xLeft)/(nnx-1.)
        fout = open(filename, 'w')
        for i in range(nnx):
            x = (xLeft+dx*i, 0.0, 0.0)
            u = self.uOfXT(x, t)
            fout.write('%12.5e  %12.5e \n' % (x[0], u))
        #
        fout.close()


class PlaneBase(SteadyState):
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0, 1.0, 1.0]):
        self.plane_n = numpy.array([cos(plane_theta)*sin(plane_phi),
                                    sin(plane_theta)*sin(plane_phi),
                                    cos(plane_phi)])
        if (plane_phi > -pi/2.0 and plane_phi < pi/2.0):
            if plane_phi == 0.0:
                logEvent("plate is in x-y plane")
                v_phi = pi/2.0
            else:
                v_phi = atan(
                    -1.0/(tan(plane_phi)*cos(v_theta-plane_theta)))
        else:
            logEvent(
                "plate is is parallel to z-axis, using angle of velocity with z-axis instead of angle with x-axis")
            v_theta = plane_theta - pi/2.0
        self.v_n = numpy.array([cos(v_theta)*sin(v_phi),
                                sin(v_theta)*sin(v_phi),
                                cos(v_phi)])
        minZ = (0.0, numpy.array([0.0, 0.0, 0.0]))
        maxZ = (0.0, numpy.array([0.0, 0.0, 0.0]))
        for i in [0, 1]:
            for j in [0, 1]:
                for k in [0, 1]:
                    x = numpy.array([i*L[0], j*L[1], k*L[2]])
                    Z = numpy.dot(x, self.plane_n)
                    if Z < minZ[0]:
                        minZ = (Z, x)
                    if Z > maxZ[0]:
                        maxZ = (Z, x)
        self.Zshift = minZ[0]
        self.Zsep = maxZ[0] - minZ[0]
        self.mu = mu
        self.grad_p = grad_p
        self.v_norm = v_norm

    def Z(self, x):
        return numpy.dot(x, self.plane_n) - self.Zshift

    def X(self, x):
        return numpy.dot(x, self.v_n)

    def U(self, x):
        Z = self.Z(x)
        return (Z*self.v_norm/self.Zsep) - self.grad_p*Z*(self.Zsep-Z)/(2.0*self.mu)


class PlanePoiseuilleFlow_u2(PlaneBase):
    """
    The exact solution for the u component of  velocity in plane Poiseuille Flow
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 plane_theta=0.0,
                 plane_phi=math.pi/2.0,
                 v_theta=math.pi/2.0,
                 v_phi=None,
                 v_norm=1.0,
                 mu=1.0,
                 grad_p=1.0,
                 L=[1.0, 1.0, 1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)

    def uOfX(self, x):
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
                 L=[1.0, 1.0, 1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)

    def uOfX(self, x):
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
                 L=[1.0, 1.0, 1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)

    def uOfX(self, x):
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
                 L=[1.0, 1.0, 1.0]):
        PlaneBase.__init__(self,
                           plane_theta,
                           plane_phi,
                           v_theta,
                           v_phi,
                           v_norm,
                           mu,
                           grad_p,
                           L)

    def uOfX(self, x):
        return self.grad_p*self.X(x)


class VortexDecay_u(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n = n
        self.Re = Re

    def uOfXT(self, x, t):
        # return -cos(x[0])*sin(x[1])*exp(-2.0*t)
        return -cos(self.n*pi*x[0])*sin(self.n*pi*x[1])*exp(-2.0 * (self.n*pi)**2 * t / self.Re)


class VortexDecay_v(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n = n
        self.Re = Re

    def uOfXT(self, x, t):
        # return sin(x[0])*cos(x[1])*exp(-2.0*t)
        return sin(self.n*pi*x[0])*cos(self.n*pi*x[1])*exp(-2.0 * (self.n*pi)**2 * t / self.Re)


class VortexDecay_p(AS_base):
    """
    The exact solution for the u component of  velocity in the vortex  decay problem
    """
    # constructor or initializer function, reads in parameters needed in the iwork and rwork arrays

    def __init__(self,
                 n=2,
                 Re=100.0):
        self.n = n
        self.Re = Re

    def uOfXT(self, x, t):
        # return -self.Re*0.25*(cos(2.0*x[0])+cos(2.0*x[1]))*exp(-4.0*t)
        return -0.25*(cos(2.0*self.n*pi*x[0]) + cos(2.0*self.n*pi*x[1]))*exp(-4.0 * (self.n*pi)**2 * t / self.Re)