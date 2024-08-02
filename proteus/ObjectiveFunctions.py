"""
A module of objective functions

.. moduleauthor:: John Chrispell
.. inheritance-diagram:: proteus.ObjectiveFunctions
   :parts: 1
"""
from .Optimizers import ObjectiveFunction_base
#John Chrispell, Summer 07
import math
import numpy

# A simple function used to test fminbound().
class SimpelFunc(ObjectiveFunction_base):
    def __init__(self,LHS_x,RHS_x):
        ObjectiveFunction_base.__init__(self,LHS_x,RHS_x)
    def getResidual(self,x):
        return (x-0.5)**2 + 0.04

# A nasty simple function used to test fminbound().
class SimpelFunc2(ObjectiveFunction_base):
    def __init__(self,LHS_x,RHS_x):
        ObjectiveFunction_base.__init__(self,LHS_x,RHS_x)
    def getResidual(self,x):
        return 0.25*math.sin(12.0*math.pi*x)+0.5*x

# A nasty simple function used to test fminbound().
class SimpelFunc3(ObjectiveFunction_base):
    def __init__(self,LHS_x,RHS_x):
        ObjectiveFunction_base.__init__(self,LHS_x,RHS_x)
    def getResidual(self,x):
        return ((x**2.0)/(x**2.0+0.5*((1.0-x)**2.0)))

# BuckleyLeverett Function with the choice of alpha=1/2
class BuckleyLeverett(object):
    def __init__(self,alpha):
        self.alpha=alpha
    def getFlux(self,s):
        return ((s**2.0)/(s**2.0+self.alpha*((1.0-s)**2.0)))

# A function used to find solution for the nonconvex scalar
# Riemann problem with arbitrary data LHS_s, and RHS_s.
class OsherFunc(ObjectiveFunction_base):
    def __init__(self,LHS_s,RHS_s,fFunc,t,x):
        ObjectiveFunction_base.__init__(self,LHS_s,RHS_s)
        self.fFunc = fFunc
        self.xi = x/t
        if(LHS_s < RHS_s):
            self.getResidual = self.Argmin
        else:
            self.getResidual = self.Argmax
    def Argmin(self,s):
        return self.fFunc.getFlux(s) - self.xi*s
    def Argmax(self,s):
        return self.xi*s - self.fFunc.getFlux(s)

# A function used to find the solution for a nonconvex scalar
# Riemann problem with arbitrary data LHS_s, and RHS_s using the
# coefficients passed to it from a numerical solver.
class OsherFuncCoef(ObjectiveFunction_base):
    def __init__(self,LHS_s,RHS_s,fFunc,t,x,useShallowCopy=True):
        ObjectiveFunction_base.__init__(self,LHS_s,RHS_s)
        self.xi = x/t
        self.t = t
        self.x = x
        # The next bit is a one dimensional dictionary for the coefficents
        # that can be passed to this Osher Function
        self.c = {('u',0):numpy.zeros((1,),'d'),
                ('m',0):numpy.zeros((1,),'d'),
                ('dm',0,0):numpy.zeros((1,),'d'),
                ('f',0):numpy.zeros((1,1),'d'),
                ('df',0,0):numpy.zeros((1,1),'d'),
                ('a',0,0):numpy.zeros((1,1,1),'d'),
                ('da',0,0,0):numpy.zeros((1,1,1),'d'),
                ('phi',0):numpy.zeros((1,),'d'),
                ('dphi',0,0):numpy.zeros((1,),'d'),
                ('r',0):numpy.zeros((1,),'d'),
                ('dr',0,0):numpy.zeros((1,),'d'),
                ('H',0):numpy.zeros((1,),'d'),
                ('dH',0,0):numpy.zeros((1,1),'d')}
        if useShallowCopy:
            self.fFunc = fFunc
        else:
            import copy
            self.fFunc  = copy.deepcopy(fFunc)
            self.fFunc.initializeElementQuadrature(self.t,self.c)
        if(LHS_s < RHS_s):
            self.getResidual = self.Argmin
        else:
            self.getResidual = self.Argmax
    def Argmin(self,s):
        self.c[('u',0)][0]=s
        self.fFunc.evaluate(self.t,self.c)
        m = self.c[('m',0)][0]
        f = self.c[('f',0)][0,0]
        return f - self.xi*m
    def Argmax(self,s):
        self.c[('u',0)][0]=s
        self.fFunc.evaluate(self.t,self.c)
        m = self.c[('m',0)][0]
        f = self.c[('f',0)][0,0]
        return self.xi*m - f