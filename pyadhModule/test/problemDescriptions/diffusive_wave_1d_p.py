from pyadh import *
from pyadh.default_p import *
import numpy
import math
"""
A module for 1D diffusive wave model verification problems.

This module contains several implementations of the diffusive wave equation coefficients and several analytical solutions.
Copyright (c) 2009 by Steve Mattis
"""

class DiffusiveWave_1D(TransportCoefficients.TC_base):
    """
    A class implementing the coefficients of the diffusive wave equation in 1D.

    This implements the regularized formulation in M. Santillana's thesis (ref).
    """
    def __init__(self,alpha=5.0/3.0,gamma=0.5,epsilon=1.0e-6,bathymetryFunc=None, reactionFunc=None, analyticalSoln=None):
        """
        Construct a coefficients object given the parameters of the Manning/Chezy formula and the regularization parameter.

        Optionally provide a function for bathymetry, right hand side (source), and an analytical solution. The
        """
        self.alpha=alpha; self.gamma=gamma; self.epsilon=epsilon
        TransportCoefficients.TC_base.__init__(self, nc=1, 
                                               variableNames=['H'],
                                               mass={0:{0:'linear'}},
                                               advection={0:{0:'constant'}},#mwf added for dg to test
                                               diffusion={0:{0:{0:'nonlinear'}}},
                                               potential={0:{0:'u'}},
                                               reaction={0:{0:'constant'}},
                                               hamiltonian={})
        self.bathymetryFunc=bathymetryFunc
        self.reactionFunc=reactionFunc
        self.analyticalSoln=analyticalSoln
    def initializeMesh(self,mesh):
        """
        Set the y component of the 1D mesh using the bathymetry function.
        """
        if self.bathymetryFunc:
            #mesh vertices are held in vector nodeArray which is nNodes_global x 3
            for nN in range(mesh.nNodes_global):
                x = mesh.nodeArray.flat[nN*3:(nN+1)*3]; z = self.bathymetryFunc(x)
                mesh.nodeArray.flat[nN*3+1]= z
    def initializeElementQuadrature(self,t,cq):
        """
        Set the y component of the element quadrature points using the bathymetry function.
        """
        if self.bathymetryFunc:
            #cq is elementQuadrature dictionary so x is nElements_global x nQuadraturePoints_element x 3 
            nPoints = cq['x'].shape[0]*cq['x'].shape[1]
            for i in range(nPoints):
                x = cq['x'].flat[i*3:(i+1)*3]; z =  self.bathymetryFunc(x)
                cq['x'].flat[i*3+1]=z
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        """
        Set the y component of the element boundary quadrature points using the bathymetry function.
        """
        if self.bathymetryFunc:
            #cebq is quadrature dictionary for local "faces" on each element so 
            #  x is nElements_global x nElementBoundaries_element x nQuadraturePoints_elementBoundary x 3 
            nPoints = cebq['x'].shape[0]*cebq['x'].shape[1]*cebq['x'].shape[2]
            for i in range(nPoints):
                x = cebq['x'].flat[i*3:(i+1)*3]; z =  self.bathymetryFunc(x)
                cebq['x'].flat[i*3+1]=z
            #cebq_global is quadrature dictionary for global "faces" so 
            #x is nElementBoundaries_global x nQuadraturePoints_elementBoundary x 3
            nPoints = cebq_global['x'].shape[0]*cebq_global['x'].shape[1]
            for i in range(nPoints):
                x = cebq_global['x'].flat[i*3:(i+1)*3]; z =  self.bathymetryFunc(x)
                cebq_global['x'].flat[i*3+1]=z

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        """
        Set the y component of the exterior element boundary quadrature points using the bathymetry function.
        """
        if self.bathymetryFunc:
            #cebqe is quadrature dictionary for global exterior "faces" 
            #  x is nExteriorElementBoundaries_global x nQuadraturePoints_elementBoundary x 3 
            nPoints = cebqe['x'].shape[0]*cebqe['x'].shape[1]
            for i in range(nPoints):
                x = cebqe['x'].flat[i*3:(i+1)*3]; z =  self.bathymetryFunc(x)
                cebqe['x'].flat[i*3+1]=z

    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        """
        Set the y component of the generlatized interpolation points using the bathymetry function.
        """
        if self.bathymetryFunc:
            #cip is dictionary for interpolation points 
            #  so x is nElements x number of interpolation points (usual dof) per element
            nPoints = cip['x'].shape[0]*cip['x'].shape[1]
            for i in range(nPoints):
                x = cip['x'].flat[i*3:(i+1)*3]; z =  self.bathymetryFunc(x)
                cip['x'].flat[i*3+1]=z

    def evaluate(self,t,c):
        """
        Evaluated the coefficients of the 1D diffusive wave model.
        """
        c[('m',0)][:]=c[('u',0)];
        c[('dm',0,0)].fill(1.0) 
        #
        for i in range(len(c[('u',0)].flat)):
            depth = max(c[('u',0)].flat[i]-c['x'].flat[i*3+1],0.0)
            #
            #c[('u',0)].flat[i]=max(c['x'].flat[i*3+1], c[('u',0)].flat[i])
            #
            c[('a',0,0)].flat[i] = (self.gamma/(self.alpha+self.gamma))**(-self.gamma)*(numpy.power(depth,self.alpha) + self.epsilon)/(numpy.power(numpy.abs(c[('grad(u)',0)].flat[i*nd:nd*(i+1)]),1.0-self.gamma) + self.epsilon)
            c[('da',0,0,0)].flat[i] = (self.gamma/(self.alpha+self.gamma))**(-self.gamma)*(self.alpha*numpy.power(depth,self.alpha-1.0) + self.epsilon)/(numpy.power(numpy.abs(c[('grad(u)',0)].flat[i*nd:nd*(i+1)]),1.0-self.gamma)+ self.epsilon)
            #c[('a',0,0)].flat[i] = (self.alpha+1.0)*numpy.power(depth,self.alpha)
            #c[('da',0,0,0)].flat[i] = self.alpha*(self.alpha+1.0)*numpy.power(depth,self.alpha-1.0)
            ####c[('a',0,0)].flat[i] = (self.gamma/(self.alpha+self.gamma))**(-self.gamma)*(numpy.power(depth,self.alpha) + self.epsilon)/(numpy.power(numpy.abs(c[('grad(u)',0)].flat[i*nd:nd*(i+1)]),1.0-self.gamma) + self.epsilon)
            ####c[('da',0,0,0)].flat[i] = (self.gamma/(self.alpha+self.gamma))**(-self.gamma)*(self.alpha*numpy.power(depth,self.alpha-1.0) + self.epsilon)/(numpy.power(numpy.abs(c[('grad(u)',0)].flat[i*nd:nd*(i+1)]),1.0-self.gamma)+ self.epsilon)
        if self.reactionFunc != None:
            for i in range(len(c[('u',0)].flat)):
                c[('r',0)].flat[i] = self.reactionFunc(c['x'].flat[i*3],t)

class DiffusiveWave_1D_nlp(DiffusiveWave_1D):
    """
    A modified diffusive wave model that uses a nonlinear potential. Only works for constant bathymetry right now!!!
    """
    def __init__(self,alpha=5.0/3.0,gamma=0.5,epsilon=1.0e-6,bathymetryFunc=None, reactionFunc=None, analyticalSoln=None):
        """
        Construct a coefficients object given the parameters of the Manning/Chezy formula and the regularization parameter. Set the potential to be nonlinear.

        Optionally provide a function for bathymetry, right hand side (source), and an analytical solution. 
        """
        self.alpha=alpha; self.gamma=gamma; self.epsilon=epsilon
        TC_base.__init__(self, 
                         variableNames=['H'],
                         mass={0:{0:'linear'}},
                         advection={},
                         diffusion={0:{0:{0:'nonlinear'}}},
                         potential={0:{0:'nonlinear'}},
                         reaction={0:{0:'constant'}},
                         hamiltonian={})
        self.bathymetryFunc=bathymetryFunc
        self.reactionFunc=reactionFunc
        self.analyticalSoln=analyticalSoln

    def evaluate(self,t,c):
        c[('m',0)][:]=c[('u',0)];
        c[('dm',0,0)].fill(1.0) 
        #
        for i in range(len(c[('u',0)].flat)):
            depth = max(c[('u',0)].flat[i]-c['x'].flat[i*3+1],0.0)
            #this only works for constant bathymetry
            if c.has_key(('grad(u)',0)):
                c[('a',0,0)].flat[i] = (self.gamma/(self.alpha+self.gamma))**(-self.gamma)/numpy.power(max(self.epsilon,numpy.abs(c[('grad(u)',0)].flat[i])),1.0-self.gamma)
                c[('da',0,0,0)].flat[i] = 0.0
            c[('phi',0)].flat[i] = numpy.power(depth,self.alpha+1.0)/(self.alpha+1.0)
            c[('dphi',0,0)].flat[i] = numpy.power(depth,self.alpha)
        if self.reactionFunc != None:
            for i in range(len(c[('u',0)].flat)):
                c[('r',0)].flat[i] = self.reactionFunc(c['x'].flat[i*3],t)

##         if self.analyticalSoln !=None:
##             for i in range(len(c[('u',0)].flat)):
##                 c[('analytical')].flat[i*3+2] = self.analyticalSoln(c['x'].flat[i*3],t)
##         print str(c[('analytical')])

##################################
nd=1
alpha=5.0/3.0
gamma=0.5
epsilon=1.0e-5

#unit interval, currently don't really support x != (0.0)
L=[1.0]; x=[0.0]
T=0.5
#diverging planes example
L=[15.0]; 
T=0.1

domain = Domain.RectangularDomain(L=L,x=x)


def artificialrightside(x,t):
    initial_time=9.5
    if x+2.0> 0.0:
        val= 2.0*(t+initial_time)*numpy.power(x+2.0,0.5)+7.0/(12.0*numpy.power(2.0,0.5))*numpy.power((100.0-numpy.power(t+initial_time,2.0)), 13.0/6.0)*numpy.power(x+2.0, -5.0/12.0)
        #val=0.0
    else:
        val =0.0
    return val

def flatBathymetry(x):
    return 0.0
def singlePeak(x):
    return -abs(x[0]-4.0) + 4.0



def BarenblattSoln(x,t,alpha=5.0/3.0, gamma=0.5, initialtime=3.0, C=10.0):

    m=1.0 +alpha/gamma
    k=(m*gamma -1.0)/(m*(gamma +1.0))*math.pow((1.0/(gamma*(m+1.0))), 1.0/gamma)

    phi=math.pow(abs(x-L[0]/2.0)*math.pow(initialtime+t, -1.0/(gamma*(m+1.0))), (gamma +1.0)/gamma)
    if (C-k*phi)>0.0:
        val= math.pow(initialtime+t,-1.0/(gamma*(m+1.0)))*math.pow((C-k*phi),gamma/(m*gamma -1.0))
    else:
        val=0.0
    return val






#simple flat domain
coefficients=DiffusiveWave_1D(alpha=alpha, gamma=gamma, epsilon=epsilon, bathymetryFunc = flatBathymetry, reactionFunc = None, analyticalSoln= BarenblattSoln)
#coefficients=DiffusiveWave_1D(alpha=alpha, gamma=gamma, epsilon=epsilon, bathymetryFunc = flatBathymetry, reactionFunc = None, analyticalSoln= BarenblattSoln)
#diverging plane example
coefficients=DiffusiveWave_1D(alpha=alpha, gamma=gamma, epsilon=epsilon, bathymetryFunc = flatBathymetry, reactionFunc = artificialrightside)

def getDBC(x,flag):
    if x[0]==0.0 or x[0]==L[0]:
        return lambda x,t:0.0
def getNoBC(x,flag):
    pass
def getrightsideDBC(x,flag):
    if x[0] in [0.0, L[0]]:
        initial_time=9.5
        if x[0]+2.0 > 0.0:
            return lambda x,t: (100.0-(initial_time+t)*(initial_time+t))*math.sqrt(x[0]+2.0)
        else:
            return lambda x,t: 0.0
        
#dirichletConditions = {0:getDBC}
#dirichletConditions = {0:getNoBC}
dirichletConditions = {0:getrightsideDBC}

class Tent:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 1.0 - abs(x[0]-0.5)*2.0

class GaussHill:
    def __init__(self,xc=0.5,bathymetryFunc=flatBathymetry):
        self.xc = xc; self.bathymetry = bathymetryFunc
    def uOfXT(self,x,t):
        z=self.bathymetry(x)
        return z + 2.0*exp(-2.0*(x[0]-self.xc)**2)

class BarenblattInitial:
    def __init__(self,alpha=5.0/3.0,gamma=0.5, initialtime=3.0,C=10.0):
        self.alpha=alpha; self.gamma=gamma; self.initialtime=initialtime; self.C=C;
        self.m=1.0 +self.alpha/self.gamma
        self.k=(self.m*self.gamma -1.0)/(self.m*(self.gamma +1.0))*math.pow((1.0/(self.gamma*(self.m+1))), 1.0/self.gamma)
        pass
    def uOfXT(self,x,t):
        #val= math.power(self.initialtime, -1.0/(self.gamma*(self.m+1.0)))*math.pow(self.C-self.k*abs(x[0]-L[0]/2.0)*math.pow(self.initialtime,-1.0/(self.gamma*(self.m+1.0))),self.gamma/(self.m*self.gamma-1.0))
        phi=math.pow(abs(x[0]-L[0]/2.0)*math.pow(self.initialtime, -1.0/(self.gamma*(self.m+1.0))), (self.gamma +1.0)/self.gamma)
        print phi
        print self.k
        print self.m
        if self.C-self.k*phi>0.0:
            val= math.pow(self.initialtime,-1.0/(self.gamma*(self.m+1.0)))*math.pow((self.C-self.k*phi),self.gamma/(self.m*self.gamma -1.0))
        else:
            val=0.0
        
        if val >= 0.0:
            return val
        else:
            return 0.0

class rightsideInitial:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        initial_time=9.5
        if x[0]+2.0 > 0.0:
            return (100.0-initial_time*initial_time)*math.sqrt(x[0]+2.0)
        else:
            return 0.0
        return 0.0

        
    

#initialConditions= {0:GaussHill(xc=L[0]/2.0,bathymetryFunc=singlePeak)}#{0:Tent()}
#initialConditions ={0:BarenblattInitial()}
initialConditions ={0:rightsideInitial()}
fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

#############################################

## class BarenblattSolution:
##     def __init__(self,alpha=5.0/3.0,gamma=0.5,C=20.0, initial_time=0.1):
##         self.alpha=alpha; self.gamma=gamma; self.C=C;self.initial_time=initial_time
##         self.m=1.0 +self.alpha/self.gamma
##         self.k=(self.m*self.gamma -1.0)/(self.m*(self.gamma +1.0))*math.pow((1.0/(self.gamma*(self.m+1))), 1.0/self.gamma)
##         pass
##     def uOfXT(self,x,t):
##         phi=math.pow(abs(x[0]-L[0]/2.0)*math.pow(t+self.initial_time, -1.0/(self.gamma*(self.m+1.0))), (self.gamma +1.0)/self.gamma)
##         if self.C-self.k*phi>0.0:
##             val= math.pow(t+self.initial_time,-1.0/(self.gamma*(self.m+1.0)))*math.pow((self.C-self.k*phi),self.gamma/(self.m*self.gamma -1.0))
##         else:
##             val=0.0
##         return val

class zeroSolution:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class BarenblattSolution:
    def __init__(self,alpha=5.0/3.0,gamma=0.5, initialtime=3.0,C=10.0):
        self.alpha=alpha; self.gamma=gamma; self.initialtime=initialtime; self.C=C;
        self.m=1.0 +self.alpha/self.gamma
        self.k=(self.m*self.gamma-1.0)/(self.m*(self.gamma+1.0))*math.pow((1.0/(self.gamma*(self.m+1.0))), 1.0/self.gamma)
        pass
    def uOfXT(self,x,t):
        phi=math.pow(abs(x[0]-L[0]/2.0)*math.pow(self.initialtime+t, -1.0/(self.gamma*(self.m+1.0))), (self.gamma+1.0)/self.gamma)
        if self.C-self.k*phi>0.0:
            val= math.pow(self.initialtime+t,-1.0/(self.gamma*(self.m+1.0)))*math.pow((self.C-self.k*phi),self.gamma/(self.m*self.gamma-1.0))
        else:
            val=0.0
        return val

class BarenblattSolution_noGrad:
    def __init__(self,alpha=5.0/3.0, initialtime=3.0e-2,C=10.0):
        self.p = alpha+1.0
        self.initialtime=initialtime
        self.C=C
    def uOfXT(self,x,t):
        T= t + self.initialtime
        xi = abs(x[0] - L[0]/2.0)
        PHI = (self.p-1.0) * (xi/(T**(1.0/(self.p+1))))**2.0 / (2.0*self.p*(self.p+1.0))
        if self.C-PHI>0.0:
            u = T**(-1.0/(self.p+1.0)) * (self.C - PHI)**(1.0/(self.p-1.0))
        else:
            u = 0.0
        return u

class rightsideSolution:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        initial_time=9.5
        if x[0]+2.0 > 0.0:
            return (100.0-(initial_time+t)*(initial_time+t))*math.sqrt(x[0]+2.0)
        else:
            return 0.0
        return 0.0
    
analyticalSolution={0:rightsideSolution()}
#analyticalSolution={0:BarenblattSolution()}
#analyticalSolution={0:BarenblattSolution_noGrad()}
initialConditions ={0:analyticalSolution[0]}
