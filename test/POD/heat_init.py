#!/usr/bin/env python
"""
Fine-scale heat equation solver
The equation is du/du - Laplace u + u + f(x,y,z,t) = 0
"""
from proteus.iproteus import *

class Heat(TransportCoefficients.TC_base):
    """
    The coefficients of the Heat equation
    The exact solution is 128*(1-t)*x*(1-x)*y*(1-y)*z*(1-z)
    """
    def __init__(self,A):
        TransportCoefficients.TC_base.__init__(self, 
                         nc=1, #number of components
                         variableNames=['u'],
			 mass = {0:{0:'linear'}},
                         diffusion = {0:{0:{0:'constant'}}},#constant means not a function of the solution
                         potential = {0:{0:'u'}}, #define the potential for the diffusion term to be the solution itself
                         reaction  = {0:{0:'linear'}})
        self.A=A;

    def rofx1(self,u,x,t):
        return -128.0*2.0*(1.0-t)*(x[...,1]*(1.0-x[...,1])*x[...,2]*(1.0-x[...,2]) + x[...,0]*(1.0-x[...,0])*x[...,2]*(1.0-x[...,2]) + x[...,0]*(1.0-x[...,0])*x[...,1]*(1.0-x[...,1])) + \
	       128.0*t*x[...,0]*(1.0-x[...,0])*x[...,1]*(1.0-x[...,1])*x[...,2]*(1.0-x[...,2]) + u
	       	
    def evaluate(self,t,c):
	c[('m',0)][:] = c[('u',0)]
	c[('dm',0,0)][:] = 1.0
        c[('a',0,0)][...,0] = self.A[0][0]
        c[('a',0,0)][...,4] = self.A[1][1]
        c[('a',0,0)][...,8] = self.A[2][2]
	c[('r',0)][:] = self.rofx1(c[('u',0)], c['x'], t)
	c[('dr',0,0)][:] = 1.0

# Below a loop may be used as an alternative to the next to last line
#	dim = c[('u',0)].shape[-1]
#	dim2= len(c[('u',0)].flat)/dim
#	for pi in range(dim2):
#	    for pi1 in range(dim):
#		u_pi=c[('u',0)][pi][pi1]
#		r_pi = self.rofx1(u_pi, c['x'][pi][pi1],t)
#		c[('r',0)].flat[pi*dim+pi1] = r_pi

#use numpy for evaluations
import numpy as np

#convenience function for holding the boundary condition
def constant_zero(x,t):
    return 0.0

# this function's job is to return another function holding the Dirichlet boundary conditions wherever they are set
def getDBC(x,flag):
        return constant_zero

class Initial(object):
    def uOfXT(self,x,t):
        return 128.0*x[0]*(1-x[0])*x[1]*(1-x[1])*x[2]*(1-x[2])

#physics

a0 = 1.0
T = 1.0

physics = default_p
physics.nd = 3; #Three dimensions
physics.L=(1.0,1.0,1.0) #spatial domain: unit cube
A  =[[a0,0.0,0.0],
     [0.0,a0,0.0],
     [0.0,0.0,a0]]
physics.coefficients=Heat(A) #the object for evaluating the coefficients   
physics.dirichletConditions = {0:getDBC}
physics.initialConditions = {0:Initial()}

#numerics

nDTout = 100
DT = T/float(nDTout)

numerics=default_n
numerics.timeIntegration=TimeIntegration.BackwardEuler
numerics.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis} #piecewise linears
numerics.elementQuadrature = Quadrature.SimplexGaussQuadrature(physics.nd,2) #Quadrature rule for elements
numerics.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(physics.nd-1,2) #Quadrature rule for element boundaries
numerics.nnx = numerics.nny = numerics.nnz = 21 #number of nodes in the x and y direction
numerics.matrix = LinearAlgebraTools.SparseMatrix #matrix type

#use petsc solvers wrapped by petsc4py
#numerics.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
#numerics.levelLinearSolver = LinearSolvers.KSP_petsc4py
#using petsc4py requires weak boundary condition enforcement
numerics.numericalFluxType = NumericalFlux.Diffusion_IIPG_exterior

#we can also use our internal wrapper for SuperLU
numerics.multilevelLinearSolver = LinearSolvers.LU
numerics.levelLinearSolver = LinearSolvers.LU
numerics.l_atol_res = 1.0e-8 #linear solver rtolerance
numerics.periodicDirichletConditions=None 

#numerics.nl_atol_res = 1.0e-4

#
# split operator options (trivial since we're not splitting)
#
so = default_so
so.sList=[default_s]
# this is a time interval from 0 to 1
so.tnList = [i*DT for i in range(nDTout+1)]

Profiling.verbose = 1
opts.logLevel = 2
