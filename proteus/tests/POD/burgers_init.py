#!/usr/bin/env python
"""
Fine-scale viscous Burgers equation solver.

The equation is

.. math: \frac{\partial u}{\partial t} + \nabla \cdot (\frac{1}{2} u^2 - \epsilon \nabla u) = 0
"""
from proteus.iproteus import *

class Burgers(TransportCoefficients.TC_base):
    """
    The coefficients of the viscout Burgers equation
    """
    def __init__(self,A,B):
        TransportCoefficients.TC_base.__init__(self,
                                               nc=1,
                                               variableNames=['u'],
                                               mass = {0:{0:'linear'}},
                                               diffusion = {0:{0:{0:'constant'}}},
                                               advection = {0:{0:'nonlinear'}},
                                               reaction = {0:{0:'nonlinear'}},
                                               potential = {0:{0:'u'}})
        self.A=A
        self.B=B
	       	
    def evaluate(self,t,c):
        u =  c[('u',0)]
        c[('m',0)][:] = u
        c[('dm',0,0)][:] = 1.0
        c[('a',0,0)][...,0] = self.A[0][0]
        c[('a',0,0)][...,4] = self.A[1][1]
        c[('a',0,0)][...,8] = self.A[2][2]
        c[('f', 0)][..., 0] = 0.5*self.B[0]*u**2
        c[('f', 0)][..., 1] = 0.5*self.B[1]*u**2
        c[('f', 0)][..., 2] = 0.5*self.B[2]*u**2
        c[('df', 0, 0)][..., 0] = self.B[0]*u
        c[('df', 0, 0)][..., 1] = self.B[1]*u
        c[('df', 0, 0)][..., 2] = self.B[2]*u
        c[('r',0)][:] = 0.0*u
        c[('dr',0,0)][:] = 0.0

#use numpy for evaluations
import numpy as np

#convenience function for holding the boundary condition
def constant_zero(x,t):
    return 0.0

def constant_one(x,t):
    return 1.0

def inflow(x):
    eps = 1.0e-8
    if x[0] <= eps or x[1] <= eps or x[2] <= eps:
        return True
    else:
        return False

def outflow(x):
    eps = 1.0e-8
    if x[0] >= 1.0-eps or x[1] >= 1.0-eps or x[2] >= 1.0-eps:
        return True
    else:
        return False

# this function's job is to return another function holding the Dirichlet boundary conditions wherever they are set
def getDBC(x,flag):
    if inflow(x):
        return constant_one
    elif outflow(x):
        return constant_zero#won't be strongly enforced if B points out

def getAFBC(x,flag):
    #return None
    return constant_zero

def getDFBC(x,flag):
    return constant_zero
#    if outflow(x):
#        return constant_zero

class Initial:
    def uOfXT(self,x,t):
        offset=np.array([0.1,0.1,0.1])
        if inflow(x-offset):
            return 1.0
        else:
            return 0.0

#physics

a0 = 1.0
T = 1.0

physics = default_p
physics.nd = 3; #Three dimensions
physics.L=(1.0,1.0,1.0) #spatial domain: unit cube
A  =[[a0,0.0,0.0],
     [0.0,a0,0.0],
     [0.0,0.0,a0]]
B = [1.0,1.0,1.0]
physics.coefficients=Burgers(A,B) #the object for evaluating the coefficients   
physics.dirichletConditions = {0:getDBC}
physics.advectiveFluxBoundaryConditions = {0:getAFBC}
physics.diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}
physics.initialConditions = {0:Initial()}

#numerics

nDTout = 100
DT = T/float(nDTout)

numerics=default_n
numerics.timeIntegration = TimeIntegration.BackwardEuler_cfl
numerics.stepController = StepControl.Min_dt_cfl_controller
numerics.runCFL=0.99
numerics.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis} #piecewise linears
numerics.elementQuadrature = Quadrature.SimplexGaussQuadrature(physics.nd,2) #Quadrature rule for elements
numerics.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(physics.nd-1,2) #Quadrature rule for element boundaries
numerics.nnx = numerics.nny = numerics.nnz = 21 #number of nodes in the x and y direction

numerics.multilevelNonlinearSolver = NonlinearSolvers.Newton
numerics.levelLinearSolver = NonlinearSolvers.Newton
numerics.tolFac = 0.0#relative tolerance
numerics.nl_atol_res = 1.0e-4 #nonlinear solver rtolerance

numerics.matrix = LinearAlgebraTools.SparseMatrix #matrix type
numerics.multilevelLinearSolver = LinearSolvers.LU
numerics.levelLinearSolver = LinearSolvers.LU
numerics.linTolFac = 0.0 #relative  tolerance 
numerics.l_atol_res = 1.0e-8 #linear solver rtolerance

numerics.periodicDirichletConditions=None 

#numerics.subgridError = SubgridError.AdvectionDiffusionReaction_ASGS(physics.coefficients,
#                                                                 physics.nd,lag=False)
#numerics.shockCapturing = ShockCapturing.ResGradQuad_SC(physics.coefficients,
#                                                            physics.nd,
#                                                            shockCapturingFactor=0.75,
#                                                             lag=True)
#numerics.numericalFluxType = NumericalFlux.Advection_DiagonalUpwind_Diffusion_SIPG_exterior
numerics.numericalFluxType = NumericalFlux.NoFlux

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

#whether or not to use deim approximation
use_deim = False

simFlagsList=None
if use_deim:
    simFlagsList=[{}]
    simFlagsList[0]['storeQuantities']=['pod_residuals']
