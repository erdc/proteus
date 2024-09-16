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
    def __init__(self,A,B,nd=3,linearize=False):
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
        self.nd = nd
        self.linearize = linearize
    def evaluate(self,t,c):
        u =  c[('u',0)]
        c[('m',0)][:] = u
        c[('dm',0,0)][:] = 1.0
        for I in range(self.nd):
            c[('a',0,0)][...,I*self.nd+I] = self.A[I][I]
            if self.linearize:
                c[('f', 0)][..., I] = self.B[I]*u
                c[('df', 0, 0)][..., I] = self.B[I]
            else:
                c[('f', 0)][..., I] = 0.5*self.B[I]*u**2
                c[('df', 0, 0)][..., I] = self.B[I]*u
        c[('r',0)][:] = 0.0*u
        c[('dr',0,0)][:] = 0.0

#use numpy for evaluations
import numpy as np

#number of spatial dimensions to use
nd = 1
#convenience function for holding the boundary condition
def constant_zero(x,t):
    return 0.0

def constant_one(x,t):
    return 1.0

def inflow(x):
    eps = 1.0e-8
    for I in range(nd):
        if x[I] <= eps:
            return True
    return False

def outflow(x):
    eps = 1.0e-8
    for I in range(nd):
        if x[I] >= 1.0-eps:
            return True
    return False

# this function's job is to return another function holding the Dirichlet boundary conditions wherever they are set
def getDBC(x,flag):
    if inflow(x):
        return constant_one
    #elif outflow(x):
    #    return constant_zero#won't be strongly enforced if B points out

def getAFBC(x,flag):
    #return None
    return constant_zero

def getDFBC(x,flag):
    return constant_zero
#    if outflow(x):
#        return constant_zero

class Initial(object):
    def uOfXT(self,x,t):
        offset=np.array([0.1,0.1,0.1])
        if inflow(x-offset):
            return 1.0
        else:
            return 0.0

#physics

a0 = 1.0e-2
T = 1.0

physics = default_p
physics.nd = nd; #
physics.name = "burgers_{0}d".format(physics.nd)
physics.L=(1.0,)*nd #spatial domain: unit cube
A  = np.zeros((nd,nd),'d')
for I in range(nd):
    A[I,I]=a0
B = np.array([1.0]*nd)
physics.coefficients=Burgers(A,B,nd=nd,linearize=False) #the object for evaluating the coefficients   
physics.dirichletConditions = {0:getDBC}
physics.advectiveFluxBoundaryConditions = {0:getAFBC}
physics.diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}
physics.initialConditions = {0:Initial()}
physics.fluxBoundaryConditions = {0:'outFlow'}
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
numerics.nn = numerics.nnx = numerics.nny = numerics.nnz = 101 #number of nodes in the x and y direction

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

numerics.subgridError = SubgridError.AdvectionDiffusionReaction_ASGS(physics.coefficients,
                                                                 physics.nd,lag=False)
numerics.shockCapturing = ShockCapturing.ResGradQuad_SC(physics.coefficients,
                                                        physics.nd,
                                                        shockCapturingFactor=0.75,
                                                        lag=True)
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
use_deim = True

simFlagsList=None
if use_deim:
    simFlagsList=[{}]
    simFlagsList[0]['storeQuantities']=['pod_residuals']


def burgers_plot3d(ns):
    #arrays for using matplotlib's unstructured plotting interface
    x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,0]
    y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,1]
    z = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,2]
    u = ns.modelList[0].levelModelList[-1].u[0].dof

    #we want to build a 3d plot of f(x,y,z0), when z0 = 0.5
    u_range = u[4410:4851]
    x_range = x[4410:4851]
    y_range = y[4410:4851]
    u_range = u_range.reshape(21,21)
    x_range = x_range.reshape(21,21)
    y_range = y_range.reshape(21,21)

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf=ax.plot_surface(x_range, y_range, u_range, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    plt.xlabel('x'); plt.ylabel('y')
    plt.title('approximate solution at t = 1');

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.05f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig("solution.png")
    plt.show()
