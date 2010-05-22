from pyadh import *
from pyadh.default_p import *
"""
Linear advection-diffusion-reaction system with five components at equilibrium in 1D.
"""

## \page Tests Test Problems 
# \ref ladr_5c_fc_ss_1d_p.py "Linear advection-diffusion-reaction, five-component, fully coupled system, at steady state"
#

##\ingroup test
#\file ladr_5c_fc_ss_1d_p.py
#
#\brief Linear advection-diffusion-reaction system with five components at equilibrium in 1D.
#\todo finish ladr_5c_fc_ss_1d_p.py doc

nd = 3

initialConditions = None

a=1.0
b=1.0

A={0:numpy.array([[a]]),
   1:numpy.array([[a]]),
   2:numpy.array([[a]]),
   3:numpy.array([[a]]),
   4:numpy.array([[a]])}
B={0:numpy.array([b]),
   1:numpy.array([-1.0*b]),
   2:numpy.array([0.5*b]),
   3:numpy.array([-0.5*b]),
   4:numpy.array([10.0*b])}
A={0:numpy.array([[a,0,0],[0,a,0],[0,0,a]]),
   1:numpy.array([[a,0,0],[0,a,0],[0,0,a]]),
   2:numpy.array([[a,0,0],[0,a,0],[0,0,a]]),
   3:numpy.array([[a,0,0],[0,a,0],[0,0,a]]),
   4:numpy.array([[a,0,0],[0,a,0],[0,0,a]])}
B={0:numpy.array([b,0,0]),
   1:numpy.array([-1.0*b,0,0]),
   2:numpy.array([0.5*b,0,0]),
   3:numpy.array([-0.5*b,0,0]),
   4:numpy.array([10.0*b,0,0])}
C={0:0.0,
   1:0.0,
   2:0.0,
   3:0.0,
   4:0.0}
M={0:0.0,
   1:0.0,
   2:0.0,
   3:0.0,
   4:0.0}        

if a > 0.0 and (b/a < 1.0/2.5e-2):
    analyticalSolution = {0:AnalyticalSolutions.LinearAD_SteadyState(b=B[0][0],a=A[0][0,0]),
                           1:AnalyticalSolutions.LinearAD_SteadyState(b=B[1][0],a=A[1][0,0]),
                           2:AnalyticalSolutions.LinearAD_SteadyState(b=B[2][0],a=A[2][0,0]),
                           3:AnalyticalSolutions.LinearAD_SteadyState(b=B[3][0],a=A[3][0,0]),
                           4:AnalyticalSolutions.LinearAD_SteadyState(b=B[4][0],a=A[4][0,0])}
else:
    analyticalSolution = None

coefficients = LinearVADR_ConstantCoefficients_full(nc=5,M=M,A=A,B=B,C=C)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    if x[0] == 0.0:
        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC,1:getDBC,2:getDBC,3:getDBC,4:getDBC}

fluxBoundaryConditions = {0:'noFlow',1:'noFlow',2:'noFlow',3:'noFlow',4:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{},3:{},4:{}}

