from pyadh import *
from pyadh.default_p import *
"""
1D, burgers equation smooth ic
"""

## \page Tests Test Problems 
# \ref burgers_sine_1d_p.py "Burgers equation Riemann problem"
#
#030508 Burgers equation with smooth ic that develops a shock. It needs
# periodic bc's to be classical test problem that is used often so
# of limited use right now

##\ingroup test
#\file burgers_sine_1d_p.py
#@{
#
# \brief Burgers equation Riemann problem.
#
# The equation is
#
#\f[
# u_j + \pd{u^2/2}{x} = 0
#\f]
#
# \todo finish burgers_sine_1d_p.py doc

nd = 1

coefficients = ViscousBurgersEqn(v=[1.0],nu=0.0,nd=1)
#now define the Dirichlet boundary conditions
from math import sin,pi
def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t:  0.25 + 0.5*sin(pi*(2.*x[0]-1.-t))#don't know if this is right for approximating periodic bcs
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

class SineIC:
    def __init__(self):
        pass
        
    def uOfXT(self,x,t):
        return 0.25 + 0.5*sin(pi*(2.*x[0]-1.))
    
analyticalSolution = {0:SineIC()}

initialConditions  = {0:SineIC()}

fluxBoundaryConditions = {0:'outFlow'}#'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.4
