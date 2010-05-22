from pyadh import *
from pyadh.default_p import *
import numpy
import math
"""
A module for 2D diffusive wave model verification problems.

This module contains several implementations of the diffusive wave equation coefficients and several analytical solutions.
Copyright (c) 2009 by Steve Mattis
"""


nd=2
alpha=5.0/3.0
gamma=0.5
epsilon=1.0e-10

#unit box
L=[10.0,10.0]; x=[0.0,0.0]
T=15.0

domain = Domain.RectangularDomain(L=L,x=x)

coefficients=DiffusiveWave_2D(alpha=alpha,gamma=gamma, epsilon=epsilon)

def getNoBC(x,flag):
    pass

dirichletConditions = {0:getNoBC}

class NoInitial:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
class PeakInitial:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if math.sqrt((x[0]-L[0]/2.0)**2 + (x[1]-L[1]/2.0)**2)<2.0:
            return 4.0*math.exp(-abs(math.sqrt((x[0]-L[0]/2.0)**2 + (x[1]-L[1]/2.0)**2)))
        else:
            return 0.0

initialConditions ={0:PeakInitial()}
fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
