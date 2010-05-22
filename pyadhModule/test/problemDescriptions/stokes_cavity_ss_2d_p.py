from pyadh import *
from pyadh.default_p import *
"""
A driven cavity problem for Stokes flow at equilibrium.
"""
##\page Tests Test Problems 
# \ref stokes_cavity_ss_2d_p.py "Stokes flow in a driven cavity at steady state"
#

##\ingroup test
#\file stokes_cavity_ss_2d_p.py
#\brief A driven cavity problem for Stokes flow at equilibrium.
#
#The boundary conditions are identical to those in stokes_cavity_2d_p.py
#
nd = 2

L=(1.0,1.0,1.0)

initialConditions = None

analyticalSolution = None

coefficients = Stokes(g=[0.0,9.8],nd=nd)#stokes_2D_tube_p.Stokes2D()
#coefficients = StokesP(g=[0.0,9.8],nd=nd)#stokes_2D_tube_p.Stokes2D()
#coefficients = Coefficients.Stokes(rho=1.0,mu=1.0,g=[0.0,0.0],nd=nd)#stokes_2D_tube_p.Stokes2D()
#coefficients = stokes_2D_tube_p.Stokes2D()
#coefficients.useC = False

def getDBC_pressure(x,flag):
    if x[1] == 0.0:
        if x[0] == 0.5:
            return lambda x,t: 0.0

def getDBC_u(x,flag):
    #set x velocity to 1 at ends and to 0 at top and bottom
    if x[0] == 0.0:
        return lambda x,t: 0.0
    elif x[0] == 1.0:
        return lambda x,t: 0.0
    elif x[1] == 0.0:
        if x[0] < 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1
    elif x[1] == 1.0:
        if x[0] < 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1

def getDBC_v(x,flag):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    elif x[1] == 1.0:
        return lambda x,t: 0.0
    elif x[0] == 0.0:
        if x[1] > 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1
    elif x[0] == 1.0:
        if x[1] > 0.5:
            return lambda x,t: 0.1
        else:
            return lambda x,t: -0.1

dirichletConditions = {0:getDBC_pressure,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getNoBC(x,flag):
    pass

advectiveFluxBoundaryConditions =  {0:getNoBC,1:getNoBC,2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getNoBC},2:{2:getNoBC}}


