from pyadh import *
from pyadh.default_p import *

nd = 2

L=(1.0,1.0,1.0)

initialConditions = None

analyticalSolution = None

coefficients = Stokes(g=[0.0,0.0],nd=nd)

#now define the Dirichlet boundary conditions
inflow = 0.1

EPS = 1.0e-8

pressureGradientBC = False #otherwise set inflow flux

def getDBC_pressure_tube(x):
    if pressureGradientBC:
        if x[0] == L[0]:
            return lambda x,t: 0.0
        elif x[0] == 0.0:
            return lambda x,t: 1.0
        else:
            pass
    else:
        if x[0] == L[0]:
            return lambda x,t: (x[1]-L[1])*coefficients.rho*coefficients.g[1]
        else:
            pass

def getDBC_u_tube(x):
    if (x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0
    elif (x[0] == 0.0 and not pressureGradientBC):
        return lambda x,t: inflow

def getDBC_v_tube(x):
    if (x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0
    elif (x[0] == 0.0 and not pressureGradientBC):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


def getAFBC_p_tube(x):
    if (x[1] == L[1] or
        x[1] == 0.0):
        return lambda x,t: 0.0
    elif (x[0] == 0 and not pressureGradientBC):
        return lambda x,t: -inflow

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

diffusiveFluxBoundaryConditions = {0:{}}

