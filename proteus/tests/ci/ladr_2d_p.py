from proteus import *
from proteus.default_p import *
from adr import *

name = "ladr_2d"
nd = 2; #Two dimensions
L=(1.0,1.0,1.0); 
T=1.0

coefficients=LAD(M=1.0,
                 A=[[0.001,0.0],
                    [0.0,0.001]],
                 B=[2.0,1.0])

def getDBC(x,flag):
    if x[0] == 0.0 or x[1] == 0.0:
        return lambda x,t: 1.0
    elif x[0] == 1.0 or x[1] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}
advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0:{}}

class IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[0] <= 0.0 or x[1] <= 0.0:
            return 1.0
        else:
            return 0.0

initialConditions  = {0:IC()}
