from pyadh import *
from pyadh.default_p import *
from adr import *
name = "ladr_2d_ldg"
nd = 2; L=(1.0,1.0,1.0); T=1.0
#mwf orig
#coefficients=LAD(M=1.0,
#                 A=[[0.001,0.0],
#                    [0.0,0.001]],
#                 B=[1.0,1.0])
coefficients=LAD(M=1.0,
                 A=[[0.01,0.0],
                    [0.0,0.01]],
                 B=[1.0,1.0])
def getDBC(x,flag):
    if x[0] == 0.0 or x[1] == 0.0:
        return lambda x,t: 1.0
    elif x[0] == 1.0 or x[1] == 1.0:
        return lambda x,t: 0.0
dirichletConditions = {0:getDBC}
diffusiveFluxBoundaryConditions = {0:{}}
sd=True
