from pyadh import *
from pyadh.default_p import *
"""
1D,Linear advection of a Gaussians
"""

## \page Tests Test Problems 
# \ref la_gauss_1d_p.py "Linear advection of a Gaussian"
#

##\ingroup test
#\file la_gauss_1d_p.py
#@{
#
# \brief Linear advecction of a guassian.
#
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_gauss_1d_p.py doc
name = "la_gauss_bdf2"

nd = 1

a0=1.0e-3
b0=1.0
A0_1c={0:numpy.array([[a0]])}
B0_1c={0:numpy.array([b0])}
C0_1c={0:0.0}
M0_1c={0:1.0}


    
coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c,
                                               useSparseDiffusion = sd)

#now define the Dirichlet boundary conditions

def getDBC(x,tag):
    if x[0] == 0.0:
        return lambda x,t:  0.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

class GaussIC:
    def __init__(self,sigma=1.0/16.,xc=0.2,b=1.0):
        self.sigma= sigma
        self.xc   = xc
        self.b    = b
    def uOfXT(self,x,t):
        xct= self.xc + self.b*t
        d2 = (x[0]-xct)**2
        return exp(-0.5*d2/self.sigma**2)
    
analyticalSolution = {0:GaussIC()}

initialConditions  = {0:GaussIC()}

fluxBoundaryConditions = {0:'outFlow'}#'noFlow'}

def getAFBC(x):
    pass
advectiveFluxBoundaryConditions =  {0:getAFBC}


def getDFBC(x):
   if x[0] == 1.0:
       return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}

T = 0.75
