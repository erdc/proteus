from pyadh import *
from pyadh.default_p import *
"""
1D,Linear advection of a slug
"""
name = "la_slug_sgi"

## \page Tests Test Problems 
# \ref la_slug_1d_p.py "Linear advection of a slug"
#

##\ingroup test
#\file la_slug_1d_p.py
#@{
#
# \brief Linear advecction of a slug
#
# The linear advection equation is
#
#\f[
# u_j + \nabla (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_gauss_1d_p.py doc

nd = 1
#use Harari heat eqn example or not
harariTestProblem = True#False#True

if harariTestProblem:
    a0 = 1.0
    b0 =  0.0
else:
    a0=0.333333#1.5e-2#6.0e-5#5.0e-5
    b0=1.0
A0_1c={0:numpy.array([[a0]])}
B0_1c={0:numpy.array([b0])}
C0_1c={0:0.0}
M0_1c={0:1.0}


coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c)

#now define the Dirichlet boundary conditions
leftBCval = 1.0
rightBCval= 0.0
def getDBC(x,tag):
    if not harariTestProblem:
        if x[0] == 0.0:
            return lambda x,t:  leftBCval
        if x[0] == 1.0:
            return lambda x,t: rightBCval
    else:
        if x[0] == 1.0:
            return lambda x,t: rightBCval

dirichletConditions = {0:getDBC}
if harariTestProblem:
    slugStart = 0.0
    slugEnd   = -1
    slugVal   = -1
    backVal   = 1.0#0.0
else:
    slugStart = -0.2#0.2
    slugEnd   = -0.4#0.4
    slugVal   = 1.0
    backVal   = 0.0

class SlugIC:
    def __init__(self,a=0.2,b=0.4,u0=backVal,u1=slugVal,v=1):
        self.a    = a
        self.b    = b
        self.u0   = u0
        self.u1   = u1
        self.v    = v
    def uOfXT(self,x,t):
        if (x[0]-self.v*t) < slugStart or (x[0]-self.v*t) > slugEnd:
            return self.u0
        return self.u1
    
analyticalSolution = {0:SlugIC()}

initialConditions  = {0:SlugIC()}

#
if harariTestProblem:
    fluxBoundaryConditions = {0:'noFlow'}
else:
    fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

if harariTestProblem:
    T = 1.0*L[0]*L[0]/a0#0.75
else:
    T = 0.5
