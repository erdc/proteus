from pyadh import *
import numpy
from math import *
details = """
Ion exchange example from Bryant Dawson etal 00
"""
#spatial domain Omega = [0,1] x [0,1]
L=(100.,1,1) #m
nd = 1
#time scale
T = 80.0

#number of components
nctrans = 2
ncflow  = 1
#porosity, dispersivities [m]
omega = 1.0
alpha_L = 1.0
alpha_T = alpha_L/10.0
K_m=1.e3
K_h=1.e4
K_w=1.e-14
Z_tot=1.e-3

C_h0 = 1.e-3
C_m0 = 1.0e-13
C_hI = 1.e-12
C_mI = 1.e-3

#hydraulic conductivities
#constantFlux
qIn= 1.0
Kx = 1.0
Ky = 1.0
#Kx = 5.0e-3# / 2.0**8 #m/s
#Ky = Kx#1.0e-6 #m/s
#sigma_x = 2.0/L[0]#3.0/L[0]
#sigma_y = 2.0/L[1]#6.0/L[1]
def hydraulicConductivity(x):
    return numpy.array([[Kx,0.0],
                        [0.0,Ky]],'d')

def hydraulicConductivity_het(x):
    return numpy.array([[Kx*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3,0.0],
                        [0.0,Ky*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3]],'d')
def fluidSourceTerms(x):
    return 0.0

#auxiliary conditions
inflowLength = 1.0
outflowLength= 1.0
outflowHead= 0.0 

transportAnalyticalSolution = None
velocityAnalyticalSolution = None
flowAnalyticalSolution = None

class constantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
    def uOfX(self,x):
        return self.cval

#inflow bc's for flow
def dummyBCs(x,tag):
    pass
def flowHeadBCs(x,tag):
    if x[0]==L[0]:
        return lambda x,t: outflowHead #m
    
def inflowFlowBCs(x,tag):
    if x[0] == 0.0:
        return lambda x,t: -qIn

#inflow bc's for concentration
def inflowConcentrationM(x,tag):
    if x[0] == 0.0:
        return lambda x,t: C_m0
def inflowConcentrationH(x,tag):
    if x[0] == 0.0:
        return lambda x,t: C_h0

def zeroInflow_tran_diff(x,tag):
    if x[0] == L[0]:
        return lambda x,t: 0.0

#numerical method stuff
nn = 101
gw_quad_order = 2

#output time steps to use
nDTout = 10

