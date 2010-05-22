from pyadh import *
import numpy
from math import *
details = """
linear transport equation in heterogeneous porous medium
"""
#spatial domain Omega = [0,1] x [0,1]
L=(15,1,1) #m
nd = 2
#time scale
T = 10

#number of components
nctrans = 3
ncflow  = 1
#porosity, dispersivities [m]
omega = 0.4*0.875
alpha_L = 0.2#1.0e-2
alpha_T = alpha_L/10.0
Kox_max=100.0
Kox_C=0.1
Kox_E=0.1
Kox_X=3.5e-3
Yield=0.05
k_d  =0.0

#hydraulic conductivities
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
inflowLengthY = 1.0 #region for inflow concentration [0,0:inflowLengthY]
outflowLengthY= 1.0 #region for outflow concentration [L[0],L[1]-inflowLengthY:L[1]]
inflowConc = [1.0,1.0,0.1]
outflowConc= 0.0 #force boundary layer possibly
backgroundConc = [0.0,1.0,3.5e-4]
inflowHead = 15.0
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
    if x[0]==0.0 and x[1] <= inflowLengthY:
        return lambda x,t: inflowHead #m
    if x[0]==L[0] and x[1] >= L[1]-outflowLengthY:
        return lambda x,t: outflowHead #m
    

def noFlowBCs(x,tag):
    if (x[1] == 0.0 or x[1] == L[1]):
        return lambda x,t: 0.0
    if (x[0] == 0.0 and x[1] > inflowLengthY):
        return lambda x,t: 0.0
    if (x[0] == L[0] and x[1] < L[1]-outflowLengthY):
        return lambda x,t: 0.0

#inflow bc's for concentration
def inflowConcentration0(x,tag):
    if x[0] == 0.0 and x[1] <= inflowLengthY:
        return lambda x,t: inflowConc[0]
def inflowConcentration1(x,tag):
    if x[0] == 0.0 and x[1] <= inflowLengthY:
        return lambda x,t: inflowConc[1]
def inflowConcentration2(x,tag):
    pass

#flux bc's for transport
def zeroInflow_tran(x,tag):
    if (x[0] == 0.0 or
        x[0] == L[0]):
        return lambda x,t: 0.0
def zeroInflow_tran_diff(x,tag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0

#numerical method stuff
nn = 41
gw_quad_order = 2

#output time steps to use
nDTout = 10

