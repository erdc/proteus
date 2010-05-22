from pyadh import *
import numpy
from math import *
details = """
linear transport equation in heterogeneous porous medium
"""

#spatial domain Omega = [0,100] x [0,100]
L=(100,100,1) #m
nd = 2
#time scale
T = 1.0e4#1.0e7 #s

#number of components
ncflow = 1
nctrans= 1
#porosity, dispersivities [m]
omega = 1.0#0.3
alpha_L = 0.0#1.0
alpha_T = alpha_L/5.0
d_mol    = 0.0#1.3e-6#0.0#1.3e-9
#hydraulic conductivities
Kx = 5.0#try m/d#5.0e-3# / 2.0**8 #m/s
Ky = Kx#1.0e-6 #m/s
sigma_x = 2.0/L[0]#3.0/L[0]
sigma_y = 2.0/L[1]#6.0/L[1]

def hydraulicConductivity(x):
    return numpy.array([[Kx*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3,0.0],
                        [0.0,Ky*(cos(2.0*pi*sigma_x*x[0]) + 1.1)**3 * (cos(2.0*pi*sigma_y*x[1]) + 1.1)**3]],'d')
def hydraulicConductivity0(x):
    return numpy.array([[Kx,0.0],
                        [0.0,Ky]],'d')
def fluidSourceTerms(x):
    return 0.0

#auxiliary conditions
inflowLengthY = 10.#10.0 #region for inflow concentration [0,0:inflowLengthY]
outflowLengthY= 10.#10.0 #region for outflow concentration [L[0],L[1]-inflowLengthY:L[1]]
inflowConc = 1.0
outflowConc= 0.0 #force boundary layer possibly
backgroundConc = 0.0
inflowHead = 1.0e1
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
def inflowConcentration(x,tag):
    if x[0] == 0.0 and x[1] <= inflowLengthY:
        return lambda x,t: inflowConc
    #mwf hack
    #if x[0] == L[0] and x[1] == 0.0:
    #    return lambda x,t: 0.0
#flux bc's for transport
def zeroInflow_tran(x,tag):
    if (x[1] == 0.0 or x[1] == L[1]):
        return lambda x,t: 0.0
    if (x[0] == 0.0 and x[1] > inflowLengthY):
        return lambda x,t: 0.0
    if (x[0] == L[0] and x[1] < L[1]-outflowLengthY):
        return lambda x,t: 0.0
def zeroInflow_tran_diff(x,tag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0


#numerical method stuf
nn = 41
gw_quad_order = 3
#output time steps to use
nDTout = 50

