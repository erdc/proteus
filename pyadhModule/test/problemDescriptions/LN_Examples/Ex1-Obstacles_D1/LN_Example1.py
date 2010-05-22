from pyadh import Quadrature,Tracking
import numpy
#parameters for Mario's Example 1

#spatial domain
nd = 2
#rlevel=1 is Mario's level 1 mesh
#if zero is used, then triangleOptions determines mesh 
rlevel=0
import sys
sys.path.append("..")
import LN_Example_Domains

domain = LN_Example_Domains.example1Domain(level=rlevel)
domain.writeAsymptote("LN_Example1_level_%s" % rlevel)
domain.writePoly("LN_Example1_level_%s" % rlevel)

#temporal domain
T = 7.#10.0

#flow equation
#left and right heads
head_left = 1.0
head_right= 0.0
initialConditions_flow = None

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0


hydraulicConductivities = {}
#no sources
sources = {}
lowConductivity = 1.0e-6
def hydraulicConductivity_0(x,t):
    return numpy.array([[1.0,0.0],[0.0,1.0]])
def hydraulicConductivity_1(x,t):
    return numpy.array([[lowConductivity,0.0],[0.0,lowConductivity]])
def nosource(x,t):
    return 0.0

#background region
hydraulicConductivities[1] = hydraulicConductivity_0
hydraulicConductivities[2] = hydraulicConductivity_1
hydraulicConductivities[3] = hydraulicConductivity_1
for i in [1,2,3]:
    sources[i]=nosource
eps_bc = 1.0e-6



##numerics
#domain
triangleOptions = "q30Dena0.0005A"#"q30Dena0.00005A"
#triangleOptions = "q30Dena0.005A"#"q30Dena0.00005A"
if rlevel >= 1:
    #uses input poly file which is Mario's mesh
    triangleOptions = "penA"

nLevels = 1

#quadrature
gw_quad_order = 6
elementQuadrature = Quadrature.SimplexGaussQuadrature(nd,gw_quad_order)#CompositeTrapezoidalTriangle(gw_quad_order)#SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(nd-1,gw_quad_order) 
#elementQuadrature =Quadrature.CompositeTrapezoidalTriangle(gw_quad_order)
#elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(nd-1,2) 

massLumping = False

from pyadh.NumericalFlux import *
numericalFlux_flow = None
#for DG
#numericalFlux_flow = Advection_DiagonalUpwind_Diffusion_IIPG
#something odd with sipg and dg-bdm
#Advection_DiagonalUpwind_Diffusion_SIPG#Diffusion_LDG#
#Advection_DiagonalUpwind_Diffusion_IIPG##


#velocity flags
vpp_flag = 'pwl'#'pwl-bdm'#'dg','pwl','p1-nc'
#ellam and tracking parameters
velocitySpaceFlag = 'rt0'#'bdm1'#'rt0'

#transport numerics
useELLAM = True
useBackwardTrackingForOldMass = False

particleTracking=None
#particle tracking info
particleTracking_params = {}
if useELLAM:
    particleTracking_params[('localVelocityRepresentationFlag',0)]=2#2 -- any RT0 except p1-nc, 1 -- p1-nc
    if velocitySpaceFlag == 'bdm1':
        particleTracking_params[('localVelocityRepresentationFlag',0)]=0

    if velocitySpaceFlag == 'rt0':
        analyticalTracking = False
        if analyticalTracking:
            particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_AnalyticalTracking_2d
            particleTracking_params[('zeroTol',0)]=1.0e-6
        else:
            particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_PT123
            particleTracking_params[('atol_tracking',0)]=1.0e-7 #RK time integration tolerances
            particleTracking_params[('rtol_tracking',0)]=0.0
            particleTracking_params[('sf_tracking',0)]  =0.9    #safety factor for RK integration
            particleTracking_params[('dn_safe_tracking',0)] = 1.0e-7 #tolerance for traking in-element tests
             
    elif velocitySpaceFlag == 'bdm1':
        particleTracking = Tracking.SteadyState_LinearAdvection_BDM1Velocity_PT123
        particleTracking_params[('atol_tracking',0)]=1.0e-7 #RK time integration tolerances
        particleTracking_params[('rtol_tracking',0)]=0.0
        particleTracking_params[('sf_tracking',0)]  =0.9    #safety factor for RK integration
        particleTracking_params[('dn_safe_tracking',0)] = 1.0e-7 #tolerance for traking in-element tests
    else:
        raise NotImplementedError


#time Integration
if useELLAM:
    runCFL = 40.5#20.5#0.2
else:
    runCFL = 0.2
nDTout = 10


#transport equation
#dispersivities
alpha_L = 1.0e-4
alpha_T = 1.0e-5
#molecular diffusion
d       = 1.0e-9
#porosity
omega   = 1.0


#boundary conditions
smoothInflow = True;
smoothDt=10.0;
inflowVal=0.0
yInflowStart = 0.8

from math import exp
class GaussIC:
    def __init__(self,sigma=1.0/16.,xc=[0.1,0.1],b=[0.0,0.0]):
        self.sigma= sigma
        self.xc   = numpy.array(xc)
        self.b    = numpy.array(b)
    def uOfXT(self,x,t):
        xct= self.xc + self.b*t
        d2 = numpy.sum((x[:nd]-xct)**2)
        return exp(-0.5*d2/self.sigma**2)

class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0
#initial conditions 
sigma=0.0316; xc=[0.1,0.9]#0.0316; xc=[0.1,0.05]

initialConditions_trans = {0:GaussIC(sigma=sigma,xc=xc)}
analyticalSolution_trans = {0:GaussIC(sigma=sigma,xc=xc)}
