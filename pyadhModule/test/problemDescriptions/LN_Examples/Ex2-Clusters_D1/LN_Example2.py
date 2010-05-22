from pyadh import Quadrature,Tracking
import numpy
#parameters for Mario's Example 2

#spatial domain
nd = 2
#rlevel=1 is Mario's level 1 mesh
#if zero is used, then triangleOptions determines mesh 
rlevel=1

import sys
sys.path.append("..")
import LN_Example_Domains
domain = LN_Example_Domains.example2Domain(level=rlevel)
domain.writeAsymptote("LN_Example2_level_%s" % rlevel)
domain.writePoly("LN_Example2_level_%s" % rlevel)

#triangleOptions = "q30Dena0.005A"
if rlevel >= 1:
    #uses input poly file which is Mario's mesh
    triangleOptions = "penA"
if rlevel == 0:
   triangleOptions = "q30Dena0.005A"

nLevels=1


#temporal domain
T = 1.0

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
def hydraulicConductivity_0(x,t):
    return numpy.array([[1.0,0.0],[0.0,1.0]])
def hydraulicConductivity_1(x,t):
    return numpy.array([[1.0e-6,0.0],[0.0,1.0e-6]])
def nosource(x,t):
    return 0.0

#background region
hydraulicConductivities[7] = hydraulicConductivity_0
for i in range(1,7):
    hydraulicConductivities[i] = hydraulicConductivity_1
for i in range(1,8):
    sources[i]=nosource

#quadrature
gw_quad_order = 3
elementQuadrature = Quadrature.SimplexGaussQuadrature(nd,gw_quad_order)#CompositeTrapezoidalTriangle(gw_quad_order)#SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(nd-1,gw_quad_order) 

massLumping = False

numericalFlux_flow = None
#for DG
#Diffusion_LDG#Advection_DiagonalUpwind_Diffusion_IIPG
#something odd with sipg and dg-bdm
#Advection_DiagonalUpwind_Diffusion_SIPG#Diffusion_LDG#
#Advection_DiagonalUpwind_Diffusion_IIPG##



#velocity flags
vpp_flag = 'pwl'#'dg','pwl','p1-nc'
#
