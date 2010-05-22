from pyadh import *
from pyadh.default_p import *
#import pome_cube_pack_gen_poly
import pome3dDomain

import pome_cube_pack_box_poly
"""
Stokes in sphere packing
"""
##\page Tests Test Problems 
# \ref stokes_pome_ss_3d_p.py "Stokes flow in sphere packing "
#

##\ingroup test
#\file stokes_pome_ss_3d_p.py
#
#\brief Stokes flow in a sphere packing
#
#The model equations are defined by the Stokes class. The
#boundary conditions are
#\f{eqnarray*}
#p(1,y,z,t) &=& (1-z) \rho  g \\
#v_x(0,y,z,t) &=& 0.1 \\
#v_x(x,y,z,t) &=& 0.0 \quad y,z \in \{0,1\}\\
#\frac{\partial v_x}{\partial x}(1,y,z,t) &=& 0 \\
#v_y(x,y,z,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#v_z(x,y,z,t) &=& 0.0 \mbox{ on } \partial \Omega\\
#\f}
polyfile = "pome_cube"

nd = 3

initialConditions = None

analyticalSolution = {}

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0,0.0],
                                             nd=3,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)

inflowU = 0.1
inflowP = 10.0
inflowPBCs = False
radius = 0.001
domain = pome3dDomain.pome3d(radius=radius,
                             nx=2,
                             ny=2,
                             nz=2,
                             nSectors = 2,
                             pad=0.0,
                             pointFacets=True)
#domain = pome3dDomain.pome3d()
domain.writePoly(polyfile)
L = domain.L
boundaryTags = domain.boundaryTags
useBoxPack = False#True
#rho=998.2; nu=1.004e-6; #single phase
Re = 0.01 #inflow * L / nu
nu = 1.0#inflowU * nspheres[0]*2.0*radius / Re #just ballpark not right esp. if using Pressure inflow
rho = 1.0


#coefficients = Stokes(rho=rho,nu=nu,g=[0.0,0.0,0.0],nd=nd)
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[0.0,0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)

if useBoxPack == False:
    L,boundaryTags = pome_cube_pack_gen_poly.genPoly(polyfile,
                                                     radius,
                                                     nspheres[0],
                                                     nspheres[1],
                                                     nspheres[2])
else:
    npad = [1,0.1,0.1] #[3,3.0,3.0] should be equivalent to pome_pack_gen_poly
    L,boundaryTags = pome_cube_pack_box_poly.genPoly(polyfile,
                                                     radius,
                                                     nspheres[0],
                                                     nspheres[1],
                                                     nspheres[2],
                                                     npad[0],
                                                     npad[1],
                                                     npad[2])
#now define the Dirichlet boundary conditions
noslipBoundaries = [boundaryTags['front'],boundaryTags['back'],boundaryTags['top'],boundaryTags['bottom'],boundaryTags['sphere']]

def getDBC_pressure_pome(x,tag):
    if tag == boundaryTags['right']:
        return lambda x,t: (L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
    #specified pressure on inflow
    if inflowPBCs and tag == boundaryTags['left']:
        return lambda x,t: inflowP+(L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
def getDBC_u_pome(x,tag):
    if tag in noslipBoundaries:
        return lambda x,t: 0.0
    if not inflowPBCs and tag == boundaryTags['left']:
        return lambda x,t: inflowU
    
def getDBC_v_pome(x,tag):
    if tag in noslipBoundaries:
        return lambda x,t: 0.0
    
def getDBC_w_pome(x,tag):
    if tag in noslipBoundaries:
        return lambda x,t: 0.0


dirichletConditions = {0:getDBC_pressure_pome,
                       1:getDBC_u_pome,
                       2:getDBC_v_pome,
                       3:getDBC_w_pome}


def getAFBC_p_pome(x,tag):
    if tag in noslipBoundaries:
        return lambda x,t: 0.0
    if not inflowPBCs and tag == boundaryTags['left']:
        return lambda x,t: -inflowU

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_pome}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{},3:{}}

