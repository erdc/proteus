from pyadh import *
from pyadh.default_p import *
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
genMesh=True
polyfile = "pome_cube"
#"pome_cube2x2x2"#"pome_cube2x2x2min"
#"pome_cube3x3x3"
nd = 3

initialConditions = None

analyticalSolution = None

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
#coefficients = Stokes(g=[0.0,0.0,0.0],nd=nd)
#coefficients = StokesP(g=[0.0,0.0,0.0],nd=nd)
inflowU = 0.1
inflowP = 10.0
inflowPBCs = True
radius = 0.001#
diameter= radius*2.
nspheres=[2,2,2]#[1,2,2]
npad=3#0 #1#3 orig

L = (0.0+(nspheres[0]+npad)*diameter,
     0.0+(nspheres[1]+npad)*diameter,
     0.0+(nspheres[2]+npad)*diameter)

#now define the Dirichlet boundary conditions
EPS = 1.0e-6
def onSphere(x,tag):
    if tag == 10:
        return True
    return False
def getDBC_pressure_pome(x,tag):
    if x[0] >= L[0]-EPS:
        return lambda x,t: (L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
    #specified pressure on inflow
    if inflowPBCs and x[0] <= EPS:
        return lambda x,t: inflowP+(L[2]-x[2])*abs(coefficients.g[2])*coefficients.rho
def getDBC_u_pome(x,tag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x,tag)):
        return lambda x,t: 0.0
    if not inflowPBCs and x[0] <= EPS:
        return lambda x,t: inflowU
    
def getDBC_v_pome(x,tag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x,tag)):
        return lambda x,t: 0.0
    
def getDBC_w_pome(x,tag):
    if (x[1] <= EPS or
        x[1] >= L[1]-EPS or
        x[2] <= EPS or
        x[2] >= L[2]-EPS or
        onSphere(x,tag)):
        return lambda x,t: 0.0


dirichletConditions = {0:getDBC_pressure_pome,
                       1:getDBC_u_pome,
                       2:getDBC_v_pome,
                       3:getDBC_w_pome}


def getAFBC_p_pome(x,tag):
    if (x[1] >= L[1]-EPS or
        x[1] <= EPS or
        x[2] >= L[2]-EPS or
        x[2] <= EPS or
        onSphere(x,tag)):
        return lambda x,t: 0.0
    if not inflowPBCs and x[0] <= EPS:
        return lambda x,t: -inflowU

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_pome}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{},3:{}}

