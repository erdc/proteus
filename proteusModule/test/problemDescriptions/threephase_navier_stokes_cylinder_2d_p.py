from pyadh import *
from pyadh.default_p import *
"""
Incompressible Navier-Stokes flow in a tube (actually between two plates) at equilibrium in 2D.
"""

##\page Tests Test Problems 
# \ref navier_stokes_tube_ss_2d_p.py "Incompressible Navier-Stokes flow between two plates at steady state"
#

##\ingroup test
#\file navier_stokes_tube_ss_2d_p.py
#\brief Incompressible Navier-Stokes flow in a tube (actually between two plates) at equilibrium in 2D.
#
#The incompressible Navier-Stokes equations are implemented in the pyadh::TransportCoefficients::NavierStokes class.
#The domain is the unit square and the initial/boundary conditions are
#\f{eqnarray*}
#p(1,y,t) &=& (1-y)\rho \|g\| \\
#u(x,y,t) &=& v(x,y,t) = 0 \mbox{ on } (0,1) \times \{0,1\} \\
#v(0,y,t) &=& 0.0001 \\
#\frac{\partial u}{\partial x}(x,y,t) &=& \frac{\partial v}{x}(x,y,t) = 0 \mbox{ on } \{0\} \times (0,1) 
#\f}
#
#\image html navier_stokes_tube_ss_p.jpg
#\image html navier_stokes_tube_ss_v.jpg

nd = 2

L = (4.0,1.0,1.0)

initialConditions = None

analyticalSolution = None

coefficients = NavierStokes(g=[0.0,-9.8],nd=nd)
coefficients = ThreephaseNavierStokes_ST_LS_SO(epsFact=1.5,
                 sigma=0.0,
                 rho_0=998.2,nu_0=1.004e-6,
                 rho_1=998.2,nu_1=1.004e-6,
                 rho_s=10000.0*998.2,nu_s=10000.0*1.004e-6,
                 g=[0.0,0.0],
                 nd=2,
                 LS_model=None,
                 KN_model=None,
                 epsFact_density=1.5,
                 stokes=False)

#now define the Dirichlet boundary conditions
inflow = 1.0e-3
T = 5.0*L[0]/inflow

def velRamp(t):
    residence_time = L[0]/inflow
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0


def getDBC_pressure_tube(x,flag):
    if x[0] == L[0]:
        return lambda x,t: -(1.0-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

def getDBC_u_tube(x,flag):
    if (x[1] == 0.0 or
        x[1] == 1.0):
        return lambda x,t: 0.0
    elif (x[0] == 0.0):
        return lambda x,t: velRamp(t)*inflow

def getDBC_v_tube(x,flag):
    if (x[1] == 0.0 or
        x[1] == 1.0):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


def getAFBC_p_tube(x,flag):
    if (x[1] == 1.0 or
          x[1] == 0.0):
        return lambda x,t: 0.0
#    elif (x[0] == 0):
#        return lambda x,t: -inflow
    #elif (x[0] == 1):
    #    return lambda x,t: inflow

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

diffusiveFluxBoundaryConditions = {0:{}}

