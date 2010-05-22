from pyadh import *
from pyadh.default_p import *
"""
Stokes flow in a duct with square cross-section.
"""
##\page Tests Test Problems 
# \ref stokes_duct_ss_3d_p.py "Stokes flow in a square duct at steady state"
#

##\ingroup test
#\file stokes_duct_ss_3d_p.py
#
#\brief Stokes flow in a duct with square cross-section.
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
nd = 3
L=(4.0,0.1,1.0)

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0,-9.8],
                                             nd=3,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)

#now define the Dirichlet boundary conditions
Re=1000.0
inflow = Re*coefficients.nu/L[0]
T = 1000.0*L[0]/inflow
tnList = [0.0,L[0]/inflow/10.0,T]

def velRamp(t):
    residence_time = L[0]/inflow
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_duct(x,flag):
    if x[0] == L[0]:#set pressure at outflow
        return lambda x,t: (x[2]-L[2])*coefficients.rho*coefficients.g[2]
    else:
        pass

def getDBC_u_duct(x,flag):
    if x[2] == 0.0:#no slip on bottom
        return lambda x,t: 0.0
    elif x[0] == 0.0:#inflow velocity
        return lambda x,t: velRamp(t)*inflow
    else:
        pass

def getDBC_v_duct(x,flag):
    if x[2] == 0.0:#no slip on bottom
        return lambda x,t: 0.0
    elif x[0] == 0.0:#inflow parallel
        return lambda x,t: 0.0
    else:
        pass

def getDBC_w_duct(x,flag):
    if x[2] == 0.0:#no slip on bottom
        return lambda x,t: 0.0
    elif x[0] == 0.0:#inflow parallel
        return lambda x,t: 0.0
    else:
        pass

dirichletConditions = {0:getDBC_pressure_duct,
                       1:getDBC_u_duct,
                       2:getDBC_v_duct,
                       3:getDBC_w_duct}

periodic = True
weakPeriodic = False

if periodic:
    def getPDBC(x,flag):
        if x[2] > 0.0 and x[2] < L[2]:
            if x[0] > 0.0 and x[0] < L[0]:
                if x[1] < 1.0e-8 or x[1] > L[1]-1.0e-8:
                    return numpy.array([x[0],0.0,x[2]])

    periodicDirichletConditions = {0:getPDBC,
                                   1:getPDBC,
                                   2:getPDBC,
                                   3:getPDBC}
    
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

def getAFBC_p_duct(x,flag):
    if (x[2] == L[2] or
        x[2] == 0.0):#no flow top and bottom
        return lambda x,t: 0.0
    elif (x[0] == 0.0):#inflow flux
        return lambda x,t: -inflow*velRamp(t)
    elif not weakPeriodic:
        if (x[1] in [0.0,L[1]]):#no flow front and back
            return lambda x,t: 0.0
    else:
        pass

def getAFBC_u_duct(x,flag):
    if x[2] == L[2]:#no flow top (bottom no slip)
        return lambda x,t: 0.0
    elif not weakPeriodic:
        if (x[1] in [0.0,L[1]]):#no flow front and back
            return lambda x,t: 0.0
    else:
        pass
def getAFBC_v_duct(x,flag):
    if x[2] == L[2]:#no flow top (bottom no slip)
        return lambda x,t: 0.0
    elif not weakPeriodic:
        if (x[1] in [0.0,L[1]]):#no flow front and back
            return lambda x,t: 0.0
    else:
        pass
def getAFBC_w_duct(x,flag):
    if x[2] == L[2]:#no flow top (bottom no slip)
        return lambda x,t: 0.0
    elif not weakPeriodic:#no flow front and back
        if (x[1] in [0.0,L[1]]):
            return lambda x,t: 0.0
    else:
        pass

advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                    1:getAFBC_u_duct,
                                    2:getAFBC_v_duct,
                                    3:getAFBC_w_duct}

def getDFBC_duct(x,flag):
    if x[2] == L[2]:#no flow top (bottom no slip)
        return lambda x,t: 0.0
    elif not weakPeriodic:#no flow front and back
        if (x[1] in [0.0,L[1]]):
            return lambda x,t: 0.0
    elif x[0] == L[0]:#outflow
        return lambda x,t: 0.0 
    else:
        pass

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_duct},
                                   2:{2:getDFBC_duct},
                                   3:{3:getDFBC_duct}}
