from pyadh import *
from pyadh.default_p import *
import step3d
from pyadh import AnalyticalSolutions
from pyadh import RANS2P

nd = 3

upstream_height=0.5
width = 0.3*upstream_height
downstream_height=1.0
upstream_length = 1.0
downstream_length = 1.0
length = upstream_length+downstream_length
polyfile = "step3d"
genMesh=True#False
boundaryTags = step3d.genPoly(fileprefix=polyfile,
                              width=width,
                              upstream_height=upstream_height,
                              downstream_height=downstream_height,
                              upstream_length=upstream_length,
                              downstream_length=downstream_length,
                              step_fun=step3d.linear_profile,
                              n_points_on_step=2,
                              step_length=0.0)


coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=998.2,nu_0=1.004e-6,
                                             rho_1=998.2,nu_1=1.004e-6,
                                             g=[0.0,0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
mu = coefficients.rho_0*coefficients.nu_0
Re = 10.0

#coefficients = NavierStokes(g=[0.0,0.0],nd=nd)
#coefficients = Stokes(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions
inflow = Re*coefficients.nu/upstream_height
residence_time = length/inflow
T=10.0*length/inflow
tnList = [0.0,0.1*residence_time,T]
tnList = [i*T for i in range(100)]#[0.0,0.5*T,T]
grad_p = -inflow/(upstream_height**2/(8.0*mu))
upstream_start_z = downstream_height - upstream_height

#triangleOptions="q30Dena%f" % (0.5*(0.05*tubeTop)**2,)
#triangleOptions="q30Den"
uProfile = AnalyticalSolutions.PlanePoiseuilleFlow_u2(plane_theta=0.0,
                                                      plane_phi=0.0,
                                                      v_theta=0.0,
                                                      v_phi=math.pi/2.0,
                                                      v_norm=0.0,
                                                      mu=mu,
                                                      grad_p=grad_p,
                                                      L=[1.0,1.0,upstream_height])
periodic = False

def velRamp(t):
    if t < residence_time:
        return 1.0-exp(-25.0*t/residence_time)
    else:
        return 1.0

def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:#set pressure on outflow to hydrostatic
        return lambda x,t: -(downstream_height-x[2])*coefficients.rho*coefficients.g[2]
    else:
        pass

bottom = [boundaryTags['upstream_bottom'],boundaryTags['step_bottom'],boundaryTags['downstream_bottom']]

def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: uProfile.uOfX(x-[0.0,0.0,upstream_start_z])*velRamp(t)
    elif (flag == boundaryTags['top'] or
          flag in bottom):
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if (flag == boundaryTags['top'] or
        flag in bottom):#no slip
        return lambda x,t: 0.0

def getDBC_w_tube(x,flag):
    if (flag == boundaryTags['top'] or
        flag in bottom):#no slip
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube,
                       3:getDBC_w_tube}

if periodic:
    def getPDBC(x,flag):
        if flag == boundaryTags['front'] or flag == boundaryTags['back']:
            return numpy.array([x[0],0.0,x[2]])

    periodicDirichletConditions = {0:getPDBC,1:getPDBC,2:getPDBC}

def getAFBC_p_tube(x,flag):
    if (flag == boundaryTags['top'] or
        (flag in bottom)):
        return lambda x,t: 0.0
    if not periodic:
        if (flag == boundaryTags['front'] or
            flag == boundaryTags['back']):
            return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

def getAFBC_p_duct(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: -inflow*velRamp(t)
    elif not periodic:
        if flag == boundaryTags['front'] or flag == boundaryTags['back']:
            return lambda x,t: 0.0

def getAFBC_u_duct(x,flag):
    if not periodic:
        if flag == boundaryTags['front'] or flag == boundaryTags['back']:
            return lambda x,t: 0.0

def getAFBC_v_duct(x,flag):
    if not periodic:
        if flag == boundaryTags['front'] or flag == boundaryTags['back']:
            return lambda x,t: 0.0

def getAFBC_w_duct(x,flag):
    if not periodic:
        if flag == boundaryTags['front'] or flag == boundaryTags['back']:
            return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                    1:getAFBC_u_duct,
                                    2:getAFBC_v_duct,
                                    3:getAFBC_w_duct}
def getDFBC_duct_u(x,flag):
    if not (flag == boundaryTags['upstream'] or
            flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0
def getDFBC_duct_v(x,flag):
    if not (flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0
def getDFBC_duct_w(x,flag):
    if not (flag == boundaryTags['top'] or
            flag in bottom):
        return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_duct_u},
                                   2:{2:getDFBC_duct_v},
                                   3:{3:getDFBC_duct_w}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_w:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}

## @}
