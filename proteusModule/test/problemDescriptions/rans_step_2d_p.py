from pyadh import *
from pyadh.default_p import *
import step2d
from rans_step_2d import *

from pyadh import AnalyticalSolutions
"""
Reynolds Averaged Incompressible Navier-Stokes flow over back step.
"""




#now try using inflow, Re = inflow*(downstream_height - upstream_height)/nu ?
#ReInflow = 1.0#250.#100.0
#inflow   = 1.0
#nu = inflow*(downstream_height - upstream_height)/ReInflow
#rho= 998.2
print "trying ReInflow = %s inflow= %s nu=%s rho=%s " % (ReInflow,inflow,nu,rho)


#mu = coefficients.rho_0*coefficients.nu_0
if useSmagorinsky:
    coefficients = ReynoldsAveragedNavierStokes_AlgebraicClosure(epsFact=0.0,
                                                                 sigma=0.0,
                                                                 rho_0=998.2,nu_0=1.004e-6,
                                                                 rho_1=998.2,nu_1=1.004e-6,
                                                                 g=[0.0,0.0],
                                                                 nd=nd,
                                                                 LS_model=None,
                                                                 KN_model=None,
                                                                 epsFact_density=None,
                                                                 stokes=False,
                                                                 turbulenceClosureFlag=0,
                                                                 smagorinskyConstant=smagorinskyConstant)
else:
    coefficients = ReynoldsAveragedNavierStokes_kEpsilon(epsFact=0.0,
                                                         sigma=0.0,
                                                         rho_0=998.2,nu_0=1.004e-6,
                                                         rho_1=998.2,nu_1=1.004e-6,
                                                         g=[0.0,0.0],
                                                         nd=nd,
                                                         LS_model=None,
                                                         KN_model=None,
                                                         epsFact_density=None,
                                                         stokes=False,
                                                         c_mu=c_mu,
                                                         KEmodelID=[1,2])#try splitting k and epsilon


mu = coefficients.rho*coefficients.nu

#now define the Dirichlet boundary conditions
#using inflow, Re = inflow*(downstream_height - upstream_height)/nu ?
#ReInflow = 100.0
#inflow = coefficients.nu*ReInflow/(downstream_height - upstream_height)
#print "trying ReInflow = %s inflow= %s " % (ReInflow,inflow)

grad_p = -inflow/(upstream_height**2/(8.0*mu))
upstream_start_z = downstream_height - upstream_height
tubeEnd = downstream_length

uProfile = AnalyticalSolutions.PlanePoiseuilleFlow_u2(plane_theta=math.pi/2.0,
                                                      plane_phi=math.pi/2.0,
                                                      v_theta=0.0,
                                                      v_phi=math.pi/2.0,
                                                      v_norm=0.0,
                                                      mu=mu,
                                                      grad_p=grad_p,
                                                      L=[1.0,upstream_height,1.0])

#mwf test flat profile
class uFlat:
    def __init__(self,val=inflow):
        self.val= val
    def uOfX(self,x):
        return self.val
class uFlat2:
    def __init__(self,val=inflow,ztop=upstream_height,zbot=0.0,delta_y=0.1):
        self.val= val; self.ztop = ztop; self.zbot=zbot
        self.delta_y = delta_y
    def uOfX(self,x):
        fact = exp(-(x[1]-self.zbot)*(self.ztop-x[1])/self.delta_y)
        #print "x[1]= %s fact=%s" % (x[1],fact)
        return self.val*(1.0-fact)
uProfile = uFlat2(val=inflow)
#print boundaryTags
def velRamp(t):
#    return 1.0
     if t < 25.0/(tubeEnd/inflow):
         return 1.0-exp(-t*tubeEnd/inflow)#1.0-exp(-t)#1.0-exp(-t*tubeEnd/inflow)
     else:
         return 1.0
def getDBC_pressure_tube(x,flag):
    if flag == boundaryTags['downstream']:
        return lambda x,t: -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]
    else:
        pass

bottom = [boundaryTags['upstream_bottom'],boundaryTags['step_bottom'],
          boundaryTags['downstream_bottom']]

def getDBC_u_tube(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t: uProfile.uOfX(x-[0.0,upstream_start_z,0.0])*velRamp(t)
    elif (flag == boundaryTags['top'] or
          (flag in bottom)):
        return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if (flag == boundaryTags['top'] or
        (flag in bottom)):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}



def getAFBC_p_tube(x,flag):
    #top and bottom
    if (flag == boundaryTags['top'] or
        (flag in bottom)):
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -(downstream_height-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self,steadyval=inflow):
        self.steadyval=steadyval
    def uOfXT(self,x,t):
        return self.steadyval

class Split_u:
    def __init__(self,steadyvalTop=inflow,
                 midpoint=downstream_height-upstream_height):
        self.steadyval=steadyvalTop
        self.midpoint=midpoint
    def uOfXT(self,x,t):
        if x[1] < self.midpoint:
            return 0.0
        return self.steadyval

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(inflow),1:Steady_u(0.0),#1:Split_u(inflow),#1:Steady_u(0.0),#0.0
                     2:Steady_v()}

## @}
