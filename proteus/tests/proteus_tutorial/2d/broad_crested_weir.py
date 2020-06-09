"""
broad crested weir 2-D
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context)
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.9,"Desired CFL restriction"),
    ("he",0.075,"Maximum element edge length"),
    ("inflow_vel",0.139,"inflow velocity for left boundary"),
    ])

waterLine_y = 0.54
waterLine_x = 1.45
outflow_level=0.04
# Water
rho_0 = 998.2
nu_0 = 1.004e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

g = [0., -9.81, 0.]
# ****************** #
# ***** GAUGES ***** #
# ****************** #


# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #
boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'airvent': np.array([-1., 0.,0.]),
                           }
boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'airvent':5,
               }

top = 1.0
vertices=[[0.0, 0.0],#0
          [1.0, 0.0],#1
          [1.0, 0.401], #2
          [1.5, 0.401],#3
          [1.5, 0.35],#4 airvent
          [1.5, 0.25],#5 airvent
          [1.5, 0.0 ],#6
          [2.5,  0.0 ],#7
          [2.5,  top ],#8
          [0.0,  top ],#9
]

vertexFlags=np.array([1, 1, 1, 1,
                      1, 1,
                      1, 1,
                      3, 3,])

segments=[[0,1],#0
          [1,2],#1
          [2,3],#2
          [3,4],#3
          [4,5],#4 airvent
          [5,6],#5
          [6,7],#6
          [7,8],#7
          [8,9],#8
          [9,0],#9
             ]

segmentFlags=np.array([1, 1, 1, 1,
                       5,
                       1, 1,
                       2,
                       3,
                       4,
                          ])

regions = [ [ 0.1 , 0.1] ]

regionFlags=np.array([1])


tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


# ----- EXTRA BOUNDARY CONDITIONS ----- #
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x-'].setTwoPhaseVelocityInlet(U=[opts.inflow_vel,0.,0.],
                                       waterLevel=waterLine_y,
                                       smoothing=1.5*opts.he,
)
tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=outflow_level,
                                                    rhoUp=rho_1,
                                                    rhoDown=rho_0,
                                                    g=g,
                                                    refLevel=top,
                                                    smoothing=1.5*opts.he,
)

tank.BC['airvent'].reset()
tank.BC['airvent'].p_dirichlet.uOfXT = lambda x, t: (1.0 - x[1])*rho_1*abs(g[1])
tank.BC['airvent'].v_dirichlet.uOfXT = lambda x, t: 0.0
tank.BC['airvent'].phi_dirichlet.uOfXT = lambda x, t: 0.0
tank.BC['airvent'].vof_dirichlet.uOfXT = lambda x, t: 1.0
tank.BC['airvent'].u_diffusive.uOfXT = lambda x, t: 0.0
tank.BC['airvent'].v_diffusive.uOfXT = lambda x, t: 0.0

domain.MeshOptions.he = opts.he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((opts.he ** 2)/2.0,)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.0

class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        if x[0] < waterLine_x and x[1] < waterLine_y:
            return -1.0
        elif x[0] > waterLine_x or x[1] > waterLine_y:
            return 1.0
        else:
            return 0.0



class PHI_IC:
    def uOfXT(self, x, t):
        phi_x = x[0] - waterLine_x
        phi_y = x[1] - waterLine_y
        phi_y_outflow = x[1] - outflow_level
        if phi_x <= 0.0:
            if phi_y < 0.0:
                return max(phi_x, phi_y)
            else:
                return phi_y
        else:
            if phi_y_outflow < 0.0:
                return phi_y_outflow
            else:
                if phi_y<0.0:
                    return min(phi_x, phi_y_outflow)
                else:
                    return min((phi_x ** 2 + phi_y ** 2)**0.5, phi_y_outflow)

class VF_IC:
    def __init__(self):
        self.phi=PHI_IC()
    def uOfXT(self, x, t):
        from proteus.ctransportCoefficients import smoothedHeaviside
        return smoothedHeaviside(1.5*opts.he,self.phi.uOfXT(x,t))

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC(),
                     'clsvof': clsvof_init_cond()}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=opts.he,
                                             domain=domain,
                                             initialConditions=initialConditions)
# copts=myTpFlowProblem.Parameters.Models.rans2p.p.CoefficientsOptions
# def getPhiDBC(x, flag):
#     if flag == boundaryTags['x-']:
#         return lambda x,t: x[1] - waterLine_y
#     elif flag == boundaryTags['x+']:
#         return lambda x,t: x[1] - outflow_level
# myTpFlowProblem.Parameters.Models.ncls.p.dirichletConditions = {0: getPhiDBC}

