"""
Wave_runup
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
from proteus import WaveTools as wt
from proteus.ctransportCoefficients import smoothedHeaviside
import math
# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",10.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.25,"Desired CFL restriction"),
    ("he",0.1,"he relative to Length of domain in x"),
    ("refinement",3,"level of refinement"),
    ("slope",(1/19.85),"Beta, slope of incline (y/x)"),
    ("slope_length",50.0,"right extent of domain x(m)"),
    ("wl",1.0,"water level"),
    ("waves",True, "wave on/off"),
    ("wave_height", 0.28, "Height of the waves in s"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("wave_type", 'solitaryWave', "type of wave"),
    ("fast", False, "switch for fast cosh calculations in WaveTools"),
    ("g", [0, -9.81, 0], "Gravity vector in m/s^2"),
    ])

toe=(-opts.wl/opts.slope)+35#shifted
structured=False
slope=opts.slope
slope_x=opts.slope_length-toe
slope_y=slope*slope_x
top=(slope_y)*1.35
he=opts.he
nnx=nny=nnz=None
wl=opts.wl
a=opts.wave_height

k=((3*a)/(4*(opts.wl**3)))**0.5
L=(2/k)*np.arccosh(0.05**(-0.5))
x0=(-opts.wl/opts.slope)-(L/2)+35 #shifted

if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.wl
    direction = opts.wave_dir
    wave = wt.SolitaryWave(	waveHeight = height, 
				mwl = wl, 
				depth = depth,
               			g = np.array(opts.g), 
				waveDir = direction,
				trans = np.array([x0, 0., 0.]),
                       		fast = opts.fast
			  )

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
domain = Domain.PlanarStraightLineGraphDomain()
nLevels = 1
 
nLayersOfOverlapForParallel = 0

boundaries=['left','right','bottom','slope','top']
boundaryOrientations = {'bottom': np.array([0., -1.,0.]),
                        'right': np.array([+1., 0.,0.]),
                        'top': np.array([0., +1.,0.]),
                        'left': np.array([-1., 0.,0.]),
                        'slope': np.array([-1., 0.,0.]),}

boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

vertices=[[0.0,0.0],#0                                                           
          [toe,0.0],#1                                                                        
          [opts.slope_length,slope_y],#2                                    
          [opts.slope_length,top],#3  
          [0.0,top]]#4

vertexFlags=[boundaryTags['left'],
             boundaryTags['bottom'],
             boundaryTags['slope'],
             boundaryTags['right'],
             boundaryTags['top']]

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,0]]
 
segmentFlags=[boundaryTags['bottom'],
              boundaryTags['slope'],
              boundaryTags['right'],
              boundaryTags['top'],
              boundaryTags['left']]

regions=[[0.1, 0.01]]

regionFlags=[1]

tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.0

class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        return x[1] - (wave.eta(x,0) + opts.wl)    

epsFact_consrv_heaviside=3.0
wavec =  np.sqrt(9.81 * (depth+opts.wave_height))

def weight(x,t):
    return 1.0-smoothedHeaviside(epsFact_consrv_heaviside*opts.he,
                                 (x[1] - (max(wave.eta(x, t%(toe/wavec)),
                                 wave.eta(x+toe, t%(toe/wavec)))
                                             +opts.wl)))
        
class vel_u(object):
    def uOfXT(self, x, t):
        if x[1] <= wave.eta(x,t) + opts.wl:
            return weight(x,t)*wave.u(x,t)[0]
        else:
            return 0.0

class vel_v(object):
    def uOfXT(self, x, t):
        if x[1] <= wave.eta(x,t) + opts.wl:
            return weight(x,t)*wave.u(x,t)[1]
        else:
            return 0.0

# ****************************** #
# ***** Boundary CONDITIONS***** #
# ****************************** #                                                                  

tank.BC['top'].setAtmosphere()
tank.BC['bottom'].setFreeSlip()
tank.BC['left'].setFreeSlip()
tank.BC['right'].setFreeSlip()
tank.BC['slope'].setFreeSlip()

he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################

outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)

initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': vel_u(),
                     'vel_v': vel_v(),
                     'clsvof': clsvof_init_cond()}

boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
    'pressure_increment_DBC': lambda x, flag: domain.bc[flag].pInc_dirichlet.init_cython(),
    'vel_u_DBC': lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
    'vel_v_DBC': lambda x, flag: domain.bc[flag].v_dirichlet.init_cython(),
    'vel_w_DBC': lambda x, flag: domain.bc[flag].w_dirichlet.init_cython(),
    'clsvof_DBC': lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython(),
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': lambda x, flag: domain.bc[flag].p_advective.init_cython(),
    'pressure_increment_AFBC': lambda x, flag: domain.bc[flag].pInc_advective.init_cython(),
    'vel_u_AFBC': lambda x, flag: domain.bc[flag].u_advective.init_cython(),
    'vel_v_AFBC': lambda x, flag: domain.bc[flag].v_advective.init_cython(),
    'vel_w_AFBC': lambda x, flag: domain.bc[flag].w_advective.init_cython(),
    'clsvof_AFBC': lambda x, flag: domain.bc[flag].vof_advective.init_cython(),
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython(),
    'vel_u_DFBC': lambda x, flag: domain.bc[flag].u_diffusive.init_cython(),
    'vel_v_DFBC': lambda x, flag: domain.bc[flag].v_diffusive.init_cython(),
    'vel_w_DFBC': lambda x, flag: domain.bc[flag].w_diffusive.init_cython(),
    'clsvof_DFBC': lambda x, flag: None}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=1,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=structured,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=False)
physical_parameters = myTpFlowProblem.physical_parameters
myTpFlowProblem.clsvof_parameters['disc_ICs']=False
