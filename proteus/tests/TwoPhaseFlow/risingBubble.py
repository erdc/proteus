"""
Rising bubble test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context)
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("test_case",1,"Rising bubble test cases"),
    ('ns_model',1,"ns_model = {rans2p,rans3p}"),
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",3,"level of refinement")
    ],mutable=True)

# Change parameters for automated testing #
opts.refinement=2

assert opts.ns_model==1, "Surface tension is only implemented with rans3p. use ns_model=1"
assert opts.test_case == 1 or opts.test_case==2, "test_case must be 1 or 2"
# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
tank_dim = (1.0,2.0)
refinement = opts.refinement
structured=True
if structured:
    nnx = 4 * refinement**2 +2
    nny = 2*nnx
    domain = Domain.RectangularDomain(tank_dim)
    boundaryTags = domain.boundaryTags
    triangleFlag=1
else:
    nnx = nny = None
    domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #
tank = Tank2D(domain, tank_dim)

# ----- EXTRA BOUNDARY CONDITIONS ----- #
tank.BC['y+'].setNoSlip()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()

he = old_div(tank_dim[0], float(4 * refinement - 1))
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.
class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        xB = 0.5
        yB = 0.5
        rB = 0.25
        zB = 0.0
        # dist to center of bubble
        r = np.sqrt((x[0]-xB)**2 + (x[1]-yB)**2)
        # dist to surface of bubble
        dB = rB - r
        return dB

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vel_w': zero(),
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
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=structured,
                                             he=he,
                                             nnx=nnx,
                                             nny=nny,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=False)
<<<<<<< HEAD
physical_parameters = myTpFlowProblem.physical_parameters
=======
physical_parameters = myTpFlowProblem.Parameters.physical
>>>>>>> TwoPhaseFlow
physical_parameters['gravity'] = [0.0, -0.98, 0.0]
if opts.test_case==1:
    physical_parameters['densityA'] = 1000.0
    physical_parameters['viscosityA'] = 10.0/physical_parameters['densityA']
    physical_parameters['densityB'] = 100.0
    physical_parameters['viscosityB'] = 1.0/physical_parameters['densityB']
    physical_parameters['surf_tension_coeff'] = 24.5
    physical_parameters['gravity'] = [0.0, -0.98, 0.0]
else: #test_case=2
    physical_parameters['densityA'] = 1000.0
    physical_parameters['viscosityA'] = 10.0/physical_parameters['densityA']
    physical_parameters['densityB'] = 1.0
    physical_parameters['viscosityB'] = 0.1/physical_parameters['densityB']
<<<<<<< HEAD
    physical_parameters['surf_tension_coeff'] = 1.96    
=======
    physical_parameters['surf_tension_coeff'] = 1.96

myTpFlowProblem.useBoundaryConditionsModule = False
myTpFlowProblem.Parameters.Models.rans3p.epsFact_viscosity = 3.
myTpFlowProblem.Parameters.Models.rans3p.epsFact_density = 3.
myTpFlowProblem.Parameters.Models.rans3p.ns_shockCapturingFactor = 0.5
myTpFlowProblem.Parameters.Models.rans3p.timeDiscretization = 'vbdf'

myTpFlowProblem.outputStepping.systemStepExact = True
>>>>>>> TwoPhaseFlow
