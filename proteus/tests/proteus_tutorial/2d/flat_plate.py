import numpy as np
from proteus import (Domain, Context,
                     WaveTools as wt)
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import Gauges as ga
from proteus.mprans import BoundaryConditions as bc 


opts=Context.Options([
    # Geometry
    ('Lx', 4., 'Domain length'),
    ('Ly', 0.4, 'Domain height'),
    # Flow
    ('U', [0.2, 0.0, 0.0], 'Set inlet velocity'),
    ('ramp',1.,"ramp time"),
    ('nu', 1e-6, 'Kinematic viscosity'),
    ('rho', 1000., 'Density'),
    ('g',np.array([0,-9.81,0]),"gravitational acceleration"),
    # Turbulence and parameters
    ("useRANS", 1, "Switch ON turbulence models: 0-None, 1-K-Epsilon, 2-K-Omega1998, 3-K-Omega1988"), # ns_closure: 1-classic smagorinsky, 2-dynamic smagorinsky, 3-k-epsilon, 4-k-omega
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("K", 0.41, "von Karman coefficient for the turbulence model"),
    ("B", 5.57, "Wall coefficient for the turbulence model"),
    ("Cmu", 0.09, "Cmu coefficient for the turbulence model"),
    # simulation options
    ('duration', 10., 'Simulation duration'),
    ('dt_init', 0.001, 'Initial timestep'),
    ('dt_output', 0.1, 'time output interval'),
    ("he", 0.03,"Mesh size"),
    ("cfl", 0.5 ,"Target cfl")
    ])
  

#########################################
#domain
#########################################
domain = Domain.PlanarStraightLineGraphDomain()
tank = st.Tank2D(domain, dim=[opts.Lx,opts.Ly])
##################################
#turbulence calculations
##################################
# Reynodls
Re0 = opts.U[0]*opts.Ly/opts.nu
# Skin friction and friction velocity for defining initial shear stress at the wall
cf = 0.045*(Re0**(-1./4.))
Ut = opts.U[0]*np.sqrt(cf/2.)
kappaP = (Ut**2)/np.sqrt(opts.Cmu)
Y_ = opts.he 
Yplus = Y_*Ut/opts.nu
dissipationP = (Ut**3)/(0.41*Y_)

# ke or kw
useRANS = opts.useRANS  # 0 -- None
                        # 1 -- K-Epsilon
                        # 2 -- K-Omega, 1998
                        # 3 -- K-Omega, 1988

model = 'ke'

if opts.useRANS >= 2:
    # k-omega in kw model w = e/k
    model = 'kw'
    dissipationP = np.sqrt(kappaP)/(opts.K*Y_*(opts.Cmu**0.25)) # dissipationP/kappaP

# inlet values 
kInflow = kappaP 
dissipationInflow = dissipationP 

#####################################################
# Boundaries
#####################################################
boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1, 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                           }
boundaryTags = {'y-': 1,
                'x+': 2,
                'y+': 3,
                'x-': 4,
}


# Attached to 'kappa' in auxiliary variables
kWallTop = bc.kWall(Y=Y_, Yplus=Yplus, nu=opts.nu)
kWallBottom = bc.kWall(Y=Y_, Yplus=Yplus, nu=opts.nu)
kWalls = [kWallTop, kWallBottom]
# Attached to 'twp' in auxiliary variables
wallTop = bc.WallFunctions(turbModel=model, kWall=kWallTop, Y=Y_, Yplus=Yplus, U0=opts.U, nu=opts.nu, Cmu=opts.Cmu, K=opts.K, B=opts.B)
wallBottom = bc.WallFunctions(turbModel=model, kWall=kWallBottom, Y=Y_, Yplus=Yplus, U0=opts.U, nu=opts.nu, Cmu=opts.Cmu, K=opts.K, B=opts.B)
walls = [wallTop, wallBottom]


tank.BC['x-'].setConstantInletVelocity(U=opts.U,ramp= opts.ramp,kk= kInflow, dd=dissipationP ,b_or=boundaryOrientations['y+'] )

tank.BC['x+'].setConstantOutletPressure(p = 0, g = opts.g, rho=opts.rho, kk=kInflow, dd= dissipationP ,b_or=boundaryOrientations['y+'])

tank.setTurbulentWall(walls)
tank.setTurbulentKWall(kWalls)
tank.BC['y+'].setWallFunction(walls[0])
tank.BC['y-'].setWallFunction(walls[1])

tank.BC['x-'].setConstantInletVelocity(opts.U,opts.ramp,kInflow,dissipationP,boundaryOrientations['x-'])
tank.BC['x+'].setConstantOutletPressure(0.,opts.rho,opts.g,kInflow,dissipationP,boundaryOrientations['x+'])
class AtRest:
    def uOfXT(self, x, t):
        return 0.0
class kIn:
    def uOfXT(self, x, t):
        return kInflow 

class dIn:
    def uOfXT(self, x, t):
        return dissipationP 

initialConditions = {'pressure':AtRest(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest(),
                     'k':kIn(),
                     'dissipation':dIn()}
                     
outputStepping = TpFlow.OutputStepping(final_time=opts.duration,
                                       dt_init=opts.dt_init,
                                       dt_output=opts.dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

                     

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=opts.he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None # set with SpatialTools,
                                               )
                     
params = myTpFlowProblem.Parameters

params.physical.densityA = opts.rho  # water
params.physical.densityB = opts.rho # air
params.physical.kinematicViscosityA = opts.nu  # water
params.physical.kinematicViscosityB = opts.nu  # air
params.physical.surf_tension_coeff = 0.
params.physical.gravity = opts.g
params.physical.useRANS = opts.useRANS

# Indices
m = params.Models
m.rans2p.index = 0
m.kappa.index = 1
m.dissipation.index = 2

########################
# Assemble domain

##########################

domain.MeshOptions.he = opts.he
st.assembleDomain(domain)
