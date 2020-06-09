from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
from proteus.WaveTools import SteadyCurrent
import numpy as np
from proteus.mprans import BodyDynamics as bd
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.mprans.SpatialTools import Tank2D
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)


################################        Context      ################################

opts=Context.Options([
    #physical parameters
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0, -9.805, 0]), "Gravity"),

    # Geometry
    ("tank_dim", (4.75,2.5), "horizontal and vertical tank dimentions"),
    # Sponge
    ("left_sponge", 2.,"length of generation/absorption zone in left side of the tank"),
    ("right_sponge",2.,"length of absorption zone in right side of the tank"),

    #initial conditions
    ("mwl", 2.0, "Height of free surface above bottom"),
    ("waterLine_x", 1000, "x- dimention of water surface used in the signed distance function"),
    ("waterLine_z", 2.0, "water depth used in the signed distance function"),

    # waves
    ("current", True, "apply current True/False"),
    ("U", np.array([0.8,0.,0.]) ,"current velocity in vector form"),
    ("rampTime", 10., "ramp time for current"),
    
    # circle2D
    ("circle2D", True, "Switch on/off circle2D"),
    ("radius", 0.125, "Radius of the circle2D"),
    ("width", 1.0, "Z-dimension of the circle2D"),
    ("mass", 46.41, "Mass of the caisson2D [kg]"),    # density of polypropylene is 946 kg/m3
    ("rotation", False, "Initial position for free oscillation"),
    ("circleBC", 'NoSlip', "circle2D boundaries: NoSlip or FreeSlip"),
    ("InputMotion", True, "If True, set a motion as input rather than calculate it"),
    ("At", [0.0, 0.075, 0.0], "Amplitude of imposed sinusoidal translational motion"),
    ("Tt", [0.0, 1.30, 0.0], "Period of imposed sinusoidal translational motion"), # fD/U 
    ("scheme", None, "Numerical scheme applied to solve motion calculation (Runge_Kutta or Central_Difference)"),
    ("springs", False, "Switch on/off soil module"),
    ("Kx", 0.0, "Horizontal stiffness in N/m"),
    ("Ky", 0.0, "Vertical stiffness in N/m"),
    ("Krot", 0.0, "Rotational stiffness in N*m"),
    ("Cx", 0.0, "Damping factor in N/m*s"),
    ("Cy", 0.0, "Damping factor in N/m*s"),
    ("Crot", 0.0, "Rotational damping factor in N*m*s"),
    ("free_x",np.array([0.0, 0.0, 0.0]), "Translational DOFs"),
    ("free_r",np.array([0.0, 0.0, 0.0]), "Rotational DOFs"),
("nPoints",28, "points used for circle definition"),

    # Turbulence
    ("useRANS", 3, "Switch ON turbulence models"),# 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega, 1998 # 3 -- K-Omega, 1988
    ("c_mu", 0.09, "mu coefficient for the turbulence model"),
    ("turbLength", 'different', "Switch for omega calculation [1 or 2]. If 1 then calculated from turbulent scale formulation. If 2 then calculated from mixing length formulation."),
    ("scaleLength", 0.05, "Only if turbLength isn't either 1 or 2. Turbulent length in terms of waterDepth"),
    ("sigma_k", 1.0, "sigma_k coefficient for the turbulence model"),
    ("sigma_e", 1.0, "sigma_e coefficient for the turbulence model"),
    ("p1_turbLength1",0.038 ," parameter multiplied for turbLength1"),
    ("p2_turbLength2",0.070 ," parameter multiplied for turbLength2"),

    # numerical options
    ("ecH", 3, "smoothing=ecH*he, smoohing coefficient"),
    ("GenZone", True, "Turn on generation zone at left side"),
    ("AbsZone", True, "Turn on absorption zone at right side"),
    ("duration", 60., "Simulation duration"),
    ("cfl", 0.5 ,"Target cfl"),
    ("freezeLevelSet", True, "No motion to the levelset"),
    ("movingDomain", True, "Moving domain and mesh option"),
    ("nd", 2 ,"used in signed distance function"),
    ("Tend",60., "total simulation time"),
    ("fract", 1., "total simulation time/ chosen simulation time"),
    ("Np", 10. ," Output points per period Tp/Np" ),
    ("dt_init", 0.001 , "initial time step"),
    ("refinement_level",2.5,"to define characteristic element size he, he=radius/refinement_level")
    
    ])

#############################        Domain & Waves      ################################
# --- Domain
domain = Domain.PlanarStraightLineGraphDomain()
tank = Tank2D(domain, opts.tank_dim)

# --- WAVE input
current=wt.SteadyCurrent(U=opts.U,
                        mwl=opts.mwl,
                        rampTime=opts.rampTime)

# --- Script on wave length
# wave_length=wave.wavelength

# --- Sponge
tank.setSponge(x_n=opts.left_sponge, x_p=None)

# --- Refinement
he=opts.radius/opts.refinement_level
smoothing=opts.ecH*he

dragAlpha = 5*(2*np.pi/opts.Tt[1])/1e-6

tank.setGenerationZones(x_n=True, 
 			waves=current,
 			dragAlpha=dragAlpha,
			smoothing = smoothing)

#tank.setAbsorptionZones(x_p=True,
#			dragAlpha=dragAlpha,
#	                )

#####################     Shape Geometry & Characteristics    ###############################

# --- Circle
b=2.0*opts.radius
h=2.0*opts.radius

# --- Circle Setup
circle = st.Circle(domain=domain,
                   radius=opts.radius,
                   coords=(0.5, opts.tank_dim[1]/2),
                   barycenter=(0.5, opts.tank_dim[1]/2),
                   nPoints=opts.nPoints)

# --- Body properties setup
circle2D = bd.RigidBody(shape=circle)

circle2D.setMass(mass=opts.mass)
I = (3.14*(opts.radius**4)/2.)*opts.mass

circle2D.It= I/opts.mass/opts.width
circle2D.setConstraints(free_x=opts.free_x, free_r=opts.free_r)
circle2D.setNumericalScheme(scheme=opts.scheme)
circle2D.inputMotion(InputMotion=opts.InputMotion, At=opts.At, Tt=opts.Tt)
circle2D.setRecordValues(filename='circle2D', all_values=True)

# --- Spring setup
  
#circle2D.setSprings(springs=opts.springs,  K=[opts.Kx, opts.Ky, 0.0, Krot], C=[opts.Cx, opts.Cy, 0.0, opts.Crot])


################################       Turbulance      ################################

turbLength1 = opts.p1_turbLength1*(2*opts.radius)  
#Usually in pipeline cases, this is a good estimator. In other case can be use(0.038*hydrDiameter), 
#where hydrDiameter=4*crossSectionArea/wetPerimeter

turbLength2 = opts.p2_turbLength2*(2*opts.radius) 
# Based on mixing length formulation. Then different omega equation.

Re=opts.U[0]*(2*opts.radius)/opts.nu_0
turbI = 0.16*( Re**(-1./8.) )
kInflow = ((turbI*opts.U[0])**2) * 3./2. 

if opts.turbLength == 1:
    dissipationInflow1 = (kInflow**0.5) / (turbLength1) # omega formulation based on turbulent length formulation
    dissipationInflow = dissipationInflow1
elif opts.turbLength == 2:
    dissipationInflow2 = (kInflow**0.5) / ( turbLength2*(opts.c_mu**0.25) ) # omega formulation based on mixing length formulation
    dissipationInflow = dissipationInflow2
else:
    dissipationInflow = (kInflow**0.5) / (opts.scaleLength*opts.mwl)

# pipeline conditions
# skin-friction calculation (see Pope, pages 279, 301)
cf = 0.664*(Re**-0.5)
Ut = opts.U[0]*sqrt(cf/2.)
kappaP = (Ut**2)/sqrt(opts.c_mu)
dissipationP = sqrt(kappaP)/((opts.c_mu**0.25)*0.41*he)

# inlet values 
kInflow = 0.0001*kappaP #None
dissipationInflow =0.0001*dissipationP # None


ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

useRANS=opts.useRANS

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

########################        Boundary Conditions      ################################

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1, 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'circle': None,
                        }
boundaryTags = {'y-': 1,
                'x+': 2,
                'y+': 3,
                'x-': 4,
                'sponge': 5,
                'circle':6,
                }

if opts.circle2D:
    for bc in circle.BC_list:
        if opts.circleBC == 'FreeSlip':
            bc.setFreeSlip()
        if opts.circleBC == 'NoSlip':
            bc.setNoSlip()
            bc.setTurbulentDirichlet(kVal=kappaP, 
                                    dissipationVal=dissipationP)

#tank.BC['x-'].setFixedNodes()

tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=current,
                                               smoothing=opts.ecH*he, 
                                               vert_axis=None)



tank.BC['x+'].setHydrostaticPressureOutletWithDepth(seaLevel=opts.mwl, 
                                                        rhoUp=opts.rho_1,
                                                        rhoDown=opts.rho_0, 
                                                        g=opts.g,
                                                        refLevel=2.5,
                                                        smoothing=3.0*he,
                                                        U=None, 
					                Uwind=None)

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()


for tb in [tank.BC['x+'], tank.BC['x-'], tank.BC['y+'], tank.BC['y-'], tank.BC['sponge']]:
        tb.setFixedNodes()


##########################        Initial Conditions      ################################

def signedDistance(x):
    phi_x = x[0]- opts.waterLine_x
    phi_y = x[1] - opts.mwl
    if phi_x < 0.0:
        if phi_y < 0.0:
            return max(phi_x, phi_y)
        else:
            return phi_y
    else:
        if phi_y < 0.0:
            return phi_x
        else:
            return np.sqrt(phi_x ** 2 + phi_y ** 2)


class P_IC:
    def __init__(self):
        self.mwl=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1] - (opts.mwl - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_dim[1] - opts.mwl)*opts.rho_1*opts.g[1]

class AtRest:
    def uOfXT(self, x, t):
        return 0.0

class kIn:
    def uOfXT(self, x, t):
        return 0.0

class dIn:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest(),
                     'k': kIn(),
                     'dissipation':dIn()}
class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions['vof'] = VOF_IC()
initialConditions['ncls'] = LS_IC()
initialConditions['rdls'] = LS_IC()

"""
################################        Gauges      ################################

gauge_dx=0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)
PG=[]
zProbes=opts.mwl*0.5
for i in probes:
    PG.append((i, zProbes, 0.),)

levelset_output=ga.PointGauges(gauges=((('phi',),PG),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='levelset_gauges.csv')
"""

################################        Numerics      ################################


Duration= opts.Tend/opts.fract
dt_output = opts.Tt[1]/opts.Np

outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                       dt_init=opts.dt_init,
                                       # cfl=cfl,
                                       dt_output=dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=None,
                                             ls_model=None,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None, 
                                             )

params = myTpFlowProblem.Parameters

myTpFlowProblem.useSuperLu=False#True
params.physical.densityA = opts.rho_0  # water
params.physical.densityB = opts.rho_1  # air
params.physical.kinematicViscosityA = opts.nu_0  # water
params.physical.kinematicViscosityB = opts.nu_1  # air
params.physical.surf_tension_coeff = opts.sigma_01

# index in order of
m = params.Models
m.moveMeshElastic.index=0
m.rans2p.index = 1
m.vof.index = 2
m.ncls.index = 3
m.rdls.index = 4
m.mcorr.index = 5
m.kappa.index = 6
m.dissipation.index = 7

# Assemble domain
domain.MeshOptions.he = he
st.assembleDomain(domain)
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']


