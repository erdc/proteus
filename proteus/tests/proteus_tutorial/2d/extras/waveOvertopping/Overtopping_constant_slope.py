from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
import numpy as np
from proteus.mprans import BodyDynamics as bd
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import Gauges as ga



# Options (Geometry, Physical Properties, Waves, Numerical Options, Obstacle Dimentions)

opts=Context.Options([
    
    
    # Geometry
    ("tank_height",0.8,"Vertical Dimention of the tank"),
    ("Lback",4.0,"Horizontal Dimention of overtopping collection tank"),
    ("tank_depth", 0.5, "depth of the tank below zero level"),
    ("obs_depth", 0.4,"depth of the structure/obstacle below zero level"),
    ("tube", 0.1,"tube dimention"),
    ("deposit_width",5.0, "width of the tank used to collect the overtopped water"),

    # Physical Properties
    ("rho_0", 998.2, "Water density"),
    ("nu_0", 1.004e-6,"Water viscosity"),
    ("rho_1", 1.205, "Air density"),
    ("nu_1",1.5e-5, "Air viscosity"),
    ("sigma_01", 0.,"surface tension"),
    ("g", np.array([0., -9.805, 0.]), "gravity"),

    # Waves
    ("Tstart", 0, "Start time"),
    ("fract", 1, "fraction of duration"),
    ("Ntotalwaves",500,"totalnumber of waves"),
    ("x0", np.array([0.,0.,0.]), "Position vector for the tinme series"),
    ("Tp", 3.5, "Peak wave period"),
    ("Hs", 0.096, "Significant wave height"),
    ("mwl", 0.4, "Mean Water level"),
    ("depth", 0.4 , "Water depth"),
    ("waveDir", np.array([1.,0.,0.]),"Direction of the waves"),
    ("N", 2000, "Number of frequency components"),
    ("bandFactor", 2.0 ,"Spectal Band Factor"),
    ("spectName", "JONSWAP","Name of Spectral Distribution"),
    ("spectral_params",{"gamma": 3.3, "TMA":False,"depth": 0.4} ,"Spectral Distribution Characteristics"),
    ("seed", 420,"Seed for random phases"),
    ("Lgen",None , "Length of the generation zone"),
    ("Nwaves", 15, "Number of waves per window"),
    ("Nfreq",32 , "Number of fourier components per window"),

    # gauges
    #("gauge_output", True, "Places Gauges in tank (10 per wavelength)"),
    ("point_gauge_output", True, "Produce point gauge output"),
    ("column_gauge_output", True, "Produce column gauge output"),
    ("gauge_dx", 0.3, "Horizontal spacing of point gauges/column gauges before structure [m]"),
    ("gauge_dx_2", 0.1,"Horizontal spacing of point/column gauges after structure [m]"),

    # Numerical Options
    ("refinement_level", 100.,"he=wavelength/refinement_level"),
    ("cfl", 0.5,"Target cfl"),
    ("ecH", 1.5,"Smoothing Coefficient"),
    ("Np", 10 ," Output points per period Tp/Np" ),
    ("dt_init", 0.001 , "initial time step" ),
    
    # Obstacle Dimensions 
    ("structure_slope", 4, "1/slope"),
    ("structureCrestLevel", 0.5, "elevation of structure crest. Equal to Water depth + Rc (crest freeboard)")
    ])

# ---- DOMAIN ----#

domain = Domain.PlanarStraightLineGraphDomain()

# Wave Input

np.random.seed(opts.seed)
phi = 2*np.pi*np.random.rand(opts.N)
Tend=opts.Ntotalwaves*opts.Tp/1.1
wave = wt.RandomWavesFast(Tstart=opts.Tstart,
                         Tend=Tend,
                         x0=opts.x0,
                         Tp=opts.Tp,
                         Hs=opts.Hs,
                         mwl=opts.mwl,
                         depth=opts.depth,
                         waveDir=opts.waveDir,
                         g=opts.g,
                         N=opts.N,
                         bandFactor=opts.bandFactor,
                         spectName=opts.spectName,
                         spectral_params=opts.spectral_params,
                         phi=phi,
                         Lgen=opts.Lgen,
                         Nwaves=opts.Nwaves,
                         Nfreq=opts.Nfreq,
                         checkAcc=True,
                         fast=True)


# Script on wave length
wave_length=wave.wavelength           
print(wave_length)
# Domain

tank_dim = (3*wave_length+opts.structureCrestLevel*opts.structure_slope+opts.deposit_width,opts.structureCrestLevel+0.2)

#shapes

tank_sponge = (wave_length, 0)

boundaryOrientations = {'y-': np.array([0., -1.,0.]),
                        'x+': np.array([+1., 0.,0.]),
                        'y+': np.array([0., +1.,0.]),
                        'x-': np.array([-1., 0.,0.]),
                        'sponge': None,
                        'porousLayer': None,
                       }

boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'sponge' : 5,
               }

################ TANK OUTLINE GEOMETRY ###########################                        

vertices=[[0.0,0.0], #0
            [wave_length,0],#1
            [wave_length,-opts.tank_depth],#2
            [3*wave_length+2*opts.tube+opts.structureCrestLevel*opts.structure_slope+1,-opts.tank_depth], #3
            [3*wave_length+2*opts.tube+opts.structureCrestLevel*opts.structure_slope+1,0.0], #4 after that abs zone
            [3*wave_length+2*opts.tube+opts.structureCrestLevel*opts.structure_slope+1+opts.Lback,0.0], #5
            [3*wave_length+2*opts.tube+opts.structureCrestLevel*opts.structure_slope+1+opts.Lback,opts.tank_height], #6
            [3*wave_length+2*opts.tube+opts.structureCrestLevel*opts.structure_slope+1,opts.tank_height],#7 abs zone upper boundary
            [wave_length,opts.tank_height], #8
            [0.0,opts.tank_height], #9
            [-wave_length,opts.tank_height], #10
            [-wave_length,0.], #11
            ]
         

vertexFlags=np.array([1, #0 
                        1, #1 
                        1, #2
                        1, #3
                        1, #4
                        2, #5
                        3, #6
                        3, #7 upper boundary abs zone
                        4, #8
                        4, #9
                        4, #10
                        4, #11
                        ])           

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,8],
          [8,9],
          [9,10],
          [10,11],
          [11,0],
          [1,8],
          [4,7],
              ]

segmentFlags=np.array([ 1, #[0,1] 
                        1,#4, #[1,2] 
                        1, #[2,3]
                        2, #[3,4]
                        1, #[4,5]
                        2, #[5,6]
                        3, #[6,7]
                        3, #[7,8]
                        3, #[8,9]
                        3, #[9,10]
                        4, #[10,11]
                        1, #[11,0]
                        5, #[1,8]
                        5, #[4,7]
                      ])


regions=[[5,0.3],[-0.5*wave_length,0.3],[3*wave_length+0.2+1+opts.structureCrestLevel*opts.structure_slope+2,0.4]]         
regionFlags =np.array([1,2,3])        

tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)

################# OBSTACLE DEFINITION  ##############################
                      
obs_boundaryOrientations = {'obstacle': None}

obs_boundaryTags = {'obstacle' : 1,}

obs_vertices=[
            [wave_length+opts.tube,0],#0
            [wave_length+opts.tube,-opts.obs_depth],#1
            [3*wave_length+opts.tube+opts.structureCrestLevel*opts.structure_slope+1,-opts.obs_depth],#2
            [3*wave_length+opts.tube+opts.structureCrestLevel*opts.structure_slope+1,0.],#3
            [3*wave_length+opts.tube+opts.structureCrestLevel*opts.structure_slope,0.],#4
            [3*wave_length+opts.tube+opts.structureCrestLevel*opts.structure_slope,opts.structureCrestLevel],#5
            [3*wave_length+opts.tube,0],#6
            ]  

obs_vertexFlags=np.array([1, #11 
                        1, #12
                        1, #13
                        1, #14
                        1, #15
                        1, #16
                        1, #17
                        ])           

obs_segments=[[0,1],
              [1,2],
              [2,3],
              [3,4],
              [4,5],
              [5,6],
              [6,0],
              ]

obs_segmentFlags=np.array([ 1, #[11,12] 
                            1, #[12,13] 
                            1, #[13,14]
                            1, #[14,15]
                            1, #[15,16]
                            1, #[16,17]
                            1, #[17,11]
                            ])

obs_regions=[[2*wave_length, -0.2]]         
obs_regionFlags =np.array([1])        

obstacle = st.CustomShape(domain, vertices=obs_vertices, vertexFlags=obs_vertexFlags,
                      segments=obs_segments, segmentFlags=obs_segmentFlags,
                      regions=obs_regions, regionFlags=obs_regionFlags,
                      boundaryTags=obs_boundaryTags, 
                      boundaryOrientations=obs_boundaryOrientations)


obstacle.setHoles([[2*wave_length, -0.2]])

# Mesh Refinement
he=wave_length/opts.refinement_level
ecH=opts.ecH
smoothing=ecH*he

# Boundary Conditions
# --- Tank
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

tank.BC['sponge'].setNonMaterial()
for bc in obstacle.BC_list:
    bc.setFreeSlip()

# -------- SET GENERATION AND ABSORPTION ZONES ------- #

dragAlpha = 5*(2*np.pi/opts.Tp)/1e-6
left = True
he=wave_length/opts.refinement_level
smoothing =he*3.
tank.setGenerationZones(flags=2,
                        epsFact_solid=wave_length/2.,
                        center=(-wave_length/2,0.35),
                        orientation=(1.,0.,0.),
                        waves=wave,
                        dragAlpha=dragAlpha)

tank.setAbsorptionZones(flags=3,
                        epsFact_solid=wave_length/2.,
                        center=(3*wave_length+0.2+1+opts.structureCrestLevel*opts.structure_slope+2,0.4),
                        orientation=(-1.,0.,0.),
                        dragAlpha=dragAlpha)

domain.MeshOptions.he = he
st.assembleDomain(domain)

# ------ INITIAL CONDITIONS -----#

def signedDistance(x):
    return x[1] - opts.mwl

class P_IC:
    def __init__(self):
        self.waterLevel=opts.mwl
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(opts.tank_height - opts.mwl)*opts.rho_1*opts.g[1] - (opts.mwl - x[1])*opts.rho_0*opts.g[1]
        else:
            return -(opts.tank_height - opts.mwl)*opts.rho_1*opts.g[1]
class AtRest:
    def uOfXT(self, x, t):
        return 0.0

class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(opts.ecH * he, signedDistance(x))

class LS_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions = {'pressure': P_IC(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest(),
                     'vel_w': AtRest(),
                     'vof': VOF_IC(),
                     'ncls': LS_IC(),
                     'rdls': LS_IC(),}


# ------ SET TIME CONTROLS ----- #

Duration= Tend/opts.fract
dt_output = opts.Tp/opts.Np

outputStepping = TpFlow.OutputStepping(final_time=Duration,
                                       dt_init=opts.dt_init,
                                       dt_output=dt_output,
                                       nDTout=None,
                                       dt_fixed=None)

# ----- ASSEMBLE THE TWO PHASE FLOW OBJECT ----- #

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=domain.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             )

# ----- ADDITIONS MODEL SETTINGS ----- #

params = myTpFlowProblem.Parameters

myTpFlowProblem.useSuperLu=False#True
params.physical.surf_tension_coeff = opts.sigma_01

# index in order of
m = params.Models
m.rans2p.index = 0
m.vof.index = 1
m.ncls.index = 2
m.rdls.index = 3
m.mcorr.index = 4
m.rdls.n.maxLineSearches=0
m.rdls.n.maxNonlinearIts=50

#Gauges 
"""
column_gauge_locations = []
point_gauge_locations = []

gauge_y = opts.mwl - 0.5 * opts.depth
number_of_gauges = (2*wave_length+opts.structureCrestLevel*opts.structure_slope) / opts.gauge_dx + 1
number_of_gauges_2 = opts.Lback / opts.gauge_dx_2 + 1

for gauge_x in np.linspace(0, 2*wave_length+opts.structureCrestLevel*opts.structure_slope, number_of_gauges):
    point_gauge_locations.append((gauge_x, gauge_y, 0))
    column_gauge_locations.append(((gauge_x, 0., 0.),
                                       (gauge_x, tank_dim[1], 0.)))

for gauge_x in np.linspace(2*wave_length+opts.structureCrestLevel*opts.structure_slope, tank_dim[0], number_of_gauges_2):
    point_gauge_locations.append((gauge_x, gauge_y, 0))
    column_gauge_locations.append(((gauge_x, 0., 0.),
                                       (gauge_x, tank_dim[1], 0.)))




tank.attachPointGauges('twp',
                           gauges=((('p',), point_gauge_locations),),
                           fileName='pressure_gaugeArray.csv')


tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),

                                  fileName='column_gauges.csv') 

"""
myTpFlowProblem.Parameters.Models.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
