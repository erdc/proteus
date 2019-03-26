"""
Optimized  Two-Phase Reynolds Averaged Navier-Stokes
"""
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import math
import proteus
import sys
from proteus.mprans.cRANS2P import *
from proteus.mprans.cRANS2P2D import *
from proteus import Profiling
from proteus import LinearAlgebraTools as LAT
from proteus.Comm import (globalSum,
                          globalMax)

class SubgridError(proteus.SubgridError.SGE_base):
    """
    Create a SubgridError  object for two-phase incompressible flow

    The VMS subgrid error approximation

    Parameters
    ----------
        coefficients  : proteus.TransportCoefficients.TC_base
            The coefficients object
        nd            : int
            Number of space  dimensions
        lag           : bool
            Use prior time step to calculate
        nStepsToDelay : int
            Lag only after nSteps
        hFactor       : float
            scaling factor based on order
        noPressureStabilization : bool
            turn off pressure stab
    """

    def __init__(self, coefficients, nd, lag=False, nStepsToDelay=1, hFactor=1.0, noPressureStabilization=False):
        self.noPressureStabilization = noPressureStabilization
        proteus.SubgridError.SGE_base.__init__(self, coefficients, nd, lag)
        coefficients.stencil[0].add(0)
        self.hFactor = hFactor
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            logEvent("RANS2P.SubgridError: lagging requested but must lag the first step; switching lagging off and delaying")
            #Ensure nStepsToDelay value makes sense by throwing an error
            if(self.nStepsToDelay is None or self.nStepsToDelay < 1):
                sys.exit("RANS2P.SubgridError: nStepsToDelay cannot be None or < 1 with lagging, please specify an integer: nStepsToDelay>=1")

    def initializeElementQuadrature(self, mesh, t, cq):
        """
        Allocated or set additional arrays for values at element quadrature

        Parameters
        ----------
        mesh : proteus.MeshTools.Mesh
           The mesh for the domain
        t    : float
           The current time
        cq   : dict
           The dictionary of element quadrature arrrays
        """
        import copy
        self.cq = cq
        #default behavior is to set v_last to point to cq[('velocity')]
        self.v_last = self.cq[('velocity', 0)]

    def updateSubgridErrorHistory(self, initializationPhase=False):
        if self.lag:
            if self.nSteps > self.nStepsToDelay:
                #at this point, v_last must be a separate object
                #update the values of v_last based on cq[('velocity')]
                self.v_last[:] = self.cq[('velocity', 0)]
            elif self.nSteps == self.nStepsToDelay:
                logEvent("RANS2P.SubgridError: switched to lagged subgrid error")
                #create a separate object identical to cq[('velocity')] 
                self.v_last = self.cq[('velocity', 0)].copy()
            else:
                pass
            self.nSteps += 1


    def calculateSubgridError(self, q):
        pass


class NumericalFlux(proteus.NumericalFlux.NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior):
    hasInterior = False

    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        proteus.NumericalFlux.NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior.__init__(self, vt, getPointwiseBoundaryConditions,
                                                                                                     getAdvectiveFluxBoundaryConditions,
                                                                                                     getDiffusiveFluxBoundaryConditions, getPeriodicBoundaryConditions)
        self.penalty_constant = 2.0
        self.includeBoundaryAdjoint = True
        self.boundaryAdjoint_sigma = 1.0
        self.hasInterior = False


class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):
    def __init__(self, coefficients, nd, shockCapturingFactor=0.25, lag=False, nStepsToDelay=1):
        proteus.ShockCapturing.ShockCapturing_base.__init__(self, coefficients, nd, shockCapturingFactor, lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            logEvent("RANS2P.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            #Ensure nStepsToDelay value makes sense by throwing an error
            if(self.nStepsToDelay is None or self.nStepsToDelay < 1):
                sys.exit("RANS2P.ShockCapturing: nStepsToDelay cannot be None or < 1 with lagging, please specify an integer: nStepsToDelay>=1")

    def initializeElementQuadrature(self, mesh, t, cq):
        self.mesh = mesh
        self.numDiff = {}
        self.numDiff_last = {}
        #default behavior is to set both numDiff_last and numDiff to point to cq[('numDiff')]
        for ci in range(1, 4):
            self.numDiff[ci] = cq[('numDiff', ci, ci)]
            self.numDiff_last[ci] = cq[('numDiff', ci, ci)]

    def updateShockCapturingHistory(self):
        if self.lag:
            if self.nSteps > self.nStepsToDelay:
                #numDiff_last is a different object
                #update the values of numDiff_last based on numDiff, 
                for ci in range(1, 4):
                    self.numDiff_last[ci][:] = self.numDiff[ci]
            elif self.nSteps == self.nStepsToDelay:
                logEvent("RANS2P.ShockCapturing: switched to lagged shock capturing")
                #if lagging, then create a separate object identical to numDiff 
                for ci in range(1, 4):
                    self.numDiff_last[ci] = self.numDiff[ci].copy()
            else:
                pass
            self.nSteps += 1

            logEvent("RANS2P: max numDiff_1 %e numDiff_2 %e numDiff_3 %e" % (globalMax(self.numDiff_last[1].max()),
                                                                         globalMax(self.numDiff_last[2].max()),
                                                                         globalMax(self.numDiff_last[3].max())))
class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for two incompresslble fluids governed by the Navier-Stokes equations and separated by a sharp interface represented by a level set function

    Notes
    -----
    The PRESSURE_PROJECTION_STABILIZATION flag allows the user to use
    the Bochev-Dohrmann-Gunzburger stabilization introduced in
    Stabilization of Low-Order Mixed Finite Elements for the  Stokes
    Equation.  This option should be turned off for most problems,
    but in some instances it may produced better preconditioning
    results than the full SGE approach.
    """
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd

    def __init__(self,
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2, nu_0=1.004e-6,
                 rho_1=1.205, nu_1=1.500e-5,
                 g=[0.0, 0.0, -9.8],
                 nd=3,
                 ME_model=0,
                 CLSVOF_model=None,
                 LS_model=None,
                 VF_model=None,
                 KN_model=None,
                 Closure_0_model=None,  # Turbulence closure model
                 Closure_1_model=None,  # Second possible Turbulence closure model
                 epsFact_density=None,
                 stokes=False,
                 sd=True,
                 movingDomain=False,
                 useVF=0.0,
                 useRBLES=0.0,
                 useMetrics=0.0,
                 useConstant_he=False,
                 dragAlpha=0.0,
                 dragBeta=0.0,
                 setParamsFunc=None,  # uses setParamsFunc if given
                 dragAlphaTypes=None,  # otherwise can use element constant values
                 dragBetaTypes=None,  # otherwise can use element constant values
                 porosityTypes=None,
                 killNonlinearDrag=False,
                 epsFact_solid=None,
                 eb_adjoint_sigma=1.0,
                 eb_penalty_constant=10.0,
                 forceStrongDirichlet=False,
                 turbulenceClosureModel=0,  # 0=No Model, 1=Smagorinksy, 2=Dynamic Smagorinsky, 3=K-Epsilon, 4=K-Omega
                 smagorinskyConstant=0.1,
                 barycenters=None,
                 NONCONSERVATIVE_FORM=0.0,
                 MOMENTUM_SGE=1.0,
                 PRESSURE_SGE=1.0,
                 VELOCITY_SGE=1.0,
                 PRESSURE_PROJECTION_STABILIZATION=0.0,
                 phaseFunction=None,
                 LAG_LES=1.0,
                 use_ball_as_particle=1,
                 ball_center=None,
                 ball_radius=None,
                 ball_velocity=None,
                 ball_angular_velocity=None,
                 ball_center_acceleration=None,
                 ball_angular_acceleration=None,
                 ball_density=None,
                 particle_velocities=None,
                 particle_centroids=None,
                 particle_sdfList=[],
                 particle_velocityList=[],
                 nParticles = 0,
                 particle_epsFact=3.0,
                 particle_alpha=1000.0,
                 particle_beta=1000.0,
                 particle_penalty_constant=1000.0,
                 particle_nitsche=1.0,
                 nullSpace='NoNullSpace'):
        self.use_pseudo_penalty = 0
        self.use_ball_as_particle = use_ball_as_particle
        self.nParticles = nParticles
        self.particle_nitsche = particle_nitsche
        self.particle_epsFact = particle_epsFact
        self.particle_alpha = particle_alpha
        self.particle_beta = particle_beta
        self.particle_penalty_constant = particle_penalty_constant
        self.particle_netForces = np.zeros((3*self.nParticles, 3), 'd')#####[total_force_1,total_force_2,...,stress_1,stress_2,...,pressure_1,pressure_2,...]  
        self.particle_netMoments = np.zeros((self.nParticles, 3), 'd')
        self.particle_surfaceArea = np.zeros((self.nParticles,), 'd')
        if ball_center is None:
            self.ball_center = 1e10*numpy.ones((self.nParticles,3),'d')
        else:
            self.ball_center = ball_center
        
        if ball_radius is None:
            self.ball_radius = 1e10*numpy.ones((self.nParticles,1),'d')
        else:
            self.ball_radius = ball_radius

        if ball_velocity is None:
            self.ball_velocity = numpy.zeros((self.nParticles,3),'d')
        else:
            self.ball_velocity = ball_velocity

        if ball_angular_velocity is None:
            self.ball_angular_velocity = numpy.zeros((self.nParticles,3),'d')
        else:
            self.ball_angular_velocity = ball_angular_velocity

        if ball_center_acceleration is None:
            self.ball_center_acceleration = numpy.zeros((self.nParticles,3),'d')
        else:
            self.ball_center_acceleration = ball_center_acceleration

        if ball_angular_acceleration is None:
            self.ball_angular_acceleration = numpy.zeros((self.nParticles,3),'d')
        else:
            self.ball_angular_acceleration = ball_angular_acceleration

        if ball_density is None:
            self.ball_density = rho_0*numpy.ones((self.nParticles,1),'d')
        else:
            self.ball_density = ball_density
        if particle_centroids is None:
            self.particle_centroids = 1e10*numpy.zeros((self.nParticles,3),'d')
        else:
            self.particle_centroids = particle_centroids
        self.particle_sdfList = particle_sdfList
        self.particle_velocityList = particle_velocityList
        
        self.LAG_LES=LAG_LES
        self.phaseFunction=phaseFunction
        self.NONCONSERVATIVE_FORM=NONCONSERVATIVE_FORM
        self.MOMENTUM_SGE=MOMENTUM_SGE
        self.PRESSURE_SGE=PRESSURE_SGE
        self.VELOCITY_SGE=VELOCITY_SGE
        self.PRESSURE_PROJECTION_STABILIZATION=PRESSURE_PROJECTION_STABILIZATION
        self.barycenters=barycenters
        self.smagorinskyConstant = smagorinskyConstant
        self.turbulenceClosureModel = turbulenceClosureModel
        self.forceStrongDirichlet = forceStrongDirichlet
        self.eb_adjoint_sigma = eb_adjoint_sigma
        self.eb_penalty_constant = eb_penalty_constant
        self.movingDomain = movingDomain
        self.epsFact_solid = epsFact_solid
        self.useConstant_he = useConstant_he
        self.useVF = useVF
        self.useRBLES = useRBLES
        self.useMetrics = useMetrics
        self.sd = sd
        if epsFact_density is not None:
            self.epsFact_density = epsFact_density
        else:
            self.epsFact_density = epsFact
        self.stokes = stokes
        self.ME_model = ME_model
        self.CLSVOF_model = CLSVOF_model
        self.LS_model = LS_model
        self.VF_model = VF_model
        self.KN_model = KN_model
        self.Closure_0_model = Closure_0_model
        self.Closure_1_model = Closure_1_model
        self.epsFact = epsFact
        self.eps = None
        self.sigma = sigma
        self.rho_0 = rho_0
        self.nu_0 = nu_0
        # cek for debugging using single phase test problems
        self.rho = rho_0
        self.nu = nu_0
        self.rho_1 = rho_1
        self.nu_1 = nu_1
        self.g = numpy.array(g)
        self.nd = nd
        #
        self.dragAlpha = dragAlpha
        self.dragBeta = dragBeta
        self.setParamsFunc = setParamsFunc
        self.dragAlphaTypes = dragAlphaTypes
        self.dragBetaTypes = dragBetaTypes
        self.porosityTypes = porosityTypes
        self.killNonlinearDrag = int(killNonlinearDrag)
        self.linearDragFactor = 1.0
        self.nonlinearDragFactor = 1.0
        if self.killNonlinearDrag:
            self.nonlinearDragFactor = 0.0
        self.nullSpace = nullSpace
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        if nd == 2:
            variableNames = ['p', 'u', 'v']
            mass = {1: {1: 'linear'},
                    2: {2: 'linear'}}
            advection = {0: {0: 'linear',
                             1: 'linear',
                             2: 'linear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'},
                         2: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'}}
            diffusion = {1: {1: {1: 'constant'}, 2: {2: 'constant'}},
                         2: {2: {2: 'constant'}, 1: {1: 'constant'}}}
            sdInfo = {(1, 1): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (1, 2): (numpy.array([0, 0, 1], dtype='i'),
                               numpy.array([0], dtype='i')),
                      (2, 2): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (2, 1): (numpy.array([0, 1, 1], dtype='i'),
                               numpy.array([1], dtype='i'))}
            potential = {1: {1: 'u'},
                         2: {2: 'u'}}
            reaction = {0: {0: 'constant'},
                        1: {1: 'nonlinear', 2: 'nonlinear'},
                        2: {1: 'nonlinear', 2: 'nonlinear'}}
            hamiltonian = {1: {0: 'linear'},
                           2: {0: 'linear'}}
            TC_base.__init__(self,
                             3,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion=sd,
                             movingDomain=movingDomain)
            self.vectorComponents = [1, 2]
        elif nd == 3:
            variableNames = ['p', 'u', 'v', 'w']
            mass = {1: {1: 'linear'},
                    2: {2: 'linear'},
                    3: {3: 'linear'}}
            advection = {0: {1: 'linear',
                             2: 'linear',
                             3: 'linear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear'},
                         2: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear'},
                         3: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear'}}
            diffusion = {1: {1: {1: 'constant'}, 2: {2: 'constant'}, 3: {3: 'constant'}},
                         2: {1: {1: 'constant'}, 2: {2: 'constant'}, 3: {3: 'constant'}},
                         3: {1: {1: 'constant'}, 2: {2: 'constant'}, 3: {3: 'constant'}}}
            sdInfo = {}
            sdInfo = {(1, 1): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i')),
                      (1, 2): (numpy.array([0, 0, 1, 1], dtype='i'), numpy.array([0], dtype='i')),
                      (1, 3): (numpy.array([0, 0, 0, 1], dtype='i'), numpy.array([0], dtype='i')),
                      (2, 1): (numpy.array([0, 1, 1, 1], dtype='i'), numpy.array([1], dtype='i')),
                      (2, 2): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i')),
                      (2, 3): (numpy.array([0, 0, 0, 1], dtype='i'), numpy.array([1], dtype='i')),
                      (3, 1): (numpy.array([0, 1, 1, 1], dtype='i'), numpy.array([2], dtype='i')),
                      (3, 2): (numpy.array([0, 0, 1, 1], dtype='i'), numpy.array([2], dtype='i')),
                      (3, 3): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i'))}
            potential = {1: {1: 'u'},
                         2: {2: 'u'},
                         3: {3: 'u'}}
            reaction = {0: {0: 'constant'},
                        1: {1: 'nonlinear', 2: 'nonlinear', 3: 'nonlinear'},
                        2: {1: 'nonlinear', 2: 'nonlinear', 3: 'nonlinear'},
                        3: {1: 'nonlinear', 2: 'nonlinear', 3: 'nonlinear'}}
            hamiltonian = {1: {0: 'linear'},
                           2: {0: 'linear'},
                           3: {0: 'linear'}}
            TC_base.__init__(self,
                             4,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion=sd,
                             movingDomain=movingDomain)
            self.vectorComponents = [1, 2, 3]

    def attachModels(self, modelList):
        # level set
        self.model = modelList[self.ME_model]
        self.model.q['phi_solid'] = self.q_phi_solid
        self.model.q['velocity_solid'] = self.q_velocity_solid
        if self.CLSVOF_model is not None: # use CLSVOF
            # LS part #
            self.q_phi = modelList[self.CLSVOF_model].q[('u', 0)]
            self.ebq_phi = None # Not used. What is this for?
            self.ebqe_phi = modelList[self.CLSVOF_model].ebqe[('u', 0)]
            self.bc_ebqe_phi = modelList[self.CLSVOF_model].ebqe[('u', 0)] #Dirichlet BCs for level set. I don't have it since I impose 1 or -1. Therefore I attach the soln at boundary
            self.q_n = modelList[self.CLSVOF_model].q[('grad(u)', 0)]
            self.ebq_n = None # Not used. What is this for?
            self.ebqe_n = modelList[self.CLSVOF_model].ebqe[('grad(u)', 0)]
            # VOF part #
            self.q_vf = modelList[self.CLSVOF_model].q[('H(u)', 0)]
            self.ebq_vf = None# Not used. What is this for?
            self.ebqe_vf = modelList[self.CLSVOF_model].ebqe[('H(u)', 0)]
            self.bc_ebqe_vf = 0.5*(1.0+modelList[self.CLSVOF_model].numericalFlux.ebqe[('u',0)]) # Dirichlet BCs for VOF. What I have is BCs for Signed function
        else: # use NCLS-RDLS-VOF-MCorr instead
            if self.LS_model is not None:
                self.q_phi = modelList[self.LS_model].q[('u', 0)]
                if ('u', 0) in modelList[self.LS_model].ebq:
                    self.ebq_phi = modelList[self.LS_model].ebq[('u', 0)]
                else:
                    self.ebq_phi = None
                self.ebqe_phi = modelList[self.LS_model].ebqe[('u', 0)]
                self.bc_ebqe_phi = modelList[self.LS_model].numericalFlux.ebqe[('u', 0)]
                # normal
                self.q_n = modelList[self.LS_model].q[('grad(u)', 0)]
                if ('grad(u)', 0) in modelList[self.LS_model].ebq:
                    self.ebq_n = modelList[self.LS_model].ebq[('grad(u)', 0)]
                else:
                    self.ebq_n = None
                self.ebqe_n = modelList[self.LS_model].ebqe[('grad(u)', 0)]
            else:
                self.q_phi = 10.0 * numpy.ones(self.model.q[('u', 1)].shape, 'd')
                self.ebqe_phi = 10.0 * numpy.ones(self.model.ebqe[('u', 1)].shape, 'd')
                self.bc_ebqe_phi = 10.0 * numpy.ones(self.model.ebqe[('u', 1)].shape, 'd')
                self.q_n = numpy.ones(self.model.q[('velocity', 0)].shape, 'd')
                self.ebqe_n = numpy.ones(self.model.ebqe[('velocity', 0)].shape, 'd')
            if self.VF_model is not None:
                self.q_vf = modelList[self.VF_model].q[('u', 0)]
                if ('u', 0) in modelList[self.VF_model].ebq:
                    self.ebq_vf = modelList[self.VF_model].ebq[('u', 0)]
                else:
                    self.ebq_vf = None
                self.ebqe_vf = modelList[self.VF_model].ebqe[('u', 0)]
                self.bc_ebqe_vf = modelList[self.VF_model].numericalFlux.ebqe[('u', 0)]
            else:
                self.q_vf = numpy.zeros(self.model.q[('u', 1)].shape, 'd')
                self.ebqe_vf = numpy.zeros(self.model.ebqe[('u', 1)].shape, 'd')
                self.bc_ebqe_vf = numpy.zeros(self.model.ebqe[('u', 1)].shape, 'd')
        # curvature
        if self.KN_model is not None:
            self.q_kappa = modelList[self.KN_model].q[('u', 0)]
            self.ebqe_kappa = modelList[self.KN_model].ebqe[('u', 0)]
            if ('u', 0) in modelList[self.KN_model].ebq:
                self.ebq_kappa = modelList[self.KN_model].ebq[('u', 0)]
            else:
                self.ebq_kappa = None
        else:
            self.q_kappa = -numpy.ones(self.model.q[('u', 1)].shape, 'd')
            self.ebqe_kappa = -numpy.ones(self.model.ebqe[('u', 1)].shape, 'd')
        # Turbulence Closures
        # only option for now is k-epsilon
        self.q_turb_var = {}
        self.q_turb_var_grad = {}
        self.ebqe_turb_var = {}
        if self.Closure_0_model is not None:
            self.q_turb_var[0] = modelList[self.Closure_0_model].q[('u', 0)]
            self.q_turb_var_grad[0] = modelList[self.Closure_0_model].q[('grad(u)', 0)]
            self.ebqe_turb_var[0] = modelList[self.Closure_0_model].ebqe[('u', 0)]
        else:
            self.q_turb_var[0] = numpy.ones(self.model.q[('u', 1)].shape, 'd')
            self.q_turb_var_grad[0] = numpy.ones(self.model.q[('grad(u)', 1)].shape, 'd')
            self.ebqe_turb_var[0] = numpy.ones(self.model.ebqe[('u', 1)].shape, 'd')
        if self.Closure_1_model is not None:
            self.q_turb_var[1] = modelList[self.Closure_1_model].q[('u', 0)]
            self.ebqe_turb_var[1] = modelList[self.Closure_1_model].ebqe[('u', 0)]
        else:
            self.q_turb_var[1] = numpy.ones(self.model.q[('u', 1)].shape, 'd')
            self.ebqe_turb_var[1] = numpy.ones(self.model.ebqe[('u', 1)].shape, 'd')
        if self.epsFact_solid is None:
            self.epsFact_solid = numpy.ones(self.model.mesh.elementMaterialTypes.max() + 1)
        assert len(self.epsFact_solid) > self.model.mesh.elementMaterialTypes.max(
        ), "epsFact_solid  array is not large  enough for the materials  in this mesh; length must be greater  than largest  material type ID"
        if self.phaseFunction != None:
            from proteus.ctransportCoefficients import smoothedHeaviside
            if self.useConstant_he:
                self.elementDiameter = self.mesh.elementDiametersArray.copy()
                self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
            else:
                self.elementDiameter = self.mesh.elementDiametersArray
            for i, quad_pts in enumerate(self.model.q['x']):
                for j, pt in enumerate(quad_pts):
                    he = self.elementDiameter[i]
                    self.q_phi[i, j] = self.phaseFunction(pt)
                    self.q_vf[i, j] = smoothedHeaviside(self.epsFact * he, self.phaseFunction(pt))
            for i, quad_pts in enumerate(self.model.ebqe['x']):
                for j, pt in enumerate(quad_pts):
                    he = self.elementDiameter[i]
                    self.ebqe_phi[i, j] = self.phaseFunction(pt)
                    self.ebqe_vf[i, j] = smoothedHeaviside(self.epsFact * he, self.phaseFunction(pt))
        self.particle_signed_distances        = 1e10*numpy.ones((self.nParticles,self.model.q['x'].shape[0],self.model.q['x'].shape[1]),'d')
        self.particle_signed_distance_normals = 1e10*numpy.ones((self.nParticles,self.model.q['x'].shape[0],self.model.q['x'].shape[1], 3),'d')
        self.particle_velocities              = 1e10*numpy.ones((self.nParticles,self.model.q['x'].shape[0],self.model.q['x'].shape[1], 3),'d')
        self.phisField                        = 1e10*numpy.ones((self.model.q['x'].shape[0],self.model.q['x'].shape[1]), 'd')
        self.ebq_global_phi_s        = numpy.ones((self.model.ebq_global['x'].shape[0],self.model.ebq_global['x'].shape[1]),'d') * 1e10
        self.ebq_global_grad_phi_s   = numpy.ones((self.model.ebq_global['x'].shape[0],self.model.ebq_global['x'].shape[1],3),'d') * 1e10
        self.ebq_particle_velocity_s = numpy.ones((self.model.ebq_global['x'].shape[0],self.model.ebq_global['x'].shape[1],3),'d') * 1e10
        self.p_old_dof = self.model.u[0].dof.copy()
        self.u_old_dof = self.model.u[1].dof.copy()
        self.v_old_dof = self.model.u[2].dof.copy()
        self.w_old_dof = self.model.u[3].dof.copy()

    def initializeMesh(self, mesh):
        
        self.phi_s = numpy.ones(mesh.nodeArray.shape[0], 'd')*1e10#
        # cek we eventually need to use the local element diameter
        self.eps_density = self.epsFact_density * mesh.h
        self.eps_viscosity = self.epsFact * mesh.h
        self.mesh = mesh
        self.elementMaterialTypes = mesh.elementMaterialTypes
        nBoundariesMax = int(globalMax(max(self.mesh.elementBoundaryMaterialTypes))) + 1
        self.wettedAreas = numpy.zeros((nBoundariesMax,), 'd')
        self.netForces_p = numpy.zeros((nBoundariesMax, 3), 'd')
        self.netForces_v = numpy.zeros((nBoundariesMax, 3), 'd')
        self.netMoments = numpy.zeros((nBoundariesMax, 3), 'd')
        if self.barycenters is None:
            self.barycenters = numpy.zeros((nBoundariesMax, 3), 'd')
        comm = Comm.get()
        import os
        if comm.isMaster():
            self.timeHistory = open(os.path.join(proteus.Profiling.logDir, "timeHistory.txt"), "w")
            self.wettedAreaHistory = open(os.path.join(proteus.Profiling.logDir, "wettedAreaHistory.txt"), "w")
            self.forceHistory_p = open(os.path.join(proteus.Profiling.logDir, "forceHistory_p.txt"), "w")
            self.forceHistory_v = open(os.path.join(proteus.Profiling.logDir, "forceHistory_v.txt"), "w")
            self.momentHistory = open(os.path.join(proteus.Profiling.logDir, "momentHistory.txt"), "w")
            if self.nParticles:
                self.particle_forceHistory = open(os.path.join(proteus.Profiling.logDir, "particle_forceHistory.txt"), "w")
                self.particle_pforceHistory = open(os.path.join(proteus.Profiling.logDir, "particle_pforceHistory.txt"), "w")
                self.particle_vforceHistory = open(os.path.join(proteus.Profiling.logDir, "particle_vforceHistory.txt"), "w")
                self.particle_momentHistory = open(os.path.join(proteus.Profiling.logDir, "particle_momentHistory.txt"), "w")
        self.comm = comm
    # initialize so it can run as single phase

    def initializeElementQuadrature(self, t, cq):
        # VRANS
        self.numerical_viscosity = numpy.zeros(cq[('u', 1)].shape, 'd')
        self.q_phi_solid = numpy.ones(cq[('u', 1)].shape, 'd')
        self.q_velocity_solid = numpy.zeros(cq[('velocity', 0)].shape, 'd')
        self.q_porosity = numpy.ones(cq[('u', 1)].shape, 'd')
        self.q_dragAlpha = numpy.ones(cq[('u', 1)].shape, 'd')
        self.q_dragAlpha.fill(self.dragAlpha)
        self.q_dragBeta = numpy.ones(cq[('u', 1)].shape, 'd')
        self.q_dragBeta.fill(self.dragBeta)
        if self.setParamsFunc is not None:
            self.setParamsFunc(cq['x'], self.q_porosity, self.q_dragAlpha, self.q_dragBeta)
        else:
            # TODO make loops faster
            if self.porosityTypes is not None:
                for eN in range(self.q_porosity.shape[0]):
                    self.q_porosity[eN, :] = self.porosityTypes[self.elementMaterialTypes[eN]]
            if self.dragAlphaTypes is not None:
                for eN in range(self.q_dragAlpha.shape[0]):
                    self.q_dragAlpha[eN, :] = self.dragAlphaTypes[self.elementMaterialTypes[eN]]
            if self.dragBetaTypes is not None:
                for eN in range(self.q_dragBeta.shape[0]):
                    self.q_dragBeta[eN, :] = self.dragBetaTypes[self.elementMaterialTypes[eN]]
        cq['velocityError'] = cq[('velocity', 0)].copy()
        #

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        # VRANS
        self.ebq_porosity = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragAlpha = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragAlpha.fill(self.dragAlpha)
        self.ebq_dragBeta = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragBeta.fill(self.dragBeta)
        if self.setParamsFunc is not None:
            self.setParamsFunc(cebq['x'], self.ebq_porosity, self.ebq_dragAlpha, self.ebq_dragBeta)
        # TODO which mean to use or leave discontinuous
        # TODO make loops faster
        if self.porosityTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
                ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
                ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 1]
                avg = 0.5 * (self.porosityTypes[self.elementMaterialTypes[eN_left]] +
                             self.porosityTypes[self.elementMaterialTypes[eN_right]])
                self.ebq_porosity[eN_left, ebN_element_left, :] = self.porosityTypes[self.elementMaterialTypes[eN_left]]
                self.ebq_porosity[eN_right, ebN_element_right, :] = self.porosityTypes[self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
                self.ebq_porosity[eN, ebN_element, :] = self.porosityTypes[self.elementMaterialTypes[eN]]
        if self.dragAlphaTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
            ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
            ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 1]
            avg = 0.5 * (self.dragAlphaTypes[self.elementMaterialTypes[eN_left]] +
                         self.dragAlphaTypes[self.elementMaterialTypes[eN_right]])
            self.ebq_dragAlpha[eN_left, ebN_element_left, :] = self.dragAlphaTypes[self.elementMaterialTypes[eN_left]]
            self.ebq_dragAlpha[eN_right, ebN_element_right, :] = self.dragAlphaTypes[self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
                self.ebq_dragAlpha[eN, ebN_element, :] = self.dragAlphaTypes[self.elementMaterialTypes[eN]]
        if self.dragBetaTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
            ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
            ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 1]
            avg = 0.5 * (self.dragBetaTypes[self.elementMaterialTypes[eN_left]] +
                         self.dragBetaTypes[self.elementMaterialTypes[eN_right]])
            self.ebq_dragBeta[eN_left, ebN_element_left, :] = self.dragBetaTypes[self.elementMaterialTypes[eN_left]]
            self.ebq_dragBeta[eN_right, ebN_element_right, :] = self.dragBetaTypes[self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
                self.ebq_dragBeta[eN, ebN_element, :] = self.dragBetaTypes[self.elementMaterialTypes[eN]]
         #

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        # VRANS
        logEvent("ebqe_global allocations in coefficients")
        self.ebqe_porosity = numpy.ones(cebqe[('u', 1)].shape, 'd')
        self.ebqe_dragAlpha = numpy.ones(cebqe[('u', 1)].shape, 'd')
        self.ebqe_dragAlpha.fill(self.dragAlpha)
        self.ebqe_dragBeta = numpy.ones(cebqe[('u', 1)].shape, 'd')
        self.ebqe_dragBeta.fill(self.dragBeta)
        logEvent("porosity and drag")
        # TODO make loops faster
        if self.setParamsFunc is not None:
            self.setParamsFunc(cebqe['x'], self.ebqe_porosity, self.ebqe_dragAlpha, self.ebqe_dragBeta)
        else:
            if self.porosityTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_porosity[ebNE, :] = self.porosityTypes[self.elementMaterialTypes[eN]]
            if self.dragAlphaTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_dragAlpha[ebNE, :] = self.dragAlphaTypes[self.elementMaterialTypes[eN]]
            if self.dragBetaTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_dragBeta[ebNE, :] = self.dragBetaTypes[self.elementMaterialTypes[eN]]
        #

    def updateToMovingDomain(self, t, c):
        pass

    def evaluate(self, t, c):
        pass

    def preStep(self, t, firstStep=False):
        self.model.dt_last = self.model.timeIntegration.dt
        if self.nParticles > 0 and self.use_ball_as_particle == 0:
            self.phi_s[:] = 1e10
            self.phisField[:] = 1e10
            self.ebq_global_phi_s[:] = 1e10
            self.ebq_global_grad_phi_s[:] = 1e10
            self.ebq_particle_velocity_s[:] = 1e10

            for i in range(self.nParticles):
                vel = lambda x: self.particle_velocityList[i](t, x)
                sdf = lambda x: self.particle_sdfList[i](t, x)

                for j in range(self.mesh.nodeArray.shape[0]):
                    sdf_at_node, sdNormals = sdf(self.mesh.nodeArray[j, :])
                    if (sdf_at_node < self.phi_s[j]):
                        self.phi_s[j] = sdf_at_node
                for eN in range(self.model.q['x'].shape[0]):
                    for k in range(self.model.q['x'].shape[1]):
                        self.particle_signed_distances[i, eN, k], self.particle_signed_distance_normals[i, eN, k] = sdf(self.model.q['x'][eN, k])
                        self.particle_velocities[i, eN, k] = vel(self.model.q['x'][eN, k])
                        if (abs(self.particle_signed_distances[i, eN, k]) < abs(self.phisField[eN, k])):
                            self.phisField[eN, k] = self.particle_signed_distances[i, eN, k]
                for ebN in range(self.model.ebq_global['x'].shape[0]):
                    for kb in range(self.model.ebq_global['x'].shape[1]):
                        sdf_ebN_kb,sdNormals = sdf(self.model.ebq_global['x'][ebN,kb])
                        if ( sdf_ebN_kb < self.ebq_global_phi_s[ebN,kb]):
                            self.ebq_global_phi_s[ebN,kb]=sdf_ebN_kb
                            self.ebq_global_grad_phi_s[ebN,kb,:]=sdNormals
                            self.ebq_particle_velocity_s[ebN,kb,:] = vel(self.model.ebq_global['x'][ebN,kb])
        # if self.comm.isMaster():
        # print "wettedAreas"
        # print self.wettedAreas[:]
        # print "Forces_p"
        # print self.netForces_p[:,:]
        # print "Forces_v"
        # print self.netForces_v[:,:]

    def postStep(self, t, firstStep=False):
        self.model.dt_last = self.model.timeIntegration.dt
        self.model.q['dV_last'][:] = self.model.q['dV']
        if self.comm.isMaster():
            logEvent("wettedAreas\n"+
                     `self.wettedAreas[:]` +
                     "\nForces_p\n" +
                     `self.netForces_p[:,:]` +
                     "\nForces_v\n" +
                     `self.netForces_v[:,:]`)
            self.timeHistory.write("%21.16e\n" % (t,))
            self.timeHistory.flush()
            self.wettedAreaHistory.write("%21.16e\n" % (self.wettedAreas[-1],))
            self.wettedAreaHistory.flush()
            self.forceHistory_p.write("%21.16e %21.16e %21.16e\n" % tuple(self.netForces_p[-1, :]))
            self.forceHistory_p.flush()
            self.forceHistory_v.write("%21.16e %21.16e %21.16e\n" % tuple(self.netForces_v[-1, :]))
            self.forceHistory_v.flush()
            self.momentHistory.write("%21.15e %21.16e %21.16e\n" % tuple(self.netMoments[-1, :]))
            self.momentHistory.flush()
            if self.nParticles:
                self.particle_forceHistory.write("%21.16e %21.16e %21.16e\n" % tuple(self.particle_netForces[0, :]))
                self.particle_forceHistory.flush()
                self.particle_pforceHistory.write("%21.16e %21.16e %21.16e\n" % tuple(self.particle_netForces[0+self.nParticles, :]))
                self.particle_pforceHistory.flush()
                self.particle_vforceHistory.write("%21.16e %21.16e %21.16e\n" % tuple(self.particle_netForces[0+2*self.nParticles, :]))
                self.particle_vforceHistory.flush()
                self.particle_momentHistory.write("%21.15e %21.16e %21.16e\n" % tuple(self.particle_netMoments[0, :]))
                self.particle_momentHistory.flush()

class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0

    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='RANS2P',
                 reuse_trial_and_test_quadrature=False,
                 sd=True,
                 movingDomain=False):
        self.eb_adjoint_sigma = coefficients.eb_adjoint_sigma
        useConstant_he = coefficients.useConstant_he  # this is a hack to test the effect of using a constant smoothing width
        self.postProcessing = True
        #
        # set the objects describing the method and boundary conditions
        #
        self.movingDomain = coefficients.movingDomain
        self.tLast_mesh = None
        #
        # cek todo clean up these flags in the optimized version
        self.bcsTimeDependent = options.bcsTimeDependent
        self.bcsSet = False
        self.name = name
        self.sd = sd
        self.lowmem = True
        self.timeTerm = True  # allow turning off  the  time derivative
        self.testIsTrial = True
        self.phiTrialIsTrial = True
        self.u = uDict
        self.Hess = False
        if isinstance(self.u[0].femSpace, C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess = True
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh  # assume the same mesh for  all components for now
        self.par_info = LinearAlgebraTools.ParInfo_petsc4py()
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList = None  # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict  # no velocity post-processing for now
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if ci in coefficients.mass:
                    for flag in list(coefficients.mass[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.advection:
                    for flag in list(coefficients.advection[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.diffusion:
                    for diffusionDict in list(coefficients.diffusion[ci].values()):
                        for flag in list(diffusionDict.values()):
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if ci in coefficients.potential:
                    for flag in list(coefficients.potential[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.reaction:
                    for flag in list(coefficients.reaction[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.hamiltonian:
                    for flag in list(coefficients.hamiltonian[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux is not None) or
                                                 (numericalFluxType is not None) or
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        self.nSpace_global = self.u[0].femSpace.nSpace_global  # assume same space dim for all variables
        self.nDOF_trial_element = [u_j.femSpace.max_nDOF_element for u_j in list(self.u.values())]
        self.nDOF_phi_trial_element = [phi_k.femSpace.max_nDOF_element for phi_k in list(self.phi.values())]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in list(self.phi.values())]
        self.nDOF_test_element = [femSpace.max_nDOF_element for femSpace in list(self.testSpace.values())]
        self.nFreeDOF_global = [dc.nFreeDOF_global for dc in list(self.dirichletConditions.values())]
        self.nVDOF_element = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self, self.nFreeVDOF_global)
        #
        # build the quadrature point dictionaries from the input (this
        # is just for convenience so that the input doesn't have to be
        # complete)
        #
        elementQuadratureDict = {}
        elemQuadIsDict = isinstance(elementQuadrature, dict)
        if elemQuadIsDict:  # set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if I in elementQuadrature:
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if I in elementQuadrature:
                        elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if ('numDiff', ci, ci) in elementQuadrature:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if I in elementBoundaryQuadrature:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        (self.elementQuadraturePoints, self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global *
                                                        self.mesh.nElementBoundaries_element *
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)
        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary, 3), 'd')
        self.ebq_global[('totalFlux', 0)] = numpy.zeros((self.mesh.nElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebq_global[('velocityAverage', 0)] = numpy.zeros((self.mesh.nElementBoundaries_global,
                                                               self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.q[('u', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m', 1)] = self.q[('u', 1)]
        self.q[('m', 2)] = self.q[('u', 2)]
        self.q[('m', 3)] = self.q[('u', 3)]
        self.q['rho'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        #self.q[('dV_u',1)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #self.q[('dV_u',2)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #self.q[('dV_u',3)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['dV'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dV_last'] = -1000 * numpy.ones((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('f', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('velocity', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q['velocity_solid'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q['phi_solid'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['x'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, 3), 'd')
        self.q[('cfl', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 1, 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 2, 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 3, 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[('u', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe['eddy_viscosity'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe['eddy_viscosity_last'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc_flag', 0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 3, 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe['penalty'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 3, 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('velocity', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        # VRANS start, defaults to RANS
        self.q[('r', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['eddy_viscosity'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['eddy_viscosity_last'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        # VRANS end
        # RANS 2eq Models start
        self.q[('grad(u)', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('grad(u)', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('grad(u)', 3)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        # probably don't need ebqe gradients
        self.ebqe[('grad(u)', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                 self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('grad(u)', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                 self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('grad(u)', 3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                 self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        # RANS 2eq Models end
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set([('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
        # use post processing tools to get conservative fluxes, None by default
        if self.postProcessing:
            self.q[('v', 0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nDOF_trial_element[0]),
                'd')
            self.q['J'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.q['det(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element),
                'd')
            self.q['inverse(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebq[('v', 0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_trial_element[0]),
                'd')
            self.ebq[('w', 0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_trial_element[0]),
                'd')
            self.ebq['x'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
            self.ebq['hat(x)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
            self.ebq['inverse(J)'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebq['g'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global - 1,
                 self.nSpace_global - 1),
                'd')
            self.ebq['sqrt(det(g))'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebq['n'] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebq[('dS_u', 0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebqe['dS'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebqe[('dS_u', 0)] = self.ebqe['dS']
            self.ebqe['n'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebqe['inverse(J)'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global,
                 self.nSpace_global),
                'd')
            self.ebqe['g'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global - 1,
                 self.nSpace_global - 1),
                'd')
            self.ebqe['sqrt(det(g))'] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
            self.ebq_global['n'] = numpy.zeros(
                (self.mesh.nElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                'd')
            self.ebq_global['x'] = numpy.zeros(
                (self.mesh.nElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 3),
                'd')
        #
        # show quadrature
        #
        logEvent("Dumping quadrature shapes for model %s" % self.name, level=9)
        logEvent("Element quadrature array (q)", level=9)
        for (k, v) in list(self.q.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Element boundary quadrature (ebq)", level=9)
        for (k, v) in list(self.ebq.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Global element boundary quadrature (ebq_global)", level=9)
        for (k, v) in list(self.ebq_global.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Exterior element boundary quadrature (ebqe)", level=9)
        for (k, v) in list(self.ebqe.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Interpolation points for nonlinear diffusion potential (phi_ip)", level=9)
        for (k, v) in list(self.phi_ip.items()):
            logEvent(str((k, v.shape)), level=9)
        #
        # allocate residual and Jacobian storage
        #
        #
        # allocate residual and Jacobian storage
        #
        self.elementResidual = [numpy.zeros(
            (self.mesh.nElements_global,
             self.nDOF_test_element[ci]),
            'd')]
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global, i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray = numpy.zeros((self.nNodes_internal,), 'i')
        for nI, n in enumerate(self.internalNodes):
            self.internalNodesArray[nI] = n
        #
        del self.internalNodes
        self.internalNodes = None
        logEvent("Updating local to global mappings", 2)
        self.updateLocal2Global()
        logEvent("Building time integration object", 2)
        logEvent(memory("inflowBC, internalNodes,updateLocal2Global", "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        logEvent(memory("TimeIntegration", "OneLevelTransport"), level=4)
        logEvent("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        if numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions:
            interleave_DOF=True
            for nDOF_trial_element_ci in self.nDOF_trial_element:
                if nDOF_trial_element_ci != self.nDOF_trial_element[0]:
                    interleave_DOF=False
        else:
            interleave_DOF=False
        self.setupFieldStrides(interleave_DOF)
        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        logEvent("initalizing numerical flux")
        logEvent(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        logEvent("initializing numerical flux penalty")
        self.numericalFlux.penalty_constant = self.coefficients.eb_penalty_constant
        # cek todo move into numerical flux initialization
        if 'penalty' in self.ebq_global:
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = old_div(self.numericalFlux.penalty_constant, \
                        (self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power))
        # penalty term
        # cek move  to Numerical flux initialization
        if 'penalty' in self.ebqe:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = old_div(self.numericalFlux.penalty_constant, \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        logEvent(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        logEvent("setting up post-processing")
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)
        logEvent(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        logEvent("initializing archiver")
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        logEvent(memory("XdmfWriters", "OneLevelTransport"), level=4)
        logEvent("flux bc objects")
        for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                if ci in self.coefficients.advection:
                    self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
            for ck, diffusiveFluxBoundaryConditionsDict in list(fbcObject.diffusiveFluxBoundaryConditionsDictDict.items()):
                self.ebqe[('diffusiveFlux_bc_flag', ck, ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc', ck, ci)].shape, 'i')
                for t, g in list(diffusiveFluxBoundaryConditionsDict.items()):
                    self.ebqe[('diffusiveFlux_bc', ck, ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag', ck, ci)][t[0], t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN = 1.0
        else:
            self.MOVING_DOMAIN = 0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape, 'd')
        # cek/ido todo replace python loops in modules with optimized code if possible/necessary
        logEvent("dirichlet conditions")
        self.forceStrongConditions = coefficients.forceStrongDirichlet
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(
                    self.u[cj].femSpace, dofBoundaryConditionsSetterDict[cj], weakDirichletConditions=False)
        logEvent("final allocations")
        compKernelFlag = 0
        if self.coefficients.useConstant_he:
            self.elementDiameter = self.mesh.elementDiametersArray.copy()
            self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        else:
            self.elementDiameter = self.mesh.elementDiametersArray
        if self.nSpace_global == 2:
            import copy
            self.u[3] = self.u[2].copy()
            self.u[3].name="w"
            self.timeIntegration.m_tmp[3] = self.timeIntegration.m_tmp[2].copy()
            self.timeIntegration.beta_bdf[3] = self.timeIntegration.beta_bdf[2].copy()
            self.coefficients.sdInfo[(1, 3)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(2, 3)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(3, 0)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(3, 1)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(3, 2)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(3, 3)] = (numpy.array([0, 1, 2], dtype='i'),
                                                numpy.array([0, 1], dtype='i'))
            self.offset.append(self.offset[2])
            self.stride.append(self.stride[2])
            self.numericalFlux.isDOFBoundary[3] = self.numericalFlux.isDOFBoundary[2].copy()
            self.numericalFlux.ebqe[('u', 3)] = self.numericalFlux.ebqe[('u', 2)].copy()
            logEvent("calling cRANS2P2D_base ctor")
            self.rans2p = cRANS2P2D_base(self.nSpace_global,
                                         self.nQuadraturePoints_element,
                                         self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                         self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                         self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                         self.u[1].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                         self.testSpace[1].referenceFiniteElement.localFunctionSpace.dim,
                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                                         compKernelFlag)
        else:
            logEvent("calling  cRANS2P_base ctor")
            self.rans2p = cRANS2P_base(self.nSpace_global,
                                       self.nQuadraturePoints_element,
                                       self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                       self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                       self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                       self.u[1].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                       self.testSpace[1].referenceFiniteElement.localFunctionSpace.dim,
                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                       compKernelFlag)
        self.velocityErrorNodal = self.u[0].dof.copy()
        logEvent('WARNING: The boundary fluxes at interpart boundaries are skipped if elementBoundaryMaterialType is 0 for RANS2P-based models. This means that DG methods are currently incompatible with RANS2P.')

        
        self.q[('force', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('force', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('force', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')

        if options is not None:
            try:
                self.chrono_model = options.chrono_model
            except AttributeError:
                logEvent('WARNING: did not find chrono model')
                pass

    def getResidual(self, u, r):
        """
        Calculate the element residuals and add in to the global residual

        Parameters
        ----------
        u : :class:`numpy.ndarray`
        r : :class:`numpy.ndarray`
            Stores the calculated residual vector.
        """

        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        # cek todo put in logic to skip if BC's don't depend on t or u
        # hack
        if self.bcsTimeDependent or not self.bcsSet:
            self.bcsSet = True
            # Dirichlet boundary conditions
            self.numericalFlux.setDirichletValues(self.ebqe)
            # Flux boundary conditions
            for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
                for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                    if ci in self.coefficients.advection:
                        self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
                for ck, diffusiveFluxBoundaryConditionsDict in list(fbcObject.diffusiveFluxBoundaryConditionsDictDict.items()):
                    for t, g in list(diffusiveFluxBoundaryConditionsDict.items()):
                        self.ebqe[('diffusiveFlux_bc', ck, ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[('diffusiveFlux_bc_flag', ck, ci)][t[0], t[1]] = 1
        r.fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 36.0
        self.coefficients.wettedAreas[:] = 0.0
        self.coefficients.netForces_p[:, :] = 0.0
        self.coefficients.netForces_v[:, :] = 0.0
        self.coefficients.netMoments[:, :] = 0.0
        self.coefficients.particle_netForces[:, :] = 0.0
        self.coefficients.particle_netMoments[:, :] = 0.0
        self.coefficients.particle_surfaceArea[:] = 0.0
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    if cj == 0:
                        self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN], self.timeIntegration.t)
                    else:
                        self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],
                                                 self.timeIntegration.t)
                        if self.MOVING_DOMAIN == 1.0:
                            self.u[cj].dof[dofN] += self.mesh.nodeVelocityArray[dofN, cj - 1]
        self.rans2p.calculateResidual(self.coefficients.NONCONSERVATIVE_FORM,
                                      self.coefficients.MOMENTUM_SGE,
                                      self.coefficients.PRESSURE_SGE,
                                      self.coefficients.VELOCITY_SGE,
                                      self.coefficients.PRESSURE_PROJECTION_STABILIZATION,
                                      self.coefficients.numerical_viscosity,
                                      # element
                                      self.u[0].femSpace.elementMaps.psi,
                                      self.u[0].femSpace.elementMaps.grad_psi,
                                      self.mesh.nodeArray,
                                      self.mesh.nodeVelocityArray,
                                      self.MOVING_DOMAIN,
                                      self.mesh.elementNodesArray,
                                      self.elementQuadratureWeights[('u', 0)],
                                      self.u[0].femSpace.psi,
                                      self.u[0].femSpace.grad_psi,
                                      self.u[0].femSpace.psi,
                                      self.u[0].femSpace.grad_psi,
                                      self.u[1].femSpace.psi,
                                      self.u[1].femSpace.grad_psi,
                                      self.u[1].femSpace.psi,
                                      self.u[1].femSpace.grad_psi,
                                      # element boundary
                                      self.u[0].femSpace.elementMaps.psi_trace,
                                      self.u[0].femSpace.elementMaps.grad_psi_trace,
                                      self.elementBoundaryQuadratureWeights[('u', 0)],
                                      self.u[0].femSpace.psi_trace,
                                      self.u[0].femSpace.grad_psi_trace,
                                      self.u[0].femSpace.psi_trace,
                                      self.u[0].femSpace.grad_psi_trace,
                                      self.u[1].femSpace.psi_trace,
                                      self.u[1].femSpace.grad_psi_trace,
                                      self.u[1].femSpace.psi_trace,
                                      self.u[1].femSpace.grad_psi_trace,
                                      self.u[0].femSpace.elementMaps.boundaryNormals,
                                      self.u[0].femSpace.elementMaps.boundaryJacobians,
                                      # physics
                                      self.eb_adjoint_sigma,
                                      self.elementDiameter,  # mesh.elementDiametersArray,
                                      self.mesh.nodeDiametersArray,
                                      self.stabilization.hFactor,
                                      self.mesh.nElements_global,
                                      self.mesh.nElementBoundaries_owned,
                                      self.coefficients.useRBLES,
                                      self.coefficients.useMetrics,
                                      self.timeIntegration.alpha_bdf,
                                      self.coefficients.epsFact_density,
                                      self.coefficients.epsFact,
                                      self.coefficients.sigma,
                                      self.coefficients.rho_0,
                                      self.coefficients.nu_0,
                                      self.coefficients.rho_1,
                                      self.coefficients.nu_1,
                                      self.coefficients.smagorinskyConstant,
                                      self.coefficients.turbulenceClosureModel,
                                      self.Ct_sge,
                                      self.Cd_sge,
                                      self.shockCapturing.shockCapturingFactor,
                                      self.numericalFlux.penalty_constant,
                                      # VRANS start
                                      self.coefficients.epsFact_solid,
                                      self.coefficients.q_phi_solid,
                                      self.coefficients.q_velocity_solid,
                                      self.coefficients.q_porosity,
                                      self.coefficients.q_dragAlpha,
                                      self.coefficients.q_dragBeta,
                                      self.q[('r', 0)],
                                      self.coefficients.q_turb_var[0],
                                      self.coefficients.q_turb_var[1],
                                      self.coefficients.q_turb_var_grad[0],
                                      self.coefficients.LAG_LES,
                                      self.q['eddy_viscosity'],
                                      self.q['eddy_viscosity_last'],
                                      self.ebqe['eddy_viscosity'],
                                      self.ebqe['eddy_viscosity_last'],
                                      # VRANS end
                                      self.u[0].femSpace.dofMap.l2g,
                                      self.u[1].femSpace.dofMap.l2g,
                                      self.l2g[0]['freeGlobal'],
                                      self.l2g[1]['freeGlobal'],
                                      self.u[0].dof,
                                      self.u[1].dof,
                                      self.u[2].dof,
                                      self.u[3].dof,
                                      self.coefficients.p_old_dof,
                                      self.coefficients.u_old_dof,
                                      self.coefficients.v_old_dof,
                                      self.coefficients.w_old_dof,
                                      self.coefficients.g,
                                      self.coefficients.useVF,
                                      self.q['rho'],
                                      self.coefficients.q_vf,
                                      self.coefficients.q_phi,
                                      self.coefficients.q_n,
                                      self.coefficients.q_kappa,
                                      self.timeIntegration.m_tmp[1],
                                      self.timeIntegration.m_tmp[2],
                                      self.timeIntegration.m_tmp[3],
                                      self.q[('f', 0)],
                                      self.timeIntegration.beta_bdf[1],
                                      self.timeIntegration.beta_bdf[2],
                                      self.timeIntegration.beta_bdf[3],
                                      self.q['dV'],
                                      self.q['dV_last'],
                                      self.stabilization.v_last,
                                      self.q[('cfl', 0)],
                                      self.q[('numDiff', 1, 1)],
                                      self.q[('numDiff', 2, 2)],
                                      self.q[('numDiff', 3, 3)],
                                      self.shockCapturing.numDiff_last[1],
                                      self.shockCapturing.numDiff_last[2],
                                      self.shockCapturing.numDiff_last[3],
                                      self.coefficients.sdInfo[(1, 1)][0], self.coefficients.sdInfo[(1, 1)][1],
                                      self.coefficients.sdInfo[(1, 2)][0], self.coefficients.sdInfo[(1, 2)][1],
                                      self.coefficients.sdInfo[(1, 3)][0], self.coefficients.sdInfo[(1, 3)][1],
                                      self.coefficients.sdInfo[(2, 2)][0], self.coefficients.sdInfo[(2, 2)][1],
                                      self.coefficients.sdInfo[(2, 1)][0], self.coefficients.sdInfo[(2, 1)][1],
                                      self.coefficients.sdInfo[(2, 3)][0], self.coefficients.sdInfo[(2, 3)][1],
                                      self.coefficients.sdInfo[(3, 3)][0], self.coefficients.sdInfo[(3, 3)][1],
                                      self.coefficients.sdInfo[(3, 1)][0], self.coefficients.sdInfo[(3, 1)][1],
                                      self.coefficients.sdInfo[(3, 2)][0], self.coefficients.sdInfo[(3, 2)][1],
                                      self.offset[0], self.offset[1], self.offset[2], self.offset[3],
                                      self.stride[0], self.stride[1], self.stride[2], self.stride[3],
                                      r,
                                      self.mesh.nExteriorElementBoundaries_global,
                                      self.mesh.exteriorElementBoundariesArray,
                                      self.mesh.elementBoundaryElementsArray,
                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                      self.coefficients.ebqe_vf,
                                      self.coefficients.bc_ebqe_vf,
                                      self.coefficients.ebqe_phi,
                                      self.coefficients.bc_ebqe_phi,
                                      self.coefficients.ebqe_n,
                                      self.coefficients.ebqe_kappa,
                                      # VRANS start
                                      self.coefficients.ebqe_porosity,
                                      self.coefficients.ebqe_turb_var[0],
                                      self.coefficients.ebqe_turb_var[1],
                                      # VRANS end
                                      self.numericalFlux.isDOFBoundary[0],
                                      self.numericalFlux.isDOFBoundary[1],
                                      self.numericalFlux.isDOFBoundary[2],
                                      self.numericalFlux.isDOFBoundary[3],
                                      self.ebqe[('advectiveFlux_bc_flag', 0)],
                                      self.ebqe[('advectiveFlux_bc_flag', 1)],
                                      self.ebqe[('advectiveFlux_bc_flag', 2)],
                                      self.ebqe[('advectiveFlux_bc_flag', 3)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 3, 3)],
                                      self.numericalFlux.ebqe[('u', 0)],
                                      self.ebqe[('advectiveFlux_bc', 0)],
                                      self.ebqe[('advectiveFlux_bc', 1)],
                                      self.ebqe[('advectiveFlux_bc', 2)],
                                      self.ebqe[('advectiveFlux_bc', 3)],
                                      self.numericalFlux.ebqe[('u', 1)],
                                      self.ebqe[('diffusiveFlux_bc', 1, 1)],
                                      self.ebqe['penalty'],
                                      self.numericalFlux.ebqe[('u', 2)],
                                      self.ebqe[('diffusiveFlux_bc', 2, 2)],
                                      self.numericalFlux.ebqe[('u', 3)],
                                      self.ebqe[('diffusiveFlux_bc', 3, 3)],
                                      self.q['x'],
                                      self.q[('velocity', 0)],
                                      self.ebqe[('velocity', 0)],
                                      self.ebq_global[('totalFlux', 0)],
                                      self.elementResidual[0],
                                      self.mesh.elementMaterialTypes,
                                      self.mesh.elementBoundaryMaterialTypes,
                                      self.coefficients.barycenters,
                                      self.coefficients.wettedAreas,
                                      self.coefficients.netForces_p,
                                      self.coefficients.netForces_v,
                                      self.coefficients.netMoments,
                                      self.q['velocityError'],
                                      self.velocityErrorNodal,
                                      self.q[('force', 0)],
                                      self.q[('force', 1)],
                                      self.q[('force', 2)],
                                      self.coefficients.use_ball_as_particle,
                                      self.coefficients.ball_center,
                                      self.coefficients.ball_radius,
                                      self.coefficients.ball_velocity,
                                      self.coefficients.ball_angular_velocity,
                                      self.coefficients.ball_center_acceleration,
                                      self.coefficients.ball_angular_acceleration,
                                      self.coefficients.ball_density,
                                      self.coefficients.particle_signed_distances,
                                      self.coefficients.particle_signed_distance_normals,
                                      self.coefficients.particle_velocities,
                                      self.coefficients.particle_centroids,
                                      self.coefficients.ebq_global_phi_s,
                                      self.coefficients.ebq_global_grad_phi_s,
                                      self.coefficients.ebq_particle_velocity_s,
                                      self.coefficients.nParticles,
                                      self.coefficients.particle_netForces,
                                      self.coefficients.particle_netMoments,
                                      self.coefficients.particle_surfaceArea,
                                      self.mesh.nElements_owned,
                                      self.coefficients.particle_nitsche,
                                      self.coefficients.particle_epsFact,
                                      self.coefficients.particle_alpha,
                                      self.coefficients.particle_beta,
                                      self.coefficients.particle_penalty_constant,
                                      self.coefficients.phi_s,
                                      self.coefficients.phisField,
                                      self.coefficients.use_pseudo_penalty)
        for i in range(self.coefficients.netForces_p.shape[0]):
            self.coefficients.wettedAreas[i] = globalSum(self.coefficients.wettedAreas[i])
            for I in range(3):
                self.coefficients.netForces_p[i, I] = globalSum(self.coefficients.netForces_p[i, I])
                self.coefficients.netForces_v[i, I] = globalSum(self.coefficients.netForces_v[i, I])
                self.coefficients.netMoments[i, I] = globalSum(self.coefficients.netMoments[i, I])
                #cek hack, testing 6DOF motion
                #self.coefficients.netForces_p[i,I] = 0.0
                #self.coefficients.netForces_v[i,I] = 0.0
                #self.coefficients.netMoments[i,I] = 0.0
                #if I==0:
                #    self.coefficients.netForces_p[i,I] = (125.0* math.pi**2 * 0.125*math.cos(self.timeIntegration.t*math.pi))/4.0
                #if I==1:
                #    self.coefficients.netForces_p[i,I] = (125.0* math.pi**2 * 0.125*math.cos(self.timeIntegration.t*math.pi) + 125.0*9.81)/4.0
                #if I==2:
                #    self.coefficients.netMoments[i,I] = (4.05* math.pi**2 * (math.pi/4.0)*math.cos(self.timeIntegration.t*math.pi))/4.0
        
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        
        comm.Allreduce(self.coefficients.wettedAreas.copy(),self.coefficients.wettedAreas)
        comm.Allreduce(self.coefficients.netForces_p.copy(),self.coefficients.netForces_p)
        comm.Allreduce(self.coefficients.netForces_v.copy(),self.coefficients.netForces_v)
        comm.Allreduce(self.coefficients.netMoments.copy(),self.coefficients.netMoments)
        
        comm.Allreduce(self.coefficients.particle_netForces.copy(),self.coefficients.particle_netForces)
        comm.Allreduce(self.coefficients.particle_netMoments.copy(),self.coefficients.particle_netMoments)
        comm.Allreduce(self.coefficients.particle_surfaceArea.copy(),self.coefficients.particle_surfaceArea)
        
        for i in range(self.coefficients.nParticles):
            logEvent("particle i=" + repr(i)+ " force " + repr(self.coefficients.particle_netForces[i]))
            logEvent("particle i=" + repr(i)+ " moment " + repr(self.coefficients.particle_netMoments[i]))
            logEvent("particle i=" + repr(i)+ " surfaceArea " + repr(self.coefficients.particle_surfaceArea[i]))
            logEvent("particle i=" + repr(i)+ " stress force " + repr(self.coefficients.particle_netForces[i+self.coefficients.nParticles]))
            logEvent("particle i=" + repr(i)+ " pressure force " + repr(self.coefficients.particle_netForces[i+2*self.coefficients.nParticles]))

        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    if cj == 0:
                        r[self.offset[cj] + self.stride[cj] * dofN] = self.u[cj].dof[dofN] - \
                            g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN], self.timeIntegration.t)
                    else:
                        r[self.offset[cj] + self.stride[cj] * dofN] = self.u[cj].dof[dofN] - \
                            g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN], self.timeIntegration.t) 
                        if self.MOVING_DOMAIN == 1.0:
                            r[self.offset[cj] + self.stride[cj] * dofN] -= self.mesh.nodeVelocityArray[dofN, cj - 1]

        cflMax = globalMax(self.q[('cfl', 0)].max()) * self.timeIntegration.dt
        logEvent("Maximum CFL = " + str(cflMax), level=2)
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        logEvent("Global residual", level=9, data=r)
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        if self.nSpace_global == 2:
            self.csrRowIndeces[(0, 3)] = self.csrRowIndeces[(0, 2)]
            self.csrColumnOffsets[(0, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrRowIndeces[(1, 3)] = self.csrRowIndeces[(0, 2)]
            self.csrColumnOffsets[(1, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrRowIndeces[(2, 3)] = self.csrRowIndeces[(0, 2)]
            self.csrColumnOffsets[(2, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrRowIndeces[(3, 0)] = self.csrRowIndeces[(2, 0)]
            self.csrColumnOffsets[(3, 0)] = self.csrColumnOffsets[(2, 0)]
            self.csrRowIndeces[(3, 1)] = self.csrRowIndeces[(2, 0)]
            self.csrColumnOffsets[(3, 1)] = self.csrColumnOffsets[(2, 0)]
            self.csrRowIndeces[(3, 2)] = self.csrRowIndeces[(2, 0)]
            self.csrColumnOffsets[(3, 2)] = self.csrColumnOffsets[(2, 0)]
            self.csrRowIndeces[(3, 3)] = self.csrRowIndeces[(2, 0)]
            self.csrColumnOffsets[(3, 3)] = self.csrColumnOffsets[(2, 0)]
            self.csrColumnOffsets_eb[(0, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(1, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(2, 3)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(3, 0)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(3, 1)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(3, 2)] = self.csrColumnOffsets[(0, 2)]
            self.csrColumnOffsets_eb[(3, 3)] = self.csrColumnOffsets[(0, 2)]

        self.rans2p.calculateJacobian(self.coefficients.NONCONSERVATIVE_FORM,
                                      self.coefficients.MOMENTUM_SGE,
                                      self.coefficients.PRESSURE_SGE,
                                      self.coefficients.VELOCITY_SGE,
                                      self.coefficients.PRESSURE_PROJECTION_STABILIZATION,
                                      #element
                                      self.u[0].femSpace.elementMaps.psi,
                                      self.u[0].femSpace.elementMaps.grad_psi,
                                      self.mesh.nodeArray,
                                      self.mesh.nodeVelocityArray,
                                      self.MOVING_DOMAIN,
                                      self.mesh.elementNodesArray,
                                      self.elementQuadratureWeights[('u', 0)],
                                      self.u[0].femSpace.psi,
                                      self.u[0].femSpace.grad_psi,
                                      self.u[0].femSpace.psi,
                                      self.u[0].femSpace.grad_psi,
                                      self.u[1].femSpace.psi,
                                      self.u[1].femSpace.grad_psi,
                                      self.u[1].femSpace.psi,
                                      self.u[1].femSpace.grad_psi,
                                      # element boundary
                                      self.u[0].femSpace.elementMaps.psi_trace,
                                      self.u[0].femSpace.elementMaps.grad_psi_trace,
                                      self.elementBoundaryQuadratureWeights[('u', 0)],
                                      self.u[0].femSpace.psi_trace,
                                      self.u[0].femSpace.grad_psi_trace,
                                      self.u[0].femSpace.psi_trace,
                                      self.u[0].femSpace.grad_psi_trace,
                                      self.u[1].femSpace.psi_trace,
                                      self.u[1].femSpace.grad_psi_trace,
                                      self.u[1].femSpace.psi_trace,
                                      self.u[1].femSpace.grad_psi_trace,
                                      self.u[0].femSpace.elementMaps.boundaryNormals,
                                      self.u[0].femSpace.elementMaps.boundaryJacobians,
                                      self.eb_adjoint_sigma,
                                      self.elementDiameter,  # mesh.elementDiametersArray,
                                      self.mesh.nodeDiametersArray,
                                      self.stabilization.hFactor,
                                      self.mesh.nElements_global,
                                      self.coefficients.useRBLES,
                                      self.coefficients.useMetrics,
                                      self.timeIntegration.alpha_bdf,
                                      self.coefficients.epsFact_density,
                                      self.coefficients.epsFact,
                                      self.coefficients.sigma,
                                      self.coefficients.rho_0,
                                      self.coefficients.nu_0,
                                      self.coefficients.rho_1,
                                      self.coefficients.nu_1,
                                      self.coefficients.smagorinskyConstant,
                                      self.coefficients.turbulenceClosureModel,
                                      self.Ct_sge,
                                      self.Cd_sge,
                                      self.shockCapturing.shockCapturingFactor,
                                      self.numericalFlux.penalty_constant,
                                      # VRANS start
                                      self.coefficients.epsFact_solid,
                                      self.coefficients.q_phi_solid,
                                      self.coefficients.q_velocity_solid,
                                      self.coefficients.q_porosity,
                                      self.coefficients.q_dragAlpha,
                                      self.coefficients.q_dragBeta,
                                      self.q[('r', 0)],
                                      self.coefficients.q_turb_var[0],
                                      self.coefficients.q_turb_var[1],
                                      self.coefficients.q_turb_var_grad[0],
                                      self.coefficients.LAG_LES,
                                      self.q['eddy_viscosity_last'],
                                      self.ebqe['eddy_viscosity_last'],
                                      # VRANS end
                                      self.u[0].femSpace.dofMap.l2g,
                                      self.u[1].femSpace.dofMap.l2g,
                                      self.u[0].dof,
                                      self.u[1].dof,
                                      self.u[2].dof,
                                      self.u[3].dof,
                                      self.coefficients.p_old_dof,
                                      self.coefficients.u_old_dof,
                                      self.coefficients.v_old_dof,
                                      self.coefficients.w_old_dof,
                                      self.coefficients.g,
                                      self.coefficients.useVF,
                                      self.coefficients.q_vf,
                                      self.coefficients.q_phi,
                                      self.coefficients.q_n,
                                      self.coefficients.q_kappa,
                                      self.timeIntegration.beta_bdf[1],
                                      self.timeIntegration.beta_bdf[2],
                                      self.timeIntegration.beta_bdf[3],
                                      self.q['dV'],
                                      self.q['dV_last'],
                                      self.stabilization.v_last,
                                      self.q[('cfl', 0)],
                                      self.shockCapturing.numDiff_last[1],
                                      self.shockCapturing.numDiff_last[2],
                                      self.shockCapturing.numDiff_last[3],
                                      self.coefficients.sdInfo[(1, 1)][0], self.coefficients.sdInfo[(1, 1)][1],
                                      self.coefficients.sdInfo[(1, 2)][0], self.coefficients.sdInfo[(1, 2)][1],
                                      self.coefficients.sdInfo[(1, 3)][0], self.coefficients.sdInfo[(1, 3)][1],
                                      self.coefficients.sdInfo[(2, 2)][0], self.coefficients.sdInfo[(2, 2)][1],
                                      self.coefficients.sdInfo[(2, 1)][0], self.coefficients.sdInfo[(2, 1)][1],
                                      self.coefficients.sdInfo[(2, 3)][0], self.coefficients.sdInfo[(2, 3)][1],
                                      self.coefficients.sdInfo[(3, 3)][0], self.coefficients.sdInfo[(3, 3)][1],
                                      self.coefficients.sdInfo[(3, 1)][0], self.coefficients.sdInfo[(3, 1)][1],
                                      self.coefficients.sdInfo[(3, 2)][0], self.coefficients.sdInfo[(3, 2)][1],
                                      self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
                                      self.csrRowIndeces[(0, 1)], self.csrColumnOffsets[(0, 1)],
                                      self.csrRowIndeces[(0, 2)], self.csrColumnOffsets[(0, 2)],
                                      self.csrRowIndeces[(0, 3)], self.csrColumnOffsets[(0, 3)],
                                      self.csrRowIndeces[(1, 0)], self.csrColumnOffsets[(1, 0)],
                                      self.csrRowIndeces[(1, 1)], self.csrColumnOffsets[(1, 1)],
                                      self.csrRowIndeces[(1, 2)], self.csrColumnOffsets[(1, 2)],
                                      self.csrRowIndeces[(1, 3)], self.csrColumnOffsets[(1, 3)],
                                      self.csrRowIndeces[(2, 0)], self.csrColumnOffsets[(2, 0)],
                                      self.csrRowIndeces[(2, 1)], self.csrColumnOffsets[(2, 1)],
                                      self.csrRowIndeces[(2, 2)], self.csrColumnOffsets[(2, 2)],
                                      self.csrRowIndeces[(2, 3)], self.csrColumnOffsets[(2, 3)],
                                      self.csrRowIndeces[(3, 0)], self.csrColumnOffsets[(3, 0)],
                                      self.csrRowIndeces[(3, 1)], self.csrColumnOffsets[(3, 1)],
                                      self.csrRowIndeces[(3, 2)], self.csrColumnOffsets[(3, 2)],
                                      self.csrRowIndeces[(3, 3)], self.csrColumnOffsets[(3, 3)],
                                      jacobian,
                                      self.mesh.nExteriorElementBoundaries_global,
                                      self.mesh.exteriorElementBoundariesArray,
                                      self.mesh.elementBoundaryElementsArray,
                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                      self.coefficients.ebqe_vf,
                                      self.coefficients.bc_ebqe_vf,
                                      self.coefficients.ebqe_phi,
                                      self.coefficients.bc_ebqe_phi,
                                      self.coefficients.ebqe_n,
                                      self.coefficients.ebqe_kappa,
                                      # VRANS start
                                      self.coefficients.ebqe_porosity,
                                      self.coefficients.ebqe_turb_var[0],
                                      self.coefficients.ebqe_turb_var[1],
                                      # VRANS end
                                      self.numericalFlux.isDOFBoundary[0],
                                      self.numericalFlux.isDOFBoundary[1],
                                      self.numericalFlux.isDOFBoundary[2],
                                      self.numericalFlux.isDOFBoundary[3],
                                      self.ebqe[('advectiveFlux_bc_flag', 0)],
                                      self.ebqe[('advectiveFlux_bc_flag', 1)],
                                      self.ebqe[('advectiveFlux_bc_flag', 2)],
                                      self.ebqe[('advectiveFlux_bc_flag', 3)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
                                      self.ebqe[('diffusiveFlux_bc_flag', 3, 3)],
                                      self.numericalFlux.ebqe[('u', 0)],
                                      self.ebqe[('advectiveFlux_bc', 0)],
                                      self.ebqe[('advectiveFlux_bc', 1)],
                                      self.ebqe[('advectiveFlux_bc', 2)],
                                      self.ebqe[('advectiveFlux_bc', 3)],
                                      self.numericalFlux.ebqe[('u', 1)],
                                      self.ebqe[('diffusiveFlux_bc', 1, 1)],
                                      self.ebqe['penalty'],
                                      self.numericalFlux.ebqe[('u', 2)],
                                      self.ebqe[('diffusiveFlux_bc', 2, 2)],
                                      self.numericalFlux.ebqe[('u', 3)],
                                      self.ebqe[('diffusiveFlux_bc', 3, 3)],
                                      self.csrColumnOffsets_eb[(0, 0)],
                                      self.csrColumnOffsets_eb[(0, 1)],
                                      self.csrColumnOffsets_eb[(0, 2)],
                                      self.csrColumnOffsets_eb[(0, 3)],
                                      self.csrColumnOffsets_eb[(1, 0)],
                                      self.csrColumnOffsets_eb[(1, 1)],
                                      self.csrColumnOffsets_eb[(1, 2)],
                                      self.csrColumnOffsets_eb[(1, 3)],
                                      self.csrColumnOffsets_eb[(2, 0)],
                                      self.csrColumnOffsets_eb[(2, 1)],
                                      self.csrColumnOffsets_eb[(2, 2)],
                                      self.csrColumnOffsets_eb[(2, 3)],
                                      self.csrColumnOffsets_eb[(3, 0)],
                                      self.csrColumnOffsets_eb[(3, 1)],
                                      self.csrColumnOffsets_eb[(3, 2)],
                                      self.csrColumnOffsets_eb[(3, 3)],
                                      self.mesh.elementMaterialTypes,
                                      self.mesh.elementBoundaryMaterialTypes,
                                      self.coefficients.use_ball_as_particle,
                                      self.coefficients.ball_center,
                                      self.coefficients.ball_radius,
                                      self.coefficients.ball_velocity,
                                      self.coefficients.ball_angular_velocity,
                                      self.coefficients.ball_center_acceleration,
                                      self.coefficients.ball_angular_acceleration,
                                      self.coefficients.ball_density,
                                      self.coefficients.particle_signed_distances,
                                      self.coefficients.particle_signed_distance_normals,
                                      self.coefficients.particle_velocities,
                                      self.coefficients.particle_centroids,
                                      self.coefficients.ebq_global_phi_s,
                                      self.coefficients.ebq_global_grad_phi_s,
                                      self.coefficients.ebq_particle_velocity_s,
                                      self.coefficients.phi_s,
                                      self.coefficients.phisField,
                                      self.coefficients.nParticles,
                                      self.mesh.nElements_owned,
                                      self.coefficients.particle_nitsche,
                                      self.coefficients.particle_epsFact,
                                      self.coefficients.particle_alpha,
                                      self.coefficients.particle_beta,
                                      self.coefficients.particle_penalty_constant,
                                      self.coefficients.use_pseudo_penalty)
        
        if not self.forceStrongConditions and max(numpy.linalg.norm(self.u[1].dof, numpy.inf), numpy.linalg.norm(self.u[2].dof, numpy.inf), numpy.linalg.norm(self.u[3].dof, numpy.inf)) < 1.0e-8:
            self.pp_hasConstantNullSpace = True
        else:
            self.pp_hasConstantNullSpace = False
        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            for cj in range(self.nc):
                for dofN in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys()):
                    global_dofN = self.offset[cj] + self.stride[cj] * dofN
                    for i in range(self.rowptr[global_dofN], self.rowptr[global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            self.nzval[i] = 1.0
                        else:
                            self.nzval[i] = 0.0
                            # print "RBLES zeroing residual cj = %s dofN= %s global_dofN= %s " % (cj,dofN,global_dofN)
        logEvent("Jacobian ", level=10, data=jacobian)
        # mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def calculateElementQuadrature(self, domainMoved=False):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
                                                     self.q['x'])
            self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                             self.q['J'],
                                                             self.q['inverse(J)'],
                                                             self.q['det(J)'])
            self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, self.q[('v', 0)])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t, self.q)
        if self.stabilization is not None and not domainMoved:
            self.stabilization.initializeElementQuadrature(self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None and not domainMoved:
            self.shockCapturing.initializeElementQuadrature(self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self, domainMoved=False):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.u[0].femSpace.elementMaps.getValuesTrace(self.elementBoundaryQuadraturePoints,
                                                          self.ebq['x'])
            self.u[0].femSpace.elementMaps.getJacobianValuesTrace(self.elementBoundaryQuadraturePoints,
                                                                  self.ebq['inverse(J)'],
                                                                  self.ebq['g'],
                                                                  self.ebq['sqrt(det(g))'],
                                                                  self.ebq['n'])
            cfemIntegrals.copyLeftElementBoundaryInfo(self.mesh.elementBoundaryElementsArray,
                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                      self.mesh.exteriorElementBoundariesArray,
                                                      self.mesh.interiorElementBoundariesArray,
                                                      self.ebq['x'],
                                                      self.ebq['n'],
                                                      self.ebq_global['x'],
                                                      self.ebq_global['n'])
            self.u[0].femSpace.elementMaps.getInverseValuesTrace(self.ebq['inverse(J)'], self.ebq['x'], self.ebq['hat(x)'])
            self.u[0].femSpace.elementMaps.getPermutations(self.ebq['hat(x)'])
            self.testSpace[0].getBasisValuesTrace(self.u[0].femSpace.elementMaps.permutations,
                                                  self.ebq['hat(x)'],
                                                  self.ebq[('w', 0)])
            self.u[0].femSpace.getBasisValuesTrace(self.u[0].femSpace.elementMaps.permutations,
                                                   self.ebq['hat(x)'],
                                                   self.ebq[('v', 0)])
            cfemIntegrals.calculateElementBoundaryIntegrationWeights(self.ebq['sqrt(det(g))'],
                                                                     self.elementBoundaryQuadratureWeights[('u', 0)],
                                                                     self.ebq[('dS_u', 0)])

    def calculateExteriorElementBoundaryQuadrature(self, domainMoved=False):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        logEvent("initalizing ebqe vectors for post-procesing velocity")
        if self.postProcessing:
            self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                        self.ebqe['x'])
            self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                                self.ebqe['inverse(J)'],
                                                                                self.ebqe['g'],
                                                                                self.ebqe['sqrt(det(g))'],
                                                                                self.ebqe['n'])
            cfemIntegrals.calculateIntegrationWeights(self.ebqe['sqrt(det(g))'],
                                                      self.elementBoundaryQuadratureWeights[('u', 0)],
                                                      self.ebqe[('dS_u', 0)])
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        logEvent("initalizing basis info")
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        logEvent("setting flux boundary conditions")
        if not domainMoved:
            self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                       self.ebqe[('x')],
                                                                                       self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                       self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                           for cj in list(self.advectiveFluxBoundaryConditionsSetterDict.keys())])
            logEvent("initializing coefficients ebqe")
            self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t, self.ebqe)
        logEvent("done with ebqe")

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        if self.postProcessing and self.conservativeFlux:
            self.rans2p.calculateVelocityAverage(self.mesh.nExteriorElementBoundaries_global,
                                                 self.mesh.exteriorElementBoundariesArray,
                                                 self.mesh.nInteriorElementBoundaries_global,
                                                 self.mesh.interiorElementBoundariesArray,
                                                 self.mesh.elementBoundaryElementsArray,
                                                 self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                 self.mesh.nodeArray,
                                                 self.mesh.nodeVelocityArray,
                                                 self.MOVING_DOMAIN,
                                                 self.mesh.elementNodesArray,
                                                 self.u[0].femSpace.elementMaps.psi_trace,
                                                 self.u[0].femSpace.elementMaps.grad_psi_trace,
                                                 self.u[0].femSpace.elementMaps.boundaryNormals,
                                                 self.u[0].femSpace.elementMaps.boundaryJacobians,
                                                 self.u[1].femSpace.dofMap.l2g,
                                                 self.u[1].dof,
                                                 self.u[2].dof,
                                                 self.u[3].dof,
                                                 self.u[1].femSpace.psi_trace,
                                                 self.ebqe[('velocity', 0)],
                                                 self.ebq_global[('velocityAverage', 0)])
            if self.movingDomain:
                logEvent("Element Quadrature", level=3)
                self.calculateElementQuadrature(domainMoved=True)
                logEvent("Element Boundary Quadrature", level=3)
                self.calculateElementBoundaryQuadrature(domainMoved=True)
                logEvent("Global Exterior Element Boundary Quadrature", level=3)
                self.calculateExteriorElementBoundaryQuadrature(domainMoved=True)
                for ci in range(len(self.velocityPostProcessor.vpp_algorithms)):
                    for cj in list(self.velocityPostProcessor.vpp_algorithms[ci].updateConservationJacobian.keys()):
                        self.velocityPostProcessor.vpp_algorithms[ci].updateWeights()
                        self.velocityPostProcessor.vpp_algorithms[ci].computeGeometricInfo()
                        self.velocityPostProcessor.vpp_algorithms[ci].updateConservationJacobian[cj] = True
        self.q['velocityError'][:] = self.q[('velocity', 0)]
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)
        # if  self.coefficients.nd ==3:
        #     self.q[('cfl',0)][:] = np.sqrt(self.q[('velocity',0)][...,0]*self.q[('velocity',0)][...,0] +
        #                                    self.q[('velocity',0)][...,1]*self.q[('velocity',0)][...,1] +
        #                                    self.q[('velocity',0)][...,2]*self.q[('velocity',0)][...,2])/self.elementDiameter[:,np.newaxis]
        # else:
        #     self.q[('cfl',0)][:] = np.sqrt(self.q[('velocity',0)][...,0]*self.q[('velocity',0)][...,0] +
        #                                    self.q[('velocity',0)][...,1]*self.q[('velocity',0)][...,1])/self.elementDiameter[:,np.newaxis]
        self.q['velocityError'] -= self.q[('velocity', 0)]
        self.q['eddy_viscosity_last'][:] = self.q['eddy_viscosity']
        self.ebqe['eddy_viscosity_last'][:] = self.ebqe['eddy_viscosity']
        
    def updateAfterMeshMotion(self):
        pass


def getErgunDrag(porosity, meanGrainSize, viscosity):
    # cek hack, this doesn't seem right
    # cek todo look up correct Ergun model for alpha and beta
    voidFrac = 1.0 - porosity
    if voidFrac > 1.0e-6:
        dragBeta = porosity * porosity * porosity * meanGrainSize * 1.0e-2 / voidFrac
    if (porosity > epsZero and meanGrainSize > epsZero):
        dragAlpha = viscosity * 180.0 * voidFrac * voidFrac / (meanGrainSize * meanGrainSize * porosity)
