import proteus
from proteus import Profiling
import numpy
from proteus import *
from proteus.Transport import *
from proteus.Transport import OneLevelTransport
import os
from proteus import cfemIntegrals, Quadrature, Norms, Comm
from proteus.NonlinearSolvers import NonlinearEquation
from proteus.FemTools import (DOFBoundaryConditions,
                              FluxBoundaryConditions,
                              C0_AffineLinearOnSimplexWithNodalBasis)
from proteus.flcbdfWrappers import globalMax
from proteus.Profiling import memory
from proteus.Profiling import logEvent as log
from proteus.Transport import OneLevelTransport
from proteus.TransportCoefficients import TC_base
from proteus.SubgridError import SGE_base
from proteus.ShockCapturing import ShockCapturing_base
import cRANS3PF

class SubgridError(proteus.SubgridError.SGE_base):

    def __init__(
            self,
            coefficients,
            nd,
            lag=False,
            nStepsToDelay=0,
            hFactor=1.0,
            noPressureStabilization=False):
        self.noPressureStabilization = noPressureStabilization
        proteus.SubgridError.SGE_base.__init__(self, coefficients, nd, lag)
        coefficients.stencil[0].add(0)
        self.hFactor = hFactor
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            log("RANS3PF.SubgridError: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay = 1
            self.lag = False

    def initializeElementQuadrature(self, mesh, t, cq):
        import copy
        self.cq = cq
        self.v_last = self.cq[('velocity', 0)]


    def updateSubgridErrorHistory(self, initializationPhase=False):
        self.nSteps += 1
        if self.lag:
            self.v_last[:] = self.cq[('velocity', 0)]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            log("RANS3PF.SubgridError: switched to lagged subgrid error")
            self.lag = True
            self.v_last = self.cq[('velocity', 0)].copy()
    def calculateSubgridError(self, q):
        pass


class NumericalFlux(
        proteus.NumericalFlux.NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior):
    hasInterior = False

    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        proteus.NumericalFlux.NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior.__init__(
            self,
            vt,
            getPointwiseBoundaryConditions,
            getAdvectiveFluxBoundaryConditions,
            getDiffusiveFluxBoundaryConditions,
            getPeriodicBoundaryConditions)
        self.penalty_constant = 2.0
        self.includeBoundaryAdjoint = True
        self.boundaryAdjoint_sigma = 1.0
        self.hasInterior = False


class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):

    def __init__(
            self,
            coefficients,
            nd,
            shockCapturingFactor=0.25,
            lag=False,
            nStepsToDelay=3):
        proteus.ShockCapturing.ShockCapturing_base.__init__(
            self, coefficients, nd, shockCapturingFactor, lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            log("RANS3PF.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay = 1
            self.lag = False

    def initializeElementQuadrature(self, mesh, t, cq):
        self.mesh = mesh
        self.numDiff = {}
        self.numDiff_last = {}
        for ci in range(0, 3):
            self.numDiff[ci] = cq[('numDiff', ci, ci)]
            self.numDiff_last[ci] = cq[('numDiff', ci, ci)]

    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(0, 3):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            log("RANS3PF.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            for ci in range(0, 3):
                self.numDiff_last[ci] = self.numDiff[ci].copy()
        log(
            "RANS3PF: max numDiff_1 %e numDiff_2 %e numDiff_3 %e" %
            (globalMax(
                self.numDiff_last[0].max()), globalMax(
                self.numDiff_last[1].max()), globalMax(
                self.numDiff_last[2].max())))


class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for two incompresslble fluids governed by the Navier-Stokes equations and separated by a sharp interface represented by a level set function
    """
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_2D_Evaluate_sd
    from proteus.ctransportCoefficients import TwophaseNavierStokes_ST_LS_SO_3D_Evaluate_sd
    from proteus.ctransportCoefficients import calculateWaveFunction3d_ref

    def __init__(self,
                 CORRECT_VELOCITY=True,
                 cMax=1.0, # For entropy viscosity (mql) 
                 cE=1.0, # For entropy viscosity (mql)
                 epsFact=1.5,
                 sigma=72.8,
                 rho_0=998.2, nu_0=1.004e-6,
                 rho_1=1.205, nu_1=1.500e-5,
                 g=[0.0, 0.0, -9.8],
                 nd=3,
                 ME_model=5,
                 PRESSURE_model=7,
                 VOS_model=0,
                 SED_model=5,
                 LS_model=None,
                 VOF_model=None,
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
                 waveFlag=None,
                 waveHeight=0.01,
                 waveCelerity=1.0,
                 waveFrequency=1.0,
                 waveNumber=2.0,
                 waterDepth=0.5,
                 Omega_s=[[0.45, 0.55], [0.2, 0.4], [0.0, 1.0]],
                 epsFact_source=1.,
                 epsFact_solid=None,
                 eb_adjoint_sigma=1.0,
                 eb_penalty_constant=10.0,
                 forceStrongDirichlet=False,
                 turbulenceClosureModel=0,
                 # 0=No Model, 1=Smagorinksy, 2=Dynamic Smagorinsky,
                 # 3=K-Epsilon, 4=K-Omega
                 smagorinskyConstant=0.1,
                 barycenters=None,
                 PSTAB=0.0,
                 set_vos=None,
                 set_sed_velocity=None,
                 aDarcy=150.0,
		 betaForch=0.0,
		 grain=0.0102,
		 packFraction=0.2,
		 packMargin=0.01,
		 maxFraction=0.635,
		 frFraction=0.57,
                 sigmaC=1.1,
                 C3e=1.2,
                 C4e=1.0,
                 eR=0.8,
 		 fContact=0.02,
                 mContact=2.0,
                 nContact=5.0,
                 angFriction=pi/6.0,
                 nParticles=0,
                 particle_epsFact=3.0,
                 particle_alpha=1000.0,
                 particle_beta=1000.0,
                 particle_penalty_constant=1000.0,
                 particle_nitsche=1.0,
                 particle_sdfList=[],
                 particle_velocityList=[],
                 granular_sdf_Calc=None,
                 granular_vel_Calc=None
		 ):
        self.CORRECT_VELOCITY=CORRECT_VELOCITY
        self.nParticles=nParticles
        self.particle_nitsche=particle_nitsche
        self.particle_epsFact=particle_epsFact
        self.particle_alpha=particle_alpha
        self.particle_beta=particle_beta
        self.particle_penalty_constant=particle_penalty_constant
        self.particle_sdfList=particle_sdfList
        self.particle_velocityList=particle_velocityList
        self.granular_sdf_Calc=granular_sdf_Calc
        self.granular_vel_Calc=granular_vel_Calc
        self.aDarcy=aDarcy
        self.betaForch=betaForch
        self.grain=grain
        self.packFraction=packFraction
        self.packMargin=packMargin
        self.maxFraction=maxFraction
        self.frFraction=frFraction
        self.sigmaC=sigmaC
        self.C3e=C3e
        self.C4e=C4e
        self.eR=eR
        self.fContact=fContact
        self.mContact=mContact
        self.nContact=nContact
        self.angFriction=angFriction
        self.set_vos=set_vos
        self.set_sed=set_sed_velocity
        self.PSTAB=PSTAB
        self.barycenters = barycenters
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
        self.PRESSURE_model = PRESSURE_model
        self.VOS_model = VOS_model
        self.SED_model = SED_model
        self.LS_model = LS_model
        self.VOF_model = VOF_model
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
        # mql: for entropy viscosity 
        self.cMax=cMax
        self.cE=cE
        #
        self.dragAlpha = dragAlpha
        self.dragBeta = dragBeta
        self.setParamsFunc = setParamsFunc
        self.dragAlphaTypes = dragAlphaTypes
        self.dragBetaTypes = dragBetaTypes
        self.porosityTypes = porosityTypes
        self.killNonlinearDrag = int(killNonlinearDrag)
        self.waveFlag = waveFlag
        self.waveHeight = waveHeight
        self.waveCelerity = waveCelerity
        self.waveFrequency = waveFrequency
        self.waveNumber = waveNumber
        self.waterDepth = waterDepth
        self.Omega_s = Omega_s
        self.epsFact_source = epsFact_source
        self.linearDragFactor = 1.0
        self.nonlinearDragFactor = 1.0
        if self.killNonlinearDrag:
            self.nonlinearDragFactor = 0.0
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        if nd == 2:
            variableNames = ['u', 'v']
            mass = {0: {0: 'linear'},
                    1: {1: 'linear'}}
            advection = {0: {0: 'nonlinear',
                             1: 'nonlinear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear'}}
            diffusion = {0: {0: {0: 'constant'}, 1: {1: 'constant'}},
                         1: {1: {1: 'constant'}, 0: {0: 'constant'}}}
            sdInfo = {(0, 0): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (0, 1): (numpy.array([0, 0, 1], dtype='i'),
                               numpy.array([0], dtype='i')),
                      (1, 1): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (1, 0): (numpy.array([0, 1, 1], dtype='i'),
                               numpy.array([1], dtype='i'))}
            potential = {0: {0: 'u'},
                         1: {1: 'u'}}
            reaction = {0: {0: 'nonlinear', 1: 'nonlinear'},
                        1: {0: 'nonlinear', 1: 'nonlinear'}}
            hamiltonian = {0: {0: 'linear'},
                           1: {1: 'linear'}}
            TC_base.__init__(self,
                             2,
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
            self.vectorComponents = [0, 1]
            self.vectorName="velocity"
        elif nd == 3:
            variableNames = ['u', 'v', 'w']
            mass = {0: {0: 'linear'},
                    1: {1: 'linear'},
                    2: {2: 'linear'}}
            advection = {0: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'},
                         2: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'}}
            diffusion = {0: {0: {0: 'constant'},
                             1: {1: 'constant'},
                             2: {2: 'constant'}},
                         1: {0: {0: 'constant'},
                             1: {1: 'constant'},
                             2: {2: 'constant'}},
                         2: {0: {0: 'constant'},
                             1: {1: 'constant'},
                             2: {2: 'constant'}}}
            sdInfo = {}
            sdInfo = {(0, 0): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i')),
                      (0, 1): (numpy.array([0, 0, 1, 1], dtype='i'), numpy.array([0], dtype='i')),
                      (0, 2): (numpy.array([0, 0, 0, 1], dtype='i'), numpy.array([0], dtype='i')),
                      (1, 0): (numpy.array([0, 1, 1, 1], dtype='i'), numpy.array([1], dtype='i')),
                      (1, 1): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i')),
                      (1, 2): (numpy.array([0, 0, 0, 1], dtype='i'), numpy.array([1], dtype='i')),
                      (2, 0): (numpy.array([0, 1, 1, 1], dtype='i'), numpy.array([2], dtype='i')),
                      (2, 1): (numpy.array([0, 0, 1, 1], dtype='i'), numpy.array([2], dtype='i')),
                      (2, 2): (numpy.array([0, 1, 2, 3], dtype='i'), numpy.array([0, 1, 2], dtype='i'))}
            potential = {0: {0: 'u'},
                         1: {1: 'u'},
                         2: {2: 'u'}}
            reaction = {0: {0: 'nonlinear', 1: 'nonlinear', 2: 'nonlinear'},
                        1: {0: 'nonlinear', 1: 'nonlinear', 2: 'nonlinear'},
                        2: {0: 'nonlinear', 1: 'nonlinear', 2: 'nonlinear'}}
            hamiltonian = {0: {0: 'linear'},
                           1: {1: 'linear'},
                           2: {2: 'linear'}}
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
            self.vectorComponents = [0, 1, 2]
            self.vectorName="velocity"

    def attachModels(self, modelList):
        # level set
        self.model = modelList[self.ME_model]
        self.model.q['phi_solid'] = self.q_phi_solid
        self.q_rho  = self.model.q[('u',0)].copy()
        self.ebqe_rho  = self.model.ebqe[('u',0)].copy()
        self.q_nu  = self.model.q[('u',0)].copy()
        self.ebqe_nu  = self.model.ebqe[('u',0)].copy()
        #DEM particles
        self.particle_netForces = np.zeros((self.nParticles,3),'d')
        self.particle_netMoments = np.zeros((self.nParticles,3),'d')
        self.particle_surfaceArea = np.zeros((self.nParticles,),'d')
        self.particle_centroids = np.zeros((self.nParticles,3),'d')
        self.particle_signed_distances=np.zeros((self.nParticles,)+self.model.q[('u',0)].shape,'d')
        self.particle_signed_distance_normals=np.zeros((self.nParticles,)+self.model.q[('velocity',0)].shape,'d')
        self.particle_velocities=np.zeros((self.nParticles,)+self.model.q[('velocity',0)].shape,'d')

        self.phisField=np.ones(self.model.q[('u',0)].shape,'d')*1e10
        # This is making a special case for granular material simulations
        # if the user inputs a list of position/velocities then the sdf are calculated based on the "spherical" particles
        # otherwise the sdf are calculated based on the input sdf list for each body
        if self.granular_sdf_Calc is not None:
            temp_1=np.zeros(self.model.q[('u',0)].shape,'d')
            temp_2=np.zeros(self.model.q[('u',0)].shape,'d')
            temp_3=np.zeros(self.model.q[('u',0)].shape,'d')
            for i in range(self.nParticles):
                print ("Attaching particle i=", i)
                for eN in range(self.model.q['x'].shape[0]):
                    for k in range(self.model.q['x'].shape[1]):
                        self.particle_signed_distances[i,eN,k],self.particle_signed_distance_normals[i,eN,k] = self.granular_sdf_Calc(self.model.q['x'][eN,k],i)    
                        self.particle_velocities[i,eN,k] = self.granular_vel_Calc(self.model.q['x'][eN,k],i)           
                #This is important to write the a field for the sdf in the domain which distinguishes solid particles
                temp_1=np.minimum(abs(self.particle_signed_distances[i]),abs(self.phisField))
                temp_2=np.minimum(self.particle_signed_distances[i],temp_1)
                temp_3=np.minimum(temp_2,self.phisField)
                self.phisField=temp_3
                self.model.q[('phis')] = temp_3
        else:
            for i,sdf,vel in zip(range(self.nParticles),
                            self.particle_sdfList, self.particle_velocityList):
                for eN in range(self.model.q['x'].shape[0]):
                    for k in range(self.model.q['x'].shape[1]):
                        self.particle_signed_distances[i,eN,k],self.particle_signed_distance_normals[i,eN,k] = sdf(0, self.model.q['x'][eN,k])
                        self.particle_velocities[i,eN,k]=vel(0,self.model.q['x'][eN,k])
                self.model.q[('phis',i)] = self.particle_signed_distances[i]
                self.model.q[('phis_vel',i)] = self.particle_velocities[i]

	
        if self.PRESSURE_model is not None:
            self.model.pressureModel = modelList[self.PRESSURE_model]
            self.model.q_p_fluid = modelList[self.PRESSURE_model].q[('u',0)]
            self.model.ebqe_p_fluid = modelList[self.PRESSURE_model].ebqe[('u',0)]
            self.model.q_grad_p_fluid = modelList[self.PRESSURE_model].q[('grad(u)',0)]
            self.model.ebqe_grad_p_fluid = modelList[self.PRESSURE_model].ebqe[('grad(u)',0)]
        if self.VOS_model is not None:
            self.vos_dof = modelList[self.VOS_model].u[0].dof
            self.q_vos = modelList[self.VOS_model].q[('u',0)]
            self.q_dvos_dt = modelList[self.VOS_model].q[('mt',0)]
            self.ebqe_vos = modelList[self.VOS_model].ebqe[('u',0)]
        else:
            if self.VOF_model is None:
                self.vos_dof = modelList[self.ME_model].u[0].dof.copy()
                self.vos_dof[:] = 0.0
                self.q_vos = modelList[self.ME_model].q[('u',0)].copy()
                self.q_vos[:] = 0.0
                self.q_dvos_dt = self.q_vos.copy()
                self.q_dvos_dt[:] = 0.0
                self.ebqe_vos = modelList[self.ME_model].ebqe[('u',0)].copy() 
                self.ebqe_vos[:] = 0.0
            else:
                self.vos_dof = modelList[self.VOF_model].u[0].dof.copy()
                self.vos_dof[:] = 0.0
                self.q_vos = modelList[self.VOF_model].coefficients.q_vos
                self.q_dvos_dt = self.q_vos.copy()
                self.q_dvos_dt[:] = 0.0
                self.ebqe_vos = modelList[self.VOF_model].coefficients.ebqe_vos
        if self.SED_model is not None:
            self.rho_s = modelList[self.SED_model].coefficients.rho_0
            self.q_velocity_solid = modelList[self.SED_model].q[('velocity',0)]
            self.ebqe_velocity_solid = modelList[self.SED_model].ebqe[('velocity',0)]
        else:
            self.rho_s = self.rho_0
            self.q_velocity_solid = self.model.q[('velocity',0)].copy()
            self.q_velocity_solid[:] = 0.0
            self.ebqe_velocity_solid = self.model.ebqe[('velocity',0)].copy()
            self.ebqe_velocity_solid[:] = 0.0
        if self.LS_model is not None:
            self.q_phi = modelList[self.LS_model].q[('u', 0)]
            if modelList[self.LS_model].ebq.has_key(('u', 0)):
                self.ebq_phi = modelList[self.LS_model].ebq[('u', 0)]
            else:
                self.ebq_phi = None
            self.ebqe_phi = modelList[self.LS_model].ebqe[('u', 0)]
            self.bc_ebqe_phi = modelList[
                self.LS_model].numericalFlux.ebqe[
                ('u', 0)]
            # normal
            self.q_n = modelList[self.LS_model].q[('grad(u)', 0)]
            if modelList[self.LS_model].ebq.has_key(('grad(u)', 0)):
                self.ebq_n = modelList[self.LS_model].ebq[('grad(u)', 0)]
            else:
                self.ebq_n = None
            self.ebqe_n = modelList[self.LS_model].ebqe[('grad(u)', 0)]
        else:
            self.q_phi = 10.0 * numpy.ones(self.model.q[('u', 0)].shape, 'd')
            self.ebqe_phi = 10.0 * \
                numpy.ones(self.model.ebqe[('u', 0)].shape, 'd')
            self.bc_ebqe_phi = 10.0 * \
                numpy.ones(self.model.ebqe[('u', 0)].shape, 'd')
            self.q_n = numpy.ones(self.model.q[('velocity', 0)].shape, 'd')
            self.ebqe_n = numpy.ones(
                self.model.ebqe[
                    ('velocity', 0)].shape, 'd')
        if self.VOF_model is not None:
            self.q_vf = modelList[self.VOF_model].q[('u', 0)]
            if modelList[self.VOF_model].ebq.has_key(('u', 0)):
                self.ebq_vf = modelList[self.VOF_model].ebq[('u', 0)]
            else:
                self.ebq_vf = None
            self.ebqe_vf = modelList[self.VOF_model].ebqe[('u', 0)]
            self.bc_ebqe_vf = modelList[
                self.VOF_model].numericalFlux.ebqe[
                ('u', 0)]
        else:
            self.q_vf = numpy.zeros(self.model.q[('u', 0)].shape, 'd')
            self.ebqe_vf = numpy.zeros(self.model.ebqe[('u', 0)].shape, 'd')
            self.bc_ebqe_vf = numpy.zeros(self.model.ebqe[('u', 0)].shape, 'd')
        # curvature
        if self.KN_model is not None:
            self.q_kappa = modelList[self.KN_model].q[('u', 0)]
            self.ebqe_kappa = modelList[self.KN_model].ebqe[('u', 0)]
            if modelList[self.KN_model].ebq.has_key(('u', 0)):
                self.ebq_kappa = modelList[self.KN_model].ebq[('u', 0)]
            else:
                self.ebq_kappa = None
        else:
            self.q_kappa = -numpy.ones(self.model.q[('u', 0)].shape, 'd')
            self.ebqe_kappa = -numpy.ones(self.model.ebqe[('u', 0)].shape, 'd')
        # Turbulence Closures
        # only option for now is k-epsilon
        self.q_turb_var = {}
        self.q_turb_var_grad = {}
        self.ebqe_turb_var = {}
        if self.Closure_0_model is not None:
            self.q_turb_var[0] = modelList[self.Closure_0_model].q[('u', 0)]
            self.q_turb_var_grad[0] = modelList[
                self.Closure_0_model].q[
                ('grad(u)', 0)]
            self.ebqe_turb_var[0] = modelList[
                self.Closure_0_model].ebqe[
                ('u', 0)]
        else:
            self.q_turb_var[0] = numpy.ones(self.model.q[('u', 0)].shape, 'd')
            self.q_turb_var_grad[0] = numpy.ones(
                self.model.q[('grad(u)', 0)].shape, 'd')
            self.ebqe_turb_var[0] = numpy.ones(
                self.model.ebqe[('u', 0)].shape, 'd')
        if self.Closure_1_model is not None:
            self.q_turb_var[1] = modelList[self.Closure_1_model].q[('u', 0)]
            self.ebqe_turb_var[1] = modelList[
                self.Closure_1_model].ebqe[
                ('u', 0)]
        else:
            self.q_turb_var[1] = numpy.ones(self.model.q[('u', 0)].shape, 'd')
            self.ebqe_turb_var[1] = numpy.ones(
                self.model.ebqe[('u', 0)].shape, 'd')
        if self.epsFact_solid is None:
            self.epsFact_solid = numpy.ones(
                self.model.mesh.elementMaterialTypes.max() + 1)
        assert len(self.epsFact_solid) > self.model.mesh.elementMaterialTypes.max(
        ), "epsFact_solid  array is not large  enough for the materials  in this mesh; length must be greater  than largest  material type ID"

    def initializeMesh(self, mesh):
        self.phi_s = numpy.zeros(mesh.nodeArray.shape[0],'d')

        if self.granular_sdf_Calc is not None:
            print ("updating",self.nParticles," particles...")
            for i in range(self.nParticles):
                for j in range(mesh.nodeArray.shape[0]):
                    sdf,sdNormals = self.granular_sdf_Calc(mesh.nodeArray[j,:],i)
                    if ( abs(sdf) < abs(self.phi_s[j]) ):
                         self.phi_s[j]=sdf
        else:
            for i,sdf in zip(range(self.nParticles),
                            self.particle_sdfList):
                for j in range(mesh.nodeArray.shape[0]):
                     self.phi_s[j],sdNormals=sdf(0,mesh.nodeArray[j,:])


        # cek we eventually need to use the local element diameter
        self.eps_density = self.epsFact_density * mesh.h
        self.eps_viscosity = self.epsFact * mesh.h
        self.mesh = mesh
        self.elementMaterialTypes = mesh.elementMaterialTypes
        self.eps_source = self.epsFact_source * mesh.h
        nBoundariesMax = int(
            globalMax(max(self.mesh.elementBoundaryMaterialTypes))) + 1
        self.wettedAreas = numpy.zeros((nBoundariesMax,), 'd')
        self.netForces_p = numpy.zeros((nBoundariesMax, 3), 'd')
        self.netForces_v = numpy.zeros((nBoundariesMax, 3), 'd')
        self.netMoments = numpy.zeros((nBoundariesMax, 3), 'd')
        if self.barycenters is None:
            self.barycenters = numpy.zeros((nBoundariesMax, 3), 'd')
        comm = Comm.get()
        if comm.isMaster():
            self.wettedAreaHistory = open(os.path.join(proteus.Profiling.logDir,"wettedAreaHistory.txt"),"w")
            self.forceHistory_p = open(os.path.join(proteus.Profiling.logDir,"forceHistory_p.txt"),"w")
            self.forceHistory_v = open(os.path.join(proteus.Profiling.logDir,"forceHistory_v.txt"),"w")
            self.momentHistory = open(os.path.join(proteus.Profiling.logDir,"momentHistory.txt"),"w")
            if self.nParticles:
                self.particle_forceHistory = open(os.path.join(proteus.Profiling.logDir,"particle_forceHistory.txt"),"w")
                self.particle_momentHistory = open(os.path.join(proteus.Profiling.logDir,"particle_momentHistory.txt"),"w")
    # initialize so it can run as single phase

    def initializeElementQuadrature(self, t, cq):
        # VRANS
        self.q_phi_solid = numpy.ones(cq[('u', 0)].shape, 'd')
        self.q_velocity_solid = numpy.zeros(cq[('velocity', 0)].shape, 'd')
        self.q_porosity = numpy.ones(cq[('u', 0)].shape, 'd')
        self.q_dragAlpha = numpy.ones(cq[('u', 0)].shape, 'd')
        self.q_dragAlpha.fill(self.dragAlpha)
        self.q_dragBeta = numpy.ones(cq[('u', 0)].shape, 'd')
        self.q_dragBeta.fill(self.dragBeta)
        if self.setParamsFunc is not None:
            self.setParamsFunc(
                cq['x'],
                self.q_porosity,
                self.q_dragAlpha,
                self.q_dragBeta)
        else:
            # TODO make loops faster
            if self.porosityTypes is not None:
                for eN in range(self.q_porosity.shape[0]):
                    self.q_porosity[
                        eN, :] = self.porosityTypes[
                        self.elementMaterialTypes[eN]]
            if self.dragAlphaTypes is not None:
                for eN in range(self.q_dragAlpha.shape[0]):
                    self.q_dragAlpha[
                        eN, :] = self.dragAlphaTypes[
                        self.elementMaterialTypes[eN]]
            if self.dragBetaTypes is not None:
                for eN in range(self.q_dragBeta.shape[0]):
                    self.q_dragBeta[
                        eN, :] = self.dragBetaTypes[
                        self.elementMaterialTypes[eN]]
        #
        cq['phisError'] = cq[('u',0)].copy()
        cq['velocityError'] = cq[('velocity',0)].copy()

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        # VRANS
        self.ebq_porosity = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragAlpha = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragAlpha.fill(self.dragAlpha)
        self.ebq_dragBeta = numpy.ones(cebq['det(J)'].shape, 'd')
        self.ebq_dragBeta.fill(self.dragBeta)
        if self.setParamsFunc is not None:
            self.setParamsFunc(
                cebq['x'],
                self.ebq_porosity,
                self.ebq_dragAlpha,
                self.ebq_dragBeta)
        # TODO which mean to use or leave discontinuous
        # TODO make loops faster
        if self.porosityTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
                ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[
                    ebN, 0]
                ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[
                    ebN, 1]
                avg = 0.5 * (self.porosityTypes[self.elementMaterialTypes[eN_left]] +
                             self.porosityTypes[self.elementMaterialTypes[eN_right]])
                self.ebq_porosity[
                    eN_left, ebN_element_left, :] = self.porosityTypes[
                    self.elementMaterialTypes[eN_left]]
                self.ebq_porosity[
                    eN_right, ebN_element_right, :] = self.porosityTypes[
                    self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                    ebN, 0]
                self.ebq_porosity[
                    eN, ebN_element, :] = self.porosityTypes[
                    self.elementMaterialTypes[eN]]
        if self.dragAlphaTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
            ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
            ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 1]
            avg = 0.5 * (self.dragAlphaTypes[self.elementMaterialTypes[eN_left]] +
                         self.dragAlphaTypes[self.elementMaterialTypes[eN_right]])
            self.ebq_dragAlpha[
                eN_left, ebN_element_left, :] = self.dragAlphaTypes[
                self.elementMaterialTypes[eN_left]]
            self.ebq_dragAlpha[
                eN_right, ebN_element_right, :] = self.dragAlphaTypes[
                self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                    ebN, 0]
                self.ebq_dragAlpha[
                    eN, ebN_element, :] = self.dragAlphaTypes[
                    self.elementMaterialTypes[eN]]
        if self.dragBetaTypes is not None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                eN_left = self.mesh.elementBoundaryElementsArray[ebN, 0]
                eN_right = self.mesh.elementBoundaryElementsArray[ebN, 1]
            ebN_element_left = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
            ebN_element_right = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 1]
            avg = 0.5 * (self.dragBetaTypes[self.elementMaterialTypes[eN_left]] +
                         self.dragBetaTypes[self.elementMaterialTypes[eN_right]])
            self.ebq_dragBeta[
                eN_left, ebN_element_left, :] = self.dragBetaTypes[
                self.elementMaterialTypes[eN_left]]
            self.ebq_dragBeta[
                eN_right, ebN_element_right, :] = self.dragBetaTypes[
                self.elementMaterialTypes[eN_right]]
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                    ebN, 0]
                self.ebq_dragBeta[
                    eN, ebN_element, :] = self.dragBetaTypes[
                    self.elementMaterialTypes[eN]]
         #

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        # VRANS
        log("ebqe_global allocations in coefficients")
        self.ebqe_velocity_last = numpy.zeros(cebqe[('velocity',0)].shape)
        self.ebqe_porosity = numpy.ones(cebqe[('u', 0)].shape, 'd')
        self.ebqe_dragAlpha = numpy.ones(cebqe[('u', 0)].shape, 'd')
        self.ebqe_dragAlpha.fill(self.dragAlpha)
        self.ebqe_dragBeta = numpy.ones(cebqe[('u', 0)].shape, 'd')
        self.ebqe_dragBeta.fill(self.dragBeta)
        log("porosity and drag")
        # TODO make loops faster
        if self.setParamsFunc is not None:
            self.setParamsFunc(
                cebqe['x'],
                self.ebqe_porosity,
                self.ebqe_dragAlpha,
                self.ebqe_dragBeta)
        else:
            if self.porosityTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_porosity[
                        ebNE, :] = self.porosityTypes[
                        self.elementMaterialTypes[eN]]
            if self.dragAlphaTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_dragAlpha[
                        ebNE, :] = self.dragAlphaTypes[
                        self.elementMaterialTypes[eN]]
            if self.dragBetaTypes is not None:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN = self.mesh.elementBoundaryElementsArray[ebN, 0]
                    self.ebqe_dragBeta[
                        ebNE, :] = self.dragBetaTypes[
                        self.elementMaterialTypes[eN]]
        #

    def updateToMovingDomain(self, t, c):
        pass

    def evaluateForcingTerms(
            self,
            t,
            c,
            mesh=None,
            mesh_trial_ref=None,
            mesh_l2g=None):
        if c.has_key('x') and len(c['x'].shape) == 3:
            if self.nd == 2:
                # mwf debug
                #import pdb
                # pdb.set_trace()
                c[('r', 0)].fill(0.0)
                eps_source = self.eps_source
                if self.waveFlag == 1:  # secondOrderStokes:
                    waveFunctions.secondOrderStokesWave(c[('r', 0)].shape[0],
                                                        c[('r', 0)].shape[1],
                                                        self.waveHeight,
                                                        self.waveCelerity,
                                                        self.waveFrequency,
                                                        self.waveNumber,
                                                        self.waterDepth,
                                                        self.Omega_s[0][0],
                                                        self.Omega_s[0][1],
                                                        self.Omega_s[1][0],
                                                        self.Omega_s[1][1],
                                                        eps_source,
                                                        c['x'],
                                                        c[('r', 0)],
                                                        t)
                elif self.waveFlag == 2:  # solitary wave
                    waveFunctions.solitaryWave(c[('r', 0)].shape[0],
                                               c[('r', 0)].shape[1],
                                               self.waveHeight,
                                               self.waveCelerity,
                                               self.waveFrequency,
                                               self.waterDepth,
                                               self.Omega_s[0][0],
                                               self.Omega_s[0][1],
                                               self.Omega_s[1][0],
                                               self.Omega_s[1][1],
                                               eps_source,
                                               c['x'],
                                               c[('r', 0)],
                                               t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave(c[('r', 0)].shape[0],
                                                    c[('r', 0)].shape[1],
                                                    self.waveHeight,
                                                    self.waveCelerity,
                                                    self.waveFrequency,
                                                    self.Omega_s[0][0],
                                                    self.Omega_s[0][1],
                                                    self.Omega_s[1][0],
                                                    self.Omega_s[1][1],
                                                    eps_source,
                                                    c['x'],
                                                    c[('r', 0)],
                                                    t)

                # mwf debug
                if numpy.isnan(c[('r', 0)].any()):
                    import pdb
                    pdb.set_trace()
            else:
                # mwf debug
                #import pdb
                # pdb.set_trace()
                c[('r', 0)].fill(0.0)
                eps_source = self.eps_source
                if self.waveFlag == 1:  # secondOrderStokes:
                    waveFunctions.secondOrderStokesWave3d(c[('r', 0)].shape[0],
                                                          c[('r', 0)].shape[1],
                                                          self.waveHeight,
                                                          self.waveCelerity,
                                                          self.waveFrequency,
                                                          self.waveNumber,
                                                          self.waterDepth,
                                                          self.Omega_s[0][0],
                                                          self.Omega_s[0][1],
                                                          self.Omega_s[1][0],
                                                          self.Omega_s[1][1],
                                                          self.Omega_s[2][0],
                                                          self.Omega_s[2][1],
                                                          eps_source,
                                                          c['x'],
                                                          c[('r', 0)],
                                                          t)
                elif self.waveFlag == 2:  # solitary wave
                    waveFunctions.solitaryWave3d(c[('r', 0)].shape[0],
                                                 c[('r', 0)].shape[1],
                                                 self.waveHeight,
                                                 self.waveCelerity,
                                                 self.waveFrequency,
                                                 self.waterDepth,
                                                 self.Omega_s[0][0],
                                                 self.Omega_s[0][1],
                                                 self.Omega_s[1][0],
                                                 self.Omega_s[1][1],
                                                 self.Omega_s[2][0],
                                                 self.Omega_s[2][1],
                                                 eps_source,
                                                 c['x'],
                                                 c[('r', 0)],
                                                 t)

                elif self.waveFlag == 0:
                    waveFunctions.monochromaticWave3d(c[('r', 0)].shape[0],
                                                      c[('r', 0)].shape[1],
                                                      self.waveHeight,
                                                      self.waveCelerity,
                                                      self.waveFrequency,
                                                      self.Omega_s[0][0],
                                                      self.Omega_s[0][1],
                                                      self.Omega_s[1][0],
                                                      self.Omega_s[1][1],
                                                      self.Omega_s[2][0],
                                                      self.Omega_s[2][1],
                                                      eps_source,
                                                      c['x'],
                                                      c[('r', 0)],
                                                      t)

        else:
            assert mesh is not None
            assert mesh_trial_ref is not None
            assert mesh_l2g is not None
            # cek hack
            pass
        #            self.calculateWaveFunction3d_ref(mesh_trial_ref,
        #                                     mesh.nodeArray,
        #                                     mesh_l2g,
        #                                     mesh.elementDiametersArray,
        #                                     numpy.array(self.Omega_s[0]),
        #                                     numpy.array(self.Omega_s[1]),
        #                                     numpy.array(self.Omega_s[2]),
        #                                     t,
        #                                     self.waveFlag,
        #                                     self.epsFact_source,
        #                                     self.waveHeight,
        #                                     self.waveCelerity,
        #                                     self.waveFrequency,
        #                                     self.waveNumber,
        #                                     self.waterDepth,
        #                                     c[('r',0)])

    def evaluate(self, t, c):
        pass

    def preStep(self, t, firstStep=False):
        # Save old solutions
        # solution at tnm1        
        self.model.u_dof_old_old[:] = self.model.u_dof_old[:]
        self.model.v_dof_old_old[:] = self.model.v_dof_old[:]
        if (self.model.nSpace_global == 3):
            self.model.w_dof_old_old[:] = self.model.w_dof_old[:]
        # solution at tn
        self.model.u_dof_old[:] = self.model.u[0].dof
        self.model.v_dof_old[:] = self.model.u[1].dof
        if (self.model.nSpace_global == 3):
            self.model.w_dof_old[:] = self.model.u[2].dof
        # Compute 2nd order extrapolation on velocity
        if (firstStep):
            self.model.q[('velocityStar',0)][:] = self.model.q[('velocity',0)]        
        else:
            if self.model.timeIntegration.timeOrder == 1: 
                r = 1.
            else:
                r = self.model.timeIntegration.dt/self.model.timeIntegration.dt_history[0] 
            self.model.q[('velocityStar',0)][:] = (1+r)*self.model.q[('velocity',0)] - r*self.model.q[('velocityOld',0)]
        self.model.q[('velocityOld',0)][:] = self.model.q[('velocity',0)]  

        # COMPUTE MATERIAL PARAMETERS AND FORCE TERMS AS FUNCTIONS (if provided)
        if self.model.hasMaterialParametersAsFunctions:
            self.model.updateMaterialParameters()
        if self.model.hasForceTermsAsFunctions: 
            self.model.updateForceTerms()

        self.model.dt_last = self.model.timeIntegration.dt
        pass

    def postStep(self, t, firstStep=False):
        if firstStep==True:
            self.model.firstStep=False
        self.model.dt_last = self.model.timeIntegration.dt
        self.model.q['dV_last'][:] = self.model.q['dV']

        self.phi_s[:]=1e10
        self.phisField=np.ones(self.model.q[('u',0)].shape,'d')*1e10
        if self.granular_sdf_Calc is not None:
            print ("updating",self.nParticles," particles...")
            for i in range(self.nParticles):
                for j in range(self.mesh.nodeArray.shape[0]):
                    vel = self.granular_vel_Calc(self.mesh.nodeArray[j,:],i)
                    sdf,sdNormals = self.granular_sdf_Calc(self.mesh.nodeArray[j,:],i)
                    sdf+=vel* self.model.dt_last
                    if ( abs(sdf) < abs(self.phi_s[j]) ):
                        self.phi_s[j]=sdf
                for eN in range(self.model.q['x'].shape[0]):
                    for k in range(self.model.q['x'].shape[1]):
                        self.particle_signed_distances[i,eN,k],self.particle_signed_distance_normals[i,eN,k] = self.granular_sdf_Calc(self.model.q['x'][eN,k],i)
                        self.particle_velocities[i,eN,k] = self.granular_vel_Calc(self.model.q['x'][eN,k],i)
                        if ( abs(self.particle_signed_distances[i,eN,k]) < abs(self.phisField[eN,k]) ):
                            self.phisField[eN,k]=self.particle_signed_distances[i,eN,k]

            self.model.q[('phis')] = self.phisField    

        else:
            for i,sdf,vel in zip(range(self.nParticles),
                            self.particle_sdfList,self.particle_velocityList
                            ):
                for j in range(self.mesh.nodeArray.shape[0]):
                    myvel = vel(t,self.mesh.nodeArray[j,:])
                    mysdf,sdNormals=sdf(t,self.mesh.nodeArray[j,:])
                    if ( abs(mysdf) < abs(self.phi_s[j])):
                        self.phi_s[j]=mysdf
                for eN in range(self.model.q['x'].shape[0]):
                    for k in range(self.model.q['x'].shape[1]):
                        self.particle_signed_distances[i,eN,k],self.particle_signed_distance_normals[i,eN,k] = sdf(t, self.model.q['x'][eN,k])
                        self.particle_velocities[i,eN,k]=vel(t,self.model.q['x'][eN,k])


         
	if self.model.comm.isMaster():
            self.wettedAreaHistory.write("%21.16e\n" % (self.wettedAreas[-1],))
            self.forceHistory_p.write("%21.16e %21.16e %21.16e\n" %tuple(self.netForces_p[-1,:]))
            self.forceHistory_p.flush()
            self.forceHistory_v.write("%21.16e %21.16e %21.16e\n" %tuple(self.netForces_v[-1,:]))
            self.forceHistory_v.flush()
            self.momentHistory.write("%21.15e %21.16e %21.16e\n" % tuple(self.netMoments[-1,:]))
            self.momentHistory.flush()
            if self.nParticles:
                self.particle_forceHistory.write("%21.16e %21.16e %21.16e\n" %tuple(self.particle_netForces[0,:]))
                self.particle_forceHistory.flush()
                self.particle_momentHistory.write("%21.15e %21.16e %21.16e\n" % tuple(self.particle_netMoments[0,:]))
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
                 name='RANS3PF',
                 reuse_trial_and_test_quadrature=True,
                 sd=True,
                 movingDomain=False,
                 bdyNullSpace=False):
        self.bdyNullSpace=bdyNullSpace
        self.firstStep=True
        self.eb_adjoint_sigma = coefficients.eb_adjoint_sigma
        # this is a hack to test the effect of using a constant smoothing width
        useConstant_he = coefficients.useConstant_he
        self.postProcessing = True
        #
        # set the objects describing the method and boundary conditions
        #
        self.movingDomain = coefficients.movingDomain
        self.tLast_mesh = None
        # mql: Kill_pressure_term? This is for debugging and for convergence of momentum equations
        if ('KILL_PRESSURE_TERM') in dir(options):
            self.KILL_PRESSURE_TERM = options.KILL_PRESSURE_TERM
        else: 
            self.KILL_PRESSURE_TERM = False
        # mql: Use entropy viscosity?
        self.STABILIZATION_TYPE = 0 #0: SUPG, 1: EV via weak residual, 2: EV via strong residual
        if ('STABILIZATION_TYPE') in dir(options):
            self.STABILIZATION_TYPE = options.STABILIZATION_TYPE 
        # mql: Check if materialParameters are declared. This is for convergence tests
        self.hasMaterialParametersAsFunctions = False
        if ('materialParameters') in dir(options):
            self.materialParameters = options.materialParameters
            self.hasMaterialParametersAsFunctions = True
        # mql: Check if forceTerms are declared
        self.hasForceTermsAsFunctions = False
        if ('forceTerms') in dir(options):
            self.forceTerms = options.forceTerms
            self.hasForceTermsAsFunctions = True
        # mql: Check if analytical pressure function is declared
        if ('analyticalPressureSolution') in dir(options):
            self.analyticalPressureSolution = options.analyticalPressureSolution
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
        if isinstance(
                self.u[0].femSpace,
                C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess = True
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.phi_s = phiDict
        self.matType = matType
        # mwf try to reuse test and trial information across components if
        # spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature  # True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1, coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[
                    0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        # Simplicial Mesh
        # assume the same mesh for  all components for now
        self.mesh = self.u[0].femSpace.mesh
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.dirichletNodeSetList = None
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        # no velocity post-processing for now
        self.conservativeFlux = conservativeFluxDict
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.advection.has_key(ci):
                    for flag in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for flag in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if coefficients.potential.has_key(ci):
                    for flag in coefficients.potential[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = (
                (self.conservativeFlux is not None) or (
                    numericalFluxType is not None) or (
                    self.fluxBoundaryConditions[ci] == 'outFlow') or (
                    self.fluxBoundaryConditions[ci] == 'mixedFlow') or (
                    self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        # assume same space dim for all variables
        self.nSpace_global = self.u[0].femSpace.nSpace_global
        self.nDOF_trial_element = [
            u_j.femSpace.max_nDOF_element for u_j in self.u.values()]
        self.nDOF_phi_trial_element = [
            phi_k.femSpace.max_nDOF_element for phi_k in self.phi.values()]
        self.n_phi_ip_element = [
            phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in self.phi.values()]
        self.nDOF_test_element = [
            femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global = [
            dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
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
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[
                        ('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff', ci, ci)):
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            ('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            'default']
                else:
                    elementQuadratureDict[
                        ('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature['default']
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
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * \
            self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints, self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[
            0]
        self.nElementBoundaryQuadraturePoints_global = (
            self.mesh.nElements_global *
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
        # mql. Some structures needed for stabilization operators
        self.u_dof_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.v_dof_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.w_dof_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.u_dof_old_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.v_dof_old_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.w_dof_old_old = numpy.zeros(self.u[0].dof.shape,'d')
        self.ML=None #lumped mass matrix
        self.MC_global=None #consistent mass matrix
        self.cterm_global=None
        # mesh
        self.ebqe['x'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.ebq_global[
            ('totalFlux',
             0)] = numpy.zeros(
            (self.mesh.nElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebq_global[
            ('velocityAverage',
             0)] = numpy.zeros(
            (self.mesh.nElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.q['p'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m', 0)] = self.q[('u', 0)]
        self.q[('m', 1)] = self.q[('u', 1)]
        self.q[('m', 2)] = self.q[('u', 2)]
        self.q[('m_last', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('mt', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dV'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dV_last'] = -1000 * \
            numpy.ones((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('dV_u',0)] = self.q[('dV')]
        self.q[('dV_u',1)] = self.q[('dV')]
        self.q[('dV_u',2)] = self.q[('dV')]
        self.q[('f', 0)] = numpy.zeros((self.mesh.nElements_global,
                                        self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[
            ('velocity',
             0)] = numpy.zeros(
                 (self.mesh.nElements_global,
                  self.nQuadraturePoints_element,
                  self.nSpace_global),
                 'd')
        self.q['velocity_solid'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        # mql: create vectors to compute div(U)
        self.q['divU'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        # mql: create vectors to compute uStar = 2*un-unm1
        self.q[
            ('velocityOld',
            0)] = numpy.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nSpace_global),
                'd')
        self.q[
            ('velocityStar',
             0)] = numpy.zeros(
                 (self.mesh.nElements_global,
                  self.nQuadraturePoints_element,
                  self.nSpace_global),
                 'd')
        self.q['phi_solid'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['x'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element, 3), 'd')
        self.q[('cfl', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 0, 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 1, 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 2, 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[
            ('u',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('u',
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('u',
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('advectiveFlux_bc_flag',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('advectiveFlux_bc_flag',
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('advectiveFlux_bc_flag',
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('diffusiveFlux_bc_flag',
             0,
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('diffusiveFlux_bc_flag',
             1,
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('diffusiveFlux_bc_flag',
             2,
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'i')
        self.ebqe[
            ('advectiveFlux_bc',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('advectiveFlux_bc',
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('advectiveFlux_bc',
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('diffusiveFlux_bc',
             0,
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe['penalty'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('diffusiveFlux_bc',
             1,
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('diffusiveFlux_bc',
             2,
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('velocity',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe[
            ('velocity',
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe[
            ('velocity',
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        # VRANS start, defaults to RANS
        self.q[('r', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['eddy_viscosity'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        # VRANS end
        # RANS 2eq Models start
        self.q[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[
            ('grad(u)',
             1)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[
            ('grad(u)',
             2)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        # probably don't need ebqe gradients
        self.ebqe[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe[
            ('grad(u)',
             1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe[
            ('grad(u)',
             2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        # RANS 2eq Models end
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set(
            [('u', ci) for ci in range(self.nc)])
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
            self.q['abs(det(J))'] = numpy.zeros(
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
        log("Dumping quadrature shapes for model %s" % self.name, level=9)
        log("Element quadrature array (q)", level=9)
        for (k, v) in self.q.iteritems():
            log(str((k, v.shape)), level=9)
        log("Element boundary quadrature (ebq)", level=9)
        for (k, v) in self.ebq.iteritems():
            log(str((k, v.shape)), level=9)
        log("Global element boundary quadrature (ebq_global)", level=9)
        for (k, v) in self.ebq_global.iteritems():
            log(str((k, v.shape)), level=9)
        log("Exterior element boundary quadrature (ebqe)", level=9)
        for (k, v) in self.ebqe.iteritems():
            log(str((k, v.shape)), level=9)
        log("Interpolation points for nonlinear diffusion potential (phi_ip)", level=9)
        for (k, v) in self.phi_ip.iteritems():
            log(str((k, v.shape)), level=9)
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
            self.inflowBoundaryBC[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
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
        log("Updating local to global mappings", 2)
        self.updateLocal2Global()
        log("Building time integration object", 2)
        log(memory("inflowBC, internalNodes,updateLocal2Global",
                   "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(
                self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration", "OneLevelTransport"), level=4)
        log("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()

        self.setupFieldStrides()
        # Aux quantity at DOFs to be filled by optimized code (MQL)
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        self.entropyResidualAtCell = numpy.zeros(self.mesh.nElements_global,'d')
        self.maxSpeed2AtCell = numpy.zeros(self.mesh.nElements_global,'d')
        self.rhoAtCell = numpy.zeros(self.mesh.nElements_global,'d')
        self.muAtCell = numpy.zeros(self.mesh.nElements_global,'d')

        # mql: material parameters defined by a function at quad points
        self.q['density'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q['dynamic_viscosity'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe['density'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe['dynamic_viscosity'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        # Initialize material parameters 
        if self.hasMaterialParametersAsFunctions: 
            self.updateMaterialParameters()
        # mql: force terms
        self.q[('force', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('force', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('force', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')

        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        log("initalizing numerical flux")
        log(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict,
                    options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        log("initializing numerical flux penalty")
        self.numericalFlux.penalty_constant = self.coefficients.eb_penalty_constant
        # cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = self.numericalFlux.penalty_constant / (
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        # penalty term
        # cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = self.numericalFlux.penalty_constant / \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        log("setting up post-processing")
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(
            self)
        log(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        log("initializing archiver")
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        log("flux bc objects")
        for ci, fbcObject in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(
                self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[
                        ('advectiveFlux_bc', ci)][
                        t[0], t[1]] = g(
                        self.ebqe[
                            ('x')][
                            t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
            for ck, diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag', ck, ci)] = numpy.zeros(
                    self.ebqe[('diffusiveFlux_bc', ck, ci)].shape, 'i')
                for t, g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[
                        ('diffusiveFlux_bc', ck, ci)][
                        t[0], t[1]] = g(
                        self.ebqe[
                            ('x')][
                            t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[
                        ('diffusiveFlux_bc_flag', ck, ci)][
                        t[0], t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN = 1.0
        else:
            self.MOVING_DOMAIN = 0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(
                self.mesh.nodeArray.shape, 'd')
        # cek/ido todo replace python loops in modules with optimized code if
        # possible/necessary
        log("dirichlet conditions")
        self.forceStrongConditions = coefficients.forceStrongDirichlet
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(
                    self.u[cj].femSpace,
                    dofBoundaryConditionsSetterDict[cj],
                    weakDirichletConditions=False)
        log("final allocations")
        compKernelFlag = 0
        if self.coefficients.useConstant_he:
            self.elementDiameter = self.mesh.elementDiametersArray.copy()
            self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        else:
            self.elementDiameter = self.mesh.elementDiametersArray
        if self.nSpace_global == 2:
            import copy
            self.u[2] = copy.deepcopy(self.u[1])
            self.timeIntegration.m_tmp[
                2] = self.timeIntegration.m_tmp[1].copy()
            self.timeIntegration.beta_bdf[
                2] = self.timeIntegration.beta_bdf[1].copy()
            self.coefficients.sdInfo[(0, 2)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(1, 2)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(2, 0)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(2, 0)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(2, 1)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.coefficients.sdInfo[(2, 2)] = (numpy.array(
                [0, 1, 2], dtype='i'), numpy.array([0, 1], dtype='i'))
            self.offset.append(self.offset[1])
            self.stride.append(self.stride[1])
            self.numericalFlux.isDOFBoundary[
                2] = self.numericalFlux.isDOFBoundary[1].copy()
            self.numericalFlux.ebqe[
                ('u', 2)] = self.numericalFlux.ebqe[
                ('u', 1)].copy()
            log("calling RANS3PF2D ctor")
            self.rans3pf = cRANS3PF.RANS3PF2D(
                self.nSpace_global,
                self.nQuadraturePoints_element,
                self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                self.nElementBoundaryQuadraturePoints_elementBoundary,
                compKernelFlag,
                self.coefficients.aDarcy,
		self.coefficients.betaForch,
		self.coefficients.grain,
		self.coefficients.packFraction,
		self.coefficients.packMargin,
		self.coefficients.maxFraction,
		self.coefficients.frFraction,
                self.coefficients.sigmaC,
                self.coefficients.C3e,
                self.coefficients.C4e,
                self.coefficients.eR,
 		self.coefficients.fContact,
                self.coefficients.mContact,
                self.coefficients.nContact,
                self.coefficients.angFriction)
        else:
            log("calling  RANS3PF_base ctor")
            self.rans3pf = cRANS3PF.RANS3PF(
                self.nSpace_global,
                self.nQuadraturePoints_element,
                self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                self.nElementBoundaryQuadraturePoints_elementBoundary,
                compKernelFlag,
                self.coefficients.aDarcy,
		self.coefficients.betaForch,
		self.coefficients.grain,
		self.coefficients.packFraction,
		self.coefficients.packMargin,
		self.coefficients.maxFraction,
		self.coefficients.frFraction,
                self.coefficients.sigmaC,
                self.coefficients.C3e,
                self.coefficients.C4e,
                self.coefficients.eR,
 		self.coefficients.fContact,
                self.coefficients.mContact,
                self.coefficients.nContact,
                self.coefficients.angFriction)
        

        self.phisErrorNodal=self.u[0].dof.copy()
        self.velocityErrorNodal=self.u[0].dof.copy()


    def updateMaterialParameters(self):
        x = self.q[('x')][:,:,0]
        y = self.q[('x')][:,:,1]
        z = self.q[('x')][:,:,2]
        X = {0:x,
             1:y,
             2:z}
        t = self.timeIntegration.t
        self.q['density'][:] = self.materialParameters['density'](X,t)
        self.q['dynamic_viscosity'][:] = self.materialParameters['dynamic_viscosity'](X,t)
        # BOUNDARY
        ebqe_X = {0:self.ebqe['x'][:,:,0],
                  1:self.ebqe['x'][:,:,1],
                  2:self.ebqe['x'][:,:,2]}
        self.ebqe['density'][:] = self.materialParameters['density'](ebqe_X,t)
        self.ebqe['dynamic_viscosity'][:] = self.materialParameters['dynamic_viscosity'](ebqe_X,t)

    def updateForceTerms(self): 
        x = self.q[('x')][:,:,0]
        y = self.q[('x')][:,:,1]
        z = self.q[('x')][:,:,2]
        X = {0:x,
             1:y,
             2:z}
        t = self.timeIntegration.t
        self.q[('force', 0)][:] = self.forceTerms[0](X,t)
        self.q[('force', 1)][:] = self.forceTerms[1](X,t)
        try: 
            self.q[('force', 2)][:] = self.forceTerms[2](X,t)
        except: 
            pass

    def getResidual(self, u, r):
        """
        Calculate the element residuals and add in to the global residual
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
            for ci, fbcObject in self.fluxBoundaryConditionsObjectsDict.iteritems():
                for t, g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                    if self.coefficients.advection.has_key(ci):
                        self.ebqe[
                            ('advectiveFlux_bc', ci)][
                            t[0], t[1]] = g(
                            self.ebqe[
                                ('x')][
                                t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[
                            ('advectiveFlux_bc_flag', ci)][
                            t[0], t[1]] = 1
                for ck, diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                    for t, g in diffusiveFluxBoundaryConditionsDict.iteritems():
                        self.ebqe[
                            ('diffusiveFlux_bc', ck, ci)][
                            t[0], t[1]] = g(
                            self.ebqe[
                                ('x')][
                                t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[
                            ('diffusiveFlux_bc_flag', ck, ci)][
                            t[0], t[1]] = 1
        r.fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 36.0
        # TODO how to request problem specific evaluations from coefficient
        # class
        if 'evaluateForcingTerms' in dir(self.coefficients):
            self.coefficients.evaluateForcingTerms(
                self.timeIntegration.t,
                self.q,
                self.mesh,
                self.u[0].femSpace.elementMaps.psi,
                self.mesh.elementNodesArray)
        self.coefficients.wettedAreas[:] = 0.0
        self.coefficients.netForces_p[:, :] = 0.0
        self.coefficients.netForces_v[:, :] = 0.0
        self.coefficients.netMoments[:, :] = 0.0
        self.coefficients.particle_netForces[:,:]=0.0
        self.coefficients.particle_netMoments[:,:]=0.0
        self.coefficients.particle_surfaceArea[:]=0.0

        if self.forceStrongConditions and self.firstStep == False:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in self.dirichletConditionsForceDOF[
                        cj].DOFBoundaryConditionsDict.iteritems():
                        self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[
                            dofN], self.timeIntegration.t)# + self.MOVING_DOMAIN * self.mesh.nodeVelocityArray[dofN, cj - 1]
                        
        if self.coefficients.set_vos:
            self.coefficients.set_vos(self.q['x'],self.coefficients.q_vos)
            self.coefficients.set_vos(self.ebqe['x'],self.coefficients.ebqe_vos)

        # mql: get sparsity pattern corresponding to a 1D system. This is to loop on DOFs
        if self.cterm_global is None:
            self.cterm={}
            self.cterm_a={}
            self.cterm_global={}
            rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
            nnz = nzval.shape[-1] #number of non-zero entries in sparse matrix
            nnz_cMatrix = nnz/self.nc/self.nc
            nzval_cMatrix = numpy.zeros(nnz_cMatrix)
            rowptr_cMatrix = numpy.zeros(self.nFreeDOF_global[0]+1,'i')
            colind_cMatrix = numpy.zeros(nnz_cMatrix,'i')
            #fill vector rowptr_cMatrix
            for i in range(1,rowptr_cMatrix.size):
                rowptr_cMatrix[i] = rowptr_cMatrix[i-1]+(rowptr[self.nc*(i-1)+1]-rowptr[self.nc*(i-1)])/self.nc 

            #fill vector colind_cMatrix
            i_cMatrix=0 #ith row of cMatrix
            for i in range(rowptr.size-1): #0 to total num of DOFs (i.e. num of rows of jacobian)
                if (i%self.nc == 0): # Just consider the rows related to the 1st variable
                    for j,offset in enumerate(range(rowptr[i],rowptr[i+1])):
                        offset_cMatrix = range(rowptr_cMatrix[i_cMatrix],rowptr_cMatrix[i_cMatrix+1])
                        if (j%self.nc == 0):
                            colind_cMatrix[offset_cMatrix[j/self.nc]] = colind[offset]/self.nc
                    i_cMatrix+=1
            # END OF SPARSITY PATTERN FOR C MATRICES 
            di = numpy.zeros((self.mesh.nElements_global,
                              self.nQuadraturePoints_element,
                              self.nSpace_global),
                             'd') #direction of derivative
            # JACOBIANS (FOR ELEMENT TRANSFORMATION)
            self.q[('J')] = numpy.zeros((self.mesh.nElements_global,
                                         self.nQuadraturePoints_element,
                                         self.nSpace_global,
                                         self.nSpace_global),
                                        'd')
            self.q[('inverse(J)')] = numpy.zeros((self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.nSpace_global,
                                                  self.nSpace_global),
                                                 'd')
            self.q[('det(J)')] = numpy.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element),
                                             'd')
            self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                             self.q['J'],
                                                             self.q['inverse(J)'],
                                                             self.q['det(J)'])
            self.q['abs(det(J))'] = numpy.abs(self.q['det(J)'])
            # SHAPE FUNCTIONS
            self.q[('w',0)] = numpy.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nDOF_test_element[0]),
                                          'd')
            self.q[('w*dV_m',0)] = self.q[('w',0)].copy()
            self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, self.q[('w',0)])
            cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u',0)],
                                                 self.q['abs(det(J))'],
                                                 self.q[('w',0)],
                                                 self.q[('w*dV_m',0)])
            #### GRADIENT OF TEST FUNCTIONS 
            self.q[('grad(w)',0)] = numpy.zeros((self.mesh.nElements_global,
                                                 self.nQuadraturePoints_element,
                                                 self.nDOF_test_element[0],
                                                 self.nSpace_global),
                                                'd')
            self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                      self.q['inverse(J)'],
                                                      self.q[('grad(w)',0)])
            self.q[('grad(w)*dV_f',0)] = numpy.zeros((self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nDOF_test_element[0],
                                                      self.nSpace_global),
                                                     'd')
            cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u',0)],
                                                          self.q['abs(det(J))'],
                                                          self.q[('grad(w)',0)],
                                                          self.q[('grad(w)*dV_f',0)])
            # LUMPED MASS MATRIX #
            #assume a linear mass term
            dm = np.ones(self.q[('u',0)].shape,'d')
            elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                          self.nDOF_test_element[0],
                                          self.nDOF_trial_element[0]),'d')
            cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                         self.q[('w',0)],
                                                         self.q[('w*dV_m',0)],
                                                         elementMassMatrix)
            self.MC_a = nzval_cMatrix.copy()
            self.MC_global = SparseMat(self.nFreeDOF_global[0],
                                       self.nFreeDOF_global[0],
                                       nnz_cMatrix,
                                       self.MC_a,
                                       colind_cMatrix,
                                       rowptr_cMatrix)
            cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.MC_global)
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0,0)]/self.nc/self.nc,
                                                                      self.csrColumnOffsets[(0,0)]/self.nc,
                                                                      elementMassMatrix,
                                                                      self.MC_global)
            diamD2 = numpy.sum(self.q['abs(det(J))'][:]*self.elementQuadratureWeights[('u',0)])  
            self.ML = np.zeros((self.nFreeDOF_global[0],),'d')
            for i in range(self.nFreeDOF_global[0]):
                self.ML[i] = self.MC_a[rowptr_cMatrix[i]:rowptr_cMatrix[i+1]].sum()
            np.testing.assert_almost_equal(self.ML.sum(), diamD2, 
                                           err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)

            # COMPUTE C-MATRICES #
            for d in range(self.nSpace_global): #spatial dimensions
                #C matrices
                self.cterm[d] = numpy.zeros((self.mesh.nElements_global,
                                             self.nDOF_test_element[0],
                                             self.nDOF_trial_element[0]),'d')
                self.cterm_a[d] = nzval_cMatrix.copy()
                self.cterm_global[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                    self.nFreeDOF_global[0],
                                                                    nnz_cMatrix,
                                                                    self.cterm_a[d],
                                                                    colind_cMatrix,
                                                                    rowptr_cMatrix)
                cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.cterm_global[d])
                di[:] = 0.0
                di[...,d] = 1.0
                cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(di,
                                                                    self.q[('grad(w)*dV_f',0)],
                                                                    self.q[('w',0)],
                                                                    self.cterm[d]) # int[(di*grad(wj))*wi*dV]
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.csrRowIndeces[(0,0)]/self.nc/self.nc,
                                                                          self.csrColumnOffsets[(0,0)]/self.nc,
                                                                          self.cterm[d],
                                                                          self.cterm_global[d])
                print "****************... ", self.nFreeDOF_global[0]

        rowptr_cMatrix, colind_cMatrix, Cx = self.cterm_global[0].getCSRrepresentation()
        rowptr_cMatrix, colind_cMatrix, Cy = self.cterm_global[1].getCSRrepresentation()
        if (self.nSpace_global == 3):
            rowptr_cMatrix, colind_cMatrix, Cz = self.cterm_global[2].getCSRrepresentation()
        else:
            Cz = numpy.zeros(Cx.shape,'d')

        # mql: select appropiate functions to compute residual and jacobian
        if (self.STABILIZATION_TYPE == 1 or self.STABILIZATION_TYPE == 2):
            self.calculateResidual = self.rans3pf.calculateResidual_entropy_viscosity
            self.calculateJacobian = self.rans3pf.calculateJacobian_entropy_viscosity
        else: 
            self.calculateResidual = self.rans3pf.calculateResidual
            self.calculateJacobian = self.rans3pf.calculateJacobian
        
        self.calculateResidual(  # element
            self.pressureModel.u[0].femSpace.elementMaps.psi,
            self.pressureModel.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.coefficients.PSTAB,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.pressureModel.u[0].femSpace.psi,
            self.pressureModel.u[0].femSpace.grad_psi,
            self.pressureModel.u[0].femSpace.psi,
            self.pressureModel.u[0].femSpace.grad_psi,
            self.pressureModel.q_p_sharp,
            self.pressureModel.q_grad_p_sharp, 
            self.pressureModel.ebqe_p_sharp,
            self.pressureModel.ebqe_grad_p_sharp,
            #self.pressureModel.q[('u',0)],
            #self.pressureModel.q[('grad(u)',0)],
            #self.pressureModel.ebqe[('u',0)],
            #self.pressureModel.ebqe[('grad(u)',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.Hessian_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.pressureModel.u[0].femSpace.elementMaps.psi_trace,
            self.pressureModel.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.pressureModel.u[0].femSpace.psi_trace,
            self.pressureModel.u[0].femSpace.grad_psi_trace,
            self.pressureModel.u[0].femSpace.psi_trace,
            self.pressureModel.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.eb_adjoint_sigma,
            self.elementDiameter,  # mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.mesh.nElements_owned,
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
            self.coefficients.epsFact_solid,
            self.coefficients.q_phi_solid,
            self.coefficients.q_velocity_solid,
            self.coefficients.q_vos,
            self.coefficients.q_dvos_dt,
            self.coefficients.q_dragAlpha,
            self.coefficients.q_dragBeta,
            self.q[('r', 0)],
            self.coefficients.q_turb_var[0],
            self.coefficients.q_turb_var[1],
            self.coefficients.q_turb_var_grad[0],
            self.q['eddy_viscosity'],
            self.pressureModel.u[0].femSpace.dofMap.l2g,
            self.u[0].femSpace.dofMap.l2g,
            self.pressureModel.u[0].dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u_dof_old,
            self.v_dof_old,
            self.w_dof_old,
            self.u_dof_old_old,
            self.v_dof_old_old,
            self.w_dof_old_old,
            self.coefficients.g,
            self.coefficients.useVF,
            self.coefficients.q_vf,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.m_tmp[0],
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.q[('f', 0)],
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.q['dV'],
            self.q['dV_last'],
            self.q[('velocityStar',0)], #mql: use uStar=2*un-unm1 to achieve 2nd order accuracy             
            self.coefficients.ebqe_velocity_last,
            self.q[('cfl', 0)],
            self.q[('numDiff', 0, 0)],
            self.q[('numDiff', 1, 1)],
            self.q[('numDiff', 2, 2)],
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[(0, 0)][0],
            self.coefficients.sdInfo[(0, 0)][1],
            self.coefficients.sdInfo[(0, 1)][0],
            self.coefficients.sdInfo[(0, 1)][1],
            self.coefficients.sdInfo[(0, 2)][0],
            self.coefficients.sdInfo[(0, 2)][1],
            self.coefficients.sdInfo[(1, 1)][0],
            self.coefficients.sdInfo[(1, 1)][1],
            self.coefficients.sdInfo[(1, 0)][0],
            self.coefficients.sdInfo[(1, 0)][1],
            self.coefficients.sdInfo[(1, 2)][0],
            self.coefficients.sdInfo[(1, 2)][1],
            self.coefficients.sdInfo[(2, 2)][0],
            self.coefficients.sdInfo[(2, 2)][1],
            self.coefficients.sdInfo[(2, 0)][0],
            self.coefficients.sdInfo[(2, 0)][1],
            self.coefficients.sdInfo[(2, 1)][0],
            self.coefficients.sdInfo[(2, 1)][1],
            self.pressureModel.offset[0],
            self.offset[0],
            self.offset[1],
            self.offset[2],
            self.pressureModel.stride[0],
            self.stride[0],
            self.stride[1],
            self.stride[2],
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
            self.coefficients.ebqe_vos,
            self.coefficients.ebqe_turb_var[0],
            self.coefficients.ebqe_turb_var[1],
            self.pressureModel.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.pressureModel.numericalFlux.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 1)],
            self.ebqe[('advectiveFlux_bc_flag', 2)],
            self.ebqe[('diffusiveFlux_bc_flag', 0, 0)],
            self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
            self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
            self.pressureModel.numericalFlux.ebqe[('u', 0)],
            self.pressureModel.numericalFlux.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 1)],
            self.ebqe[('advectiveFlux_bc', 2)],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('diffusiveFlux_bc', 0, 0)],
            self.ebqe['penalty'],
            self.numericalFlux.ebqe[('u', 1)],
            self.ebqe[('diffusiveFlux_bc', 1, 1)],
            self.numericalFlux.ebqe[('u', 2)],
            self.ebqe[('diffusiveFlux_bc', 2, 2)],
            self.q['x'],
            self.q[('velocity', 0)],
            self.ebqe[('velocity', 0)],
            self.q[('grad(u)',0)],
            self.q[('grad(u)',1)],
            self.q[('grad(u)',2)],
            self.q['divU'],
            self.ebqe[('grad(u)',0)],
            self.ebqe[('grad(u)',1)],
            self.ebqe[('grad(u)',2)],
            self.ebq_global[('totalFlux', 0)],
            self.elementResidual[0],
            self.mesh.elementMaterialTypes,
            self.mesh.elementBoundaryMaterialTypes,
            self.coefficients.barycenters,
            self.coefficients.wettedAreas,
            self.coefficients.netForces_p,
            self.coefficients.netForces_v,
            self.coefficients.netMoments,
            self.coefficients.q_rho,
            self.coefficients.ebqe_rho,
            self.coefficients.q_nu,
            self.coefficients.ebqe_nu,
            self.coefficients.nParticles,
	    self.coefficients.particle_epsFact,
	    self.coefficients.particle_alpha,
	    self.coefficients.particle_beta,
	    self.coefficients.particle_penalty_constant,
	    self.coefficients.particle_signed_distances,
	    self.coefficients.particle_signed_distance_normals,
	    self.coefficients.particle_velocities,
	    self.coefficients.particle_centroids,
            self.coefficients.particle_netForces,
            self.coefficients.particle_netMoments,
            self.coefficients.particle_surfaceArea,
            self.coefficients.particle_nitsche, 
            self.STABILIZATION_TYPE,
            self.elementQuadratureWeights[('u',0)].sum(),
            self.coefficients.cMax, 
            self.coefficients.cE, 
            self.q[('force', 0)],
            self.q[('force', 1)],
            self.q[('force', 2)], 
            self.KILL_PRESSURE_TERM,
            self.timeIntegration.dt, 
            self.entropyResidualAtCell,
            self.maxSpeed2AtCell,
            self.maxSpeed2AtCell.max(),
            self.rhoAtCell, 
            self.muAtCell, 
            self.quantDOFs, 
            self.nFreeDOF_global[0],
            rowptr_cMatrix,
            colind_cMatrix, 
            self.ML, 
            Cx, 
            Cy, 
            Cz, 
            self.hasMaterialParametersAsFunctions, 
            self.q['density'], 
            self.q['dynamic_viscosity'], 
            self.ebqe['density'],
            self.ebqe['dynamic_viscosity'],             
            self.u[0].femSpace.order)

        # mql: Save the solution in 'u' to allow SimTools.py to compute the errors
        for dim in range(self.nSpace_global):
            self.q[('u',dim)][:] = self.q[('velocity',0)][:,:,dim]

        from proteus.flcbdfWrappers import globalSum
        for i in range(self.coefficients.netForces_p.shape[0]):
            self.coefficients.wettedAreas[i] = globalSum(
                self.coefficients.wettedAreas[i])
            for I in range(3):
                self.coefficients.netForces_p[i, I] = globalSum(
                    self.coefficients.netForces_p[i, I])
                self.coefficients.netForces_v[i, I] = globalSum(
                    self.coefficients.netForces_v[i, I])
                self.coefficients.netMoments[i, I] = globalSum(
                    self.coefficients.netMoments[i, I])
        for i in range(self.coefficients.particle_netForces.shape[0]):
            for I in range(3):
                self.coefficients.particle_netForces[i, I] = globalSum(
                    self.coefficients.particle_netForces[i, I])
                self.coefficients.particle_netMoments[i, I] = globalSum(
                    self.coefficients.particle_netMoments[i, I])
            self.coefficients.particle_surfaceArea[i] = globalSum(
                self.coefficients.particle_surfaceArea[i])
            logEvent("particle i="+`i`+ " force "+`self.coefficients.particle_netForces[i]`)
            logEvent("particle i="+`i`+ " moment "+`self.coefficients.particle_netMoments[i]`)
            logEvent("particle i="+`i`+ " surfaceArea "+`self.coefficients.particle_surfaceArea[i]`)

        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    r[self.offset[cj] + self.stride[cj] * dofN] = self.u[cj].dof[dofN] - g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[
                        dofN], self.timeIntegration.t)# - self.MOVING_DOMAIN * self.mesh.nodeVelocityArray[dofN, cj - 1]

        cflMax = globalMax(self.q[('cfl', 0)].max()) * self.timeIntegration.dt
        log("Maximum CFL = " + str(cflMax), level=2)
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual", level=9, data=r)
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        if self.nSpace_global == 2:
            self.csrRowIndeces[(0, 2)] = self.csrRowIndeces[(0, 1)]
            self.csrColumnOffsets[(0, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrRowIndeces[(0, 2)] = self.csrRowIndeces[(0, 1)]
            self.csrColumnOffsets[(0, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrRowIndeces[(1, 2)] = self.csrRowIndeces[(0, 1)]
            self.csrColumnOffsets[(1, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrRowIndeces[(2, 0)] = self.csrRowIndeces[(1, 0)]
            self.csrColumnOffsets[(2, 0)] = self.csrColumnOffsets[(1, 0)]
            self.csrRowIndeces[(2, 0)] = self.csrRowIndeces[(1, 0)]
            self.csrColumnOffsets[(2, 0)] = self.csrColumnOffsets[(1, 0)]
            self.csrRowIndeces[(2, 1)] = self.csrRowIndeces[(1, 0)]
            self.csrColumnOffsets[(2, 1)] = self.csrColumnOffsets[(1, 0)]
            self.csrRowIndeces[(2, 2)] = self.csrRowIndeces[(1, 0)]
            self.csrColumnOffsets[(2, 2)] = self.csrColumnOffsets[(1, 0)]
            self.csrColumnOffsets_eb[(0, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(0, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(1, 2)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(2, 0)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(2, 0)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(2, 1)] = self.csrColumnOffsets[(0, 1)]
            self.csrColumnOffsets_eb[(2, 2)] = self.csrColumnOffsets[(0, 1)]

        self.calculateJacobian(  # element
            self.pressureModel.u[0].femSpace.elementMaps.psi,
            self.pressureModel.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.coefficients.PSTAB,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.pressureModel.u[0].femSpace.psi,
            self.pressureModel.u[0].femSpace.grad_psi,
            self.pressureModel.u[0].femSpace.psi,
            self.pressureModel.u[0].femSpace.grad_psi,
            self.pressureModel.q_p_sharp,
            self.pressureModel.q_grad_p_sharp,
            self.pressureModel.ebqe_p_sharp,
            self.pressureModel.ebqe_grad_p_sharp,
            #self.pressureModel.q[('u',0)],
            #self.pressureModel.q[('grad(u)',0)],
            #self.pressureModel.ebqe[('u',0)],
            #self.pressureModel.ebqe[('grad(u)',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.Hessian_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.pressureModel.u[0].femSpace.elementMaps.psi_trace,
            self.pressureModel.u[0].femSpace.elementMaps.grad_psi_trace,
            self.pressureModel.elementBoundaryQuadratureWeights[('u', 0)],
            self.pressureModel.u[0].femSpace.psi_trace,
            self.pressureModel.u[0].femSpace.grad_psi_trace,
            self.pressureModel.u[0].femSpace.psi_trace,
            self.pressureModel.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.eb_adjoint_sigma,
            self.elementDiameter,  # mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.stabilization.hFactor,
            self.mesh.nElements_global,
            self.mesh.nElements_owned,
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
            self.coefficients.q_vos,
            self.coefficients.q_dvos_dt,
            self.coefficients.q_dragAlpha,
            self.coefficients.q_dragBeta,
            self.pressureModel.q[('r', 0)],
            self.coefficients.q_turb_var[0],
            self.coefficients.q_turb_var[1],
            self.coefficients.q_turb_var_grad[0],
            # VRANS end
            self.pressureModel.u[0].femSpace.dofMap.l2g,
            self.u[0].femSpace.dofMap.l2g,
            self.pressureModel.u[0].dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.coefficients.g,
            self.coefficients.useVF,
            self.coefficients.q_vf,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.q['dV'],
            self.q['dV_last'],
            self.q[('velocityStar',0)], #mql: use uStar=2*un-unm1 to achieve 2nd order accuracy 
            self.coefficients.ebqe_velocity_last,
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[
                (0, 0)][0], self.coefficients.sdInfo[
                (0, 0)][1],
            self.coefficients.sdInfo[
                (0, 1)][0], self.coefficients.sdInfo[
                (0, 1)][1],
            self.coefficients.sdInfo[
                (0, 2)][0], self.coefficients.sdInfo[
                (0, 2)][1],
            self.coefficients.sdInfo[
                (1, 1)][0], self.coefficients.sdInfo[
                (1, 1)][1],
            self.coefficients.sdInfo[
                (1, 0)][0], self.coefficients.sdInfo[
                (1, 0)][1],
            self.coefficients.sdInfo[
                (1, 2)][0], self.coefficients.sdInfo[
                (1, 2)][1],
            self.coefficients.sdInfo[
                (2, 2)][0], self.coefficients.sdInfo[
                (2, 2)][1],
            self.coefficients.sdInfo[
                (2, 0)][0], self.coefficients.sdInfo[
                (2, 0)][1],
            self.coefficients.sdInfo[
                (2, 1)][0], self.coefficients.sdInfo[
                (2, 1)][1],
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 1)], self.csrColumnOffsets[(0, 1)],
            self.csrRowIndeces[(0, 2)], self.csrColumnOffsets[(0, 2)],
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 1)], self.csrColumnOffsets[(0, 1)],
            self.csrRowIndeces[(0, 2)], self.csrColumnOffsets[(0, 2)],
            self.csrRowIndeces[(1, 0)], self.csrColumnOffsets[(1, 0)],
            self.csrRowIndeces[(1, 0)], self.csrColumnOffsets[(1, 0)],
            self.csrRowIndeces[(1, 1)], self.csrColumnOffsets[(1, 1)],
            self.csrRowIndeces[(1, 2)], self.csrColumnOffsets[(1, 2)],
            self.csrRowIndeces[(2, 0)], self.csrColumnOffsets[(2, 0)],
            self.csrRowIndeces[(2, 0)], self.csrColumnOffsets[(2, 0)],
            self.csrRowIndeces[(2, 1)], self.csrColumnOffsets[(2, 1)],
            self.csrRowIndeces[(2, 2)], self.csrColumnOffsets[(2, 2)],
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
            self.coefficients.ebqe_vos,
            self.coefficients.ebqe_turb_var[0],
            self.coefficients.ebqe_turb_var[1],
            # VRANS end
            self.pressureModel.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.pressureModel.numericalFlux.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 1)],
            self.ebqe[('advectiveFlux_bc_flag', 2)],
            self.ebqe[('diffusiveFlux_bc_flag', 0, 0)],
            self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
            self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
            self.pressureModel.numericalFlux.ebqe[('u', 0)],
            self.pressureModel.numericalFlux.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 1)],
            self.ebqe[('advectiveFlux_bc', 2)],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('diffusiveFlux_bc', 0, 0)],
            self.ebqe['penalty'],
            self.numericalFlux.ebqe[('u', 1)],
            self.ebqe[('diffusiveFlux_bc', 1, 1)],
            self.numericalFlux.ebqe[('u', 2)],
            self.ebqe[('diffusiveFlux_bc', 2, 2)],
            self.csrColumnOffsets_eb[(0, 0)],
            self.csrColumnOffsets_eb[(0, 0)],
            self.csrColumnOffsets_eb[(0, 1)],
            self.csrColumnOffsets_eb[(0, 2)],
            self.csrColumnOffsets_eb[(0, 0)],
            self.csrColumnOffsets_eb[(0, 0)],
            self.csrColumnOffsets_eb[(0, 1)],
            self.csrColumnOffsets_eb[(0, 2)],
            self.csrColumnOffsets_eb[(1, 0)],
            self.csrColumnOffsets_eb[(1, 0)],
            self.csrColumnOffsets_eb[(1, 1)],
            self.csrColumnOffsets_eb[(1, 2)],
            self.csrColumnOffsets_eb[(2, 0)],
            self.csrColumnOffsets_eb[(2, 0)],
            self.csrColumnOffsets_eb[(2, 1)],
            self.csrColumnOffsets_eb[(2, 2)],
            self.mesh.elementMaterialTypes,
            self.coefficients.nParticles,
	    self.coefficients.particle_epsFact,
	    self.coefficients.particle_alpha,
	    self.coefficients.particle_beta,
	    self.coefficients.particle_penalty_constant,
	    self.coefficients.particle_signed_distances,
	    self.coefficients.particle_signed_distance_normals,
	    self.coefficients.particle_velocities,
	    self.coefficients.particle_centroids,
            self.coefficients.particle_nitsche, 
            self.KILL_PRESSURE_TERM,
            self.hasMaterialParametersAsFunctions, 
            self.q['density'], 
            self.q['dynamic_viscosity'], 
            self.ebqe['density'],
            self.ebqe['dynamic_viscosity'])

        if not self.forceStrongConditions and max(
            numpy.linalg.norm(
                self.u[0].dof, numpy.inf), numpy.linalg.norm(
                self.u[1].dof, numpy.inf), numpy.linalg.norm(
                self.u[2].dof, numpy.inf)) < 1.0e-8:
            self.pp_hasConstantNullSpace = True
        else:
            self.pp_hasConstantNullSpace = False
        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            for cj in range(self.nc):
                for dofN in self.dirichletConditionsForceDOF[
                        cj].DOFBoundaryConditionsDict.keys():
                    global_dofN = self.offset[cj] + self.stride[cj] * dofN
                    for i in range(
                        self.rowptr[global_dofN], self.rowptr[
                            global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            self.nzval[i] = 1.0
                        else:
                            self.nzval[i] = 0.0
                            # print "RBLES zeroing residual cj = %s dofN= %s
                            # global_dofN= %s " % (cj,dofN,global_dofN)
        log("Jacobian ", level=10, data=jacobian)
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
            self.u[0].femSpace.elementMaps.getValues(
                self.elementQuadraturePoints, self.q['x'])
            self.u[0].femSpace.elementMaps.getJacobianValues(
                self.elementQuadraturePoints,
                self.q['J'],
                self.q['inverse(J)'],
                self.q['det(J)'])
            self.q['abs(det(J))'][:] = numpy.abs(self.q['det(J)'])
            self.u[0].femSpace.getBasisValues(
                self.elementQuadraturePoints, self.q[('v', 0)])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisHessianValuesRef(
            self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(
            self.timeIntegration.t, self.q)
        if self.stabilization is not None and not domainMoved:
            self.stabilization.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None and not domainMoved:
            self.shockCapturing.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self, domainMoved=False):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.u[0].femSpace.elementMaps.getValuesTrace(
                self.elementBoundaryQuadraturePoints, self.ebq['x'])
            self.u[0].femSpace.elementMaps.getJacobianValuesTrace(
                self.elementBoundaryQuadraturePoints,
                self.ebq['inverse(J)'],
                self.ebq['g'],
                self.ebq['sqrt(det(g))'],
                self.ebq['n'])
            cfemIntegrals.copyLeftElementBoundaryInfo(
                self.mesh.elementBoundaryElementsArray,
                self.mesh.elementBoundaryLocalElementBoundariesArray,
                self.mesh.exteriorElementBoundariesArray,
                self.mesh.interiorElementBoundariesArray,
                self.ebq['x'],
                self.ebq['n'],
                self.ebq_global['x'],
                self.ebq_global['n'])
            self.u[0].femSpace.elementMaps.getInverseValuesTrace(
                self.ebq['inverse(J)'], self.ebq['x'], self.ebq['hat(x)'])
            self.u[0].femSpace.elementMaps.getPermutations(self.ebq['hat(x)'])
            self.testSpace[0].getBasisValuesTrace(
                self.u[0].femSpace.elementMaps.permutations, self.ebq['hat(x)'], self.ebq[
                    ('w', 0)])
            self.u[0].femSpace.getBasisValuesTrace(
                self.u[0].femSpace.elementMaps.permutations, self.ebq['hat(x)'], self.ebq[
                    ('v', 0)])
            cfemIntegrals.calculateElementBoundaryIntegrationWeights(
                self.ebq['sqrt(det(g))'], self.elementBoundaryQuadratureWeights[
                    ('u', 0)], self.ebq[
                    ('dS_u', 0)])

    def calculateExteriorElementBoundaryQuadrature(self, domainMoved=False):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        log("initalizing ebqe vectors for post-procesing velocity")
        if self.postProcessing:
            self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(
                self.elementBoundaryQuadraturePoints, self.ebqe['x'])
            self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace(
                self.elementBoundaryQuadraturePoints,
                self.ebqe['inverse(J)'],
                self.ebqe['g'],
                self.ebqe['sqrt(det(g))'],
                self.ebqe['n'])
            cfemIntegrals.calculateIntegrationWeights(
                self.ebqe['sqrt(det(g))'], self.elementBoundaryQuadratureWeights[
                    ('u', 0)], self.ebqe[
                    ('dS_u', 0)])
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        log("initalizing basis info")
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(
            self.elementBoundaryQuadraturePoints, self.ebqe['x'])
        log("setting flux boundary conditions")
        if not domainMoved:
            self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                       self.ebqe[('x')],
                                                                                       self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                       self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                           for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
            log("initializing coefficients ebqe")
            self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(
                self.timeIntegration.t, self.ebqe)
        log("done with ebqe")

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        if self.postProcessing and self.conservativeFlux:
            self.rans3pf.calculateVelocityAverage(
                self.mesh.nExteriorElementBoundaries_global,
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
                self.u[0].femSpace.dofMap.l2g,
                self.u[0].dof,
                self.u[1].dof,
                self.u[2].dof,
                self.coefficients.vos_dof,
                self.u[0].femSpace.psi_trace,
                self.ebqe[
                    ('velocity',
                     0)],
                self.ebq_global[
                    ('velocityAverage',
                     0)])
            if self.movingDomain:
                log("Element Quadrature", level=3)
                self.calculateElementQuadrature(domainMoved=True)
                log("Element Boundary Quadrature", level=3)
                self.calculateElementBoundaryQuadrature(domainMoved=True)
                log("Global Exterior Element Boundary Quadrature", level=3)
                self.calculateExteriorElementBoundaryQuadrature(
                    domainMoved=True)
                for ci in range(
                        len(self.velocityPostProcessor.vpp_algorithms)):
                    for cj in self.velocityPostProcessor.vpp_algorithms[
                            ci].updateConservationJacobian.keys():
                        self.velocityPostProcessor.vpp_algorithms[
                            ci].updateWeights()
                        self.velocityPostProcessor.vpp_algorithms[
                            ci].computeGeometricInfo()
                        self.velocityPostProcessor.vpp_algorithms[
                            ci].updateConservationJacobian[cj] = True


        self.q['velocityError'][:]=self.q[('velocity',0)]                 
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)
        self.q['velocityError']-=self.q[('velocity',0)]
        if self.q.has_key('phis') and self.coefficients.granular_sdf_Calc is not None:
            self.q['phisError']=self.q[('phis')]
        else:#this needs to be fixed for the case that multiple bodies are present
            if self.q.has_key(('phis',0)):
                self.q['phisError'][:]=self.q[('phis',0)]

    def updateAfterMeshMotion(self):
        pass


def getErgunDrag(porosity, meanGrainSize, viscosity):
    # cek hack, this d0oesn't seem right
    # cek todo look up correct Ergun model for alpha and beta
    voidFrac = 1.0 - porosity
    if voidFrac > 1.0e-6:
        dragBeta = porosity * porosity * porosity * meanGrainSize * 1.0e-2 / voidFrac
    if (porosity > epsZero and meanGrainSize > epsZero):
        dragAlpha = viscosity * 180.0 * voidFrac * voidFrac / \
            (meanGrainSize * meanGrainSize * porosity)
