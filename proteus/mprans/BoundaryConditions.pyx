"""
Module for creating boundary conditions. Imported in mprans.SpatialTools.py
"""
import sys
import numpy as np
from proteus import AuxiliaryVariables
from proteus.BoundaryConditions import (BC_Base,
                                        constantBC,
                                        linearBC)
from proteus.Profiling import logEvent as log
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)


class BC_RANS(BC_Base):
    """
    class defining boundary conditions for a facet or a segment
    """
    def __init__(self, shape=None, name=None, b_or=None, b_i=None):
        super(BC_RANS, self).__init__(shape, name, b_or, b_i)
        # _dirichlet
        self.p_dirichlet = None  # pressure
        self.u_dirichlet = None  # velocity u
        self.v_dirichlet = None  # velocity v
        self.w_dirichlet = None  # velocity w
        self.vof_dirichlet = None  # VOF
        self.k_dirichlet = None  # kappa
        self.dissipation_dirichlet = None  # dissipation
        # _advective
        self.p_advective = None
        self.u_advective = None
        self.v_advective = None
        self.w_advective = None
        self.vof_advective = None
        self.k_advective = None
        self.dissipation_advective = None
        # _diffusive
        self.u_diffusive = None
        self.v_diffusive = None
        self.w_diffusive = None
        self.k_diffusive = None
        self.dissipation_diffusive = None
        # moveMesh boundary conditions
        self.hx_dirichlet = None
        self.hy_dirichlet = None
        self.hz_dirichlet = None
        self.u_stress = 0.
        self.v_stress = 0.
        self.w_stress = 0.

    def reset(self):
        """
        Resets all boundary conditions to None, apart from:
        - moveMesh BCs
        """
        self.BC_type = 'None'
        self.p_dirichlet = None
        self.u_dirichlet = None
        self.v_dirichlet = None
        self.w_dirichlet = None
        self.vof_dirichlet = None
        self.k_dirichlet = None
        self.dissipation_dirichlet = None
        self.p_advective = None
        self.u_advective = None
        self.v_advective = None
        self.w_advective = None
        self.vof_advective = None
        self.k_advective = None
        self.dissipation_advective = None
        self.u_diffusive = None
        self.v_diffusive = None
        self.w_diffusive = None
        self.k_diffusive = None
        self.dissipation_diffusive = None

    def setNonMaterial(self):
        """
        For non-material boundaries.
        Sets diffusive flux and advective vof to zero
        """
        self.reset()
        self.BC_type = 'NonMaterial'
        self.vof_advective = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    def setTank(self):
        b_or = self._b_or[self._b_i]
        if b_or[0] == 1 or b_or[0] == -1:
            self.hx_dirichlet = constantBC(0.)
            self.u_stress = None
        elif b_or[1] == 1 or b_or[1] == -1:
            self.hy_dirichlet = constantBC(0.)
            self.v_stress = None
        elif len(b_or) > 2 and (b_or[2] == 1 or b_or[2] == -1):
            self.hz_dirichlet = constantBC(0.)
            self.w_stress = None

    def setFixedNodes(self):
        """
        For moving domains: fixes nodes/boundary
        """
        self.hx_dirichlet = constantBC(0.)
        self.hy_dirichlet = constantBC(0.)
        self.hz_dirichlet = constantBC(0.)
        self.u_stress = None
        self.v_stress = None
        self.w_stress = None

    def setNoSlip(self):
        """
        sets no slip conditions at the boundary
        """
        self.reset()
        self.BC_type = 'NoSlip'
        self.u_dirichlet = constantBC(0.)
        self.v_dirichlet = constantBC(0.)
        self.w_dirichlet = constantBC(0.)
        self.p_advective = constantBC(0.)
        self.vof_advective = constantBC(0.)
        self.k_dirichlet = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setFreeSlip(self):
        """
        sets free slip conditions at the boundary
        """
        self.reset()
        self.BC_type = 'FreeSlip'
        self.p_advective = constantBC(0.)
        self.u_advective = constantBC(0.)
        self.v_advective = constantBC(0.)
        self.w_advective = constantBC(0.)
        self.vof_advective = constantBC(0.)
        self.k_dirichlet = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setOpenAir(self, orientation=None):
        """
        sets open boundary conditions (water can come out)
        """
        self.BC_type = 'OpenAir'

        def get_ux_dirichlet(i):
            if b_or[i] == 1. or b_or[i] == -1.:
                return None
            else:
                return constantBC(0.)

        if orientation is None and self._b_or[self._b_i] is not None:
            b_or = self._b_or[self._b_i]
        elif orientation is not None:
            b_or = orientation
        else:
            print('Boundary orientation needs to be defined')
        self.reset()
        self.p_dirichlet = constantBC(0.)
        self.u_dirichlet = get_ux_dirichlet(0)
        self.v_dirichlet = get_ux_dirichlet(1)
        if len(b_or) > 2:
            self.w_dirichlet = get_ux_dirichlet(2)
        self.vof_dirichlet = constantBC(1.)  # air
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setMoveMesh(self, last_pos, h=(0., 0., 0.), rot_matrix=None):
        """
        sets rigid body boundary conditions for moving the mesh
        :param last_pos: position (barycentre) of body pre-calculation
        :param h: displacement during time step
        :param rot_matrix: rotation matrix of body pre-calculation
        (!) should not be set manually
        """
        if rot_matrix is None:
            rot_matrix = np.array([[1., 0., 0.],
                                   [0., 1., 0.],
                                   [0., 0., 1.]])
        def get_DBC_h(i):
            def DBC_h(x, t):
                x_0 = x-last_pos
                new_x_0 = np.dot(x_0, rot_matrix)
                hx = new_x_0-x_0+h
                return hx[i]
            return DBC_h
        self.hx_dirichlet = get_DBC_h(i=0)
        self.hy_dirichlet = get_DBC_h(i=1)
        if len(last_pos) > 2:
            self.hz_dirichlet = get_DBC_h(i=2)

    def setUnsteadyTwoPhaseVelocityInlet(self, wave, vert_axis=None,
                                         windSpeed=(0., 0., 0.), air=1.,
                                         water=0.):
        """
        Imposes a velocity profile lower than the sea level and an open
        boundary for higher than the sealevel.
        :param U: Velocity vector at the global system.
        :param vert_axis: index of vertical in position vector, must always be
                        aligned with gravity, by default set to 1].
        :param air: Volume fraction for air (1.0 by default).
        :param water: Volume fraction for water (0.0 by default).
        Below the seawater level, the condition returns the _dirichlet and
        p_advective condition according to the inflow velocity.
        Above the sea water level, the condition returns the gravity as zero,
        and sets _dirichlet condition to zero, only if there is a zero inflow
        velocity component.
        (!) This condition is best used for boundaries and gravity aligned with
            one of the main axes.
        (!) Boundary condition relies on specific variables defined in Context
        """
        self.reset()

        windSpeed=np.array(windSpeed)
        if vert_axis is None:
            vert_axis = self.Shape.Domain.nd-1

        def get_inlet_ux_dirichlet(i):
            def ux_dirichlet(x, t):
                from proteus import Context
                ct = Context.get()
                waveHeight = wave.mwl+wave.eta(x, t)
                wavePhi = x[vert_axis]-waveHeight
                if wavePhi <= 0:
                    waterSpeed = wave.u(x, t)
                else:
                    x_max = list(x)
                    x_max[vert_axis] = waveHeight
                    waterSpeed = wave.u(x_max, t)
                he, ecH = ct.domain.MeshOptions.he, ct.epsFact_consrv_heaviside
                # smoothing only above wave, only on half the VOF smoothing length
                H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
                ux = H*windSpeed + (1-H)*waterSpeed
                return ux[i]
            return ux_dirichlet

        def inlet_vof_dirichlet(x, t):
            from proteus import Context
            ct = Context.get()
            level = wave.mwl + wave.eta(x,t)
            mesh = ct.domain.MeshOptions
            he, ecH = ct.domain.MeshOptions.he, ct.epsFact_consrv_heaviside
            H = smoothedHeaviside(ecH*he,x[vert_axis]-level)
            return H

        def inlet_p_advective(x, t):
            from proteus import Context
            ct = Context.get()
            # This is the normal velocity, based on the outwards boundary
            # orientation b_or
            # needs to be equal to -ux_dirichlet
            b_or = self._b_or[self._b_i]
            nd = len(b_or)
            waterSpeed = np.array(wave.u(x, t))
            waveHeight = wave.mwl+wave.eta(x, t)
            wavePhi = x[vert_axis]-waveHeight
            he = ct.domain.MeshOptions.he
            if wavePhi <= 0:
                waterSpeed = wave.u(x, t)
            else:
                x_max = list(x)
                x_max[vert_axis] = waveHeight
                waterSpeed = wave.u(x_max, t)
            he, ecH = ct.domain.MeshOptions.he, ct.epsFact_consrv_heaviside
            # smoothing only above wave, only on half the VOF smoothing length
            H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
            U = H*windSpeed + (1-H)*waterSpeed
            u_p = np.sum(U[:nd]*b_or)
            return u_p

        self.u_dirichlet = get_inlet_ux_dirichlet(0)
        self.v_dirichlet = get_inlet_ux_dirichlet(1)
        self.w_dirichlet = get_inlet_ux_dirichlet(2)
        self.vof_dirichlet = inlet_vof_dirichlet
        self.p_advective = inlet_p_advective

    # FOLLOWING BOUNDARY CONDITION IS UNTESTED #
    def setTwoPhaseVelocityInlet(self, U, waterLevel, vert_axis=None, air=1.,
                                 water=0.):
        """
        Imposes a velocity profile lower than the sea level and an open
        boundary for higher than the sealevel.
        :param U: Velocity vector at the global system.
        :param waterLevel: water level at global coordinate system.
        :param vert_axis: index of vertical in position vector, must always be
                        aligned with gravity, by default set to 1].
        :param air: Volume fraction for air (1.0 by default).
        :param water: Volume fraction for water (0.0 by default).
        Below the seawater level, the condition returns the _dirichlet and
        p_advective condition according to the inflow velocity.
        Above the sea water level, the condition returns the gravity as zero,
        and sets _dirichlet condition to zero, only if there is a zero inflow
        velocity component.
        (!) This condition is best used for boundaries and gravity aligned with
            one of the main axes.
        """
        self.reset()
        self.BC_type = 'TwoPhaseVelocityInlet'

        U = np.array(U)
        if vert_axis is None:
            vert_axis = self.Shape.Domain.nd - 1

        def get_inlet_ux_dirichlet(ux):
            def ux_dirichlet(x, t):
                if x[vert_axis] < waterLevel:
                    return ux
                elif x[vert_axis] >= waterLevel and ux == 0:
                    return 0.
            return ux_dirichlet

        def inlet_vof_dirichlet(x, t):
            if x[vert_axis] < waterLevel:
                return water
            elif x[vert_axis] >= waterLevel:
                return air

        def inlet_p_advective(x, t, u=U):
            b_or = self._b_or[self._b_i]
            u_p = np.sum(U*b_or)
            # This is the normal velocity, based on the inwards boundary
            # orientation -b_or
            u_p = -u_p
            if x[vert_axis] < waterLevel:
                return u_p
            elif x[vert_axis] >= waterLevel:
                return None

        self.u_dirichlet = get_inlet_ux_dirichlet(U[0])
        self.v_dirichlet = get_inlet_ux_dirichlet(U[1])
        if len(U) == 3:
                self.w_dirichlet = get_inlet_ux_dirichlet(U[2])
        self.vof_dirichlet = inlet_vof_dirichlet
        self.p_advective = inlet_p_advective

    def setHydrostaticPressureOutlet(self, rho, g, refLevel, vof, pRef=0.0,
                                     vert_axis=-1):
        self.reset()
        a0 = pRef - rho*g[vert_axis]*refLevel
        a1 = rho*g[vert_axis]
        # This is the normal velocity, based on the boundary orientation

        def get_outlet_ux_dirichlet(i):
            def ux_dirichlet(x, t):
                b_or = self._b_or[self._b_i]
                if b_or[i] == 0:
                    return 0.
            return ux_dirichlet
        self.u_dirichlet = get_outlet_ux_dirichlet(0)
        self.v_dirichlet = get_outlet_ux_dirichlet(1)
        if len(g) == 3:
            self.w_dirichlet = get_outlet_ux_dirichlet(2)
        self.p_dirichlet = linearBC(a0, a1, vert_axis)
        self.vof_dirichlet = constantBC(vof)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    # FOLLOWING BOUNDARY CONDITION IS UNTESTED #
    def hydrostaticPressureOutletWithDepth(self, seaLevel, rhoUp, rhoDown, g,
                                           refLevel, pRef=0.0, vert_axis=None,
                                           air=1.0, water=0.0):
        """Imposes a hydrostatic pressure profile and open boundary conditions
        with a known otuflow depth
        :param rhoUp: Phase density of the upper part.
        :param rhoDown: Phase density of the lower part.
        :param g: Gravitational acceleration vector.
        :param refLevel: Level at which pressure = pRef.
        :param pRef: Reference value for the pressure at x[vert_axis]=refLevel,
                   be default set to 0.
        :param vert_axis: index of vertical in position vector, must always be
                        aligned with gravity, by default set to 1.
        :return: hydrostaticPressureOutlet except when the pressure and the
                 vof are defined. Then it returns the pressure and vof profile
                 based on the known depth.
        If the boundary is aligned with one of the main axes, sets the
        tangential velocity components to zero as well.
        (!) This condition is best used for boundaries and gravity aligned with
            one of the main axes.
        """
        self.reset()

        if vert_axis is None:
            vert_axis = self.Shape.Domain.nd - 1

        def hydrostaticPressureOutletWithDepth_p_dirichlet(x, t):
            if x[vert_axis] < seaLevel:
                a0 = pRef-rhoUp*g[vert_axis]*(refLevel-seaLevel)-rhoDown*g[vert_axis]*seaLevel
                a1 = rhoDown*g[vert_axis]
                return a0 + a1*x[vert_axis]

        def hydrostaticPressureOutletWithDepth_vof_dirichlet(x, t):
            if x[vert_axis] < seaLevel:
                return water

        self.hydrostaticPressureOutlet(rhoUp, g, refLevel, pRef, vert_axis,
                                       air)
        self.p_dirichlet = hydrostaticPressureOutletWithDepth_p_dirichlet
        self.vof_dirichlet = hydrostaticPressureOutletWithDepth_vof_dirichlet



# for regions

class RelaxationZone:
    def __init__(self, shape, zone_type, center, orientation, waves, windSpeed,
                 epsFact_solid, dragAlpha, dragBeta, porosity):
        self.Shape = shape
        self.zone_type = zone_type
        self.center = center
        self.orientation = orientation
        self.waves = waves
        self.windSpeed = windSpeed
        if zone_type == 'absorption' or zone_type == 'porous':
            self.u = self.v = self.w = lambda x, t: 0.
        elif zone_type == 'generation':
            self.u = self.setGenerationFunctions(0)
            self.v = self.setGenerationFunctions(1)
            self.w = self.setGenerationFunctions(2)
        else:
            log('Wrong zone type: ' + self.zone_type)
            sys.exit()
        self.epsFact_solid = epsFact_solid
        self.dragAlpha = dragAlpha
        self.dragBeta = dragBeta
        self.porosity = porosity

    def setGenerationFunctions(self, i):
        """
        Sets the functions necessary for generation zones
        """
        def twp_flowVelocity(x, t):
            from proteus import Context
            ct = Context.get()
            vert_axis = self.Shape.Domain.nd-1
            waveHeight = self.waves.mwl+self.waves.eta(x, t)
            wavePhi = x[vert_axis]-waveHeight
            if wavePhi <= 0:
                waterSpeed = self.waves.u(x, t)
            else:
                x_max = np.copy(x)
                x_max[vert_axis] = waveHeight
                waterSpeed = self.waves.u(x_max, t)
            he, ech = ct.domain.MeshOptions.he, ct.epsFact_consrv_heaviside
            H = smoothedHeaviside(0.5*ech*he, wavePhi-0.5*ech*he)
            return H*self.windSpeed[i] + (1-H)*waterSpeed[i]
        return twp_flowVelocity


class RelaxationZoneWaveGenerator(AuxiliaryVariables.AV_base):
    """
    Prescribe a velocity penalty scaling in a material zone via a
    Darcy-Forchheimer penalty
    :param zones: A dictionary mapping integer material types to Zones, where a
    Zone is a named tuple specifying the x coordinate of the zone center and
    the velocity components
    """
    def __init__(self, zones, nd):
        assert isinstance(zones, dict)
        self.zones = zones
        self.nd = nd

    def calculate(self):
        for l, m in enumerate(self.model.levelModelList):
            for eN in range(m.coefficients.q_phi.shape[0]):
                mType = m.mesh.elementMaterialTypes[eN]
                if mType in self.zones:
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN, k]
                        zone = self.zones[mType]
                        coeff = m.coefficients
                        ori = zone.orientation
                        nd = zone.Shape.Domain.nd
                        if zone.zone_type == 'porous':
                            coeff.q_phi_solid[eN, k] = zone.epsFact_solid
                        else:
                            coeff.q_phi_solid[eN, k] = np.dot(ori, zone.center[:nd]-x[:nd])
                        coeff.q_velocity_solid[eN, k, 0] = zone.u(x, t)
                        coeff.q_velocity_solid[eN, k, 1] = zone.v(x, t)
                        if self.nd > 2:
                            coeff.q_velocity_solid[eN, k, 2] = zone.w(x, t)
            m.q['phi_solid'] = m.coefficients.q_phi_solid
            m.q['velocity_solid'] = m.coefficients.q_velocity_solid

