#!python
#cython: profile=True

"""
Module for creating boundary conditions. Imported in mprans.SpatialTools.py
"""
import sys
import numpy as np
from proteus import AuxiliaryVariables
from proteus.BoundaryConditions import (BC_Base,
                                        BoundaryCondition,
                                        # constantBC,
                                        # linearBC
                                        )
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)


class BC_RANS(BC_Base):
    """
    Class regrouping boundary conditions for two-phase flows
    """
    def __init__(self, shape=None, name=None, b_or=None, b_i=None):
        super(BC_RANS, self).__init__(shape, name, b_or, b_i)
        # _dirichlet
        self.p_dirichlet = BoundaryCondition()  # pressure
        self.u_dirichlet = BoundaryCondition()  # velocity u
        self.v_dirichlet = BoundaryCondition()  # velocity v
        self.w_dirichlet = BoundaryCondition()  # velocity w
        self.vof_dirichlet = BoundaryCondition()  # VOF
        self.k_dirichlet = BoundaryCondition()  # kappa
        self.dissipation_dirichlet = BoundaryCondition()  # dissipation
        # _advective
        self.p_advective = BoundaryCondition()
        self.u_advective = BoundaryCondition()
        self.v_advective = BoundaryCondition()
        self.w_advective = BoundaryCondition()
        self.vof_advective = BoundaryCondition()
        self.k_advective = BoundaryCondition()
        self.dissipation_advective = BoundaryCondition()
        # _diffusive
        self.u_diffusive = BoundaryCondition()
        self.v_diffusive = BoundaryCondition()
        self.w_diffusive = BoundaryCondition()
        self.k_diffusive = BoundaryCondition()
        self.dissipation_diffusive = BoundaryCondition()
        # moveMesh boundary conditions
        self.hx_dirichlet = BoundaryCondition()
        self.hy_dirichlet = BoundaryCondition()
        self.hz_dirichlet = BoundaryCondition()
        self.u_stress = BoundaryCondition()
        self.u_stress.uOfXT = 0.
        self.v_stress = BoundaryCondition()
        self.v_stress.uOfXT = 0.
        self.w_stress = BoundaryCondition()
        self.w_stress.uOfXT = 0.

    def reset(self):
        """
        Resets all BoundaryCondtion functions to None, apart from the BCs
        affecting: moving mesh
        """
        self.BC_type = 'None'
        self.p_dirichlet.uOfXT = None
        self.u_dirichlet.uOfXT = None
        self.v_dirichlet.uOfXT = None
        self.w_dirichlet.uOfXT = None
        self.vof_dirichlet.uOfXT = None
        self.k_dirichlet.uOfXT = None
        self.dissipation_dirichlet.uOfXT = None
        self.p_advective.uOfXT = None
        self.u_advective.uOfXT = None
        self.v_advective.uOfXT = None
        self.w_advective.uOfXT = None
        self.vof_advective.uOfXT = None
        self.k_advective.uOfXT = None
        self.dissipation_advective.uOfXT = None
        self.u_diffusive.uOfXT = None
        self.v_diffusive.uOfXT = None
        self.w_diffusive.uOfXT = None
        self.k_diffusive.uOfXT = None
        self.dissipation_diffusive.uOfXT = None

    def setNonMaterial(self):
        """
        Sets non-material boundary conditions (diffusive flux and advective vof
        to 0.).
        """
        self.reset()
        self.BC_type = 'NonMaterial'
        self.vof_advective.setConstantBC(0.)
        self.u_diffusive.setConstantBC(0.)
        self.v_diffusive.setConstantBC(0.)
        self.w_diffusive.setConstantBC(0.)

    def setTank(self):
        b_or = self._b_or[self._b_i]
        if b_or[0] == 1 or b_or[0] == -1:
            self.hx_dirichlet.setConstantBC(0.)
            self.u_stress.uOfXT = None
        elif b_or[1] == 1 or b_or[1] == -1:
            self.hy_dirichlet.setConstantBC(0.)
            self.v_stress.uOfXT = None
        elif len(b_or) > 2 and (b_or[2] == 1 or b_or[2] == -1):
            self.hz_dirichlet.setConstantBC(0.)
            self.w_stress.uOfXT = None

    def setFixedNodes(self):
        """
        For moving domains: fixes nodes/boundary
        """
        self.hx_dirichlet.setConstantBC(0.)
        self.hy_dirichlet.setConstantBC(0.)
        self.hz_dirichlet.setConstantBC(0.)
        self.u_stress.uOfXT = None
        self.v_stress.uOfXT = None
        self.w_stress.uOfXT = None

    def setNoSlip(self):
        """
        Sets no slip conditions at the boundary
        """
        self.reset()
        self.BC_type = 'NoSlip'
        self.u_dirichlet.setConstantBC(0.)
        self.v_dirichlet.setConstantBC(0.)
        self.w_dirichlet.setConstantBC(0.)
        self.p_advective.setConstantBC(0.)
        self.vof_advective.setConstantBC(0.)
        self.k_dirichlet.setConstantBC(0.)
        self.dissipation_diffusive.setConstantBC(0.)

    def setFreeSlip(self):
        """
        Sets free slip conditions at the boundary
        """
        self.reset()
        self.BC_type = 'FreeSlip'
        self.p_advective.setConstantBC(0.)
        self.u_advective.setConstantBC(0.)
        self.v_advective.setConstantBC(0.)
        self.w_advective.setConstantBC(0.)
        self.vof_advective.setConstantBC(0.)
        self.k_dirichlet.setConstantBC(0.)
        self.u_diffusive.setConstantBC(0.)
        self.v_diffusive.setConstantBC(0.)
        self.w_diffusive.setConstantBC(0.)
        self.dissipation_diffusive.setConstantBC(0.)

    def setAtmosphere(self, orientation=None, vof_air=1.):
        """
        Sets atmosphere boundary conditions (water can come out)
        (!) pressure dirichlet set to 0 for this BC

        Parameters
        ----------
        orientation: Optional[array_like]
            orientation of the boundary. Optional if orientation was already
            passed when creating the BC_RANS class instance.
        vof_air: Optional[float]
            VOF value of air (default is 1.)
        """
        self.BC_type = 'OpenAir'

        def get_ux_dirichlet(i):
            if b_or[i] == 1. or b_or[i] == -1.:
                return None
            else:
                return lambda x, t: 0.

        if orientation is None and self._b_or[self._b_i] is not None:
            b_or = self._b_or[self._b_i]
        elif orientation is not None:
            b_or = orientation
        else:
            raise ValueError('Boundary orientation needs to be defined')
        self.reset()
        self.p_dirichlet.setConstantBC(0.)
        self.u_dirichlet.uOfXT = get_ux_dirichlet(0)
        self.v_dirichlet.uOfXT = get_ux_dirichlet(1)
        if len(b_or) > 2:
            self.w_dirichlet.uOfXT = get_ux_dirichlet(2)
        self.vof_dirichlet.setConstantBC(vof_air)  # air
        self.u_diffusive.setConstantBC(0.)
        self.v_diffusive.setConstantBC(0.)
        self.w_diffusive.setConstantBC(0.)
        self.k_diffusive.setConstantBC(0.)
        self.dissipation_diffusive.setConstantBC(0.)

    def setMoveMesh(self, last_pos, h=(0., 0., 0.), rot_matrix=None):
        """
        Sets boundary conditions for moving the mesh with a rigid body

        Parameters
        ----------
        last_pos: array_like
            last position of rigig body
        h: array_like
            displacement of the body
        rot_matrix:
            rotation matrix describing displament due to rotation between last
            position and new position (3x3 array)

        (!) if set manually, the input arrays should be updated externally
            without loosing their memory address
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
        self.hx_dirichlet.uOfXT = get_DBC_h(0)
        self.hy_dirichlet.uOfXT = get_DBC_h(1)
        if len(last_pos) > 2:
            self.hz_dirichlet.uOfXT = get_DBC_h(2)

    def setUnsteadyTwoPhaseVelocityInlet(self, wave, vert_axis=None,
                                         wind_speed=(0., 0., 0.), vof_air=1.,
                                         vof_water=0.):
        """
        Imposes a velocity profile on the fluid with input wave and wind
        conditions.

        Parameters
        ----------
        wave: proteus.WaveTools
            class describing a wave (from proteus.WaveTools)
        vert_axis: Optional[int]
            index of vertical position vector (x:0, y:1, z:2), must always be
            aligned with gravity. If not set, will be 1 in 2D (y), 2 in 3D (z).
        wind_speed: Optional[array_like]
        vof_air: Optional[float]
            VOF value of air (default is 1.)
        vof_water: Optional[float]
            VOF value of water (default is 0.)

        Below the sea water level: fluid velocity to wave speed.
        Above the sea water level: fluid velocity set to wind speed.
        (!) Boundary condition relies on specific variables defined in Context:
            he (mesh element size) and ecH (number of elements for smoothing)
        """
        self.reset()

        if vert_axis is None:
            vert_axis = self.Shape.Domain.nd-1

        def get_ux_dirichlet_cython(i):
            def ux_dirichlet_cython():
                wave_mwl = wave.mwl
                wave_eta = wave.eta
                wave_u = wave.u
                wind_speed_arr = np.array(wind_speed)
                he = self.ct.he
                ecH = self.ct.ecH
                def ux_dirichlet(x, t):
                    waveHeight = wave_mwl+wave_eta(x, t)
                    wavePhi = x[vert_axis]-waveHeight
                    if wavePhi <= 0:
                        water_speed = wave_u(x, t)
                    elif wavePhi > 0 and wavePhi < 0.5*ecH*he:
                        x_max = list(x)
                        x_max[vert_axis] = waveHeight
                        water_speed = wave_u(x_max, t)
                    else:
                        water_speed = np.array([0., 0., 0.])
                        # smoothing only above wave, only on half the VOF smoothing length
                    H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
                    ux = H*wind_speed_arr + (1-H)*water_speed
                    return ux[i]
                return ux_dirichlet
            return ux_dirichlet_cython

        def vof_dirichlet_cython():
            wave_mwl = wave.mwl
            wave_eta = wave.eta
            he = self.ct.he
            ecH = self.ct.ecH
            def vof_dirichlet(x, t):
                level = wave_mwl + wave_eta(x,t)
                H = smoothedHeaviside(ecH*he,x[vert_axis]-level)
                return H*vof_air+(1-H)*vof_water
            return vof_dirichlet

        def p_advective_cython():
            # This is the normal velocity, based on the outwards boundary
            # orientation b_or
            # needs to be equal to -ux_dirichlet
            b_or = self._b_or[self._b_i]
            cdef int nd = len(b_or)
            wave_mwl = wave.mwl
            wave_eta = wave.eta
            wave_u = wave.u
            wind_speed_arr = np.array(wind_speed)
            he = self.ct.he
            ecH = self.ct.ecH
            def p_advective(x, t):
                waveHeight = wave_mwl+wave_eta(x, t)
                wavePhi = x[vert_axis]-waveHeight
                if wavePhi <= 0:
                    water_speed = wave_u(x, t)
                elif wavePhi > 0 and wavePhi < 0.5*ecH*he:
                    x_max = list(x)
                    x_max[vert_axis] = waveHeight
                    water_speed = wave_u(x_max, t)
                else:
                    water_speed = np.array([0., 0., 0.])
                H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
                U = H*wind_speed_arr + (1-H)*water_speed
                u_p = np.sum(U[:nd]*b_or)
                return u_p
            return p_advective

        self.u_dirichlet.init_cython = get_ux_dirichlet_cython(0)
        self.v_dirichlet.init_cython = get_ux_dirichlet_cython(1)
        self.w_dirichlet.init_cython = get_ux_dirichlet_cython(2)
        self.vof_dirichlet.init_cython = vof_dirichlet_cython
        self.p_advective.init_cython = p_advective_cython

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

        def get_inlet_ux_dirichlet_cython(i):
            def get_inlet_ux_dirichlet():
                def ux_dirichlet(x, t):
                    if x[vert_axis] < waterLevel:
                        return U[i]
                    elif x[vert_axis] >= waterLevel and U[i] == 0:
                        return 0.
                return ux_dirichlet
            return get_inlet_ux_dirichlet

        def inlet_vof_dirichlet_cython():
            def inlet_vof_dirichlet(x, t):
                if x[vert_axis] < waterLevel:
                    return water
                elif x[vert_axis] >= waterLevel:
                    return air
            return inlet_vof_dirichlet

        def inlet_p_advective_cython():
            def inlet_p_advective(x, t):
                b_or = self._b_or[self._b_i]
                u_p = np.sum(U * b_or)
                # This is the normal velocity, based on the inwards boundary
                # orientation -b_or
                u_p = -u_p
                if x[vert_axis] < waterLevel:
                    return u_p
                elif x[vert_axis] >= waterLevel:
                    return None
            return inlet_p_advective

        self.u_dirichlet.init_cython = get_inlet_ux_dirichlet_cython(0)
        self.v_dirichlet.init_cython = get_inlet_ux_dirichlet_cython(1)
        if len(U) == 3:
                self.w_dirichlet.init_cython = get_inlet_ux_dirichlet_cython(2)
        self.vof_dirichlet.init_cython = inlet_vof_dirichlet_cython
        self.p_advective.init_cython = inlet_p_advective_cython

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
        self.u_dirichlet.uOfXT = get_outlet_ux_dirichlet(0)
        self.v_dirichlet.uOfXT = get_outlet_ux_dirichlet(1)
        if len(g) == 3:
            self.w_dirichlet.uOfXT = get_outlet_ux_dirichlet(2)
        self.p_dirichlet.setLinearBC(a0, a1, vert_axis)
        self.vof_dirichlet.setConstantBC(vof)
        self.u_diffusive.setConstantBC(0.)
        self.v_diffusive.setConstantBC(0.)
        self.w_diffusive.setConstantBC(0.)

    # FOLLOWING BOUNDARY CONDITION IS UNTESTED #
    def setHydrostaticPressureOutletWithDepth(self, seaLevel, rhoUp, rhoDown, g,
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

        self.setHydrostaticPressureOutlet(rhoUp, g, refLevel, pRef, vert_axis)
        self.p_dirichlet.uOfXT = hydrostaticPressureOutletWithDepth_p_dirichlet
        self.vof_dirichlet.uOfXT = hydrostaticPressureOutletWithDepth_vof_dirichlet



# for regions

class RelaxationZone:
    """
    Holds information about a relaxation zone (wave generation/absorption
    or porous zone)

    Parameters
    ----------
    zone_type: string
        type of zone, can be set to 'absorption', 'generation', or 'porous'
    center: array_like
        coordinates of center of the zone
    orientation: array_like
        orientation for absorption/generation zones: from boundary to tank
    epsFact_solid: float
        half the zone length
    waves: Optional[proteus.WaveTools]
        class instance of a wave from proteus.WaveTools (must be set for
        generation zone)
    shape: Optional[proteus.SpatialTools.Shape]
        shape class instance containing region
    dragAlpha: Optional[float]
        parameter for porous zones (default: 0.5/1.005e-6)
    dragBeta: Optional[float]
        parameter for porous zones (default: 0.)
    porosity: Optional[float]
        parameter for porous zone (default: 1.)
    """
    def __init__(self, zone_type, center, orientation, epsFact_solid,
                 waves=None, shape=None, wind_speed=(0.,0.,0.),
                 dragAlpha=0.5/1.005e-6, dragBeta=0., porosity=1.):
        self.Shape = shape
        self.zone_type = zone_type
        self.center = center
        self.orientation = orientation
        self.waves = waves
        self.wind_speed = wind_speed
        self.epsFact_solid = epsFact_solid
        self.dragAlpha = dragAlpha
        self.dragBeta = dragBeta
        self.porosity = porosity


class RelaxationZoneWaveGenerator(AuxiliaryVariables.AV_base):
    """
    Prescribe a velocity penalty scaling in a material zone via a
    Darcy-Forchheimer penalty

    Parameters
    ----------
    zones: dict
        dictionary with key as the region flag and values as a RelaxationZone
        class
    nd: int
        number of dimensions of domain
    """
    def __init__(self, zones, nd):
        assert isinstance(zones, dict)
        self.zones = zones
        self.nd = nd

    def calculate_init(self):
        from proteus import Context
        ct = Context.get()
        for key, zone in self.zones.iteritems():
            #print zone #[temp] a loose print statement - may be useful, but supressed for now
            if zone.zone_type == 'absorption' or zone.zone_type == 'porous':
                zone.u = zone.v = zone.w = lambda x, t: 0.
            elif zone.zone_type == 'generation':
                waves_u = zone.waves.u
                waves_mwl = zone.waves.mwl
                waves_eta = zone.waves.eta
                wind_speed = zone.wind_speed
                he = ct.he
                ecH = ct.ecH
                vert_axis = zone.Shape.Domain.nd-1
                def get_twp_flowVelocity(i):
                    def twp_flowVelocity(x, t):
                        waveHeight = waves_mwl+waves_eta(x, t)
                        wavePhi = x[vert_axis]-waveHeight
                        if wavePhi <= 0:
                            waterSpeed = waves_u(x, t)
                            H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
                        elif wavePhi > 0 and wavePhi < 0.5*ecH*he:
                            x_max = np.copy(x)
                            x_max[vert_axis] = waveHeight
                            waterSpeed = waves_u(x_max, t)
                        else:
                            waterSpeed = (0., 0., 0.)
                        H = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
                        return H*wind_speed[i] + (1-H)*waterSpeed[i]
                    return twp_flowVelocity
                zone.u = get_twp_flowVelocity(0)
                zone.v = get_twp_flowVelocity(1)
                zone.w = get_twp_flowVelocity(2)



    def calculate(self):
        cdef int nl = len(self.model.levelModelList)
        def iterate():
            m = self.model.levelModelList[l]
            cdef int nE = m.coefficients.q_phi.shape[0]
            cdef int nk = m.coefficients.q_phi.shape[1]
            for eN in range(nE):
                mType = m.mesh.elementMaterialTypes[eN]
                if mType in self.zones:
                    for k in range(nk):
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
        for l in range(nl):
            iterate()

