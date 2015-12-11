"""
Module for creating boundary conditions. Imported in Shape.py
"""
import numpy as np
from proteus import AuxiliaryVariables
from proteus.Profiling import logEvent as log
from proteus.ctransportCoefficients import (smoothedHeaviside,
                                            smoothedHeaviside_integral)



def constantBC(value):
    """
    function returning constant BC
    """
    return lambda x, t: value


def linearBC(a0, a1, i):
    """
    function returning linear BC
    :arg a0:
    :arg a1:
    :arg i:
    """
    return lambda x, t: a0 + a1*x[i]


class BoundaryConditions:
    """
    class defining boundary conditions for a facet or a segment
    """
    def __init__(self, b_or=None, b_i=None):
        self._b_or = b_or  # array of orientation of all boundaries of shape
        self._b_i = b_i  # indice for this boundary in list of boundaries
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

    def setParallelFlag0(self):
        """
        This function should be used only on BC flag 0 (index 0 in the list
        domain.bc)
        """
        self.vof_advective = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    def setTank(self):
        b_or = self._b_or[self._b_i].tolist()
        if b_or[0] == 1 or b_or[0] == -1:
            self.hx_dirichlet = constantBC(0.)
            self.u_stress = None
        elif b_or[1] == 1 or b_or[1] == -1:
            self.hy_dirichlet = constantBC(0.)
            self.v_stress = None
        elif len(b_or) > 2 and (b_or[2] == 1 or b_or[2] == -1):
            self.hz_dirichlet = constantBC(0.)
            self.w_stress = None

    def setNoSlip(self):
        """
        sets no slip conditions at the boundary
        """
        self.reset()
        self.u_dirichlet = constantBC(0.)
        self.v_dirichlet = constantBC(0.)
        self.w_dirichlet = constantBC(0.)
        self.k_dirichlet = constantBC(0.)
        self.dissipation_dirichlet = constantBC(0.)
        self.p_advective = constantBC(0.)
        self.u_advective = constantBC(0.)
        self.v_advective = constantBC(0.)
        self.w_advective = constantBC(0.)
        self.vof_advective = constantBC(0.)
        self.k_advective = constantBC(0.)
        self.dissipation_advective = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setFreeSlip(self):
        """
        sets free slip conditions at the boundary
        """
        self.reset()
        self.p_advective = constantBC(0.)
        self.u_advective = constantBC(0.)
        self.v_advective = constantBC(0.)
        self.w_advective = constantBC(0.)
        self.k_advective = constantBC(0.)
        self.dissipation_advective = constantBC(0.)
        self.vof_advective = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setClosed(self):
        self.k_advective = constantBC(0.)
        self.k_diffusive = constantBC(0.)

    def setOpenAir(self, orientation=None):
        """
        sets open boundary conditions (water can come out)
        """

        def get_ux_dirichlet(i):
            if b_or[i] == 1. or b_or[i] == -1.:
                return None
            else:
                return constantBC(0.)

        if orientation is None and self._b_or is not None:
            b_or = self._b_or[self._b_i]
        elif orientation is not None:
            b_or = orientation
        else:
            print('Boundary orientation needs to be defined')

        self.p_dirichlet = constantBC(0.)
        self.u_dirichlet = get_ux_dirichlet(0)
        self.v_dirichlet = get_ux_dirichlet(1)
        if len(b_or) > 2:
            self.w_dirichlet = get_ux_dirichlet(2)
        self.vof_dirichlet = constantBC(1.)  # air
        self.k_dirichlet = constantBC(0.)
        self.dissipation_dirichlet = constantBC(0.)
        self.p_advective = None
        self.u_advective = None
        self.v_advective = None
        self.w_advective = None
        self.vof_advective = None
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setMoveMesh(self, body):
        """
        sets rigid body boundary conditions for moving the mesh
        :arg last_position: position (barycentre) of body pre-calculation
        :arg position: position (barycentre) of body post-calculation
        :arg last_rotation: rotation matrix of body pre-calculation
        :arg rotation: rotation matrix of body post-calculation
        (!) should not be set manually
        """
        def get_DBC_h(i):
            def DBC_h(x, t):
                x_0 = x-body.last_position
                new_x_0 = np.dot(x_0, body.rotation_matrix)
                hx = new_x_0-x_0+body.h
                return hx[i]
            return DBC_h
        self.hx_dirichlet = get_DBC_h(i=0)
        self.hy_dirichlet = get_DBC_h(i=1)
        if len(body.last_position) > 2:
            self.hz_dirichlet = get_DBC_h(i=2)

    def setTwoPhaseVelocityInlet(self, U, waterLevel, vert_axis=-1, air=1.,
                                 water=0.):
        """
        Imposes a velocity profile lower than the sea level and an open
        boundary for higher than the sealevel.
        :arg U: Velocity vector at the global system.
        :arg waterLevel: water level at global coordinate system.
        :arg vert_axis: index of vertical in position vector, must always be
                        aligned with gravity, by default set to 1].
        :arg air: Volume fraction for air (1.0 by default).
        :arg water: Volume fraction for water (0.0 by default).
        Below the seawater level, the condition returns the _dirichlet and
        p_advective condition according to the inflow velocity.
        Above the sea water level, the condition returns the gravity as zero,
        and sets _dirichlet condition to zero, only if there is a zero inflow
        velocity component.
        (!) This condition is best used for boundaries and gravity aligned with
            one of the main axes.
        """
        self.reset()
        U = np.array(U)

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
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    def setHydrostaticPressureOutlet(self, rho, g, refLevel, pRef=0.0,
                                     vert_axis=-1, air=1.0):
        self.reset()
        a0 = pRef - rho*g[vert_axis]*refLevel
        a1 = rho*g[vert_axis]
        # This is the normal velocity, based on the boundary orientation

        def get_outlet_v_dirichletel(i):
            def u_dirichletx(x, t):
                b_or = self._b_or[self._b_i]
                if b_or[i] == 0:
                    return 0.
        self.u_dirichlet = get_outlet_v_dirichletel(0)
        self.v_dirichlet = get_outlet_v_dirichletel(1)
        if len(g) == 3:
            self.w_dirichlet = get_outlet_v_dirichletel(2)
        self.p_dirichlet = linearBC(a0, a1, vert_axis)
        self.vof_dirichlet = constantBC(air)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    def hydrostaticPressureOutletWithDepth(self, seaLevel, rhoUp, rhoDown, g,
                                           refLevel, pRef=0.0, vert_axis=-1,
                                           air=1.0, water=0.0):
        """Imposes a hydrostatic pressure profile and open boundary conditions
        with a known otuflow depth
        :arg rhoUp: Phase density of the upper part.
        :arg rhoDown: Phase density of the lower part.
        :arg g: Gravitational acceleration vector.
        :arg refLevel: Level at which pressure = pRef.
        :arg pRef: Reference value for the pressure at x[vert_axis]=refLevel,
                   be default set to 0.
        :arg vert_axis: index of vertical in position vector, must always be
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
    def __init__(self, center_x, sign, u, v, w):
        self.center_x = center_x
        self.sign = sign
        self.u = u
        self.v = v
        self.w = w


class RelaxationZoneWaveGenerator(AuxiliaryVariables.AV_base):
    """
    Prescribe a velocity penalty scaling in a material zone via a
    Darcy-Forchheimer penalty
    :param zones: A dictionary mapping integer material types to Zones, where a
    Zone is a named tuple specifying the x coordinate of the zone center and
    the velocity components
    """
    def __init__(self, zones, shape, waves=None, wind=None):
        assert isinstance(zones, dict)
        self.zones = zones
        self.shape = shape
        self.waves = waves
        self.windVel = wind
        shape.auxiliaryVariables += [self]
                
    def setGenerationFunctions(self, i):
        """
        Sets the functions necessary for generation zones
        """
        if i == 0:
            s = 'x'
        elif i == 1:
            s = 'y'
        elif i == 2:
            s = 'z'
        def twp_flowVelocity(x, t):
            vert_axis = self.shape.domain.nd - 1
            waterSpeed = self.waves.u(x[0], x[1], x[2], t, s)[i]
            waveHeight = self.waves.eta(x[0], x[1], x[2], t)
            wavePhi = x[vert_axis] - waveHeight
            he = self.shape.domain.Mesh.he
            # !!!!!!!!!!!!!!!!!!!!!!!
            # epsFact_consrv_heaviside should be called from context!
            # !!!!!!!!!!!!!!!!!!!!!!!
            epsFact_consrv_heaviside = 3.0
            H = smoothedHeaviside(epsFact_consrv_heaviside*he,
                                wavePhi-epsFact_consrv_heaviside*he)
            return H*self.windVel[0] + (1-H)*waterSpeed
        return twp_flowVelocity

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
                        sign = zone.sign
                        coeff.q_phi_solid[eN, k] = sign*(zone.center_x-x[0])
                        coeff.q_velocity_solid[eN, k, 0] = zone.u(x, t)
                        coeff.q_velocity_solid[eN, k, 1] = zone.v(x, t)
                        if self.shape.domain.nd > 2:
                            coeff.q_velocity_solid[eN, k, 2] = zone.w(x, t)
            m.q['phi_solid'] = m.coefficients.q_phi_solid
            m.q['velocity_solid'] = m.coefficients.q_velocity_solid

# --------------------------------------------------------------------------- #
# --------------------------- INITIAL CONDITIONS ---------------------------- #
# --------------------------------------------------------------------------- #

# from proteus.ctransportCoefficients import (smoothedHeaviside,
#                                             smoothedHeaviside_integral)

# class InitialConditions:

#     def __init__(self, domain):
#         self.domain = domain
#         self.nd = domain.nd
#         self.initial_condition_func = 0.

#     def twpflowPressure_init(self, x, t, waterLevel, ceiling=None,
#                              axis=None, rho_0=998.2, rho_1=1.205, p_L=0.):
#         if axis is None:
#             axis = self.nd-1
#         if ceiling is None:
#             ceiling = self.domain[axis]
#         p_L = p_L
#         phi_L = ceiling - waterLevel
#         phi = x[axis] - waterLevel
#         return p_L -g[axis]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
#                                                             -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

#     def uOfXT(self, x, t):
#         return self.initial_condition_func




# class P_IC:
#     def uOfXT(self, x, t):
#         return ct.twpflowPressure_init(x, t)

# class U_IC:
#     def uOfXT(self, x, t):
#         return 0.0

# class V_IC:
#     def uOfXT(self, x, t):
#         return 0.0

# class W_IC:
#     def uOfXT(self, x, t):
#         return 0.0
