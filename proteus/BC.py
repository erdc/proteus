"""
Module for creating boundary conditions. Imported in Shape.py
"""
import numpy as np


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
        self._b_or = b_or  # array of orientation of all facets/segments (boundaries) of shape
        self._b_i = b_i  # indice for this boundary in list of boundaries in shape
        self.reset()
        # moveMesh boundary conditions
        self.hx_dirichlet = None
        self.hy_dirichlet = None
        self.hz_dirichlet = None
        self.u_stress = 0.
        self.v_stress = 0.
        self.w_stress = 0.

    def reset(self):
        # _dirichlet
        self.p_dirichlet = None  # pressure
        self.u_dirichlet = None  # velocity u
        self.v_dirichlet = None  # velocity v
        self.w_dirichlet = None  # velocity w
        self.vof_dirichlet = None  # VOF
        self.k_dirichlet = None  # kappa
        self.dissipation_dirichlet = None  # dissipation
        # _advective
        self.p_advective = constantBC(0.)
        self.u_advective = constantBC(0.)
        self.v_advective = constantBC(0.)
        self.w_advective = constantBC(0.)
        self.vof_advective = constantBC(0.)
        self.k_advective = None
        self.dissipation_advective = None
        # _diffusive
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = None
        self.dissipation_diffusive = None

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

    def setFreeSlip(self):
        """
        sets free slip conditions at the boundary
        """
        self.reset()

    def setClosed(self):
        self.k_advective = constantBC(0.)
        self.k_diffusive = constantBC(0.)

    def setOpenAir(self):
        """
        sets open boundary conditions (water can come out)
        """
        self.p_dirichlet = constantBC(0.)
        self.u_dirichlet = constantBC(0.)
        self.v_dirichlet = constantBC(0.)
        self.w_dirichlet = constantBC(0.)
        self.vof_dirichlet = constantBC(1.)
        self.p_advective = None
        self.u_advective = None
        self.v_advective = None
        self.w_advective = None
        self.k_advective = constantBC(0.)
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)
        self.k_diffusive = constantBC(0.)
        self.dissipation_diffusive = constantBC(0.)

    def setObstacle(self):
        """
        sets rigid body boundary conditions
        """
        self.u_dirichlet = constantBC(0.)
        self.v_dirichlet = constantBC(0.)
        self.w_dirichlet = constantBC(0.)
        self.k_dirichlet = constantBC(0.)
        self.u_advective = None
        self.v_advective = None
        self.w_advective = None
        self.vof_advective = constantBC(0.)
        self.u_diffusive = None
        self.v_diffusive = None
        self.w_diffusive = None
        self.DFBC_d = constantBC(0.)

    def setMoveMesh(self, body):
        """
        sets rigid body boundary conditions for moving the mesh
        :arg last_position: position of (barycentre of) body before last calculation step
        :arg position: position of (barycentre of) body after last calculation step
        :arg last_rotation: rotation matrix of body before last calculation step
        :arg rotation: rotation matrix of body after last calculation step
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

    def setTwoPhaseVelocityInlet(self, U, waterLevel, vert_axis=-1, air=1., water=0.):
        """
        Imposes a velocity profile lower than the sea level and an open boundary for higher than the sealevel
        :arg U: Velocity vector at the global system
        :arg waterLevel: water level at global coordinate system
        :arg vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1]
        :arg air: Volume fraction for air (1.0 by default)
        :arg water: Volume fraction for water (0.0 by default)
        Below the seawater level, the condition returns the _dirichlet and p_advective condition according to the inflow velocity
        Above the sea water level, the condition returns the gravity as zero, and sets _dirichlet condition to zero, only if there is an 
        zero inflow velocity component
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.reset()
        U = np.array(U)

        def get_inlet_v_dirichletel(ux):
            def u_dirichletx(x, t):
                if x[vert_axis] < waterLevel:
                    return ux
                elif x[vert_axis] >= waterLevel and ux==0:
                    return 0.
            return u_dirichletx

        def inlet_vof_dirichlet(x, t):
            if x[vert_axis] < waterLevel:
                return water
            elif x[vert_axis] >= waterLevel:
                return air

        def inlet_p_advective(x, t, u=U):
            b_or = self._b_or[self._b_i]
            u_p = np.sum(U*b_or)
            # This is the normal velocity, based on the inwards boundary orientation -b_or
            u_p = -u_p
            if x[vert_axis] < waterLevel:
                return u_p
            elif x[vert_axis] >= waterLevel:
                return None

        self.u_dirichlet = get_inlet_v_dirichletel(U[0])
        self.v_dirichlet = get_inlet_v_dirichletel(U[1])
        if len(U) == 3:
                self.w_dirichlet = get_inlet_v_dirichletel(U[2])
        self.vof_dirichlet = inlet_vof_dirichlet
        self.p_advective = inlet_p_advective
        self.u_diffusive = constantBC(0.)
        self.v_diffusive = constantBC(0.)
        self.w_diffusive = constantBC(0.)

    def setHydrostaticPressureOutlet(self, rho, g, refLevel, pRef=0.0, vert_axis=-1, air=1.0):
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

    def hydrostaticPressureOutletWithDepth(self, seaLevel, rhoUp, rhoDown, g, refLevel, pRef=0.0, vert_axis=-1, air=1.0, water=0.0):
        """Imposes a hydrostatic pressure profile and open boundary conditions with a known otuflow depth
        :arg rhoUp: Phase density of the upper part
        :arg rhoDown: Phase density of the lower part
        :arg g: Gravitational acceleration vector
        :arg refLevel: Level at which pressure = pRef
        :arg pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
        :arg vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
        :return: hydrostaticPressureOutlet except when the pressure and the vof are defined. Then it returns the pressure and vof profile based on the known depth
        If the boundary is aligned with one of the main axes, sets the tangential velocity components to zero as well
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.reset()

        def hydrostaticPressureOutletWithDepth_p_dirichlet(x, t):
            if x[vert_axis] < seaLevel:
                a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
                a1 = rhoDown*g[vert_axis]
                return a0 + a1*x[vert_axis]

        def hydrostaticPressureOutletWithDepth_vof_dirichlet(x, t):
            if x[vert_axis] < seaLevel:
                a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
                a1 = rhoDown*g[vert_axis]
                return water

        def hydrostaticPressureOutletWithDepth_vof_dirichlet(x, t):
            if x[vert_axis] < seaLevel:
                return water

        self.hydrostaticPressureOutlet(rhoUp, g, refLevel, pRef, vert_axis, air)
        self.p_dirichlet = hydrostaticPressureOutletWithDepth_p_dirichlet
        self.vof_dirichlet = hydrostaticPressureOutletWithDepth_vof_dirichlet
