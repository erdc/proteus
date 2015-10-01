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

def noBC(x, t):
    """
    function returning None
    :arg x: coordinates of point
    :return: None
    """
    return None


def zeroBC(x, t):
    """
    function returning 0
    :arg x: coordinates of point
    :return: 0
    """
    return 0.0


class BoundaryConditions:
    """
    class defining boundary conditions for a facet or a segment
    """
    def __init__(self, b_or=None, b_i=None):
        self._b_or = b_or  # array of orientation of all facets/segments (boundaries) of shape
        self._b_i = b_i  # indice for this boundary in list of boundaries in shape
        self.reset()
        # moveMesh boundary conditions
        self.DBC_hx = noBC
        self.DBC_hy = noBC
        self.DBC_hz = noBC
        self.stress_u = zeroBC
        self.stress_v = zeroBC
        self.stress_w = zeroBC

    def reset(self):
        # Dirichlet
        self.DBC_p = noBC  # pressure
        self.DBC_u = noBC  # velocity u
        self.DBC_v = noBC  # velocity v
        self.DBC_w = noBC  # velocity w
        self.DBC_vof = noBC  # VOF
        self.DBC_k = noBC  # kappa
        self.DBC_dissipation = noBC  # dissipation
        # Advective
        self.AFBC_p = noBC
        self.AFBC_u = noBC
        self.AFBC_v = noBC
        self.AFBC_w = noBC
        self.AFBC_vof = noBC
        self.AFBC_k = noBC
        self.AFBC_dissipation = noBC
        # Diffusive
        self.DFBC_u = noBC
        self.DFBC_v = noBC
        self.DFBC_w = noBC
        self.DFBC_k = noBC
        self.DFBC_dissipation = noBC

    def setTank(self):
        b_or = self._b_or[self._b_i].tolist()
        if b_or[0] == 1 or b_or[0] == -1:
            self.DBC_hx = constantBC(0.)
            self.stress_u = noBC
        elif b_or[1] == 1 or b_or[1] == -1:
            self.DBC_hy = constantBC(0.)
            self.stress_v = noBC
        elif len(b_or) > 2 and (b_or[2] == 1 or b_or[2] == -1):
            self.DBC_hz = constantBC(0.)
            self.stress_w = noBC

    def setNoSlip(self):
        """
        sets no slip conditions at the boundary
        """
        self.reset()
        self.DBC_u = constantBC(0.)
        self.DBC_v = constantBC(0.)
        self.DBC_w = constantBC(0.)
        self.AFBC_p = constantBC(0.)
        self.AFBC_u = constantBC(0.)
        self.AFBC_v = constantBC(0.)
        self.AFBC_w = constantBC(0.)
        self.AFBC_vof = constantBC(0.)

    def setFreeSlip(self):
        """
        sets free slip conditions at the boundary
        """
        self.reset()
        self.AFBC_p = constantBC(0.)
        self.AFBC_u = constantBC(0.)
        self.AFBC_v = constantBC(0.)
        self.AFBC_w = constantBC(0.)
        self.AFBC_vof = constantBC(0.)
        self.DFBC_u = constantBC(0.)
        self.DFBC_v = constantBC(0.)
        self.DFBC_w = constantBC(0.)

    def setClosed(self):
        self.AFBC_k = constantBC(0.)
        self.DFBC_k = constantBC(0.)

    def setOpenAir(self):
        """
        sets open boundary conditions (water can come out)
        """
        self.DBC_p = constantBC(0.)
        self.DBC_vof = constantBC(1.)
        self.AFBC_k = constantBC(0.)
        self.DFBC_k = constantBC(0.)
        self.DFBC_dissipation = constantBC(0.)

    def setObstacle(self):
        """
        sets rigid body boundary conditions
        """
        self.reset()
        self.DBC_hx = constantBC(0.)  # initial mesh conditions
        self.DBC_hy = constantBC(0.)
        self.DBC_hz = constantBC(0.)
        self.DBC_u = constantBC(0.)
        self.DBC_v = constantBC(0.)
        self.DBC_w = constantBC(0.)
        self.DBC_k = constantBC(0.)
        self.AFBC_p = constantBC(0.)
        self.AFBC_vof = constantBC(0.)
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
        self.DBC_hx = get_DBC_h(i=0)
        self.DBC_hy = get_DBC_h(i=1)
        if len(body.last_position) > 2:
            self.DBC_hz = get_DBC_h(i=2)

    def setTwoPhaseVelocityInlet(self, U, waterLevel, vert_axis=-1, air=1., water=0.):
        """
        Imposes a velocity profile lower than the sea level and an open boundary for higher than the sealevel
        :arg U: Velocity vector at the global system
        :arg waterLevel: water level at global coordinate system
        :arg vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1]
        :arg air: Volume fraction for air (1.0 by default)
        :arg water: Volume fraction for water (0.0 by default)
        Below the seawater level, the condition returns the Dirichlet and pAdvective condition according to the inflow velocity
        Above the sea water level, the condition returns the gravity as zero, and sets Dirichlet condition to zero, only if there is an 
        zero inflow velocity component
        THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
        """
        self.reset()
        U = np.array(U)

        def get_inlet_DBC_vel(ux):
            def DBC_ux(x, t):
                if x[vert_axis] < waterLevel:
                    return ux
                elif x[vert_axis] >= waterLevel and ux==0:
                    return 0.
            return DBC_ux

        def inlet_DBC_vof(x, t):
            if x[vert_axis] < waterLevel:
                return water
            elif x[vert_axis] >= waterLevel:
                return air

        def inlet_AFBC_p(x, t, u=U):
            b_or = self._b_or[self._b_i]
            u_p = np.sum(U*b_or)
            # This is the normal velocity, based on the inwards boundary orientation -b_or
            u_p = -u_p
            if x[vert_axis] < waterLevel:
                return u_p
            elif x[vert_axis] >= waterLevel:
                return None

        self.DBC_u = get_inlet_DBC_vel(U[0])
        self.DBC_v = get_inlet_DBC_vel(U[1])
        if len(U) == 3:
                self.DBC_w = get_inlet_DBC_vel(U[2])
        self.DBC_vof = inlet_DBC_vof
        self.AFBC_p = inlet_AFBC_p
        self.DFBC_u = constantBC(0.)
        self.DFBC_v = constantBC(0.)
        self.DFBC_w = constantBC(0.)

    def setHydrostaticPressureOutlet(self, rho, g, refLevel, pRef=0.0, vert_axis=-1, air=1.0):
        self.reset()
        a0 = pRef - rho*g[vert_axis]*refLevel
        a1 = rho*g[vert_axis]
        # This is the normal velocity, based on the boundary orientation

        def get_outlet_DBC_vel(i):
            def DBC_ux(x, t):
                b_or = self._b_or[self._b_i]
                if b_or[i] == 0:
                    return 0.
        self.DBC_u = get_outlet_DBC_vel(0)
        self.DBC_v = get_outlet_DBC_vel(1)
        if len(g) == 3:
            self.DBC_w = get_outlet_DBC_vel(2)
        self.DBC_p = linearBC(a0, a1, vert_axis)
        self.DBC_vof = constantBC(air)
        self.DFBC_u = constantBC(0.)
        self.DFBC_v = constantBC(0.)
        self.DFBC_w = constantBC(0.)

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

        def hydrostaticPressureOutletWithDepth_DBC_p(x, t):
            if x[vert_axis] < seaLevel:
                a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
                a1 = rhoDown*g[vert_axis]
                return a0 + a1*x[vert_axis]

        def hydrostaticPressureOutletWithDepth_DBC_vof(x, t):
            if x[vert_axis] < seaLevel:
                a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
                a1 = rhoDown*g[vert_axis]
                return water

        def hydrostaticPressureOutletWithDepth_DBC_vof(x, t):
            if x[vert_axis] < seaLevel:
                return water

        self.hydrostaticPressureOutlet(rhoUp, g, refLevel, pRef, vert_axis, air)
        self.DBC_p = hydrostaticPressureOutletWithDepth_DBC_p
        self.DBC_vof = hydrostaticPressureOutletWithDepth_DBC_vof
