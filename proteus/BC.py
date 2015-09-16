"""
Module for creating boundary conditions. Imported in Shape.py
"""



def constantBC(value):
    if value is None:
        return None
    else:
        return lambda x, t: value

def linearBC(a0, a1, i):
    return lambda x, t: a0 + a1*x[i]

class BoundaryConditions:
    def __init__(self):
        # Dirichlet
        self.DBC_p = None  # pressure
        self.DBC_u = None  # velocity u
        self.DBC_v = None  # velocity v
        self.DBC_w = None  # velocity w
        self.DBC_vof = None  # VOF
        self.DBC_k = None  # kappa
        self.DBC_dissipation = None
        # Advective
        self.AFBC_p = None
        self.AFBC_u = None
        self.AFBC_v = None
        self.AFBC_w = None
        self.AFBC_vof = None
        self.AFBC_k = None  # weak dirichlet / outflow
        self.AFBC_dissipation = None
        # Diffusive
        self.DFBC_u = None
        self.DFBC_v = None
        self.DFBC_w = None
        self.DFBC_k = None  # weak dirichlet
        self.DFBC_dissipation = constantBC(0.)

    def setNoSlip(self):
        self.DBC_u = constantBC(0.)
        self.DBC_v = constantBC(0.)
        self.DBC_w = constantBC(0.)
        self.DBC_vof = constantBC(1.)
        self.AFBC_p = constantBC(0.)  ############
        self.AFBC_u = constantBC(0.)
        self.AFBC_v = constantBC(0.)
        self.AFBC_w = constantBC(0.)

    def setFreeSlip(self):
        self.AFBC_p = constantBC(0.)
        self.AFBC_u = constantBC(0.)
        self.AFBC_v = constantBC(0.)
        self.AFBC_w = constantBC(0.)
        self.DBC_vof = constantBC(1.)

    def setClosed(self):
        self.AFBC_k = constantBC(0.)
        self.DFBC_k = constantBC(0.)
        
    def setOpen(self):
        self.DBC_k = 'inflow'
        self.DFBC_dissipation = 'inflow'
        self.DFBC_k = constantBC(0.)
        
    def setInflow(self):
        self.DBC_p = None
        self.DBC_k = 'inflow'
        self.DBC_dissipation = 'inflow'  ########## to modify CHECK dissipation_p.py
        self.DFBC_u = constantBC(0.)
        self.DFBC_v = constantBC(0.)
        self.DFBC_w = constantBC(0.)
        self.DFBC_dissipation = None  # weak dirichlet

    def setObstacle(self):
        self.DBC_u = constantBC(0.)
        self.DBC_v = constantBC(0.)
        self.DBC_w = constantBC(0.)
        self.DBC_k = constantBC(0.)


#########################################################################
# following functions taken and modified from other script, not tested yet
#########################################################################

    # def twoPhaseVelocityInlet(self, x, U, seaLevel, b_or, vert_axis=1, air=1.0, water=0.0):
    #     """Imposes a velocity profile lower than the sea level and an open boundary for higher than the sealevel
    #     Takes the following arguments
    #     BCType: Type of boundary condition 
    #     x: Position vector
    #     U: Velocity vector at the global system
    #     seaLevel: water level at global coordinate system
    #     b_or is the orientation of the boundary (always positive outwards vector)
    #     vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1]
    #     air: Volume fraction for air (1.0 by default)
    #     water: Volume fraction for water (
    #     Below the seawater level, the condition returns the Dirichlet and pAdvective condition according to the inflow velocity
    #     Above the sea water level, the condition returns the gravity as zero, and sets Dirichlet condition to zero, only if there is an 
    #     zero inflow velocity component
    #     THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
    #     """
    #     if len(U) == 3:
    #         u, v, w = map(float, U)
    #         u_p = u*b_or[0]+v*b_or[1]+w*b_or[2]
    #     elif len(U) == 2:
    #         u, v = map(float, U)
    #         w = None
    #         u_p = u*b_or[0]+v*b_or[1]
    #     # This is the normal velocity, based on the inwards boundary orientation -b_or
    #     u_p = -u_p
    #     if x[vert_axis] < seaLevel:
    #         self.DBC_u = constantBC(u)
    #         self.DBC_v = constantBC(v)
    #         self.DBC_w = constantBC(w)
    #         self.DBC_vof = constantBC(water)
    #         self.AFBC_p = constantBC(u_p)
    #     elif x[vert_axis] >= seaLevel:
    #         if u == 0:
    #             self.DBC_u = constantBC(0.)
    #         if v == 0:
    #             self.DBC_v = constantBC(0.)
    #         if w == 0:
    #             self.DBC_w = constantBC(0.)
    #         self.DBC_vof = constantBC(air)


#     def hydrostaticPressureOutlet(self,BCType,rho,g,refLevel,b_or,pRef=0.0,vert_axis=1,air = 1.0):
#         """Imposes a hydrostatic pressure profile and  open boundary conditions
#         Takes the following arguments
#         BCType: Type of boundary condition 
#         x: Position vector
#         rho: Phase density
#         g: Gravitational acceleration vector
#         refLevel: Level at which pressure = pRef
#         pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
#         b_or is the orientation of the boundary (always positive outwards vector)
#         vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
#         Returns a hydrostatic profile for pDirichlet and zero Diffusive conditions. 
#         If the boundary is aligned with one of the main axes, sets the tangential velocity components to zero as well
#         THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
#         """
#         a0 = pRef - rho*g[vert_axis]*refLevel
#         a1 = rho*g[vert_axis]
#         # This is the normal velocity, based on the boundary orientation
#         if b_or[0] == 0:
#             self.DBC_u = constantBC(0.)
#         if b_or[1] == 0:
#             self.DBC_v = constantBC(0.)
#         if len(b_or) > 2 and b_or[2] == 0:
#             self.DBC_w = constantBC(0.)
#         self.DBC_p = self.linear(a0, a1, vert_axis)
#         self.DBC_vof = constantBC(air)
#         self.DFBC_u = constantBC(0.)
#         self.DFBC_v = constantBC(0.)
#         if len(b_or) > 2:
#             self.DFBC_w = constantBC(0.)

#     def hydrostaticPressureOutletWithDepth(self,BCType,x,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef=0.0,vert_axis=1,air=1.0,water=0.0):
#         """Imposes a hydrostatic pressure profile and open boundary conditions with a known otuflow depth
#         Takes the following arguments
#         BCType: Type of boundary condition 
#         x: Position vector
#         rhoUp: Phase density of the upper part
#         rhoDown: Phase density of the lower part
#         g: Gravitational acceleration vector
#         refLevel: Level at which pressure = pRef
#         pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
#         b_or is the orientation of the boundary (always positive outwards vector)
#         vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
#         Returns the hydrostaticPressureOutlet except when the pressure and the vof are defined. Then it returns the pressure 
#         and vof profile based on the known depth
#         If the boundary is aligned with one of the main axes, sets the tangential velocity components to zero as well
#         THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
#         """
#         self.hydrostaticPressureOutlet(BCType,rhoUp,g,refLevel,b_or,pRef,vert_axis,air)
#         if x[vert_axis] < seaLevel:
#             a0 = pRef - rhoUp*g[vert_axis]*(refLevel - seaLevel) - rhoDown*g[vert_axis]*seaLevel
#             a1 = rhoDown*g[vert_axis]
#             self.DBC_p = self.linear(a0, a1, vert_axis)
#             self.DBC_vof = constantBC(water)

#     def forcedOutlet(self,BCType,x,U,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef=0.0,vert_axis=1,air=1.,water=0.):
#         """Imposes a known velocit & pressure profile at the outflow
#         Takes the following arguments
#         BCType: Type of boundary condition 
#         BCType: Type of boundary condition 
#         x: Position vector
#         rhoUp: Phase density of the upper part
#         rhoDown: Phase density of the lower part
#         g: Gravitational acceleration vector
#         refLevel: Level at which pressure = pRef
#         pRef: Reference value for the pressure at x[vert_axis]=refLevel, be default set to 0
#         b_or is the orientation of the boundary (always positive outwards vector)
#         vert_axis: index of vertical in position vector, must always be aligned with gravity, by default set to 1
#         Returns only the dirichlet condition usinghydrostaticPressureOutletWithDepth for pressure and twoPhaseVelocityInlet 
#         for velocity, wuth inversed velocity
#         THIS CONDITION IS BEST USED FOR BOUNDARIES AND GRAVITY ALIGNED WITH ONE OF THE MAIN AXES
#         """
#         self.hydrostaticPressureOutletWithDepth(BCType,x,seaLevel,rhoUp,rhoDown,g,refLevel,b_or,pRef,vert_axis,air,water)
#         self.twoPhaseVelocityInlet(BCType,x,U,seaLevel,b_or,vert_axis,air,water)  # AFBC_p defined here!

######################################
#     For twp_navier_stokes_p.py     #
######################################

# from proteus import Context
# ct = Context.get()

# def getDBC_p(x, flag):
#     return ct.domain.bc[flag].DBC_p

# def getDBC_u(x, flag):
#     return ct.domain.bc[flag].DBC_u

# def getDBC_v(x, flag):
#     return ct.domain.bc[flag].DBC_v

# def getDBC_w(x, flag):
#     return ct.domain.bc[flag].DBC_w

# def getAFBC_p(x, flag):
#     return ct.domain.bc[flag].AFBC_p

# def getAFBC_u(x, flag):
#     return ct.domain.bc[flag].AFBC_u

# def getAFBC_v(x, flag):
#     return ct.domain.bc[flag].AFBC_v

# def getAFBC_w(x, flag):
#     return ct.domain.bc[flag].AFBC_w

# def getDFBC_u(x, flag):
#     return ct.domain.bc[flag].DFBC_u

# def getDFBC_v(x, flag):
#     return ct.domain.bc[flag].DFBC_v

# def getDFBC_w(x, flag):
#     return ct.domain.bc[flag].DFBC_w

# dirichletConditions = {0:getDBC_p,
#                        1:getDBC_u,
#                        2:getDBC_v,
#                        3:getDBC_w}

# advectiveFluxBoundaryConditions =  {0:getAFBC_p,
#                                     1:getAFBC_u,
#                                     2:getAFBC_v
#                                     3:getAFBC_w}

# diffusiveFluxBoundaryConditions = {0:{},
#                                    1:{1:getDFBC_u},
#                                    2:{2:getDFBC_v}}
#                                    3:{3:getDFBC_w}}

#####################################################
