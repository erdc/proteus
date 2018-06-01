import numpy as np
from numpy import linalg as LA

# return the force exerted on i-th particle by the j-th particle
def get_repulsive_force_glowinski(Xi, Ri, Xj, Rj, force_range, stiffness_parameter):
    
    dij = LA.norm(Xi-Xj, 2)
    Fij = np.zeros_like(Xi)
    if dij <= Ri+Rj+force_range:
        Fij[:] = Xi - Xj
        Fij *= (Ri+Rj+force_range-dij)**2/stiffness_parameter

    return Fij

def get_repulsive_force_from_vetical_wall_glowinski(Xi, Ri, xl, force_range, stiffness_parameter):
    
    Xj = np.copy(Xi)
    Xj[0] = xl - (Xj[0]-xl)
    Rj = Ri
    return get_repulsive_force_glowinski(Xi, Ri, Xj, Rj, force_range, stiffness_parameter)

def get_repulsive_force_from_horizontal_wall_glowinski(Xi, Ri, xd, force_range, stiffness_parameter):
    
    Xj = np.copy(Xi)
    Xj[1] = xd - (Xj[1]-xd)
    Rj = Ri
    return get_repulsive_force_glowinski(Xi, Ri, Xj, Rj, force_range, stiffness_parameter)
