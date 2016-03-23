
"""
Testing module for proteus.BoundaryConditions
Work in progress

TO ADD TO TESTS:
setTwoPhaseVelocityInlet()
setHydrostaticPressureOutlet()
hydrostaticPressureOutletWithDepth()
"""


import os, sys
import random
import unittest
import numpy.testing as npt
import numpy as np
from proteus import Comm, Profiling
from proteus.Profiling import logEvent as log
from proteus.BoundaryConditions import (BC_Base,
                                        constantBC,
                                        linearBC)
from proteus.mprans.BoundaryConditions import (BC_RANS)


comm = Comm.init()
Profiling.procID = comm.rank()
log("Testing BoundaryConditions")

# for creating a generic context used in some test cases
from proteus import Domain
epsFact_consrv_heaviside = 3.
domain = Domain.PlanarStraightLineGraphDomain()
domain.MeshOptions.he = 0.2

def create_BC(folder=None, b_or=None, b_i=None):
    if folder is None:
        return BC_Base(b_or=b_or, b_i=b_i)
    if folder == 'mprans':
        return BC_RANS(b_or=b_or, b_i=b_i)

def get_random_x(start=0., stop=10.):
    x1 = random.uniform(start, stop)
    x2 = random.uniform(start, stop)
    x3 = random.uniform(start, stop)
    return [x1, x2, x3]

def get_time_array(start=0, stop=5, steps=100):
    return np.linspace(0, 5, 100)

class TestBC(unittest.TestCase):

    def test_bc_base(self):
        BC = create_BC()
        BC.newGlobalBC('test1', 4)

        npt.assert_equal(BC.test1, 4)

    def test_constantBC(self):
        t_list = get_time_array()
        constant = 4
        constants = []
        values = []
        BC_func = constantBC(constant)
        for t in t_list:
            x = get_random_x()
            constants += [constant]
            values += [BC_func(x, t)]
        npt.assert_equal(values, constants)

    def test_non_material(self):
        BC = create_BC(folder='mprans')
        BC.setNonMaterial()
        vof_adv, u_diff, v_diff, w_diff = [], [], [], []
        t_list = get_time_array()
        for t in t_list:
            x = get_random_x()
            vof_adv += [BC.vof_advective(x, t)]
            u_diff += [BC.u_diffusive(x, t)]
            v_diff += [BC.v_diffusive(x, t)]
            w_diff += [BC.w_diffusive(x, t)]
        zeros = np.zeros(len(t_list))
        npt.assert_equal(vof_adv, zeros)
        npt.assert_equal(u_diff, zeros)
        npt.assert_equal(v_diff, zeros)
        npt.assert_equal(w_diff, zeros)

    def test_mprans_no_slip(self):
        BC = create_BC(folder='mprans')
        BC.setNoSlip()
        u_dir, v_dir, w_dir, p_adv, k_dir, d_diff, vof_adv = [], [], [], [], [], [], []
        t_list = get_time_array()
        for t in t_list:
            x = get_random_x()
            u_dir += [BC.u_dirichlet(x, t)]
            v_dir += [BC.v_dirichlet(x, t)]
            w_dir += [BC.w_dirichlet(x, t)]
            p_adv += [BC.p_advective(x, t)]
            k_dir += [BC.k_dirichlet(x, t)]
            d_diff += [BC.dissipation_diffusive(x, t)]
            vof_adv += [BC.vof_advective(x, t)]
        zeros = np.zeros(len(t_list))
        npt.assert_equal(BC.p_dirichlet, None)
        npt.assert_equal(u_dir, zeros)
        npt.assert_equal(v_dir, zeros)
        npt.assert_equal(w_dir, zeros)
        npt.assert_equal(BC.vof_dirichlet, None)
        npt.assert_equal(k_dir, zeros)
        npt.assert_equal(BC.dissipation_dirichlet, None)
        npt.assert_equal(p_adv, zeros)
        npt.assert_equal(BC.u_advective, None)
        npt.assert_equal(BC.v_advective, None)
        npt.assert_equal(BC.w_advective, None)
        npt.assert_equal(vof_adv, zeros)
        npt.assert_equal(BC.k_advective, None)
        npt.assert_equal(BC.dissipation_advective, None)
        npt.assert_equal(BC.u_diffusive, None)
        npt.assert_equal(BC.v_diffusive, None)
        npt.assert_equal(BC.w_diffusive, None)
        npt.assert_equal(BC.k_diffusive, None)
        npt.assert_equal(d_diff, zeros)

    def test_mprans_free_slip(self):
        BC = create_BC(folder='mprans')
        BC.setFreeSlip()
        u_adv, v_adv, w_adv,p_adv, u_diff, v_diff, w_diff, k_dir, d_diff, vof_adv = [], [], [], [], [], [], [], [], [], []
        t_list = get_time_array()
        for t in t_list:
            x = get_random_x()
            u_adv += [BC.u_advective(x, t)]
            v_adv += [BC.v_advective(x, t)]
            w_adv += [BC.w_advective(x, t)]
            p_adv += [BC.p_advective(x, t)]
            u_diff += [BC.u_diffusive(x, t)]
            v_diff += [BC.v_diffusive(x, t)]
            w_diff += [BC.w_diffusive(x, t)]
            k_dir += [BC.k_dirichlet(x, t)]
            d_diff += [BC.dissipation_diffusive(x, t)]
            vof_adv += [BC.vof_advective(x, t)]
        zeros = np.zeros(len(t_list))
        npt.assert_equal(BC.p_dirichlet, None)
        npt.assert_equal(BC.u_dirichlet, None)
        npt.assert_equal(BC.v_dirichlet, None)
        npt.assert_equal(BC.w_dirichlet, None)
        npt.assert_equal(BC.vof_dirichlet, None)
        npt.assert_equal(k_dir, zeros)
        npt.assert_equal(BC.dissipation_dirichlet, None)
        npt.assert_equal(p_adv, zeros)
        npt.assert_equal(u_adv, zeros)
        npt.assert_equal(v_adv, zeros)
        npt.assert_equal(w_adv, zeros)
        npt.assert_equal(vof_adv, zeros)
        npt.assert_equal(BC.k_advective, None)
        npt.assert_equal(BC.dissipation_advective, None)
        npt.assert_equal(u_diff, zeros)
        npt.assert_equal(v_diff, zeros)
        npt.assert_equal(w_diff, zeros)
        npt.assert_equal(d_diff, zeros)
        # check if other BC are None

    def test_open_air(self):
        BC = create_BC(folder='mprans')
        BC.setFreeSlip()
        u_adv, v_adv, w_adv,p_adv, u_diff, v_diff, w_diff, k_dir, d_diff, vof_adv = [], [], [], [], [], [], [], [], [], []
        t_list = get_time_array()
        for t in t_list:
            x = get_random_x()
            u_adv += [BC.u_advective(x, t)]
            v_adv += [BC.v_advective(x, t)]
            w_adv += [BC.w_advective(x, t)]
            p_adv += [BC.p_advective(x, t)]
            u_diff += [BC.u_diffusive(x, t)]
            v_diff += [BC.v_diffusive(x, t)]
            w_diff += [BC.w_diffusive(x, t)]
            k_dir += [BC.k_dirichlet(x, t)]
            d_diff += [BC.dissipation_diffusive(x, t)]
            vof_adv += [BC.vof_advective(x, t)]
        zeros = np.zeros(len(t_list))
        npt.assert_equal(BC.p_dirichlet, None)
        npt.assert_equal(BC.u_dirichlet, None)
        npt.assert_equal(BC.v_dirichlet, None)
        npt.assert_equal(BC.w_dirichlet, None)
        npt.assert_equal(BC.vof_dirichlet, None)
        npt.assert_equal(k_dir, zeros)
        npt.assert_equal(BC.dissipation_dirichlet, None)
        npt.assert_equal(p_adv, zeros)
        npt.assert_equal(u_adv, zeros)
        npt.assert_equal(v_adv, zeros)
        npt.assert_equal(w_adv, zeros)
        npt.assert_equal(vof_adv, zeros)
        npt.assert_equal(BC.k_advective, None)
        npt.assert_equal(BC.dissipation_advective, None)
        npt.assert_equal(u_diff, zeros)
        npt.assert_equal(v_diff, zeros)
        npt.assert_equal(w_diff, zeros)
        npt.assert_equal(d_diff, zeros)
        # check if other BC are None
    # def test_unsteady_two_phase_velocity_inlet(self):

    def test_set_tank(self):
        BC = create_BC(folder='mprans', b_or=[[0., 1., 0.]], b_i=0)
        BC.setTank()
        # checking if other BC leaves setTank BC as it should be
        BC.setFreeSlip()
        BC.reset()
        # ----
        t_list = get_time_array()
        hy_dir, u_stress, w_stress = [], [], []
        for t in t_list:
            x = get_random_x()
            hy_dir += [BC.hy_dirichlet(x, t)]
            u_stress += [BC.u_stress]
            w_stress += [BC.w_stress]
        zeros = np.zeros(len(t_list))
        npt.assert_equal(BC.hx_dirichlet, None)
        npt.assert_equal(hy_dir, zeros)
        npt.assert_equal(BC.hz_dirichlet, None)
        npt.assert_equal(hy_dir, zeros)
        npt.assert_equal(u_stress, zeros)
        npt.assert_equal(BC.v_stress, None)
        npt.assert_equal(w_stress, zeros)

    def test_move_mesh(self):
        BC = create_BC(folder='mprans')
        last_pos = [1., 1., 1.]
        h = [0., 0., 0.]
        rot_matrix = np.eye(3)
        BC.setMoveMesh(last_pos=last_pos, h=h, rot_matrix=rot_matrix)
        # checking if other BC leaves setTank BC as it should be
        BC.setFreeSlip()
        BC.reset()
        # ----
        t_list = get_time_array()
        hx_dir, hy_dir, hz_dir = [], [], []
        displacement = []
        from proteus.SpatialTools import rotation3D
        for t in t_list:
            new_pos = np.array(get_random_x())
            h[:] = new_pos-last_pos
            rot_angle = get_random_x()[0]
            rot_axis = get_random_x()
            rot_matrix[:] = rotation3D(rot_matrix, rot=rot_angle, axis=rot_axis)
            x = np.array(get_random_x())
            hx_dir += [BC.hx_dirichlet(x, t)]
            hy_dir += [BC.hy_dirichlet(x, t)]
            hz_dir += [BC.hz_dirichlet(x, t)]
            x0 = x-last_pos
            displacement += [(np.dot(x0, rot_matrix)-x0)+h]
            last_pos[:] = x
        zeros = np.zeros(len(t_list))
        displacement = np.array(displacement)
        npt.assert_equal(hx_dir, displacement[:, 0])
        npt.assert_equal(hy_dir, displacement[:, 1])
        npt.assert_equal(hz_dir, displacement[:, 2])

    def test_unsteady_two_phase_velocity_inlet(self):
        from proteus.WaveTools import MonochromaticWaves
        b_or = [[0., -1., 0.]]
        b_i = 0
        BC = create_BC(folder='mprans', b_or=b_or, b_i=b_i)
        # creating a wave
        period = 0.8
        height = 0.029
        mwl = depth = 0.9
        direction = np.array([1., 0., 0.])
        g = np.array([0., -9.81, 0.])
        waves = MonochromaticWaves(period, height, mwl, depth, g, direction)
        # need to set epsFact and he with context as they are called in BC...
        from proteus import Context
        case = sys.modules[__name__]
        Context.setFromModule(case)
        ct = Context.get()
        ecH = ct.epsFact_consrv_heaviside
        he = ct.domain.MeshOptions.he
        from proteus.ctransportCoefficients import smoothedHeaviside
        #-----
        # set BC
        windSpeed=np.array([1., 2., 3.4])
        BC.setUnsteadyTwoPhaseVelocityInlet(waves, vert_axis=1,
                                            windSpeed=windSpeed)
        u_dir, v_dir, w_dir, vof_dir, p_adv = [], [], [], [], []
        u_calc, vof_calc, p_calc = [], [], []
        t_list = get_time_array()
        for t in t_list:
            x = np.array(get_random_x())
            u_dir += [BC.u_dirichlet(x, t)]
            v_dir += [BC.v_dirichlet(x, t)]
            w_dir += [BC.w_dirichlet(x, t)]
            vof_dir += [BC.vof_dirichlet(x, t)]
            p_adv += [BC.p_advective(x, t)]
            # calculations
            waveHeight = waves.mwl+waves.eta(x, t)
            wavePhi = x[1]-waveHeight
            if wavePhi <= 0:
                wave_u = waves.u(x, t)
            else:
                x_max = list(x)
                x_max[1] = waveHeight
                wave_u = waves.u(x_max, t)
            Hu = smoothedHeaviside(0.5*ecH*he, wavePhi-0.5*ecH*he)
            U = Hu*windSpeed + (1-Hu)*wave_u
            u_calc += [U]
            p_calc += [np.sum(U*b_or[b_i])]
            Hvof = smoothedHeaviside(ecH*he, wavePhi)
            vof_calc += [Hvof]
        u_calc = np.array(u_calc)
        vof_calc = np.array(vof_calc)
        npt.assert_equal(BC.p_dirichlet, None)
        npt.assert_equal(u_dir, u_calc[:, 0])
        npt.assert_equal(v_dir, u_calc[:, 1])
        npt.assert_equal(w_dir, u_calc[:, 2])
        npt.assert_equal(vof_dir, vof_calc)
        npt.assert_equal(BC.k_dirichlet, None)
        npt.assert_equal(BC.dissipation_dirichlet, None)
        npt.assert_equal(p_adv, p_calc)
        npt.assert_equal(BC.u_advective, None)
        npt.assert_equal(BC.v_advective, None)
        npt.assert_equal(BC.w_advective, None)
        npt.assert_equal(BC.vof_advective, None)
        npt.assert_equal(BC.k_advective, None)
        npt.assert_equal(BC.dissipation_advective, None)
        npt.assert_equal(BC.u_diffusive, None)
        npt.assert_equal(BC.v_diffusive, None)
        npt.assert_equal(BC.w_diffusive, None)
        npt.assert_equal(BC.k_diffusive, None)
        npt.assert_equal(BC.dissipation_diffusive, None)





if __name__ == '__main__':

    unittest.main(verbosity=2)

    from proteus import Domain
    epsFact_consrv_heaviside = ecH = 3.
    domain = Domain.PlanarStraightLineGraphDomain()
    domain.MeshOptions.he = he = 0.2
