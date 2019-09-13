import numpy as np
import numpy.testing as npt
import unittest
from proteus.mbd import CouplingFSI as fsi
import pytest
import pychrono as chrono

def getSystem():
    system = fsi.ProtChSystem()
    system.ChSystem.Set_G_acc(chrono.ChVectorD(0., 0., -9.81))
    system.setTimeStep(1e-1)
    system.ChSystem.SetSolverType(chrono.ChSolver.Type_MINRES)
    system.ChSystem.SetSolverWarmStarting(True)
    system.ChSystem.SetMaxItersSolverSpeed(100)
    system.ChSystem.SetMaxItersSolverStab(100)
    system.ChSystem.SetTolForce(1e-10)
    return system


class TestCable(unittest.TestCase):
    def testHangingCableANCF(self):
        system = getSystem()
        mesh = fsi.ProtChMesh(system)
        L = np.array([5.])
        nb_elems = np.array([3])
        d = np.array([1e-3])
        rho = np.array([1000.])
        E = np.array([1e10])
        beam_type = b"CableANCF"
        mooring = fsi.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=beam_type)
        mooring.external_forces_manual = True # tri: should work without this line
        # vertical cable
        fairlead_body = fsi.ProtChBody(system)
        fairlead_body.ChBody.SetBodyFixed(True)
        mooring.setNodesPositionFunction(lambda s: np.array([0., 0., s]), lambda s: np.array([0., 0., 1.]))
        mooring.setNodesPosition()
        mooring.buildNodes()
        mooring.attachBackNodeToBody(fairlead_body)
        system.calculate_init()
        system.calculate(0.5)
        T = mooring.getTensionBack()
        g = np.array([0., 0., -9.81])
        strain = mooring.getNodesTension(eta=1.)[-1]*np.pi*d**2/4*E
        T_sol = -np.ones(3)*g*rho*(np.pi*d**2/4.*L)
        npt.assert_almost_equal(-T, T_sol)

    def testHydrodynamicForces(self):
        g = np.array([0.,0.,-9.81])
        system = getSystem()
        mesh = fsi.ProtChMesh(system)
        L = np.array([5.])
        nb_elems = np.array([2])
        d = np.array([1e-3])
        rho = np.array([1000.])
        E = np.array([1e10])
        beam_type = b"CableANCF"
        mooring = fsi.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=beam_type)
        Cdt = 1.
        Cdn = 1.
        Cat = 1.
        Can = 1.
        mooring.external_forces_manual = True # tri: should work without this line
        mooring.setDragCoefficients(tangential=Cdt, normal=Cdn, segment_nb=0)
        mooring.setAddedMassCoefficients(tangential=Cat, normal=Can, segment_nb=0)
        # vertical cable
        fairlead_body = fsi.ProtChBody(system)
        fairlead_body.ChBody.SetBodyFixed(True)
        mooring.setNodesPositionFunction(lambda s: np.array([0., 0., s]), lambda s: np.array([0., 0., 1.]))
        mooring.setNodesPosition()
        mooring.buildNodes()
        mooring.attachBackNodeToBody(fairlead_body)
        fluid_density = np.array([1000., 1000., 1000.])
        fluid_velocity = np.array([[0.7, 0.5, 1.], [-0.8, 0.7, 0.3], [2.2, 0.4, 1.3]])
        mooring.setFluidDensityAtNodes(fluid_density)
        mooring.setFluidVelocityAtNodes(fluid_velocity)
        # mooring.setFluidAccelerationAtNodes(fluid_acceleration)
        system.calculate_init()
        mooring.fluid_velocity_array_previous = np.zeros((nb_elems[0]+1, 3))
        system.calculate(1.)

        fluid_velocity_chrono = mooring.getFluidVelocity()
        fluid_acceleration_chrono = mooring.getFluidAcceleration()
        fluid_acceleration = mooring.getFluidAcceleration()
        forces_drag_chrono = mooring.getDragForces()
        forces_addedmass_chrono = mooring.getAddedMassForces()

        # calculate forces
        forces_drag = np.zeros((nb_elems[0]+1, 3))
        forces_addedmass = np.zeros((nb_elems[0]+1, 3))
        t_dir = np.array([0., 0., 1.])
        for ii in range(nb_elems[0]+1):
            urel = fluid_velocity[ii]-0  # 0 because no velocity of cable
            arel = fluid_acceleration[ii]-0  # 0 because no velocity of cable
            Vt = np.dot(urel, t_dir)*t_dir  # tangential
            Vn = urel-Vt  # normal
            Fdt = 0.5*fluid_density[ii]*Cdt*np.pi*d[0]*np.linalg.norm(Vt)*Vt
            Fdn = 0.5*fluid_density[ii]*Cdn*d[0]*np.linalg.norm(Vn)*Vn
            forces_drag[ii] = Fdt+Fdn
            Vt = np.dot(arel, t_dir)*t_dir  # tangential
            Vn = arel-Vt  # normal
            Fat = fluid_density[ii]*Cat*np.pi*d[0]**2/4.*Vt
            Fan = fluid_density[ii]*Can*np.pi*d[0]**2/4.*Vn
            Faf = fluid_density[ii]*np.pi*d[0]**2/4.*fluid_acceleration[ii]
            forces_addedmass[ii] = Fat+Fan+Faf

        for ii in range(nb_elems[0]):
            npt.assert_almost_equal(fluid_velocity_chrono[ii], fluid_velocity[ii])
            npt.assert_almost_equal(fluid_acceleration_chrono[ii], fluid_acceleration[ii])
            npt.assert_almost_equal(forces_drag_chrono[ii], forces_drag[ii])
            npt.assert_almost_equal(forces_addedmass_chrono[ii], forces_addedmass[ii])

    # # CURRENTLY NOT FULLY WORKING
    # def testHangingCableEuler(self):
    #     g = np.array([0.,0.,-9.81])
    #     system = getSystem()
    #     mesh = fsi.ProtChMesh(system)
    #     L = np.array([5.])
    #     nb_elems = np.array([3])
    #     d = np.array([1e-3])
    #     rho = np.array([1000.])
    #     E = np.array([1e10])
    #     cable_type = b"BeamEuler"
    #     mooring = fsi.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=cable_type)
    #     mooring.external_forces_manual = True # tri: should work without this line
    #     # vertical cable
    #     fairlead_body = fsi.ProtChBody(system)
    #     fairlead_body.ChBody.SetBodyFixed(True)
    #     mooring.setNodesPositionFunction(lambda s: np.array([0., 0., s]), lambda s: np.array([0., 0., 1.]))
    #     mooring.setNodesPosition()
    #     mooring.buildNodes()
    #     mooring.attachBackNodeToBody(fairlead_body)
    #     system.calculate_init()
    #     system.calculate(0.5)
    #     T = mooring.getTensionBack()
    #     T_sol = -np.ones(3)*g*rho*(np.pi*d**2/4.*L)
    #     npt.assert_almost_equal(T, T_sol)


