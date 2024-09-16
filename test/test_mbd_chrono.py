import numpy as np
import numpy.testing as npt
import unittest
from proteus.mbd import CouplingFSI as fsi
import pytest
import pychrono as chrono

class TestCable(unittest.TestCase):
    @pytest.mark.skip()
    def testHangingCableANCF(self):
        g = np.array([0.,0.,-9.81])
        system = fsi.ProtChSystem()
        system.ChSystem.SetGravitationalAcceleration(chrono.ChVector3d(g[0], g[1], g[2]))
        system.setTimeStep(1e-1)
        timestepper = chrono.ChTimestepperEulerImplicitLinearized()
        system.ChSystem.SetTimestepper(timestepper)
        solver = chrono.ChSolverMINRES()
        system.ChSystem.SetSolver(solver)
        solver.SetMaxIterations(100)
        solver.EnableWarmStart(True)
        solver.EnableDiagonalPreconditioner(True)
        solver.SetVerbose(True)
        system.ChSystem.GetSolver().AsIterative().SetTolerance(1e-10)
        #system.ChSystem.GetSolver().AsIterative().SetTolerancePrimal(1e-10)
        system.ChSystem.GetSolver().AsIterative().SetMaxIterations(100)
        mesh = fsi.ProtChMesh(system)
        L = np.array([5.])
        nb_elems = np.array([3])
        d = np.array([1e-3])
        rho = np.array([1000.])
        E = np.array([1e10])
        cable_type = b"CableANCF"
        fairlead_body = fsi.ProtChBody(system)
        fairlead_body.ChBody.SetFixed(True)
        mooring = fsi.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=cable_type)
        # vertical cable
        mooring.setNodesPositionFunction(lambda s: np.array([0., 0., s]), lambda s: np.array([0., 0., 1.]))
        mooring.setNodesPosition()
        mooring.buildNodes()
        mooring.attachBackNodeToBody(fairlead_body)
        system.calculate_init()
        system.calculate(0.5)
        T = mooring.getTensionBack()
        strain = mooring.getNodesTension(eta=1.)[-1]*np.pi*d**2/4*E
        T_sol = -np.ones(3)*g*rho*(np.pi*d**2/4.*L)
        npt.assert_almost_equal(-T, T_sol)

    def testSetterGetter(self):
        g = np.array([0.,0.,-9.81])
        system = fsi.ProtChSystem()
        system.setGravitationalAcceleration(g)
        g0 = system.getGravitationalAcceleration()
        npt.assert_almost_equal(g0, g)

        mesh = fsi.ProtChMesh(system)
        mesh.getChronoObject()

        body = fsi.ProtChBody(system)
        # position
        pos = np.array([1.,2.,3.])
        body.setPosition(pos)
        pos0 = body.getPosition()
        npt.assert_almost_equal(pos0, pos)
        # inertia
        inertia = np.array([[1., 4., 5.],
                            [4., 2., 6.],
                            [5., 6., 3.]])
        inertiaXX = np.array([inertia[0, 0], inertia[1, 1], inertia[2, 2]])
        inertiaXY = np.array([inertia[0, 1], inertia[0, 2], inertia[1, 2]])
        body.setInertiaXX(inertiaXX)
        # inertiaXX0 = body.getInertiaXX()
        body.setInertiaXY(inertiaXY)
        # inertiaXY0 = body.getInertiaXY()
        inertia0 = body.getInertia()
        # npt.assert_almost_equal(inertiaXX0, inertiaXX)
        # npt.assert_almost_equal(inertiaXY0, inertiaXY)
        npt.assert_almost_equal(inertia0, inertia)
        # mass
        mass = 10.2
        body.setMass(mass)
        mass0 = body.getMass()
        npt.assert_almost_equal(mass0, mass)

    #CURRENTLY NOT FULLY WORKING
    @pytest.mark.skip()
    def testHangingCableEuler(self):
        g = np.array([0.,0.,-9.81])
        system = fsi.ProtChSystem()
        system.setGravitationalAcceleration(g)
        system.setTimeStep(1e-1)
        mesh = fsi.ProtChMesh(system)
        mesh.getChronoObject()
        L = np.array([5.])
        nb_elems = np.array([3])
        d = np.array([1e-3])
        rho = np.array([1000.])
        E = np.array([1e10])
        cable_type = b"BeamEuler"
        mooring = fsi.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=cable_type)
        mooring.external_forces_manual = True # tri: should work without this line
        # vertical cable
        fairlead_body = fsi.ProtChBody(system)
        fairlead_body.ChBody.SetFixed(True)
        mooring.setNodesPositionFunction(lambda s: np.array([0., 0., s]), lambda s: np.array([0., 0., 1.]))
        mooring.setNodesPosition()
        mooring.buildNodes()
        mooring.attachBackNodeToBody(fairlead_body)
        system.calculate_init()
        system.calculate(0.5)
        T = mooring.getTensionBack()
        T_sol = -np.ones(3)*g*rho*(np.pi*d**2/4.*L)
        npt.assert_almost_equal(T, T_sol)
