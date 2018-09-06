import numpy as np
import numpy.testing as npt
import unittest
from proteus.mbd import ChRigidBody as crb
import pytest

class TestCable(unittest.TestCase):
    # @pytest.mark.skip(reason="unpredictable segfaults")
    def testHangingCableANCF(self):
        g = np.array([0.,0.,-9.81])
        system = crb.ProtChSystem(gravity=g)
        system.setTimeStep(1e-1)
        mesh = crb.Mesh(system)
        L = np.array([5.])
        nb_elems = np.array([3])
        d = np.array([1e-3])
        rho = np.array([1000.])
        E = np.array([1e10])
        cable_type = "CableANCF"
        mooring = crb.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=cable_type)
        mooring.external_forces_manual = True # tri: should work without this line
        # vertical cable
        fairlead_body = crb.ProtChBody(system)
        fairlead_body.ChBody.SetBodyFixed(True)
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

    # CURRENTLY NOT FULLY WORKING
    # def testHangingCableEuler(self):
    #     g = np.array([0.,0.,-9.81])
    #     system = crb.ProtChSystem(gravity=g)
    #     system.setTimeStep(1e-1)
    #     mesh = crb.Mesh(system)
    #     L = np.array([5.])
    #     nb_elems = np.array([3])
    #     d = np.array([1e-3])
    #     rho = np.array([1000.])
    #     E = np.array([1e10])
    #     cable_type = "BeamEuler"
    #     mooring = crb.ProtChMoorings(system=system, mesh=mesh, length=L, nb_elems=nb_elems, d=d, rho=rho, E=E, beam_type=cable_type)
    #     mooring.external_forces_manual = True # tri: should work without this line
    #     # vertical cable
    #     fairlead_body = crb.ProtChBody(system)
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


