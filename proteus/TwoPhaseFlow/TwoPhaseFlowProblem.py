from __future__ import division
from past.utils import old_div
from proteus import FemTools as ft
from proteus import MeshTools as mt
from proteus.Profiling import logEvent
from builtins import object
from proteus.TwoPhaseFlow.utils import Parameters

class TwoPhaseFlowProblem:
    """ TwoPhaseFlowProblem """

    def __init__(self,
                 ns_model=0, #0: rans2p, 1: rans3p
                 ls_model=1, #0: vof+ncls+rdls+mcorr, 1: clsvof
                 nd=2,
                 # TIME STEPPING #
                 cfl=0.33,
                 outputStepping=None,
                 # DOMAIN AND MESH #
                 structured=False,
                 he=None,
                 nnx=None,
                 nny=None,
                 nnz=None,
                 domain=None,
                 triangleFlag=1,
                 # INITIAL CONDITIONS #
                 initialConditions=None,
                 # BOUNDARY CONDITIONS #
                 boundaryConditions=None,
                 # OTHERS #
                 useSuperlu=False):
        # ***** SET OF ASSERTS ***** #
        if ns_model is not None:
            assert ns_model in [0,1], "ns_model={0,1} for rans2p or rans3p respectively"
        if ls_model is not None:
            assert ls_model in [0,1], "ls_model={0,1} for vof+ncls+rdls+mcorr or clsvof respectively"
        assert nd in [2,3], "nd={2,3}"
        assert cfl <= 1, "Choose cfl <= 1"
        assert isinstance (outputStepping,OutputStepping), "Provide an object from the OutputStepping class"
        assert type(he)==float , "Provide (float) he (characteristic mesh size)"
        assert domain is not None, "Provide a domain"
        if structured:
            assert type(nnx)==int and type(nny)==int, "Provide (int) nnx and nny"
            if nd==3:
                assert type(nnz)==int, "Provide (int) nnz"
        else:
            assert domain.MeshOptions.triangleOptions != 'q30DenA', "Set domain.MeshOptions.triangleOptions"
        assert triangleFlag in [0,1,2], "triangleFlag must be 1, 2 or 3"
        if initialConditions is not None:
            assert type(initialConditions)==dict, "Provide dict of initial conditions"
            # assertion now done in TwoPhaseFlow_so.py
        if boundaryConditions is not None:
            assert type(boundaryConditions)==dict, "Provide dict of boundary conditions"
            # assertion now done in TwoPhaseFlow_so.py

        # ***** SAVE PARAMETERS ***** #
        self.domain=domain
        self.Parameters = Parameters.ParametersHolder(ProblemInstance=self)
        self.ns_model=ns_model
        self.ls_model = ls_model
        self.nd=nd
        self.cfl=cfl
        self.outputStepping=outputStepping
        self.outputStepping.setOutputStepping()
        self.Parameters.mesh.he = he
        self.Parameters.mesh.nnx = nnx
        self.Parameters.mesh.nny = nny
        self.Parameters.mesh.nnz = nnz
        self.Parameters.mesh.triangleFlag = triangleFlag
        self.Parameters.mesh.setTriangleOptions()
        self.initialConditions=initialConditions
        self.boundaryConditions=boundaryConditions
        self.useSuperlu = useSuperlu
        self.movingDomain = False
        self.archiveAllSteps = False
        # to use proteus.mprans.BoundaryConditions
        # but only if SpatialTools was used to make the domain
        self.useBoundaryConditionsModule = True


        # ***** CREATE FINITE ELEMENT SPACES ***** #
        self.FESpace = FESpace(self.ns_model, self.nd)
        self.FESpace.setFESpace()

        # ***** DEFINE PHYSICAL AND NUMERICAL PARAMETERS ***** #
        self.physical_parameters = default_physical_parameters
        self.rans2p_parameters = default_rans2p_parameters
        self.rans3p_parameters = default_rans3p_parameters
        self.clsvof_parameters = default_clsvof_parameters

        # set indice of models if ns_model and ls_model is set
        ind = 0
        if self.ls_model == 0:
            if self.ns_model == 0:
                self.Parameters.Models.rans2p.index = ind
                ind += 1
            else:
                self.Parameters.Models.rans3p.index = ind
                ind += 1
                self.Parameters.Models.pressureIncrement.index = ind
                ind += 1
                self.Parameters.Models.pressure.index = ind
                ind += 1
                self.Parameters.Models.pressureInitial.index = ind
                ind += 1
            self.Parameters.Models.vof.index = ind
            ind += 1
            self.Parameters.Models.ncls.index = ind
            ind += 1
            self.Parameters.Models.rdls.index = ind
            ind += 1
            self.Parameters.Models.mcorr.index = ind
            ind += 1
        elif self.ls_model == 1:
            if self.ns_model == 0:
                self.Parameters.Models.rans2p.index = ind
                ind += 1
                self.Parameters.Models.clsvof.index = ind
                ind += 1
            else:
                self.Parameters.Models.clsvof.index = ind
                ind += 1
                self.Parameters.Models.rans3p.index = ind
                ind += 1
                self.Parameters.Models.pressureIncrement.index = ind
                ind += 1
                self.Parameters.Models.pressure.index = ind
                ind += 1
                self.Parameters.Models.pressureInitial.index = ind
                ind += 1
        # ***** DEFINE OTHER GENERAL NEEDED STUFF ***** #
        self.general = default_general

    def assert_initialConditions(self):
        initialConditions = self.initialConditions
        nd = self.nd
        ns_model = self.ns_model
        ls_model = self.ls_model
        if ns_model is not None:
            assert 'pressure' in initialConditions, 'Provide pressure in ICs'
            assert 'vel_u' in initialConditions, 'Provide vel_u in ICs'
            assert 'vel_v' in initialConditions, 'Provide vel_v in ICs'
            if nd==3:
                assert 'vel_w' in initialConditions, 'Provide vel_w in ICs'
            if self.ns_model == 1: #rans3p
                assert 'pressure_increment' in initialConditions, 'Provide pressure_increment in ICs'
        if ls_model == 0:
            assert 'vof' in initialConditions, 'Provide vof in ICs'
            assert 'ncls' in initialConditions, 'Provide ncls in ICs'
        elif self.ls_model == 1:
            assert 'clsvof' in initialConditions or ('ncls' in initialConditions and 'vof' in initialConditions), 'Provide clsvof or ncls and vof in ICs'
    #
    def assert_boundaryConditions(self):
        boundaryConditions = self.boundaryConditions
        nd = self.nd
        ns_model = self.ns_model
        ls_model = self.ls_model
        if boundaryConditions is not None:
            # check dirichlet BCs
            if ns_model is not None:
                assert 'pressure_DBC' in boundaryConditions, "Provide pressure_DBC"
                assert 'vel_u_DBC' in boundaryConditions, "Provide vel_u_DBC"
                assert 'vel_v_DBC' in boundaryConditions, "Provide vel_v_DBC"
                if nd==3:
                    assert 'vel_w_DBC' in boundaryConditions, "Provide vel_w_DBC"
            if ls_model == 0:
                assert 'vof_DBC' in boundaryConditions, "Provide vof_DBC"
                assert 'ncls_DBC' in boundaryConditions, "Provide ncls_DBC"
            elif ls_model == 1:
                assert 'clsvof_DBC' in boundaryConditions, "Provide clsvof_DBC"
            # check advective flux BCs
            if ns_model is not None:
                assert 'pressure_AFBC' in boundaryConditions, "Provide pressure_AFBC"
                assert 'vel_u_AFBC' in boundaryConditions, "Provide vel_u_AFBC"
                assert 'vel_v_AFBC' in boundaryConditions, "Provide vel_v_AFBC"
                if nd==3:
                    assert 'vel_w_AFBC' in boundaryConditions, "Provide vel_w_AFBC"
            if ls_model == 0:
                assert 'clsvof_AFBC' in boundaryConditions, "Provide clsvof_AFBC"
            if ls_model == 0:
                assert 'vof_AFBC' in boundaryConditions, "Provide vof_AFBC"
            # check diffusive flux BCs
            if ns_model is not None:
                assert 'vel_u_DFBC' in boundaryConditions, "provide vel_u_DFBC"
                assert 'vel_v_DFBC' in boundaryConditions, "provide vel_v_DFBC"
                if nd==3:
                    assert 'vel_w_DFBC' in boundaryConditions, "provide vel_w_DFBC"
            if ls_model == 1:
                assert 'clsvof_DFBC' in boundaryConditions, "provide clsvof_DFBC"
            if ns_model==1: #rans3p
                # check dirichlet BCs
                assert 'pressure_increment_DBC' in boundaryConditions, "Provide pressure_increment_DBC"
                # check advective flux BCs
                assert 'pressure_increment_AFBC' in boundaryConditions,"Provide pressure_increment_AFBC"
                # check diffusive flux BCs
                assert 'pressure_increment_DFBC' in boundaryConditions,"Provide pressure_increment_DFBC"
        else:
            assert self.domain.useSpatialTools is True, 'Either define boundaryConditions dict or use proteus.mprans.SpatialTools to set Boundary Conditions and run function assembleDomain'

    def initializeAll(self):
        # initial conditions
        self.assert_initialConditions()
        # boundary conditions
        self.assert_boundaryConditions()
        # parameters
        self.Parameters.initializeParameters()
        # mesh
        if self.Parameters.mesh.outputFiles['poly'] is True:
            self.domain.writePoly(self.Parameters.mesh.outputFiles_name)
        if self.Parameters.mesh.outputFiles['ply'] is True:
            self.domain.writePLY(self.Parameters.mesh.outputFiles_name)
        if self.Parameters.mesh.outputFiles['asymptote'] is True:
            self.domain.writeAsymptote(self.Parameters.mesh.outputFiles_name)
        if self.Parameters.mesh.outputFiles['geo'] is True or self.Parameters.mesh.use_gmsh is True:
            self.domain.writeGeo(self.Parameters.mesh.outputFiles_name)
            self.domain.use_gmsh = True

class OutputStepping:
    """
    OutputStepping handles how often the solution is outputted.
    """
    def __init__(self,
                 final_time,
                 dt_init=0.001,
                 dt_output=None,
                 nDTout=None,
                 dt_fixed=None):
        self.final_time=final_time
        self.dt_init=dt_init
        assert not (dt_output is None and nDTout is None), "Provide dt_output or nDTout"
        self.dt_output=dt_output
        self.nDTout = nDTout
        self.dt_fixed = dt_fixed
        self.systemStepExact = False

    def __getitem__(self, key):
        return self.__dict__[key]

    def setOutputStepping(self):
        # COMPUTE dt_init #
        self.dt_init = min(0.1 * self.dt_output, self.dt_init)
        if self.nDTout is None:
            self.nDTout = int(round(old_div(self.final_time, self.dt_output)))
        else:
            self.dt_output = float(self.final_time)/float(self.nDTout)
#

class FESpace:
    """
    Create FE Spaces.
    """
    def __init__(self,ns_model,nd):
        if ns_model is not None:
            assert ns_model == 0 or ns_model == 1, "ns_model must be 0 (rans2p) or 1 (rans3p)"
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        self.ns_model=ns_model
        self.nd=nd
        # For now we just support rans2p with: p1-p1 and rans3p with: p2-p1
        if ns_model == 0 or ns_model is None: # rans2p or None
            self.velSpaceOrder=1
            self.pSpaceOrder=1
        else: #rans3p
            self.velSpaceOrder=2
            self.pSpaceOrder=1

    def __getitem__(self, key):
        return self.__dict__[key]

    def setFESpace(self):
        ##################
        # VELOCITY SPACE #
        ##################
        if self.velSpaceOrder == 1: # p1 space
            self.hFactor = 1.0
            self.velBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            self.hFactor = 0.5
            self.velBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ##################
        # PRESSURE SPACE #
        ##################
        if self.pSpaceOrder == 1: # p1 space
            self.pBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis
        else: # p2 space
            self.pBasis = ft.C0_AffineQuadraticOnSimplexWithNodalBasis
        ###################
        # LEVEL SET SPACE #
        ###################
        self.lsBasis = ft.C0_AffineLinearOnSimplexWithNodalBasis # p1 space
        ###################
        # QUADRATURE RULE #
        ###################
        if max(self.velSpaceOrder,self.pSpaceOrder)==1:
            self.elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 3)
            self.elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 3)
        else:
            self.elementQuadrature = ft.SimplexGaussQuadrature(self.nd, 5)
            self.elementBoundaryQuadrature = ft.SimplexGaussQuadrature(self.nd - 1, 5)

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
default_physical_parameters ={'densityA': 998.2,
                              'viscosityA': 1.004e-6,
                              'densityB': 1.205,
                              'viscosityB': 1.500e-5,
                              'surf_tension_coeff': 72.8E-3,
                              'gravity': [0.0, -9.8, 0.0]}

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
default_rans2p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 3.0,
                             'epsFact_density': 3.0,
                             'ns_forceStrongDirichlet': False,
                             'weak_bc_penalty_constant': 1.0E6,
                             'useRBLES': 0.0,
                             'useRANS': 0.0,
                             'ns_closure': 0,
                             'useVF': 1.0,
                             'ns_shockCapturingFactor': 0.25,
                             'ns_lag_shockCapturing': True,
                             'ns_lag_subgridError': True,
                             'timeDiscretization': 'vbdf',
                             'timeOrder': 2}
default_rans3p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 3.0,
                             'epsFact_density': 3.0,
                             'ns_forceStrongDirichlet': False,
                             'ns_sed_forceStrongDirichlet': False,
                             'weak_bc_penalty_constant': 1.0E6,
                             'useRBLES': 0.0,
                             'useRANS': 0.0,
                             'ns_closure': 0,
                             'useVF': 1.0,
                             'ns_shockCapturingFactor': 0.5,
                             'ns_lag_shockCapturing': True,
                             'ns_lag_subgridError': True,
                             'timeDiscretization': 'vbdf',
                             'timeOrder': 2,
                             'PSTAB': 0,
                             'USE_SUPG': True,
                             'ARTIFICIAL_VISCOSITY': 2,
                             'cE': 1.0,
                             'cMax': 1.0}
default_clsvof_parameters = {'useMetrics': 1.0,
                             'epsFactHeaviside': 1.5,
                             'epsFactDirac': 1.5,
                             'epsFactRedist': 0.33,
                             'lambdaFact': 10.0,
                             'outputQuantDOFs': True,
                             'computeMetrics': 1,
                             'eps_tolerance_clsvof': False}
default_general = {'nLevels': 1,
                   'nLayersOfOverlapForParallel': 0}
