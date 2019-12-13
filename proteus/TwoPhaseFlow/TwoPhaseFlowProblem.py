from __future__ import division
from past.utils import old_div
from proteus import FemTools as ft
from proteus import MeshTools as mt
from proteus import SplitOperator
from proteus.Archiver import ArchiveFlags
from proteus.Profiling import logEvent
from builtins import object
from proteus.TwoPhaseFlow.utils import Parameters
from proteus.defaults import System_base

class TwoPhaseFlowProblem:
    """ TwoPhaseFlowProblem """

    def __init__(self,
                 models=None, #list of models
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
                 triangleFlag=0,
                 # INITIAL CONDITIONS #
                 initialConditions=None,
                 # BOUNDARY CONDITIONS #
                 boundaryConditions=None,
                 # AUXILIARY VARIABLES #
                 auxVariables=None,
                 # OTHERS #
                 useSuperlu=True,
                 fastArchive=False):
        # ***** SET OF ASSERTS ***** #
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

        # ***** SAVE PARAMETERS ***** #
        self.domain=domain
        self.Parameters = Parameters.ParametersHolder(ProblemInstance=self)
        self.modelList=[] #list used for internal tracking of models
        self.modelDict={} #dict used to allow user access to models
        self.modelIdxDict = {} #dict used for internal tracking of model indices
        self.nd=nd
        self.cfl=cfl
        self.outputStepping=outputStepping
        self.outputStepping.setOutputStepping()
        self.so = System_base()
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

        # ***** CREATE FINITE ELEMENT  ***** #
        self.FESpace = FESpace(self.nd)

        # ***** DEFINE OTHER GENERAL NEEDED STUFF ***** #
        self.general = default_general
        self.fastArchive = fastArchive
        self.usePETScOptionsFileExternal = False

    def addModel(self,modelObject,name):
        self.modelList.append(modelObject)
        self.modelDict[name] = modelObject

    def attachModels(self):
        #attach problem object to each model, identify index, and then form dictionary
        for (idx,model) in enumerate(self.modelList):
            model._Problem = self
            model.index = idx
            model._Problem.modelIdxDict[model.name]=idx

    def assert_initialConditions(self):
        initialConditions = self.initialConditions
        nd = self.nd
        #ns_model = self.ns_model
        #ls_model = self.ls_model
        #if ns_model is not None:
        #    assert 'pressure' in initialConditions, 'Provide pressure in ICs'
        #    assert 'vel_u' in initialConditions, 'Provide vel_u in ICs'
        #    assert 'vel_v' in initialConditions, 'Provide vel_v in ICs'
        #    if nd==3:
        #        assert 'vel_w' in initialConditions, 'Provide vel_w in ICs'
        #    if self.ns_model == 1: #rans3p
        #        assert 'pressure_increment' in initialConditions, 'Provide pressure_increment in ICs'
        #if ls_model == 0:
        #    assert 'vof' in initialConditions, 'Provide vof in ICs'
        #    assert 'ncls' in initialConditions, 'Provide ncls in ICs'
        #elif self.ls_model == 1:
        #    assert 'clsvof' in initialConditions or ('ncls' in initialConditions and 'vof' in initialConditions), 'Provide clsvof or ncls and vof in ICs'
    #
    def assert_boundaryConditions(self):
        boundaryConditions = self.boundaryConditions
        nd = self.nd
        #ns_model = self.ns_model
        #ls_model = self.ls_model
        #if boundaryConditions is not None:
        #    # check dirichlet BCs
        #    if ns_model is not None:
        #        assert 'pressure_DBC' in boundaryConditions, "Provide pressure_DBC"
        #        assert 'vel_u_DBC' in boundaryConditions, "Provide vel_u_DBC"
        #        assert 'vel_v_DBC' in boundaryConditions, "Provide vel_v_DBC"
        #        if nd==3:
        #            assert 'vel_w_DBC' in boundaryConditions, "Provide vel_w_DBC"
        #    if ls_model == 0:
        #        assert 'vof_DBC' in boundaryConditions, "Provide vof_DBC"
        #        assert 'ncls_DBC' in boundaryConditions, "Provide ncls_DBC"
        #    elif ls_model == 1:
        #        assert 'clsvof_DBC' in boundaryConditions, "Provide clsvof_DBC"
        #    # check advective flux BCs
        #    if ns_model is not None:
        #        assert 'pressure_AFBC' in boundaryConditions, "Provide pressure_AFBC"
        #        assert 'vel_u_AFBC' in boundaryConditions, "Provide vel_u_AFBC"
        #        assert 'vel_v_AFBC' in boundaryConditions, "Provide vel_v_AFBC"
        #        if nd==3:
        #            assert 'vel_w_AFBC' in boundaryConditions, "Provide vel_w_AFBC"
        #    if ls_model == 1:
        #        assert 'clsvof_AFBC' in boundaryConditions, "Provide clsvof_AFBC"
        #    if ls_model == 0:
        #        assert 'vof_AFBC' in boundaryConditions, "Provide vof_AFBC"
        #    # check diffusive flux BCs
        #    if ns_model is not None:
        #        assert 'vel_u_DFBC' in boundaryConditions, "provide vel_u_DFBC"
        #        assert 'vel_v_DFBC' in boundaryConditions, "provide vel_v_DFBC"
        #        if nd==3:
        #            assert 'vel_w_DFBC' in boundaryConditions, "provide vel_w_DFBC"
        #    if ls_model == 1:
        #        assert 'clsvof_DFBC' in boundaryConditions, "provide clsvof_DFBC"
        #    if ns_model==1: #rans3p
        #        # check dirichlet BCs
        #        assert 'pressure_increment_DBC' in boundaryConditions, "Provide pressure_increment_DBC"
        #        # check advective flux BCs
        #        assert 'pressure_increment_AFBC' in boundaryConditions,"Provide pressure_increment_AFBC"
        #        # check diffusive flux BCs
        #        assert 'pressure_increment_DFBC' in boundaryConditions,"Provide pressure_increment_DFBC"
        #else:
        #    assert self.domain.useSpatialTools is True, 'Either define boundaryConditions dict or use proteus.mprans.SpatialTools to set Boundary Conditions and run function assembleDomain'

    def initializeAll(self):

        #model organization
        self.attachModels()
        self.Parameters.model_list=self.modelList

        # initial conditions
        self.assert_initialConditions()
        # boundary conditions
        self.assert_boundaryConditions()

        # ***** SET FINITE ELEMENT  ***** #
        self.FESpace.setFESpace(self.modelIdxDict)

        # parameters
        self.Parameters.initializeParameters()
        # mesh
        # if self.Parameters.mesh.outputFiles['poly'] is True:
        #     self.domain.writePoly(self.Parameters.mesh.outputFiles_name)
        # if self.Parameters.mesh.outputFiles['ply'] is True:
        #     self.domain.writePLY(self.Parameters.mesh.outputFiles_name)
        # if self.Parameters.mesh.outputFiles['asymptote'] is True:
        #     self.domain.writeAsymptote(self.Parameters.mesh.outputFiles_name)
        # if self.Parameters.mesh.outputFiles['geo'] is True or self.Parameters.mesh.use_gmsh is True:
        #     self.domain.writeGeo(self.Parameters.mesh.outputFiles_name)
        #     self.domain.use_gmsh = True
        # split operator

        self.initializeSO()

    def initializeSO(self):
        so = self.so
        params = self.Parameters
        # list of models
        so.pnList = [None for i in range(len(params.model_list))]
        for model in params.model_list:
            so.pnList[model.index] = (model.p, model.n)
        so.systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
        if self.outputStepping.dt_fixed:
            so.dt_system_fixed = self.outputStepping.dt_fixed
        # rans3p specific options
        if 'rans3p' in self.modelIdxDict:
            PINIT_model = self.modelIdxDict['pressureInitial'] 
            assert PINIT_model is not None, 'must set pressureInitial model index when using rans3p'
            so.modelSpinUpList = [PINIT_model]
            from proteus.default_so import defaultSystem
            class Sequential_MinAdaptiveModelStepPS(SplitOperator.Sequential_MinAdaptiveModelStep):
                def __init__(self,modelList,system=defaultSystem,stepExact=True):
                    SplitOperator.Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
                    self.modelList = modelList[:len(so.pnList)-1]
            #
            so.systemStepControllerType = Sequential_MinAdaptiveModelStepPS
        # others
        so.needEBQ_GLOBAL = False
        so.needEBQ = False
        so.measureSpeedOfCode = False
        so.fastArchive = self.fastArchive
        # archiving time
        outputStepping = self.outputStepping
        tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]
        if len(tnList) > 2 and tnList[1] == tnList[2]:
            del tnList[1]
        if outputStepping['dt_output'] is None:
            if outputStepping['dt_fixed'] > 0:
                if outputStepping['dt_init'] < outputStepping['dt_fixed']:
                    tnList = [0., outputStepping['dt_init'], outputStepping['dt_fixed'], outputStepping['final_time']]
                else:
                    tnList = [0., outputStepping['dt_fixed'], outputStepping['final_time']]
            else:
                tnList = [0., outputStepping['dt_init'], outputStepping['final_time']]
        so.tnList = tnList
        so.archiveFlag = ArchiveFlags.EVERY_USER_STEP
        if self.archiveAllSteps is True:
            so.archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
        so.systemStepExact = outputStepping.systemStepExact

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
        self.dt_init = min(self.dt_output, self.dt_init)
        if self.nDTout is None:
            self.nDTout = int(round(old_div(self.final_time, self.dt_output)))
        else:
            self.dt_output = float(self.final_time)/float(self.nDTout)
#

class FESpace:
    """
    Create FE Spaces.
    """
    def __init__(self,nd):
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        self.nd=nd
        self.velSpaceOrder=None
        self.pSpaceOrder=None

    def __getitem__(self, key):
        return self.__dict__[key]

    def setFESpace(self,modelIdxDict):
        """
        User is allowed to define velSpaceOrder and pSpaceOrder. If not specified, defaults will be chosen based on 
        """
        
        # For now we just support rans2p with: p1-p1 and rans3p with: p2-p1
        if(self.velSpaceOrder is None and self.pSpaceOrder is None):
            if 'rans2p' in modelIdxDict: # rans2p or None
                self.velSpaceOrder=1
                self.pSpaceOrder=1
            elif 'rans3p' in modelIdxDict: #rans3p
                self.velSpaceOrder=2
                self.pSpaceOrder=1
            else:
                assert False, "the rans2p or rans3p model needs to be added to the problem"
        assert self.velSpaceOrder is not None, "need to set problem.FESpace.velSpaceOrder"
        assert self.pSpaceOrder is not None, "need to set problem.FESpace.pSpaceOrder"


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
                              'kinematicViscosityA': 1.004e-6,
                              'densityB': 1.205,
                              'kinematicViscosityB': 1.500e-5,
                              'surf_tension_coeff': 72.8E-3,
                              'gravity': [0.0, -9.8, 0.0]}

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
default_rans2p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 1.5,
                             'epsFact_density': 1.5,
                             'ns_forceStrongDirichlet': False,
                             'weak_bc_penalty_constant': 1.0E6,
                             'useRBLES': 0.0,
                             'ns_closure': 0,
                             'useVF': 1.0,
                             'ns_shockCapturingFactor': 0.25,
                             'ns_lag_shockCapturing': True,
                             'ns_lag_subgridError': True,
                             'timeDiscretization': 'vbdf',
                             'timeOrder': 2}
default_rans3p_parameters = {'useMetrics': 1.0,
                             'epsFact_viscosity': 1.5,
                             'epsFact_density': 1.5,
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
                             'USE_SUPG': False,
                             'ARTIFICIAL_VISCOSITY': 3, #edge based with smoothness indicator
                             'INT_BY_PARTS_PRESSURE': 1,
                             'cE': 1.0,
                             'cMax': 1.0}
default_clsvof_parameters = {'useMetrics': 1.0,
                             'epsFactHeaviside': 1.5,
                             'epsFactDirac': 1.5,
                             'epsFactRedist': 0.33,
                             'lambdaFact': 10.0,
                             'outputQuantDOFs': True,
                             'computeMetrics': 1,
                             'computeMetricsForBubble': False,
                             'eps_tolerance_clsvof': False,
                             'disc_ICs': False}
default_general = {'nLevels': 1,
                   'nLayersOfOverlapForParallel': 0}
