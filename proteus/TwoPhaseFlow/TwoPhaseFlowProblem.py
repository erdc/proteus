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
import collections
import sys

class TwoPhaseFlowProblem(Parameters.FreezableClass):
    """ TwoPhaseFlowProblem """

    def __init__(self,
                 # DOMAIN AND MESH #
                 domain=None,
                 mesh=None,
                 # INITIAL CONDITIONS #
                 initialConditions=None,
                 # BOUNDARY CONDITIONS #
                 boundaryConditions=None,
                 ):

        self.domain=domain
        self.mesh = mesh

        # ***** CREATE SYSTEM PHYSICS OBJECT ***** #
        self.SystemPhysics = SystemPhysics(ProblemInstance=self)
        self.SystemPhysics.initialConditions = initialConditions
        self.SystemPhysics.boundaryConditions= boundaryConditions

        #Model tracking objects
        self.SystemPhysics._modelIdxDict = {}           
        self.SystemPhysics.modelDict = collections.OrderedDict()           

        # ***** CREATE SYSTEM NUMERICS OBJECT ***** #
        self.SystemNumerics = SystemNumerics(ProblemInstance=self)

        # ***** CREATE MODEL PARAMETERS OBJECT ***** #
        self.Parameters = Parameters.ParametersHolder(ProblemInstance=self)

        # ***** CREATE FINITE ELEMENT SPACES ***** #
        self.FESpace = FESpace(ProblemInstance=self)

        # ***** CREATING OUTPUT MANAGEMENT OBJECTS ***** #
        self.so = System_base()
        self.outputStepping=OutputStepping()

        self._freeze()

    def checkProblem(self):
        # ***** SET OF ASSERTS ***** #
        assert self.domain.nd in [2,3], "nd={2,3}"
        assert self.SystemNumerics.cfl <= 1, "Choose cfl <= 1"
        assert type(self.domain.MeshOptions.he)==float , "Provide (float) he (characteristic mesh size)"
        assert self.domain is not None, "Provide a domain"
        if self.domain.MeshOptions.structured:
            assert type(self.domain.MeshOptions.nnx)==int and type(self.domain.MeshOptions.nny)==int, "Provide (int) nnx and nny"
            if self.domain.nd==3:
                assert type(self.domain.MeshOptions.nnz)==int, "Provide (int) nnz"
        else:
            assert self.domain.MeshOptions.triangleOptions != 'q30DenA', "Set domain.MeshOptions.triangleOptions"
        assert self.domain.MeshOptions.triangleFlag in [0,1,2], "triangleFlag must be 1, 2 or 3"
        if self.SystemPhysics.initialConditions is not None:
            assert type(self.SystemPhysics.initialConditions)==dict, "Provide dict of initial conditions"
            # assertion now done in TwoPhaseFlow_so.py
        if self.SystemPhysics.boundaryConditions is not None:
            assert type(self.SystemPhysics.boundaryConditions)==dict, "Provide dict of boundary conditions"

    def assert_initialConditions(self):
        initialConditions = self.SystemPhysics.initialConditions
        nd = self.domain.nd
        for model in self.SystemPhysics.modelDict.values():
            for key in model.p.initialConditions.__dict__.keys():
                if(key == 'name'):
                    continue
                assert model.p.initialConditions[key] is not None, 'Need to provide initial conditions for variable '+key+' in model '+model.name

    def assert_boundaryConditions(self):
        boundaryConditions = self.SystemPhysics.boundaryConditions
        nd = self.domain.nd
        ns_model = None
        ls_model = None

        for model in self.SystemPhysics.modelDict.values():
            if(isinstance(model,Parameters.ParametersModelRANS3PF)):                 
                ns_model = 'rans3p'
            elif(isinstance(model,Parameters.ParametersModelRANS2P)):
                ns_model = 'rans2p'
            elif isinstance(model,Parameters.ParametersModelCLSVOF):
                ls_model = 'clsvof'
            elif isinstance(model,Parameters.ParametersModelVOF):
                ls_model = 'vof'

        if boundaryConditions is not None:
            # check dirichlet BCs
            if ns_model is not None:
                assert 'pressure_DBC' in boundaryConditions, "Provide pressure_DBC"
                assert 'vel_u_DBC' in boundaryConditions, "Provide vel_u_DBC"
                assert 'vel_v_DBC' in boundaryConditions, "Provide vel_v_DBC"
                if nd==3:
                    assert 'vel_w_DBC' in boundaryConditions, "Provide vel_w_DBC"
            if ls_model == 'vof':
                assert 'vof_DBC' in boundaryConditions, "Provide vof_DBC"
                assert 'ncls_DBC' in boundaryConditions, "Provide ncls_DBC"
            elif ls_model == 'clsvof':
                assert 'clsvof_DBC' in boundaryConditions, "Provide clsvof_DBC"
            # check advective flux BCs
            if ns_model is not None:
                assert 'pressure_AFBC' in boundaryConditions, "Provide pressure_AFBC"
                assert 'vel_u_AFBC' in boundaryConditions, "Provide vel_u_AFBC"
                assert 'vel_v_AFBC' in boundaryConditions, "Provide vel_v_AFBC"
                if nd==3:
                    assert 'vel_w_AFBC' in boundaryConditions, "Provide vel_w_AFBC"
            if ls_model == 'vof':
                assert 'vof_AFBC' in boundaryConditions, "Provide vof_AFBC"
            if ls_model == 'clsvof':
                assert 'clsvof_AFBC' in boundaryConditions, "Provide clsvof_AFBC"
            # check diffusive flux BCs
            if ns_model is not None:
                assert 'vel_u_DFBC' in boundaryConditions, "provide vel_u_DFBC"
                assert 'vel_v_DFBC' in boundaryConditions, "provide vel_v_DFBC"
                if nd==3:
                    assert 'vel_w_DFBC' in boundaryConditions, "provide vel_w_DFBC"
            if ls_model == 'clsvof':
                assert 'clsvof_DFBC' in boundaryConditions, "provide clsvof_DFBC"
            if ns_model=='rans3p':
                # check dirichlet BCs
                assert 'pressure_increment_DBC' in boundaryConditions, "Provide pressure_increment_DBC"
                # check advective flux BCs
                assert 'pressure_increment_AFBC' in boundaryConditions,"Provide pressure_increment_AFBC"
                # check diffusive flux BCs
                assert 'pressure_increment_DFBC' in boundaryConditions,"Provide pressure_increment_DFBC"
        else:
            assert self.domain.useSpatialTools is True, 'Either define boundaryConditions dict or use proteus.mprans.SpatialTools to set Boundary Conditions and run function assembleDomain'

    def initializeAll(self):

        #Set dimension
        assert self.domain is not None, "Need to define domain before proceeding"

        # ***** SET FINITE ELEMENT  ***** #
        self.FESpace.setFESpace()

        self.outputStepping.setOutputStepping()
        #self.Parameters = Parameters.ParametersHolder(ProblemInstance=self)
        
        #preparing to extricate the mesh generation process from the workflow
        self.Parameters.mesh = self.domain.MeshOptions

        #model organization
        self.SystemPhysics.attachModels()
        self.Parameters.model_list=self.SystemPhysics.modelDict.values()
        
        self.checkProblem()
        # initial conditions
        self.assert_initialConditions()
        # boundary conditions
        self.assert_boundaryConditions()

        # parameters
        self.Parameters.initializeParameters()
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
        if 'rans3p' in self.SystemPhysics._modelIdxDict:
            PINIT_model = self.SystemPhysics._modelIdxDict['pressureInitial'] 
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
        if outputStepping.archiveAllSteps is True:
            so.archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
        so.systemStepExact = outputStepping.systemStepExact

class OutputStepping:
    """
    OutputStepping handles how often the solution is outputted.
    """
    def __init__(self):
        self.final_time=None 
        self.dt_init=0.001
        self.dt_output=None
        self.nDTout = None
        self.dt_fixed = None
        self.systemStepExact = False
        self.archiveAllSteps = False

    def __getitem__(self, key):
        return self.__dict__[key]

    def setOutputStepping(self):
        assert not (self.dt_output is None and nDTout is None), "Provide dt_output or nDTout"
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
    def __init__(self,ProblemInstance):
        self.velSpaceOrder=None
        self.pSpaceOrder=None
        self.Problem = ProblemInstance

    def __getitem__(self, key):
        return self.__dict__[key]

    def setFESpace(self):
        nd = self.Problem.domain.nd
        useExact = self.Problem.SystemNumerics.useExact
        assert nd in [2,3], 'number of dimensions must be 2 or 3'
        for key,value in self.Problem.SystemPhysics.modelDict.items():
            
            if(isinstance(value,Parameters.ParametersModelRANS2P)):
                self.velSpaceOrder=1
                self.pSpaceOrder=1
            elif(isinstance(value,Parameters.ParametersModelRANS3PF)):
                self.velSpaceOrder=2
                self.pSpaceOrder=1
        assert self.velSpaceOrder is not None

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
            if useExact:
                quadOrder=6
            else:
                quadOrder=3
            self.elementQuadrature = ft.SimplexGaussQuadrature(nd, quadOrder)
            self.elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, quadOrder)
        else:
            if useExact:
                quadOrder=6
            else:
                quadOrder=5
            self.elementQuadrature = ft.SimplexGaussQuadrature(nd, quadOrder)
            self.elementBoundaryQuadrature = ft.SimplexGaussQuadrature(nd - 1, quadOrder)

class SystemPhysics(Parameters.FreezableClass):
    def __init__(self,ProblemInstance):
        super(SystemPhysics, self).__init__(name='SystemPhysics')

        self._Problem = ProblemInstance
        self.gravity = None 
        
        #phase 0
        self.nu_0 = None
        self.rho_0 = None

        #phase 1
        self.nu_1 = None
        self.rho_1 = None
        
        self.surf_tension_coeff = None
        self.movingDomain = False

        self.initialConditions = None
        self.boundaryConditions=None

        # to use proteus.mprans.BoundaryConditions
        # but only if SpatialTools was used to make the domain
        self.useBoundaryConditionsModule = True

        self.useRANS=False

    def setDefaults(self):
        self.rho_0 = 998.2
        self.nu_0 = 1.004e-6
        self.rho_1 = 1.205
        self.nu_1 = 1.500e-5
        dim = self._Problem.domain.nd
        self.surf_tension_coeff = 72.8E-3
        if(dim==2):
            self.gravity =  [0.0, -9.81, 0.0]
        elif(dim==3):
            self.gravity =  [0.0, 0.0, -9.81]

    def useDefaultModels(self,flowModel=0,interfaceModel=0):
        """
        this provides the same functionality as the previous ns_model, ls_model flags
        but using the addModels API
        """

        modelNames = {
            'RANS2P': 'flow',
            'RANS3PF': 'flow',
            'VOF': 'vof',
            'LS': 'ncls',
            'RDLS': 'rdls',
            'MCorr': 'mcorr',
            'CLSVOF': 'clsvof',
            'PressureIncrement': 'pressureInc',
            'Pressure': 'pressure',
            'PressureInitial': 'pressureInit'
        }
 
        assert flowModel in [0,1] and interfaceModel in [0,1], "flowModel and interfaceModel must either be 0 or 1"
        if(flowModel == 0):
            self.addModel(Parameters.ParametersModelRANS2P,modelNames['RANS2P'])
        elif(flowModel == 1):
            self.addModel(Parameters.ParametersModelRANS3PF,modelNames['RANS3PF'])
            self.addModel(Parameters.ParametersModelPressureIncrement,modelNames['PressureIncrement'])
            self.addModel(Parameters.ParametersModelPressure,modelNames['Pressure'])
            self.addModel(Parameters.ParametersModelPressureInitial,modelNames['PressureInitial'])

        if(interfaceModel==0):
            self.addModel(Parameters.ParametersModelVOF,modelNames['VOF'])
            self.addModel(Parameters.ParametersModelNCLS,modelNames['LS'])
            self.addModel(Parameters.ParametersModelRDLS,modelNames['RDLS'])
            self.addModel(Parameters.ParametersModelMCorr,modelNames['MCorr']) 
        else:
            self.addModel(Parameters.ParametersModelCLSVOF,modelNames['CLSVOF'])
            if(flowModel == 1):
                self.modelDict.move_to_end(modelNames['CLSVOF'],last=False)

    def addModel(self,modelObject,name):
        #attach problem to model

        modelInitialized = modelObject(ProblemInstance=self._Problem) 
        self.modelDict[name] = modelInitialized

    def attachModels(self):
        #attach problem object to each model, identify index, and then form dictionary
        for (idx,model) in enumerate(self.modelDict.values()):
            #model._Problem = self
            model.index = idx
            #model._Problem.modelIdxDict[model.name]=idx
            self._modelIdxDict[model.name]=idx


class SystemNumerics(Parameters.FreezableClass):
    def __init__(self,ProblemInstance):
        super(SystemNumerics, self).__init__(name='SystemNumerics')

        self._Problem = ProblemInstance
        self.useSuperlu = True
        self.cfl = 0.33
        self.useExact = False
     
        # ***** DEFINE OTHER GENERAL NEEDED STUFF ***** #
        self.usePETScOptionsFileExternal = False
