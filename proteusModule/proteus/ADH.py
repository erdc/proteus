"""
A module for the python interface to cADH

The "level 0" interface is a function, cadhRun, that simply allows you
to to call pre_adh and adh from Python. This provides no ability to
interact with PyADH modules or use any archiving or visualization
features, but it should run at exactly the same speed as cADH from the
command line.

The "level 1" interface is provided by the ADH_NumericalSolution
class.  This class allows creation of an object representing an adh
run as well as access to the solution and mesh, which can then be
viewed and archived using PyADH tool. ADH solves a subset of equations
from a fixed set of 8 model types on a single level with a fixed
method. For that reason it is most clostly related to PyADH's
NumericalSolution abstraction, which represents a set of fully
discrete model equations. We therefore try to access most of ADH
through the ADH_NumericalSolution class. The other PyADH abstractions
will mostly be defined through the ADH_NumericalSolution object.

The "level 2" interface is provided by the ADH_OneLevelTransport
class. This class provides a way to define a fully discrete nonlinear
system and Jacobian using cADH's discretization, which can then be
solved by PyADH's NonlinearSolver and LinearSolver objects.
"""
from NumericalSolution import NS_base
from Transport import OneLevelTransport
import Comm
import cadh
import MeshTools,FemTools
import proteusGraphical
import Archiver
from Archiver import ArchiveFlags
import default_so,default_p,default_n
import Profiling
import os
log = Profiling.logEvent
import Viewers
import numpy

def cadhRun(filename,runname=None):
    """
    Call pre_adh and adh executables on input filename
    """
    import os
    import subprocess
    ADH_HOME = os.getenv("PYADH_PACKAGES")+"/adh"
    PRE_ADH = ADH_HOME+"/bin/pre_adh"
    ADH = ADH_HOME+"/bin/adh"
    PRE_ADH_CALL = [PRE_ADH,filename]
    ADH_CALL = [ADH,filename]
    if runname != None:
        ADH_CALL.append(runname)
    #pre_adh
    try:
        retcode = subprocess.check_call(PRE_ADH_CALL)
    except OSError:
        print "Calling preadh with "+`PRE_ADH_CALL`+" indicates failure but may have been succesful"
        pass
    except subprocess.CalledProcessError:
        print "Calling preadh with "+`PRE_ADH_CALL`+"failed, with return code",retcode
        raise subprocess.CalledProcessError
    #adh
    try:
        retcode = subprocess.check_call(ADH_CALL)
    except OSError:
        print "Calling adh with "+`ADH_CALL`+" indicates failure but may have been successful"
        pass
    except subprocess.CalledProcessError:
        print "Calling adh with "+`ADH_CALL`+"failed, with return code",retcode
        raise subprocess.CalledProcessError

class ADH_InputTranslator:
    """
    A class for making cADH style input files available to proteus

    For now this will just use an ADH_NumericalSolution object to allocate the needed information
    """
    def __init__(self,adhInput):
        os.system(os.getenv('PYADH_PACKAGES')+"/adh/bin/pre_adh "+`adhInput`);
        self.cinputTranslator = cadh.cADH_InputTranslator(adhInput)

class ADH_MultilevelTransport:
    """
    A dummy class since cADH solves on a single level
    """
    def __init__(self,transport):
        self.levelModelList=[transport]

class ADH_TransportCoefficients:
    """
    A dummy class since cADH doesn't provide access to the pde coefficients seperately from the discretization
    """
    def __init__(self,nc=1,names=['u'],vectorComponents=None,vectorName='adh_velocity'):
        self.nc = nc
        self.variableNames = names
        self.vectorComponents = vectorComponents
        self.vectorName=vectorName

class ADH_NumericalSolution(NS_base):
    r"""
    A class for NumericalSolutions calculated using pure cADH.
    
    For now use cADH i/o as well, but set it up so that we can use
    translations of so/p/n files when they're compatable.

    TODO:
    add partitioning step to keep adh's partition
    """
    def __init__(self,so=None,pList=None,nList=None,sList=None,opts=None,simFlagsList=None,adhInput=None,runname=None,comm=None,petscMatrix=False,adh_oneLevelTransport=None):
        #set the adh input filename and runname
        if adhInput != None:
            self.adhInput = adhInput
        elif so != None:
            self.adhInput = so.name
        elif pList != None:
            self.adhInput = pList[0].name
        if runname == None:
            runname = "default_runname"
        #get a comm reference to make sure mpifinalize has not been called before adh wrappers go out of scope
        if comm != None:
            self.comm = comm
        else:
            self.comm = Comm.get()
        #allocate the C extensions including the cADH global variables
        self.cadh_ns = cadh.cADH_NumericalSolution(self.adhInput,runname,comm)
        #construct archiver
        if so == None:
            self.so = default_so
        else:
            self.so = so
        #build a PyADH representation of the cADH mesh
        self.adhMesh = cadh.cADH_Mesh()
        #self.nd = self.adhMesh.nSpace_global
        #self.mesh = self.adhMesh.generateMeshToolsMesh()
        #todo control partitionMesh so that it matches adh partition
        #self.mesh.partitionMesh(nLayersOfOverlap=1,parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)
        self.mesh = None
        self.buildMesh(self.adhMesh)
        #build one level transport objects some of this below could be
        #rearranged and some of the stuff that currently happens in
        #ADH_OneLevelTransport should probably happen here and have
        #the result passed in. I'm not going to mess with it right now
        self.cadh_transport = cadh.cADH_OneLevelTransport(self.cadh_ns)
        #Archiver
        self.ar = {0:Archiver.XdmfArchive(opts.dataDir,
                                          adhInput,
                                          useTextArchive=opts.useTextArchive,
                                          gatherAtClose=opts.gatherArchive,
                                          hotStart=opts.hotStart)}
        self.archiveFlag= default_so.archiveFlag
        self.p = default_p
        vectorComponents_adh = self.cadh_transport.getVectorSolutionComponentIds()
        vectorComponents = None
        if len(vectorComponents_adh) > 0:
            vectorComponents = vectorComponents_adh
        self.p.coefficients = ADH_TransportCoefficients(self.cadh_transport.nc,['adh_u%s' % i for i in range(self.cadh_transport.nc)],vectorComponents=vectorComponents)
        self.coefficients   = self.p.coefficients
        if adh_oneLevelTransport == None:
            self.adh_transport = ADH_OneLevelTransport(adh_ns = self,petscMatrix=petscMatrix)
            self.owns_adh_transport = True
            #adh_transport ctor should be completed
            self.u = self.adh_transport.u 
            self.t = self.adh_transport.t
        else:
            self.adh_transport = adh_oneLevelTransport
            self.owns_adh_transport = False
            #may not be through with adh_transport's ctor
            self.u = None 
            self.t = None 
        self.nc = self.cadh_transport.nc#adh_transport grabs this from adh_ns.cadh_transport self.adh_transport.nc
        self.mlvt = ADH_MultilevelTransport(self.adh_transport)
        self.viewer = Viewers.V_base(p=self.p)
        
    def __del__(self):
        if self.owns_adh_transport:
            del self.adh_transport
        del self.cadh_transport
        del self.cadh_ns
    def calculateSolution(self):
        #need to deal with u here if constructed numerical solution part of the way through one-level transport ctor
        if self.u == None:
            self.u = self.adh_transport.u
        self.archiveInitialSolution()
        self.initializeViewSolution()
        self.tCount = 1
        self.t = self.cadh_transport.t
        stepping = True
        while stepping:
            self.cadh_ns.step()
            while(self.cadh_transport.solve()):
                log("ADH solver failed at t = " +`self.t`)
            self.archiveSolution()
            self.viewSolution()
            self.tCount+=1 
            stepping = self.cadh_ns.stepTaken()
            self.t = self.cadh_transport.t; self.adh_transport.t = self.t
            if self.cadh_transport.meshAdaptedForNextStep():
                self.adhMesh.update()
                self.buildMesh(self.adhMesh)
                self.mlvt.levelModelList[-1].mesh = self.mesh
                self.mlvt.levelModelList[-1].buildFiniteElementFunctions()
        self.ar[0].close()
        self.viewer.postprocess(self.mlvt,tsim=self.t)        
    def preStep(self,model):
        pass
    def postStep(self,model):
        pass
    def setWeakDirichletConditions(self,model):
        pass
    def restrictFromFineMesh(self,model):
        pass
    def archiveInitialSolution(self,model=None,index=None):
        import xml.etree.ElementTree as ElementTree
        self.ar[0].domain = ElementTree.SubElement(self.ar[0].tree.getroot(),"Domain")
        self.tCount=0
        for ci in range(self.adh_transport.nc):
            self.u[ci].femSpace.writeMeshXdmf(self.ar[0],self.u[ci].name,self.adh_transport.t0,init=True,meshChanged=True,tCount=self.tCount)
            self.u[ci].femSpace.writeFunctionXdmf(self.ar[0],self.u[ci],self.tCount)
        self.ar[0].sync();   
    def archiveSolution(self,model=None,index=None,t=None):
        #archive
        for ci in range(self.nc):
            self.u[ci].femSpace.writeMeshXdmf(self.ar[0],self.u[ci].name,self.t,init=self.cadh_transport.meshAdaptedForNextStep(),meshChanged=True,tCount=self.tCount)
            self.u[ci].femSpace.writeFunctionXdmf(self.ar[0],self.u[ci],self.tCount)
        self.ar[0].sync();  
    def closeArchive(self,model,index):
        if self.archiveFlag == None:
            return
        if self.so.useOneArchive:
            if index==0:
                log("Closing solution archive for "+self.adhInput)
                self.ar[index].close()
        else:
            log("Closing solution archive for "+self.adhInput)
            self.ar[index].close()
    def initializeViewSolution(self,model=None):
        self.viewer.preprocess(self.mlvt,self.t)
    def viewSolution(self,model=None,initialCondition=False):
        self.viewer.processTimeLevel(self.mlvt,tsim=self.t)
    def finalizeViewSolution(self,model):
        pass
    def buildMesh(self,adhMesh):
        if self.mesh != None:
            del self.mesh
        self.mesh = adhMesh.generateMeshToolsMesh()
        self.nd = adhMesh.nSpace_global
        #todo control partitionMesh so that it matches adh partition
        self.mesh.partitionMesh(nLayersOfOverlap=1,parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)

class ADH_OneLevelTransport(OneLevelTransport):
    r"""
    A class for model equations discretized in pure cADH

    I'm workiing towards first allowing this to be used by an NS_base
    numerical solution with the cADH nonlinear solver and time step
    control. Next I'll try to allow PyADH nonlinear solvers and time
    step control.
    """
    def __init__(self,
#                  uDict,
#                  phiDict,
#                  testSpaceDict,
#                  matType,
#                  dofBoundaryConditionsDict,
#                  dofBoundaryConditionsSetterDict,
#                  coefficients,
#                  elementQuadrature,
#                  elementBoundaryQuadrature,
#                  fluxBoundaryConditionsDict=None,
#                  advectiveFluxBoundaryConditionsSetterDict=None,
#                  diffusiveFluxBoundaryConditionsSetterDictDict=None,
#                  stressTraceBoundaryConditionsSetterDictDict=None,
#                  stabilization=None,
#                  shockCapturing=None,
#                  conservativeFluxDict=None,
#                  numericalFluxType=None,
#                  TimeIntegrationClass=None,
#                  massLumping=False,
#                  reactionLumping=False,
                  options=None,
#                  name='defaultName',
#                  reuse_trial_and_test_quadrature=True,
#                  sd = True,
#                  movingDomain=False
                 comm=None,
                 adh_ns=None,
                 adhInput=None,
                 runname=None,
                 petscMatrix=False):#,
        #use this approach to make sure mpifinalize has not been called before adh wrappers go out of scope
        if comm != None:
            self.comm = comm
        else:
            self.comm = Comm.get()
        #make sure we have adh global variables
        if adh_ns == None:
            self.adh_ns = ADH_NumericalSolution(adhInput=adhInput,runname=runname,comm=self.comm,opts=options,petscMatrix=petscMatrix,adh_oneLevelTransport=self)
            self.owns_adh_ns = True
        else:
            self.adh_ns = adh_ns
            self.owns_adh_ns = False
        self.coefficients = self.adh_ns.coefficients
        self.mesh = self.adh_ns.mesh
        self.nd = self.adh_ns.nd
        self.nc = self.adh_ns.cadh_transport.nc
        self.t0 = self.adh_ns.cadh_transport.t0
        self.t  = self.t0
        self.nSpace_global = self.adh_ns.adhMesh.nSpace_global
        ###finite element solutions
        self.trialSpaceDict = {}; self.u = {};
        #for now need a separate array for plotting vector valued components because of memory layout issues
        #with adh arrays of x,y,z structs
        self.u_arrays_for_plotting = {}
        self.buildFiniteElementFunctions()
        self.petscMatrix = petscMatrix
        if petscMatrix:
            self.cadh_petsc_interface = cadh.cADH_PETSc_Interface(self.adh_ns.cadh_transport)
            self.cadh_petsc_interface.update()
    def __del__(self):
        if self.petscMatrix:
            del self.cadh_petsc_interface
        if self.owns_adh_ns:
            self.adh_ns.__del__()
    def buildFiniteElementFunctions(self):
        for ci in range(self.nc):
            self.trialSpaceDict[ci]= FemTools.C0_AffineLinearOnSimplexWithNodalBasis(self.mesh.subdomainMesh,self.nd)
        for ci in range(self.nc):
            ci_dim = self.adh_ns.cadh_transport.getSolutionComponentDimension(ci)
            if ci_dim == 1: #scalars are easy
                self.u[ci] =FemTools.FiniteElementFunction(self.trialSpaceDict[ci],dof=self.adh_ns.cadh_transport.getSolutionComponent(ci),
                                                            dim_dof=1,name="adh_u_%s" % ci,
                                                            isVector=False)
                self.u_arrays_for_plotting[ci] = self.u[ci].dof
            else:
                #entire vector-valued solution for which ci is a part
                ci_base = self.adh_ns.cadh_transport.getFirstSolutionComponentForVectorComponent(ci)
                #this is a record array, stored as array of tuples, i
                #contains logical entries '0',...,'ci_dim-1' access all of logical component i as u_vec['i']
                u_vec = self.adh_ns.cadh_transport.getSolutionComponent(ci)
                #relative index for component ci in the vector-valued solution
                ci_adh = '%s' % (int(ci)-int(ci_base))
                self.u[ci] = FemTools.FiniteElementFunction(self.trialSpaceDict[ci],dof=u_vec[ci_adh],
                                                            dim_dof=1,name="adh_u_%s" % ci,
                                                            isVector=False)
                self.u_arrays_for_plotting[ci] = numpy.ascontiguousarray(self.u[ci].dof,dtype='d')
        
    def getResidual(self,u,r):
        #get u and r copied into ADH storage
        #the finite element function (self.u) 
        #call the the resudal function
        #need to know which ADH model to call
        self.r=r
        self.adh.gw_resid(self.iwhich_flag,self.bc_mask,r)
    def getJacobian(self,jacobian):
        #evaluate the difference jacobian
        #copy into petsc storage
        self.adh.gw_load(self.iwhich_flag,
                         self.bc_mask,
                         jacobian,
                         self.r)
    def calculateElementQuadrature(self):
        pass
    def calculateElementBoundaryQuadrature(self):
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        pass
    def estimate_mt(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
    def viewSolutionVTK(self,plotOffSet=None,titleModifier='',dgridnx=50,dgridny=50,dgridp=16.,
                        pause=False):
        """
        special purpose viewSolutionVTK for ADH wrapper necessary because of difficulty/inability to
          represent vector components as double* 's using shallow copies
        """
        import Viewers
	from proteusGraphical import vtkViewers
        if plotOffSet != None:
            windowNumberSave = Viewers.windowNumber
            Viewers.windowNumber=plotOffSet
        else:
            windowNumberSave = None
        if self.nSpace_global == 1:
            for ci in range(self.coefficients.nc):
                title = self.coefficients.variableNames[ci]+titleModifier
                if isinstance(self.u[ci].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
                    vtkViewers.viewScalar_1D(self.mesh.nodeArray[:,0],self.u[ci].dof,"x",self.u[ci].name,title,Viewers.windowNumber,
                                                Pause=pause,sortPoints=True)
                    Viewers.newPlot()
                    Viewers.newWindow()             
        elif self.nSpace_global == 2:
            for ci in range(self.coefficients.nc):                
                title = self.coefficients.variableNames[ci]+titleModifier
                ci_dim = self.adh_ns.cadh_transport.getSolutionComponentDimension(ci)
                if ci_dim > 1:
                    #force copy
                    self.u_arrays_for_plotting[ci][:] = self.u[ci].dof
                if (isinstance(self.u[ci].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis)):
                    vtkViewers.viewScalar_tri3_2D(self.mesh, 
                                                  self.u_arrays_for_plotting[ci][:self.mesh.nodeArray.shape[0]],
                                                  self.coefficients.variableNames[ci]+titleModifier,
                                                  Viewers.windowNumber,
                                                  viewTypes=['colorMapped'],#,'contour','warp'],
                                                  IsoSurface=True, 
                                                  Pause=pause,
                                                  Adapted = self.adh_ns.cadh_transport.meshAdaptedForNextStep())                    
                    Viewers.newPlot()
                    Viewers.newWindow()
            if self.coefficients.vectorComponents != None and not self.adh_ns.cadh_transport.meshAdaptedForNextStep(): #turn off for adaption right now?
                #
                scale_x = max(numpy.absolute(self.u[self.coefficients.vectorComponents[0]].dof.flat))
                scale_y = max(numpy.absolute(self.u[self.coefficients.vectorComponents[1]].dof.flat))
                L = min((max(self.mesh.nodeArray[:,0]),max(self.mesh.nodeArray[:,1])))
                scale=10.0*max((scale_x,scale_y,1.0e-16))/L
                #assume all components the same FemSpace for now
                if isinstance(self.u[self.coefficients.vectorComponents[0]].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
                    vtkViewers.viewVector_tri3_2D(self.mesh,
                                                  self.u[self.coefficients.vectorComponents[0]].dof,
                                                  self.u[self.coefficients.vectorComponents[1]].dof,
                                                  self.coefficients.vectorName+titleModifier)
                    Viewers.newPlot()
                    Viewers.newWindow()
        elif self.nSpace_global == 3:
            (slice_x,slice_y,slice_z) = self.mesh.nodeArray[self.mesh.nodeArray.shape[0]/2,:]
            for ci in range(self.coefficients.nc):
                title = self.coefficients.variableNames[ci]+titleModifier
                ci_dim = self.adh_ns.cadh_transport.getSolutionComponentDimension(ci)
                if ci_dim > 1: 
                    #force copy
                    self.u_arrays_for_plotting[ci][:] = self.u[ci].dof
                if isinstance(self.u[ci].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis):
                    vtkViewers.viewScalar_tet4_3D(self.mesh, 
                                                  self.u_arrays_for_plotting[ci],
                                                  self.coefficients.variableNames[ci]+titleModifier,
                                                  Viewers.windowNumber, Pause=pause,
                                                  Adapted = self.adh_ns.cadh_transport.meshAdaptedForNextStep())
                    Viewers.newPlot()
                    Viewers.newWindow()
            if self.coefficients.vectorComponents != None and not self.adh_ns.cadh_transport.meshAdaptedForNextStep():
                #
                scale_x = max(numpy.absolute(self.u[self.coefficients.vectorComponents[0]].dof.flat))
                scale_y = max(numpy.absolute(self.u[self.coefficients.vectorComponents[1]].dof.flat))
                scale_z = max(numpy.absolute(self.u[self.coefficients.vectorComponents[2]].dof.flat))
                L = min((max(self.mesh.nodeArray[:,0]),max(self.mesh.nodeArray[:,1]),max(self.mesh.nodeArray[:,2])))
                scale=10.0*max((scale_x,scale_y,scale_z,1.0e-16))/L
                #assume all components the same FemSpace for now
                if (isinstance(self.u[self.coefficients.vectorComponents[0]].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis) and
                    isinstance(self.u[self.coefficients.vectorComponents[1]].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis) and
                    isinstance(self.u[self.coefficients.vectorComponents[2]].femSpace,FemTools.C0_AffineLinearOnSimplexWithNodalBasis)) :
                    vtkViewers.viewVector_tet4_3D(self.mesh,
                                                  self.u[self.coefficients.vectorComponents[0]].dof,
                                                  self.u[self.coefficients.vectorComponents[1]].dof,
                                                  self.u[self.coefficients.vectorComponents[2]].dof,
                                                  self.coefficients.vectorName+titleModifier)
                    Viewers.newPlot()
                    Viewers.newWindow()

            #vector components 
        if windowNumberSave != None:
            Viewers.windowNumber = windowNumberSave
        return Viewers.windowNumber
        
    
