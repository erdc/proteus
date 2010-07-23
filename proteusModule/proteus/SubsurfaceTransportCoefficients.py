"""
TransportCoefficients for flow and transport in porous media
"""
from math import *
from TransportCoefficients import TC_base
import numpy
import Profiling
log = Profiling.logEvent
from proteus import FemTools
from proteus import Transport
import subsurfaceTransportFunctions as stfuncs

######################################################################
#Utility classes for dealing with common aspects of flow and transport
######################################################################

class BlockHeterogeneousCoefficients:
    """
    Basic data structures and functionality for keeping track of a 
    block heterogeneity 
    """
    def __init__(self,mesh):
        self.mesh = mesh
    def initializeMaterialTypes(self):
        """
        returns material type identifiers for mesh topology
          in tuple
        element,exterior_element_boundaries,element_boundaries

        note element_boundaries is nElementBoundaries_global x 2 and gives the 
         element material property to the left and right of a global element boundary 
        """
        elementMaterialTypes = self.mesh.elementMaterialTypes
        exteriorElementBoundaryTypes =  numpy.zeros((self.mesh.nExteriorElementBoundaries_global),'i')
        stfuncs.setExteriorElementBoundaryTypes(self.mesh.nExteriorElementBoundaries_global,
                                                self.mesh.exteriorElementBoundariesArray,
                                                self.mesh.elementBoundaryElementsArray,
                                                elementMaterialTypes,
                                                exteriorElementBoundaryTypes)
        elementBoundaryTypes = numpy.zeros((self.mesh.nElementBoundaries_global,2),'i')
        stfuncs.setElementBoundariesArray(self.mesh.nElementBoundaries_global,
                                          self.mesh.elementBoundaryElementsArray,
                                          self.mesh.elementMaterialTypes,
                                          elementBoundaryTypes)
        return elementMaterialTypes,exteriorElementBoundaryTypes,elementBoundaryTypes
##################################################
#Single-phase flow
##################################################
class SinglePhaseDarcyCoefficients(TC_base):
    """
    -\deld ( K_i(x,t) \grad h_i ) + r(x,t) = 0 i=1,nc
    """
    def __init__(self,K_types,source_types,nc=1,nd=2,
                 timeVaryingCoefficients=False,
                 materialValuesLocallyConstant=False):
        self.K_types = K_types
        self.source_types = source_types
        self.nd = nd
        self.timeVaryingCoefficients=timeVaryingCoefficients
        self.materialValuesLocallyConstant = materialValuesLocallyConstant
        self.K_types_const=None
        self.source_types_const=None
        if self.materialValuesLocallyConstant:
            self.K_types_const = {}; self.source_types_const= {}
            x = numpy.zeros((3,),'d'); t0 = 0.0
            for mat in self.K_types.keys():
                self.K_types_const[mat]     =self.K_types[mat](x,t0)
            for mat in self.source_types.keys():
                self.source_types_const[mat]=self.source_types[mat](x,t0) 
        #for matching shapes of quadrature arrays
        self.q_x_shape = None;  self.ebqe_x_shape=None
        self.ebq_global_x_shape=None;  self.ebq_shape=None
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            reaction[i]  = {i : 'constant'}
            advection[i] = {i : 'constant'} #now include for gravity type terms
            potential[i] = {i : 'u'}
        #end i
        sparseDiffusionTensors = {}
        for ci in range(nc):
            sparseDiffusionTensors[(ci,ci)]=(numpy.arange(start=0,stop=nd**2+1,step=nd,dtype='i'),
                                             numpy.array([range(nd) for row in range(nd)],dtype='i'))
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)

        
    def initializeMesh(self,mesh):
        self.elementMaterialTypes,self.exteriorElementBoundaryTypes,self.elementBoundaryTypes = BlockHeterogeneousCoefficients(mesh).initializeMaterialTypes()
        self.elementBoundariesArray = mesh.elementBoundariesArray
    def initializeElementQuadrature(self,t,cq):
        self.evaluateHeterogeneity_element(t,cq)
        self.q_x_shape = cq['x'].shape
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.evaluateHeterogeneity_globalElementBoundary(t,cebq_global)
        self.ebq_global_x_shape = cebq_global['x'].shape
        self.evaluateHeterogeneity_elementBoundary(t,cebq)
        self.ebq_x_shape = cebq['x'].shape
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.evaluateHeterogeneity_exteriorElementBoundary(t,cebqe)
        self.ebqe_x_shape = cebqe['x'].shape
    def evaluateHeterogeneity_element(self,t,cq):
        if self.materialValuesLocallyConstant:
            x = numpy.zeros((3,),'d'); t0 = 0.0
            for mat in self.K_types.keys():
                self.K_types_const[mat]     =self.K_types[mat](x,t0)
            for mat in self.source_types.keys():
                self.source_types_const[mat]=self.source_types[mat](x,t0) 
            for ci in range(self.nc):
                cq[('f',ci)].flat[:] = 0.0

                stfuncs.setScalarMaterialFunctionOverElements(self.elementMaterialTypes,
                                                              cq[('r',ci)],
                                                              self.source_types_const)
                cq[('r',ci)] *= -1.0
                #sparse matrix storage for tensors
                stfuncs.setVectorMaterialFunctionOverElements(self.elementMaterialTypes,
                                                              cq[('a',ci,ci)],
                                                              self.K_types_const)
        else:
            for ci in range(self.nc):
                cq[('f',ci)].flat[:] = 0.0

                stfuncs.evaluateScalarMaterialFunctionOverElements(t,self.elementMaterialTypes,
                                                                  cq['x'],
                                                                  cq[('r',ci)],
                                                                  self.source_types)
                cq[('r',ci)] *= -1.0
                #sparse matrix storage for tensors
                stfuncs.evaluateVectorMaterialFunctionOverElements(t,self.elementMaterialTypes,
                                                                   cq['x'],
                                                                   cq[('a',ci,ci)],
                                                                   self.K_types)
    def evaluateHeterogeneity_elementBoundary(self,t,cebq):
        #use harmonic average for a, arith average for f
        if self.materialValuesLocallyConstant:
            x = numpy.zeros((3,),'d'); t0 = 0.0
            for mat in self.K_types.keys():
                self.K_types_const[mat]     =self.K_types[mat](x,t0)
            for mat in self.source_types.keys():
                self.source_types_const[mat]=self.source_types[mat](x,t0) 
            for ci in range(self.nc):
                if cebq.has_key(('f',ci)): cebq[('f',ci)].flat[:] = 0.0
                if cebq.has_key(('r',ci)):
                    stfuncs.setScalarMaterialFunctionOverElementBoundaries_arithmeticAverage(self.elementBoundariesArray,
                                                                                             self.elementBoundaryTypes,
                                                                                             cebq[('r',ci)],
                                                                                             self.source_types_const)
                    cebq[('r',ci)] *= -1.0
                if cebq.has_key(('a',ci,ci)):
                    stfuncs.setSparseTensorMaterialFunctionOverElementBoundaries_harmonicAverage(self.nd,
                                                                                                 self.elementBoundariesArray,
                                                                                                 self.elementBoundaryTypes,
                                                                                                 cebq[('a',ci,ci)],
                                                                                                 self.K_types_const)
        
        else:
            for ci in range(self.nc):
                if cebq.has_key(('f',ci)): cebq[('f',ci)].flat[:] = 0.0
                if cebq.has_key(('r',ci)):
                    stfuncs.evaluateScalarMaterialFunctionOverElementBoundaries_arithmeticAverage(t,
                                                                                                  self.elementBoundariesArray,
                                                                                                  self.elementBoundaryTypes,
                                                                                                  cebq['x'],
                                                                                                  cebq[('r',ci)],
                                                                                                  self.source_types)
                    cebq[('r',ci)] *= -1.0
                if cebq.has_key(('a',ci,ci)):
                    stfuncs.evaluateSparseTensorMaterialFunctionOverElementBoundaries_harmonicAverage(self.nd,
                                                                                                      t,
                                                                                                      self.elementBoundariesArray,
                                                                                                      self.elementBoundaryTypes,
                                                                                                      cebq['x'],
                                                                                                      cebq[('a',ci,ci)],
                                                                                                      self.K_types)
    def evaluateHeterogeneity_globalElementBoundary(self,t,cebq_global):
        #use harmonic average for a, arith average for f
        if self.materialValuesLocallyConstant:
            x = numpy.zeros((3,),'d'); t0 = 0.0
            for mat in self.K_types.keys():
                self.K_types_const[mat]     =self.K_types[mat](x,t0)
            for mat in self.source_types.keys():
                self.source_types_const[mat]=self.source_types[mat](x,t0) 
            for ci in range(self.nc):
                if cebq_global.has_key(('f',ci)): cebq_global[('f',ci)].flat[:] = 0.0
                if cebq_global.has_key(('r',ci)):
                    stfuncs.setScalarMaterialFunctionOverGlobalElementBoundaries_arithmeticAverage(self.elementBoundariesArray,
                                                                                                   self.elementBoundaryTypes,
                                                                                                   cebq_global[('r',ci)],
                                                                                                   self.source_types_const)
                    cebq_global[('r',ci)] *= -1.0
                if cebq_global.has_key(('a',ci,ci)):
                    stfuncs.setSparseTensorMaterialFunctionOverGlobalElementBoundaries_harmonicAverage(self.nd,
                                                                                                       self.elementBoundariesArray,
                                                                                                       self.elementBoundaryTypes,
                                                                                                       cebq_global[('a',ci,ci)],
                                                                                                       self.K_types_const)
        
        else:
            for ci in range(self.nc):
                if cebq_global.has_key(('f',ci)): cebq_global[('f',ci)].flat[:] = 0.0
                if cebq_global.has_key(('r',ci)):
                    stfuncs.evaluateScalarMaterialFunctionOverGlobalElementBoundaries_arithmeticAverage(t,
                                                                                                        self.elementBoundariesArray,
                                                                                                        self.elementBoundaryTypes,
                                                                                                        cebq_global['x'],
                                                                                                        cebq_global[('r',ci)],
                                                                                                        self.source_types)
                    cebq_global[('r',ci)] *= -1.0
                if cebq_global.has_key(('a',ci,ci)):
                    stfuncs.evaluateSparseTensorMaterialFunctionOverGlobalElementBoundaries_harmonicAverage(self.nd,
                                                                                                            t,
                                                                                                            self.elementBoundariesArray,
                                                                                                            self.elementBoundaryTypes,
                                                                                                            cebq_global['x'],
                                                                                                            cebq_global[('a',ci,ci)],
                                                                                                            self.K_types)

    def evaluateHeterogeneity_exteriorElementBoundary(self,t,cebqe):
        if self.materialValuesLocallyConstant:
            x = numpy.zeros((3,),'d'); t0 = 0.0
            for mat in self.K_types.keys():
                self.K_types_const[mat]     =self.K_types[mat](x,t0)
            for mat in self.source_types.keys():
                self.source_types_const[mat]=self.source_types[mat](x,t0) 
            for ci in range(self.nc):
                cebqe[('f',ci)].flat[:] = 0.0

                stfuncs.setScalarMaterialFunctionOverElements(self.exteriorElementBoundaryTypes,
                                                              cebqe[('r',ci)],
                                                              self.source_types_const)
                cebqe[('r',ci)] *= -1.0
                #sparse matrix storage for tensors
                stfuncs.setVectorMaterialFunctionOverElements(self.exteriorElementBoundaryTypes,
                                                              cebqe[('a',ci,ci)],
                                                              self.K_types_const)
        else:
            for ci in range(self.nc):
                cebqe[('f',ci)].flat[:] = 0.0

                stfuncs.evaluateScalarMaterialFunctionOverElements(t,self.exteriorElementBoundaryTypes,
                                                                  cebqe['x'],
                                                                  cebqe[('r',ci)],
                                                                  self.source_types)
                cebqe[('r',ci)] *= -1.0
                #sparse matrix storage for tensors
                stfuncs.evaluateVectorMaterialFunctionOverElements(t,self.exteriorElementBoundaryTypes,
                                                                   cebqe['x'],
                                                                   cebqe[('a',ci,ci)],
                                                                   self.K_types)
    def evaluate(self,t,c):
        if self.timeVaryingCoefficients == True:
            if c['x'].shape == self.q_x_shape:
                self.evaluateHeterogeneity_element(t,c)
            elif c['x'].shape == self.ebqe_x_shape:
                self.evaluateHeterogeneity_exteriorElementBoundary(t,c)
            elif c['x'].shape == self.ebq_global_x_shape:
                self.evaluateHeterogeneity_globalElementBoundary(t,c)
            elif c['x'].shape == self.ebq_x_shape:
                self.evaluateHeterogeneity_elementBoundary(t,c)
            else:
                raise NotImplementedError
            #mwf debug
            #print "eval c[('a',ci,ci)]= %s" % c[('a',0,0)]
    #end def
########################################
#Richards equation
########################################
class ConservativeHeadRichardsMualemVanGenuchten(TC_base):
    """
    version of Re where element material type id's used in evals
    """
    from ctransportCoefficients import conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2
    from ctransportCoefficients import conservativeHeadRichardsMualemVanGenuchten_sd_het
    def __init__(self,
                 nd,
                 Ksw_types,
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 gravity,
                 density,
                 beta,
                 diagonal_conductivity=True,
                 getSeepageFace=None):
        variableNames=['pressure_head']
        nc=1
        mass={0:{0:'nonlinear'}}
        advection={0:{0:'nonlinear'}}
        diffusion={0:{0:{0:'nonlinear'}}}
        potential={0:{0:'u'}}
        reaction={0:{0:'linear'}}
        hamiltonian={}
        self.getSeepageFace=getSeepageFace
        self.gravity=gravity
        self.rho = density
        self.beta=beta
        self.vgm_n_types = vgm_n_types
        self.vgm_alpha_types = vgm_alpha_types
        self.thetaR_types    = thetaR_types
        self.thetaSR_types   = thetaSR_types
        self.elementMaterialTypes = None
        self.exteriorElementBoundaryTypes  = None
        self.materialTypes_q    = None
        self.materialTypes_ebq  = None
        self.materialTypes_ebqe  = None
        self.nd = nd
        self.nMaterialTypes = len(thetaR_types)
        #try to allow some flexibility in input of permeability/conductivity tensor
        self.diagonal_conductivity = diagonal_conductivity
        self.Ksw_types_in = Ksw_types
        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

            assert len(Ksw_types.shape) in [1,2], "if diagonal conductivity true then Ksw_types scalar or vector of diagonal entries"
            #allow scalar input Ks
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I] = Ksw_types
            else:
                self.Ksw_types = Ksw_types
        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}
            assert len(Ksw_types.shape) in [1,2], "if full tensor conductivity true then Ksw_types scalar or 'flattened' row-major representation of entries"
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd**2),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I*self.nd+I] = Ksw_types
            else:
                assert Ksw_types.shape[1] == self.nd**2
                self.Ksw_types = Ksw_types
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
            
    def initializeMesh(self,mesh):
        #import pdb
        #pdb.set_trace()
        self.elementMaterialTypes,self.exteriorElementBoundaryTypes,self.elementBoundaryTypes = BlockHeterogeneousCoefficients(mesh).initializeMaterialTypes()
        #want element boundary material types for evaluating heterogeneity
        #not boundary conditions
        self.isSeepageFace = numpy.zeros((mesh.nExteriorElementBoundaries_global),'i')
        if self.getSeepageFace != None:
            for ebNE in range(mesh.nExteriorElementBoundaries_global):
                #mwf missing ebNE-->ebN?
                ebN = mesh.exteriorElementBoundariesArray[ebNE]
                #print "eb flag",mesh.elementBoundaryMaterialTypes[ebN]
                #print self.getSeepageFace(mesh.elementBoundaryMaterialTypes[ebN])
                self.isSeepageFace[ebNE] = self.getSeepageFace(mesh.elementBoundaryMaterialTypes[ebN])
        #print self.isSeepageFace
    def initializeElementQuadrature(self,t,cq):
        self.materialTypes_q = self.elementMaterialTypes
        self.q_shape = cq[('u',0)].shape
#        cq['Ks'] = numpy.zeros(self.q_shape,'d')
#        for k in range(self.q_shape[1]):
#            cq['Ks'][:,k] = self.Ksw_types[self.elementMaterialTypes,0]
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.materialTypes_ebq = numpy.zeros(cebq[('u',0)].shape[0:2],'i')
        self.ebq_shape = cebq[('u',0)].shape
        for ebN_local in range(self.ebq_shape[1]):
            self.materialTypes_ebq[:,ebN_local] = self.elementMaterialTypes
        #for eN in range(self.elementMaterialTypes.shape[0]):
        #    self.materialTypes_ebq[eN,:] = self.elementMaterialTypes[eN]
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.materialTypes_ebqe = self.exteriorElementBoundaryTypes
        self.ebqe_shape = cebqe[('u',0)].shape
        #
    def evaluate(self,t,c):
        if c[('u',0)].shape == self.q_shape:
            materialTypes = self.materialTypes_q
        elif c[('u',0)].shape == self.ebqe_shape:
            materialTypes = self.materialTypes_ebqe
        elif c[('u',0)].shape == self.ebq_shape:
            materialTypes = self.materialTypes_ebq
        else:
            assert False, "no materialType found to match c[('u',0)].shape= %s " % c[('u',0)].shape
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.conservativeHeadRichardsMualemVanGenuchten_sd_het(self.sdInfo[(0,0)][0],
                                                               self.sdInfo[(0,0)][1],
                                                               materialTypes,
                                                               self.rho,
                                                               self.beta,
                                                               self.gravity,
                                                               self.vgm_alpha_types,
                                                               self.vgm_n_types,
                                                               self.thetaR_types,
                                                               self.thetaSR_types,
                                                               self.Ksw_types,
                                                               c[('u',0)],
                                                               c[('m',0)],
                                                               c[('dm',0,0)],
                                                               c[('f',0)],
                                                               c[('df',0,0)],
                                                               c[('a',0,0)],
                                                               c[('da',0,0,0)])

#         self.conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(materialTypes,
#                                                                      self.rho,
#                                                                      self.beta,
#                                                                      self.gravity,
#                                                                      self.vgm_alpha_types,
#                                                                      self.vgm_n_types,
#                                                                      self.thetaR_types,
#                                                                      self.thetaSR_types,
#                                                                      self.Ksw_types,
#                                                                      c[('u',0)],
#                                                                      c[('m',0)],
#                                                                      c[('dm',0,0)],
#                                                                      c[('f',0)],
#                                                                      c[('df',0,0)],
#                                                                      c[('a',0,0)],
#                                                                      c[('da',0,0,0)])
        #mwf debug
        if (numpy.isnan(c[('da',0,0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('df',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any() or
            numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any()):
            import pdb
            pdb.set_trace()

#         #mwf debug
#         if c[('u',0)].shape == self.q_shape:
#             c[('visPerm',0)]=c[('a',0,0)][:,:,0,0]

from proteus import cfemIntegrals
class RE_NCP1_OneLevelTransport(Transport.OneLevelTransport):
    """
    OneLevelTransport designed specifically for Non-Conforming P^1 approximation to RE
    Approximation uses nodal quadrature and upwinding
    """
    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=False,
                 sd = True,
                 movingDomain=False):#,
        Transport.OneLevelTransport.__init__(self,
                                             uDict,
                                             phiDict,
                                             testSpaceDict,
                                             matType,
                                             dofBoundaryConditionsDict,
                                             dofBoundaryConditionsSetterDict,
                                             coefficients,
                                             elementQuadrature,
                                             elementBoundaryQuadrature,
                                             fluxBoundaryConditionsDict=fluxBoundaryConditionsDict,
                                             advectiveFluxBoundaryConditionsSetterDict=advectiveFluxBoundaryConditionsSetterDict,
                                             diffusiveFluxBoundaryConditionsSetterDictDict=diffusiveFluxBoundaryConditionsSetterDictDict,
                                             stressTraceBoundaryConditionsSetterDict=stressTraceBoundaryConditionsSetterDict,
                                             stabilization=stabilization,
                                             shockCapturing=shockCapturing,
                                             conservativeFluxDict=conservativeFluxDict,
                                             numericalFluxType=numericalFluxType,
                                             TimeIntegrationClass=TimeIntegrationClass,
                                             massLumping=massLumping,
                                             reactionLumping=reactionLumping,
                                             options=options,
                                             name=name,
                                             reuse_trial_and_test_quadrature=reuse_trial_and_test_quadrature,
                                             sd = sd,
                                             movingDomain=movingDomain)
        
        assert self.nQuadraturePoints_element == self.nSpace_global+1
        #particular arrays needed for Darcy flow
        self.q[('k_r',0,0)]    = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dk_r',0,0,0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('f_lin',0)]    = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('a_lin',0,0)]  = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.coefficients.sdInfo[(0,0)][0][self.nSpace_global]),'d')
        
        #need to check that points are the interpolation points too

#         #quadrature dictionary for solution interpolation points
#         self.n_u_ip_element = [u_j.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  u_j in self.u.values()]
#         self.u_ip = {}
#         #coefficients
#         self.points_u_ip =StorageSet(shape=(self.mesh.nElements_global,self.n_u_ip_element[0],3))
#         self.scalars_u_ip=StorageSet(shape=(self.mesh.nElements_global,self.n_u_ip_element[0]))
#         self.vectors_u_ip=StorageSet(shape=(self.mesh.nElements_global,self.n_u_ip_element[0],self.nSpace_global))
#         self.tensors_u_ip=StorageSet(shape={})
#         #shape
#         self.trial_shape_u_ip=StorageSet(shape={})
#         self.test_shape_u_ip = StorageSet(shape={}) 
#         self.test_shapeGradient_u_ip = StorageSet(shape={})
#         self.trial_shapeGradient_u_ip= StorageSet(shape={})
#         trial_shape_u_ip_duplicate = StorageSet(shape={})
#         trial_shape_u_ip_duplicate_map = {}
        
#         #
#         # build sets
#         #
#         #u
#         self.scalars_u_ip |= set([('u',ci) for ci  in range(self.nc)])
#         self.vectors_u_ip |= set([('grad(u)',ci) for ci in range(self.nc)]) 
#         #trial
#         self.trial_shape_u_ip |= set([('v',ci) for ci in self.unique_test_trial_range])
#         trial_shape_u_ip_duplicate |= set([('v',ci) for ci in self.duplicate_test_trial_range])
#         for ci in self.duplicate_test_trial_range:
#             trial_shape_u_ip_duplicate_map[('v',ci)] = ('v',self.unique_test_trial_range[0]) 
#         #test
# 	self.test_shape_u_ip |= set([('w',ci) for ci in range(self.nc)])
#         self.test_shapeGradient_u_ip |= set([('grad(w)',ci) for ci in range(self.nc)])
#         #m
#         self.scalars_u_ip |= set([('m',ci) for ci  in self.coefficients.mass.keys()])
#         self.scalars_u_ip |= set([('mt',ci) for ci  in self.coefficients.mass.keys()])
#         for ci,cjDict in self.coefficients.mass.iteritems():
#             self.scalars_u_ip |= set([('dm',ci,cj) for cj  in cjDict.keys()])
#         #f
#         self.vectors_u_ip |= set([('f',ci) for ci  in self.coefficients.advection.keys()])
#         #special for Darcy flow
#         self.vectors_u_ip |= set([('f_lin',ci) for ci  in self.coefficients.advection.keys()])
#         for ci,cjDict in self.coefficients.advection.iteritems():
#             self.vectors_u_ip |= set([('df',ci,cj) for cj  in cjDict.keys()])
#         #a and phi
#         self.vectors_u_ip |= set([('velocity',ci) for ci  in self.coefficients.diffusion.keys()])
#         self.scalars_u_ip |= set([('phi',ck) for ck  in self.coefficients.potential.keys()])
#         for ci,ckDict in self.coefficients.diffusion.iteritems():
#             self.tensors_u_ip |= set([('a',ci,ck) for ck  in ckDict.keys()])
#             #special for Darcy flow
#             self.tensors_u_ip |= set([('a_lin',ci,ck) for ck  in ckDict.keys()])
#             self.scalars_u_ip |= set([('k_r',ci,ck) for ck  in self.ckDict.keys()])
            
#             self.vectors_u_ip |= set([('grad(phi)',ck) for ck  in ckDict.keys()])            
#             for ck,cjDict in ckDict.iteritems():
#                 self.tensors_u_ip |= set([('da',ci,ck,cj) for cj  in cjDict.keys()])
#                 #special for Darcy flow
#                 self.scalars_u_ip |= set([('dk_r',ci,ck,cj) for cj  in cjDict.keys()])
#                 #mwf trial-dup start
# 		#mwf orig self.trial_shape_quadrature |= set([('v',cj) for cj in cjDict.keys()])
#                 unique_cj   = set.intersection(set(self.unique_test_trial_range),set(cjDict.keys()))
#                 duplicate_cj= set.intersection(set(self.duplicate_test_trial_range),set(cjDict.keys()))
# 		self.trial_shape_u_ip |= set([('v',cj) for cj in unique_cj])
# 		trial_shape_u_ip_duplicate |= set([('v',cj) for cj in duplicate_cj])
#                 #any reason to do this again
#                 for cj in self.duplicate_test_trial_range:
#                     trial_shape_u_ip_duplicate_map[('v',cj)] = ('v',self.unique_test_trial_range[0]) 
#                 #mwf trial-dup end
#                 for cj in self.coefficients.potential[ck].keys():
#                     self.scalars_u_ip |= set([('dphi',ck,cj)])
#         #mwf what if miss nonlinearities in diffusion that don't match potentials?
#         for comb in self.dphi.keys():
#             self.scalars_u_ip |= set([('dphi',comb[0],comb[1])])
#         #r
#         self.scalars_u_ip |= set([('r',ci) for ci  in self.coefficients.reaction.keys()])
#         for ci,cjDict in self.coefficients.reaction.iteritems():
#             self.scalars_u_ip |= set([('dr',ci,cj) for cj  in cjDict.keys()])
#         #mesh
#         self.points_u_ip |= set([('x',ci) for ci in self.unique_test_trial_range])
#         #
#         memory()
#         #
#         for ci in self.unique_test_trial_range:
#             self.u_ip[('x',ci)] = self.u[ci].femSpace.updateInterpolationPoints()
#         for ci in self.duplicate_test_trial_range:
#             self.u_ip[('x',ci)] = self.u_ip[('x',self.unique_test_trial_range[0])]
#         #
#         #allocations
#         #
#         self.scalars_u_ip.allocate(self.u_ip)
#         #points 
#         self.points_u_ip.allocate(self.u_ip)
#         #scalars
#         self.scalars_u_ip.allocate(self.u_ip)
#         #vectors
#         self.vectors_u_ip.allocate(self.u_ip)
#         #allocate tensors element quadrature
#         for k in self.tensors_u_ip:
#             if (self.sd
#                 and k[0] in ['a','da']
#                 and self.coefficients.sdInfo != None
#                 and (k[1],k[2]) in self.coefficients.sdInfo.keys()):                
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.n_u_ip_element[k[2]],
#                      self.coefficients.sdInfo[(k[1],k[2])][0][self.nSpace_global]),
#                     'd')
#             else:
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self..n_u_ip_element[k[2]],
#                      self.nSpace_global,
#                      self.nSpace_global),
#                     'd')
#         #
#         # shape
#         #
#         def makeAlias(sd,kw,vstring='v',wstring='w'):
#             aliased=False
#             kv0 = kw[0].replace(wstring,vstring)
#             if len(kw) > 1:
#                 kv = (kv0,)+kw[1:]
#             else:
#                 kv = kv0
#             if self.testIsTrial:
#                 if kv in sd.keys():
#                     log("Shallow copy of trial shape is being used for test shape %s " % kw[0],level=4)
#                     sd[kw] = sd[kv]
#                     aliased=True
#             return aliased
#         def makeAliasForComponent1(sd,k,entryStrings,refi=0):
#             aliased=False
#             #need to add combinatorial options
#             ks = k[0]
#             if ks not in entryStrings:
#                 return False
#             k0 = (ks,)+(refi,)
#             if self.reuse_test_trial_quadrature and refi != None:
#                 if k0 in sd.keys():
#                     log("Shallow copy of trial shape %s is being used for trial shape %s" % (k0,k),level=4)
#                     sd[k] = sd[k0]
#                     aliased=True
#             return aliased
#         #
#         # element quadrature
#         #
#         #allocate trial shape functions
#         for k in sorted(self.trial_shape_u_ip):
#             if not makeAliasForComponent1(self.u_ip,k,['v'],refi=0):#need to handle multiple component combinations 
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.n_u_ip_element,
#                      self.nDOF_trial_element[k[-1]]),
#                     'd')
#         for k in trial_shape_quadrature_duplicate:
#             self.u_ip[k] = self.u_ip[trial_shape_quadrature_duplicate_map[k]]
#         #mwf trial-dup end
#         log(memory("solution interpolation points, test/trial functions trial_shape","OneLevelTransport"),level=4)
#         #allocate test shape functions 
#         for k in self.test_shape_quadrature:
#             if not makeAlias(self.u_ip,k):
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.n_u_ip_element,
#                      self.nDOF_test_element[k[-1]]),
#                     'd')
#         log(memory("solution interpolation points, test/trial functions test_shape","OneLevelTransport"),level=4)
#         #allocate trial shape function gradients
#         for k in sorted(self.trial_shapeGradient_quadrature):
#             if not makeAliasForComponent1(self.u_ip,k,['grad(v)'],refi=0):#need to handle multiple component combinations 
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.n_u_ip_element,
#                      self.nDOF_trial_element[k[-1]],
#                      self.nSpace_global),
#                     'd')        
#         log(memory("solution interpolation points, test/trial functions trial_shapeGradient","OneLevelTransport"),level=4)
#         #allocate test shape function gradients
#         for k in self.test_shapeGradient_quadrature:
#             if not makeAlias(self.u_ip,k):
#                 self.u_ip[k]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.n_u_ip_element,
#                      self.nDOF_test_element[k[-1]],
#                      self.nSpace_global),
#                     'd')
    def calculateElementCoefficients(self):
        """
        calculate the nonlinear coefficients at the quadrature points and nodes 
        include interpolation points explicitly here now
        """
        
	for cj in range(self.nc):
	    self.u[cj].getValues(self.q[('v',cj)],
				 self.q[('u',cj)])
	    if self.q.has_key(('grad(u)',cj)):
		self.u[cj].getGradientValues(self.q[('grad(v)',cj)],
					     self.q[('grad(u)',cj)])
        
        #can skip this after first call
        stfuncs.RE_NCP1_evaluateElementCoefficients_Linear(self.coefficients.rho,
                                                                                self.coefficients.gravity,
                                                                                self.coefficients.sdInfo[(0,0)][0],
                                                                                self.coefficients.sdInfo[(0,0)][1],
                                                                                self.coefficients.Ksw_types,
                                                                                self.nSpace_global,
                                                                                self.mesh.nElements_global,
                                                                                self.mesh.nElementBoundaries_element,
                                                                                self.mesh.elementNeighborsArray,
                                                                                self.mesh.elementMaterialTypes,
                                                                                self.q[('f_lin',0)],
                                                                                self.q[('a_lin',0,0)])

        stfuncs.RE_NCP1_evaluateElementCoefficients_VGM(self.coefficients.rho,
                                                                             self.coefficients.beta,
                                                                             self.coefficients.gravity,
                                                                             self.coefficients.vgm_alpha_types,
                                                                             self.coefficients.vgm_n_types,
                                                                             self.coefficients.thetaR_types,
                                                                             self.coefficients.thetaSR_types,
                                                                             self.nSpace_global,
                                                                             self.mesh.nElements_global,
                                                                             self.mesh.nElementBoundaries_element,
                                                                             self.mesh.elementMaterialTypes,
                                                                             self.nDOF_trial_element[0],
                                                                             self.u[0].femSpace.dofMap.l2g,
                                                                             self.u[0].dof,
                                                                             self.q['x'],
                                                                             self.q[('u',0)],
                                                                             self.q[('m',0)],
                                                                             self.q[('dm',0,0)],
                                                                             self.q[('r',0)],
                                                                             self.q[('k_r',0,0)],
                                                                             self.q[('dk_r',0,0,0)])

        if self.movingDomain and self.coefficients.movingDomain:
            self.coefficients.updateToMovingDomain(self.timeIntegration.t,self.q)
        if self.timeTerm:
            self.timeIntegration.calculateElementCoefficients(self.q)
        #cek need to clean up calculation of dimless numbers, might as well do it all the time and pass to subgrid error
        #mwf figure out what to do with this
        #what happens if stabilization didn't compute cfl?
        for ci in range(self.nc):
            #for two phase flow would need this
            self.q[('dphi',ci,ci)].fill(1.0)
            if self.sd:
                cfemIntegrals.calculateDimensionlessNumbersADR_sd(self.mesh.nElements_global,
                                                                  self.nQuadraturePoints_element,
                                                                  self.nSpace_global,
                                                                  self.coefficients.sdInfo[(ci,ci)][0],self.coefficients.sdInfo[(ci,ci)][1],
                                                                  self.elementEffectiveDiametersArray, 
                                                                  self.q[('df',ci,ci)],
                                                                  self.q[('a',ci,ci)],
                                                                  self.q[('dphi',ci,ci)],
                                                                  self.q[('dr',ci,ci)],
                                                                  self.q[('dmt',ci,ci)],
                                                                  self.q[('pe',ci)],
                                                                  self.q[('cfl',ci)])
            else:
                cfemIntegrals.calculateDimensionlessNumbersADR(self.mesh.nElements_global,
                                                               self.nQuadraturePoints_element,
                                                               self.nSpace_global,
                                                               self.elementEffectiveDiametersArray, 
                                                               self.q[('df',ci,ci)],
                                                               self.q[('a',ci,ci)],
                                                               self.q[('dphi',ci,ci)],
                                                               self.q[('dr',ci,ci)],
                                                               self.q[('dmt',ci,ci)],
                                                               self.q[('pe',ci)],
                                                               self.q[('cfl',ci)])
                    
        if self.shockCapturing != None:
            self.shockCapturing.calculateNumericalDiffusion(self.q)

    def calculateElementResidual(self):
        """Calculate all the element residuals"""
        import pdb
        #mwf debug
        #Transport.OneLevelTransport.calculateElementResidual(self)
        #self.elementResidual_save = numpy.copy(self.elementResidual[0])
        for ci in range(self.nc):
            self.elementResidual[ci].fill(0.0)
        #pdb.set_trace()
        stfuncs.RE_NCP1_getElementResidual(self.coefficients.gravity,
                                           self.coefficients.sdInfo[(0,0)][0],
                                           self.coefficients.sdInfo[(0,0)][1],
                                           self.nSpace_global,
                                           self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_global,
                                           self.mesh.elementNeighborsArray,
                                           self.mesh.elementBarycentersArray,
                                           self.nDOF_test_element[0],
                                           self.q[('u',0)],
                                           self.q[('grad(u)',0)],
                                           self.q[('grad(w)',0)],
                                           self.q['abs(det(J))'],
                                           self.q[('m',0)],
                                           self.q[('mt',0)],
                                           self.q[('r',0)],
                                           self.q[('k_r',0,0)],
                                           self.q[('f_lin',0)],
                                           self.q[('a_lin',0,0)],
                                           self.elementResidual[0])

        
    def calculateElementJacobian(self):
        import pdb
        #Transport.OneLevelTransport.calculateElementJacobian(self)
        #self.elementJacobian_save = numpy.copy(self.elementJacobian[0][0])
        for ci in range(self.nc):
            for cj in self.coefficients.stencil[ci]:
                self.elementJacobian[ci][cj].fill(0.0)

        #mwf debug
        #pdb.set_trace()
        stfuncs.RE_NCP1_getElementJacobian(self.coefficients.gravity,
                                           self.coefficients.sdInfo[(0,0)][0],
                                           self.coefficients.sdInfo[(0,0)][1],
                                           self.nSpace_global,
                                           self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_global,
                                           self.mesh.elementNeighborsArray,
                                           self.mesh.elementBarycentersArray,
                                           self.nDOF_test_element[0],
                                           self.nDOF_trial_element[0],
                                           self.q[('u',0)],
                                           self.q[('grad(u)',0)],
                                           self.q[('grad(w)',0)],
                                           self.q[('grad(v)',0)],
                                           self.q['abs(det(J))'],
                                           self.q[('m',0)],
                                           self.q[('dm',0,0)],
                                           self.q[('mt',0)],
                                           self.q[('dmt',0,0)],
                                           self.q[('r',0)],
                                           self.q[('k_r',0,0)],
                                           self.q[('dk_r',0,0,0)],
                                           self.q[('f_lin',0)],
                                           self.q[('a_lin',0,0)],
                                           self.elementJacobian[0][0])

########################################
#Immiscible two-phase flow
########################################
class TwophaseDarcyFlow_base(TC_base):
    """
    base class for two-phase flow implementations
    holds information for 

        fluid properties
        eos model tag
        psk model tag
        material heterogeneity information (number of types, lookup arrays)
    default fluids are water and air assumed to be incompressible
    """
    default_density_w_parameters = {'model':'Exponential',
                                    'nParameters':3,
                                    'rho_0':998.2,#Water kg/m^3
                                    'psi_0':0.0,
                                    'beta':0.0}
    default_density_n_parameters = {'model':'Exponential',
                                    'nParameters':3,
                                    'rho_0':1.205,#Air   kg/m^3
                                    'psi_0':0.0,
                                    'beta':0.0}
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=default_density_w_parameters,
                 density_n_parameters=default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 diagonal_conductivity=True): #otherwise Full
        self.nd = nd
        self.g  = dimensionless_gravity
        ##fluid properties
        #normalized density terms
        self.rhon = density_n/density_w
        self.rhow = density_w/density_w
        self.b  = density_n/density_w
        #normalized viscosities
        self.muw= viscosity_w/viscosity_w
        self.mun= viscosity_n/viscosity_w
        #density eos options
        self.density_types={'Exponential':1,
                            'IdealGas':2}
        self.density_w_model = density_w_parameters['model']
        self.density_n_model = density_n_parameters['model']
        assert self.density_w_model in self.density_types
        assert self.density_n_model in self.density_types
        self.density_w_parameters = density_w_parameters
        self.density_n_parameters = density_n_parameters
        for params,rwork in zip([self.density_w_parameters,self.density_n_parameters],
                                ['rwork_density_w','rwork_density_n']):
            if params != None:
                if params['model'] == 'Exponential':
                    setattr(self,rwork,numpy.array([params['rho_0']/params['rho_0'],#normalize by phase density
                                                    params['psi_0'],
                                                    params['beta']],dtype='d'))
                elif params['model'] == 'IdealGas':
                    setattr(self,rwork,numpy.array([params['T'],
                                                    params['W']/params['rho_0'],#normalize by phase density
                                                    params['R'],
                                                    params['headToPressure'],
                                                    params['rho_0']/params['rho_0'],#normalize by phase density
                                                    params['psi_0']],dtype='d'))
                else:
                    assert False, 'TwophaseDarcy_base density params= %s not found ' % params
        #

        ##psk relations and heterogeneity
        self.psk_model=psk_model
        self.psk_types={'simp':0,
                        'VGM':1,
                        'VGB':2,
                        'BCM':3,
                        'BCB':4}
        self.psk_tolerances={'default':{'eps_small':1.0e-16},
                             'VGM':{'eps_small':1.0e-16,'ns_del':1.0e-8}}
        self.nPskTolerances=1
        assert(self.psk_model in self.psk_types.keys())
        self.nMaterialTypes = nMaterialTypes
        #psk rwork array lengths 
        if self.psk_model == 'simp':
            self.nPskParams=2
        elif self.psk_model in ['VGM','VGB']:
            self.nPskParams=4
            if self.psk_model == 'VGM':
                self.nPskTolerances=2
        elif self.psk_model in ['BCM','BCB']:
            self.nPskParams=4
        self.rwork_psk = None
        self.rwork_psk_tolerances = numpy.zeros((self.nPskTolerances,),'d')
        for i,tol in enumerate(self.psk_tolerances[self.psk_model].values()):
            self.rwork_psk_tolerances[i] = tol
        self.diagonal_conductivity=diagonal_conductivity
    def initializeMesh(self,mesh):
        self.elementMaterialTypes,self.exteriorElementBoundaryTypes,self.elementBoundaryTypes = BlockHeterogeneousCoefficients(mesh).initializeMaterialTypes()
        self.mesh = mesh #for debugging
    def initializeElementQuadrature(self,t,cq):
        self.materialTypes_q = self.elementMaterialTypes
        self.q_shape = cq[('u',0)].shape
        #
        cq['psi_n'] = numpy.zeros(cq[('u',0)].shape,'d')
        cq[('dpsi_n',0)] = numpy.zeros(cq[('u',0)].shape,'d')
        cq[('dpsi_n',1)] = numpy.zeros(cq[('u',0)].shape,'d')
        #for visualization
        cq['Ks'] = numpy.zeros(self.q_shape,'d')
        for k in range(self.q_shape[1]):
            cq['Ks'][:,k] = self.Ksw_types[self.elementMaterialTypes,0]
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.materialTypes_ebq = numpy.zeros(cebq[('u',0)].shape[0:2],'i')
        self.ebq_shape = cebq[('u',0)].shape
        for ebN_local in range(self.ebq_shape[1]):
            self.materialTypes_ebq[:,ebN_local] = self.elementMaterialTypes
        if cebq.has_key(('u',0)):
            cebq['psi_n'] = numpy.zeros(cebq[('u',0)].shape,'d')
            cebq[('dpsi_n',0)] = numpy.zeros(cebq[('u',0)].shape,'d')
            cebq[('dpsi_n',1)] = numpy.zeros(cebq[('u',0)].shape,'d')
        if cebq_global.has_key(('u',0)):
            cebq_global['psi_n'] = numpy.zeros(cebq_global[('u',0)].shape,'d')
            cebq_global[('dpsi_n',0)] = numpy.zeros(cebq_global[('u',0)].shape,'d')
            cebq_global[('dpsi_n',1)] = numpy.zeros(cebq_global[('u',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.materialTypes_ebqe = self.exteriorElementBoundaryTypes
        self.ebqe_shape = cebqe[('u',0)].shape
        #
        cebqe['psi_n'] = numpy.zeros(cebqe[('u',0)].shape,'d')
        cebqe[('dpsi_n',0)] = numpy.zeros(cebqe[('u',0)].shape,'d')
        cebqe[('dpsi_n',1)] = numpy.zeros(cebqe[('u',0)].shape,'d')
    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        self.materialTypes_ip = self.elementMaterialTypes
        self.ip_shape = cip[('u',0)].shape
        #
        cip['psi_n'] = numpy.zeros(cip[('u',0)].shape,'d')
        cip[('dpsi_n',0)] = numpy.zeros(cip[('u',0)].shape,'d')
        cip[('dpsi_n',1)] = numpy.zeros(cip[('u',0)].shape,'d')
    def setMaterialTypes(self,
                         Ksw_types=[1.0],
                         omega_types  = [0.4],
                         Sw_max_types = [1.0],
                         Sw_min_types = [0.0],
                         bc_lambda_types = None,
                         bc_pd_types = None,
                         vg_alpha_types = None,
                         vg_m_types = None):
        self.nMaterialTypes=len(omega_types)
        self.omega_types= omega_types
        if self.psk_model == 'simp':
            assert self.nPskParams == 2
            self.rwork_psk = numpy.zeros((self.nMaterialTypes,self.nPskParams),'d')
            for Sw_min,Sw_max,i in zip(Sw_min_types,Sw_max_types,range(self.nMaterialTypes)):
                self.rwork_psk[i,0] = Sw_min
                self.rwork_psk[i,1] = Sw_max
        elif self.psk_model in ['VGM','VGB']:
            assert(vg_alpha_types != None and vg_m_types  != None)
            assert self.nPskParams == 4
            self.rwork_psk = numpy.zeros((self.nMaterialTypes,self.nPskParams),'d')
            for Sw_min,Sw_max,vg_alpha,vg_m,i in zip(Sw_min_types,Sw_max_types,vg_alpha_types,vg_m_types,range(self.nMaterialTypes)):
                self.rwork_psk[i,0] = Sw_min
                self.rwork_psk[i,1] = Sw_max
                self.rwork_psk[i,2] = vg_alpha
                self.rwork_psk[i,3] = vg_m
        elif self.psk_model in ['BCM','BCB']:
            assert(bc_lambda_types != None and bc_lambda_types  != None)
            assert self.nPskParams == 4
            self.rwork_psk = numpy.zeros((self.nMaterialTypes,self.nPskParams),'d')
            for Sw_min,Sw_max,bc_pd,bc_lambda,i in zip(Sw_min_types,Sw_max_types,bc_pd_types,bc_lambda_types,range(self.nMaterialTypes)):
                self.rwork_psk[i,0] = Sw_min
                self.rwork_psk[i,1] = Sw_max
                self.rwork_psk[i,2] = bc_pd
                self.rwork_psk[i,3] = bc_lambda
        #try to allow some flexibility in input of permeability/conductivity tensor
        #keep track of input for debugging now
        self.Ksw_types_in = Ksw_types
        if self.diagonal_conductivity:
            assert len(Ksw_types.shape) in [1,2], "if diagonal conductivity true then Ksw_types scalar or vector of diagonal entries"
            #allow scalar input Ks
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I] = Ksw_types
            else:
                self.Ksw_types = Ksw_types
        else: #full
            assert len(Ksw_types.shape) in [1,2], "if full tensor conductivity true then Ksw_types scalar or or 'flattened' row-major representation of entries"
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd**2),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I*self.nd+I] = Ksw_types
            else:
                assert Ksw_types.shape[1] == self.nd**2
                self.Ksw_types = Ksw_types
        #not the best place for this
        self.rwork_psk_tolerances = numpy.zeros((self.nPskTolerances,),'d')
        for i,tol in enumerate(self.psk_tolerances[self.psk_model].values()):
            self.rwork_psk_tolerances[i] = tol
#         #mwf do some debugging
#         assert self.elementMaterialTypes.max() < self.nMaterialTypes
#         assert self.elementMaterialTypes.min() == 0
#         assert self.exteriorElementBoundaryTypes.max()  < self.nMaterialTypes
#         assert self.exteriorElementBoundaryTypes.min() == 0
#         assert self.elementBoundaryTypes.max() <  self.nMaterialTypes
#         assert self.elementBoundaryTypes.min() == 0
    #

class TwophaseDarcy_fc(TwophaseDarcyFlow_base):
    """
    continuity equation for each phase

      \pd{m_w}{t} - \deld (\ten{a}_w \grad \phi_w) + r_w = 0

      \pd{m_n}{t} - \deld (\ten{a}_n \grad \phi_n) + r_n = 0

    (normalized) mass for each phase 
       m_i = \theta_i \rho_i,  i=w,n
       
       \theta_i = \theta_s S_i
       \rho_i   = \varrho_i/\varrho_{i,0}, i=w,n

     where S_i is the saturation for each phase, \varrho_{i} is the density and
     \varrho_{i,0} is a reference value

    (normalized) mass flux for each phase is

       \vec \sigma_i = - \ten{a}_{i}\grad \phi_i,    i=w,n
    
       \ten{a}_i = \rho_i k_{r,i}/\hat{\mu}_i \ten{K}_s
    and the potentials are defined as
       \phi_{w}  = \psi_{w} - \rho_w \vec g\cdot \vec x
       \phi_{n}  = \psi_{n} - b \rho_n \vec g\cdot \vec x

    where
       \psi_{i} = p_{i}/|\vec g|\rho_{w,0}
       b        = \varrho_{n,0}/\varrho_{w,0}
       \hat{mu}_{i} = \mu_i/\mu_{w}

    for the pressure of each phase, p_i, and we have the capillary pressure relation
       \psi_{c} = \psi_{n} - \psi_{w}

    The dependent variables are
      S_w, and \psi_w

    Note S_n = 1-S_w

    needs to be implemented
    r_i = r_i(\vec x,t) 


    slight compressibility assumed in spatial gradients

    TODO:
   """
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_fc_sd_het_matType
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=TwophaseDarcyFlow_base.default_density_w_parameters,
                 density_n_parameters=TwophaseDarcyFlow_base.default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 diagonal_conductivity=True):
        TwophaseDarcyFlow_base.__init__(self,
                                        nd=nd,
                                        dimensionless_gravity=dimensionless_gravity,
                                        density_w=density_w,
                                        density_n=density_n,
                                        viscosity_w=viscosity_w,
                                        viscosity_n=viscosity_n,
                                        density_w_parameters=density_w_parameters,
                                        density_n_parameters=density_n_parameters,
                                        psk_model=psk_model,
                                        nMaterialTypes=nMaterialTypes,
                                        diagonal_conductivity=diagonal_conductivity)

        nc=2
        variableNames=['s_w','psi_w']
        mass = {0:{0:'linear',1:'nonlinear'},
                1:{0:'nonlinear',1:'nonlinear'}}
        advection = {}
        hamiltonian={}
        potential = {0:{1:'nonlinear',#actually this is just linear if we use slight compressibility
                        0:'nonlinear'},#don't really need saturation dependence for potential
                     1:{0:'nonlinear',
                        1:'nonlinear'}}
        diffusion = {0:{0:{0:'nonlinear',
                           1:'nonlinear'}},
                     1:{1:{0:'nonlinear',
                           1:'nonlinear'}}}
        reaction = {0:{0:'linear'},
                    1:{1:'linear'}}
        

        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i')),
                                      (1,1):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i')),
                                      (1,1):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}


        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
			 variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)


    def evaluate(self,t,c):
        if c[('u',0)].shape == self.q_shape:
            materialTypes = self.materialTypes_q
        elif c[('u',0)].shape == self.ebqe_shape:
            materialTypes = self.materialTypes_ebqe
        elif c[('u',0)].shape == self.ip_shape:
            materialTypes = self.materialTypes_ip
        elif c[('u',0)].shape == self.ebq_shape:
            materialTypes = self.materialTypes_ebq
        else:
            assert False, "no materialType found to match c[('u',0)].shape= %s " % c[('u',0)].shape
        assert self.rwork_psk != None
        #mwf do some debugging
        assert materialTypes.max() < self.nMaterialTypes
        assert materialTypes.min() == 0
        
        self.twophaseDarcy_fc_sd_het_matType(self.psk_types[self.psk_model],
                                             self.density_types[self.density_w_model],
                                             self.density_types[self.density_n_model],
                                             self.sdInfo[(0,0)][0],
                                             self.sdInfo[(0,0)][1],
                                             materialTypes,
                                             self.muw,
                                             self.mun,
                                             self.omega_types,
                                             self.Ksw_types,
                                             self.b,
                                             self.rwork_psk,
                                             self.rwork_psk_tolerances,
                                             self.rwork_density_w,
                                             self.rwork_density_n,
                                             self.g[:self.nd],#todo get consistent on dimension setting
                                             c['x'],
                                             c[('u',0)],
                                             c[('u',1)],
                                             c[('m',0)],
                                             c[('dm',0,0)],
                                             c[('dm',0,1)],
                                             c[('m',1)],
                                             c[('dm',1,0)],
                                             c[('dm',1,1)],
                                             c['psi_n'],
                                             c[('dpsi_n',0)],
                                             c[('dpsi_n',1)],
                                             c[('phi',0)],
                                             c[('dphi',0,1)],
                                             c[('phi',1)],
                                             c[('dphi',1,1)],
                                             c[('dphi',1,0)],
                                             c[('a',0,0)],
                                             c[('da',0,0,0)],
                                             c[('da',0,0,1)],
                                             c[('a',1,1)],
                                             c[('da',1,1,0)],
                                             c[('da',1,1,1)])
 

        #mwf debug
        if (numpy.isnan(c[('da',0,0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('da',1,1,0)]).any() or
            numpy.isnan(c[('a',1,1)]).any() or
            numpy.isnan(c[('u',0)]).any() or
            numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any()):
            import pdb
            pdb.set_trace()
#

class FullyCoupledMualemVanGenuchten(TwophaseDarcy_fc):
    """
    Formulation using phase continuity equations and
     Van-Genuchten Mualem psk relations

    Basically a convenience wrapper for fully coupled approximation
     with volume-fraction based inputs as in Richards' equation formulations

    """
    def __init__(self,
                 nd,
                 Ksw_types, 
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 dimensionless_gravity,
                 density_w,
                 density_n,
                 viscosity_w,
                 viscosity_n,
                 density_w_params,
                 density_n_params,
                 diagonal_conductivity=True,
                 vgm_small_eps=1.0e-16,
                 vgm_ns_del=1.0e-8):

        TwophaseDarcy_fc.__init__(self,
                                  nd=nd,
                                  dimensionless_gravity=dimensionless_gravity,
                                  density_w=density_w,
                                  density_n=density_n,
                                  viscosity_w=viscosity_w,
                                  viscosity_n=viscosity_n,
                                  density_w_parameters=density_w_params,
                                  density_n_parameters=density_n_params,
                                  psk_model='VGM',
                                  nMaterialTypes=len(thetaR_types),
                                  diagonal_conductivity=diagonal_conductivity)
        for input in [vgm_n_types,vgm_alpha_types,thetaR_types,thetaSR_types]:
            assert len(input)==self.nMaterialTypes
        
        vgm_m_types    = 1.0-1.0/vgm_n_types
        thetaS_types   = thetaSR_types + thetaR_types
        Sw_max_types   = numpy.ones((self.nMaterialTypes,),'d')
        Sw_min_types   = thetaR_types/thetaS_types
        self.psk_tolerances['VGM']['eps_small']=vgm_small_eps
        self.psk_tolerances['VGM']['ns_del']   =vgm_ns_del
        self.setMaterialTypes(Ksw_types=Ksw_types,
                              omega_types=thetaS_types,
                              Sw_max_types=Sw_max_types,
                              Sw_min_types=Sw_min_types,
                              vg_alpha_types=vgm_alpha_types,
                              vg_m_types=vgm_m_types)


#
class TwophaseDarcy_split_pressure_base(TwophaseDarcyFlow_base):
    """
    Base class for 'pressure' or total flow conservation equation in fractional flow formulations. This

    The primary functionality of the base class is to handle synchronization with a 'saturation' model to
    get the saturation, S_w, and capillary pressure (head), \psi_c, variables 

   """
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=TwophaseDarcyFlow_base.default_density_w_parameters,
                 density_n_parameters=TwophaseDarcyFlow_base.default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nSatModel=1,
                 diagonal_conductivity=True,
                 #for debugging
                 swConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcyFlow_base.__init__(self,
                                        nd=nd,
                                        dimensionless_gravity=dimensionless_gravity,
                                        density_w=density_w,
                                        density_n=density_n,
                                        viscosity_w=viscosity_w,
                                        viscosity_n=viscosity_n,
                                        density_w_parameters=density_w_parameters,
                                        density_n_parameters=density_n_parameters,
                                        psk_model=psk_model,
                                        nMaterialTypes=nMaterialTypes,
                                        diagonal_conductivity=diagonal_conductivity)

    	self.nc=1
        self.nd=nd
        #
        self.nSatModel=nSatModel
        #for debugging
        self.swConstant=swConstant
        self.capillaryDiffusionScaling=capillaryDiffusionScaling
    def attachModels(self,modelList):
        if self.nSatModel == None:
            print 'Warning TwophaseDarcy_split_pressure_base nSatModel == None returning in attachModels'
            return
        #not ideal, but need a way to force nonlinear potential to be evaluated in saturation model
        modelList[self.nSatModel].calculateElementCoefficients()
        self.q_s_w   = modelList[self.nSatModel].q[('u',0)]
        self.ebqe_s_w = modelList[self.nSatModel].ebqe[('u',0)]
        if modelList[self.nSatModel].ebq.has_key(('u',0)):
            self.ebq_s_w = modelList[self.nSatModel].ebq[('u',0)]
        #mwf need to check
        assert modelList[self.nSatModel].phi_ip.has_key(('u',0))
        assert self.ip_s_w.shape ==  modelList[self.nSatModel].phi_ip[('u',0)].shape
        self.ip_s_w = modelList[self.nSatModel].phi_ip[('u',0)]
        #mwf hack 04/03/09 skip
        #assert modelList[self.nSatModel].phi_ip.has_key(('grad(phi)',0))
        self.ip_grad_psic = None#modelList[self.nSatModel].phi_ip[('grad(phi)',0)]
        #mwf end ip stuff
        self.q_grad_psic   = modelList[self.nSatModel].q[('grad(phi)',0)]
        self.ebqe_grad_psic = modelList[self.nSatModel].ebqe[('grad(phi)',0)]
        if modelList[self.nSatModel].ebq.has_key(('grad(phi)',0)):
            self.ebq_grad_psic = modelList[self.nSatModel].ebq[('grad(phi)',0)]
        self.q_psic   = modelList[self.nSatModel].q[('phi',0)]
        self.ebqe_psic= modelList[self.nSatModel].ebqe[('phi',0)]
        if modelList[self.nSatModel].ebq.has_key(('phi',0)):
            self.ebq_psic = modelList[self.nSatModel].ebq[('phi',0)]
        assert modelList[self.nSatModel].phi_ip.has_key(('phi',0))
        assert self.ip_psic.shape ==  modelList[self.nSatModel].phi_ip[('phi',0)].shape
        self.ip_psic = modelList[self.nSatModel].phi_ip[('phi',0)]

    def initializeElementQuadrature(self,t,cq):
        TwophaseDarcyFlow_base.initializeElementQuadrature(self,t,cq)
    	#set up dummy values in case we're not running the other model
        self.q_s_w   = numpy.zeros(cq[('u',0)].shape,'d')
        self.q_s_w[:] = self.swConstant
        for i in range(len(self.q_s_w.flat)/2,len(self.q_s_w.flat)):
                self.q_s_w.flat[i] = 1.0e-4
        self.q_grad_psic   = numpy.zeros(cq[('f',0)].shape,'d')
        self.q_psic        = numpy.zeros(cq[('u',0)].shape,'d')
        #mwf not sure if this is ok
        cq['psi_n'] = numpy.zeros(cq[('u',0)].shape,'d')
        cq[('dpsi_n',0)] = numpy.ones(cq[('u',0)].shape,'d')
        
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        TwophaseDarcyFlow_base.initializeElementBoundaryQuadrature(self,t,cebq,cebq_global)
    	#set up dummy values in case we're not running the other model
        self.ebq_s_w = numpy.zeros(cebq[('u',0)].shape,'d')
        self.ebq_s_w[:]=self.swConstant
        for i in range(len(self.ebq_s_w.flat)/2,len(self.ebq_s_w.flat)):
                self.ebq_s_w.flat[i] = 1.0e-4
        self.ebq_grad_psic = numpy.zeros(cebq[('f',0)].shape,'d')
        self.ebq_psic = numpy.zeros(cebq[('u',0)].shape,'d')
        if cebq.has_key(('u',0)):
            cebq['psi_n'] = numpy.zeros(cebq[('u',0)].shape,'d')
            cebq[('dpsi_n',0)] = numpy.ones(cebq[('u',0)].shape,'d')
        if cebq_global.has_key(('u',0)):
            cebq_global['psi_n'] = numpy.zeros(cebq_global[('u',0)].shape,'d')
            cebq_global[('dpsi_n',0)] = numpy.ones(cebq_global[('u',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        TwophaseDarcyFlow_base.initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe)
    	#set up dummy values in case we're not running the other model
        self.ebqe_s_w = numpy.zeros(cebqe[('u',0)].shape,'d')
        self.ebqe_s_w[:]=self.swConstant
        for i in range(len(self.ebqe_s_w.flat)/2,len(self.ebqe_s_w.flat)):
                self.ebqe_s_w.flat[i] = 1.0e-4
        self.ebqe_grad_psic = numpy.zeros(cebqe[('f',0)].shape,'d')
        self.ebqe_psic = numpy.zeros(cebqe[('u',0)].shape,'d')
        cebqe['psi_n'] = numpy.zeros(cebqe[('u',0)].shape,'d')
        cebqe[('dpsi_n',0)] = numpy.ones(cebqe[('u',0)].shape,'d')
    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        TwophaseDarcyFlow_base.initializeGeneralizedInterpolationPointQuadrature(self,t,cip)
    	#set up dummy values in case we're not running the other model
        self.ip_s_w = numpy.zeros(cip[('u',0)].shape,'d')
        self.ip_s_w[:]=self.swConstant
        for i in range(len(self.ip_s_w.flat)/2,len(self.ip_s_w.flat)):
                self.ip_s_w.flat[i] = 1.0e-4
        self.ip_grad_psic = numpy.zeros(cip[('f',0)].shape,'d')
        self.ip_psic = numpy.zeros(cip[('u',0)].shape,'d')
        cip['psi_n'] = numpy.zeros(cip[('u',0)].shape,'d')
        cip[('dpsi_n',0)] = numpy.ones(cip[('u',0)].shape,'d')


#
class TwophaseDarcy_split_saturation_base(TwophaseDarcyFlow_base):
    """
    Base class for aqueous phase mass conservation equation (saturation equation) in
    a fractional flow formulation.

    The primary responsibility of the base class is to handle
    synchronization with the 'pressure' equation to get the total flow
    velocity variable, q_t, and aqueous phase pressure head, psi_w

   """
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=TwophaseDarcyFlow_base.default_density_w_parameters,
                 density_n_parameters=TwophaseDarcyFlow_base.default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nPressModel=0,
                 diagonal_conductivity=True,
                 #for debugging
                 qScalarConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcyFlow_base.__init__(self,
                                        nd=nd,
                                        dimensionless_gravity=dimensionless_gravity,
                                        density_w=density_w,
                                        density_n=density_n,
                                        viscosity_w=viscosity_w,
                                        viscosity_n=viscosity_n,
                                        density_w_parameters=density_w_parameters,
                                        density_n_parameters=density_n_parameters,
                                        psk_model=psk_model,
                                        nMaterialTypes=nMaterialTypes,
                                        diagonal_conductivity=diagonal_conductivity)

    	self.nc=1
        self.nd=nd
        self.nPressModel=nPressModel
        #for debugging
        self.qScalarConstant=1.0
        self.capillaryDiffusionScaling=capillaryDiffusionScaling
    def attachModels(self,modelList):
        if self.nPressModel == None:
            print 'Warning TwophaseDarcy_split_saturation_base nPressModel == None returning in attachModels'
            return
        self.flowModel = modelList[self.nPressModel]
        #
        self.q_q_t    = modelList[self.nPressModel].q[('velocity',0)]
        self.ebqe_q_t  = modelList[self.nPressModel].ebqe[('velocity',0)]
        if modelList[self.nPressModel].ebq.has_key(('velocity',0)):
            self.ebq_q_t  = modelList[self.nPressModel].ebq[('velocity',0)]
        #do we really need other model values for q_t in potential calculation?
        assert self.ip_psiw.shape == modelList[self.nPressModel].phi_ip[('u',0)].shape
        self.ip_psiw = modelList[self.nPressModel].phi_ip[('u',0)]
        self.q_psiw    = modelList[self.nPressModel].q[('u',0)]
        self.ebqe_psiw = modelList[self.nPressModel].ebqe[('u',0)]
        if modelList[self.nPressModel].ebq.has_key(('u',0)):
            self.ebq_psiw = modelList[self.nPressModel].ebq[('u',0)]
    def initializeElementQuadrature(self,t,cq):
        TwophaseDarcyFlow_base.initializeElementQuadrature(self,t,cq)
    def initializeElementQuadrature(self,t,cq):
        TwophaseDarcyFlow_base.initializeElementQuadrature(self,t,cq)
    	#set up dummy values in case we're not running the other model
        self.q_q_t   = numpy.zeros(cq[('f',0)].shape,'d')
        self.q_q_t[:] = self.qScalarConstant
        self.q_psiw   = numpy.ones(cq[('u',0)].shape,'d')
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        TwophaseDarcyFlow_base.initializeElementBoundaryQuadrature(self,t,cebq,cebq_global)
    	#set up dummy values in case we're not running the other model
        self.ebq_q_t = numpy.zeros(cebq[('f',0)].shape,'d')
        self.ebq_q_t[:] = self.qScalarConstant
        self.ebq_psiw = numpy.ones(cebq[('u',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        TwophaseDarcyFlow_base.initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe)
    	#set up dummy values in case we're not running the other model
        self.ebqe_q_t = numpy.zeros(cebqe[('f',0)].shape,'d')
        self.ebqe_q_t[:] = self.qScalarConstant
        self.ebqe_psiw = numpy.ones(cebqe[('u',0)].shape,'d')
    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        TwophaseDarcyFlow_base.initializeGeneralizedInterpolationPointQuadrature(self,t,cip)
    	#set up dummy values in case we're not running the other model
        self.ip_q_t = numpy.zeros(cip[('f',0)].shape,'d')
        self.ip_q_t[:] = self.qScalarConstant
        self.ip_psiw = numpy.ones(cip[('u',0)].shape,'d')

class TwophaseDarcy_incompressible_split_pressure(TwophaseDarcy_split_pressure_base):
    """
    Total flow conservation equation in an incompressible fractional flow formulation

Saturation equation
\begin{eqnarray}
\label{eq:2p-ff-mb-w}
\pd{m_w}{t} + \deld\left(\vec f_w - \ten{a}_{w}\grad \phi_w \right) + r_w &=& 0
\end{eqnarray}
and total flow conservation equation
\begin{eqnarray}
\label{eq:2p-ff-mb-m}
\deld\left(\vec f_m - \ten{a}_m \grad \phi_m \right) + r_m &=& 0 
\end{eqnarray}

\begin{table}[ht]
\caption{Coefficient definitions for Two-phase flow, \eqn{2p-ff-mb-w} and \eqn{2p-ff-mb-m}
\label{tab:2p-ff-coef-1}
}
\begin{tabular}{cc}
\hline
Var. & Def. \\
\hline
$u_w$ & $S_w $ \\
$u_n$ & $\psi_w $ \\
$\phi_w$ & $\psi_c$ \\
$\phi_m$ & $\psi_w$  \\
$m_w$ & $\theta_s \rho_{w}S_w$  \\
$\vec f_w$ & $\gvec{\sigma}_t F_w - \ten{K}_s\lambda_wF_n\left(b\rho_n - \rho_w\right)\vec g_u$  \\
$\vec f_m$ &$-\ten{K}_s\lambda_tF_n\grad \psi_c + \ten{K}_s\vec g\lambda_t\left[\rho_w + F_n\left(b\rho_n - \rho_w\right)\right]$\\
$\ten a_{w,w}$ & $-\lambda_wF_n \ten{K}_{s}$ \\
$\ten a_{m,m}$ & $\lambda_t \ten{K}_{s}$ \\
\hline
$F_i $ & $\lambda_i/\lambda_t $, $i=w,n$  \\
$\grad \psi_c $ & $\od{\psi_c}{S_w}\grad S_w $  \\
$\gvec {\sigma}_t$ & $\gvec \sigma_w + \gvec \sigma_n$\\
\hline 
\end{tabular}
\end{table}

     Here S_i is the saturation for each phase, \varrho_{i} is the density and
     \varrho_{i,0} is a reference value

    (normalized) mass flux for each phase is

       \vec \sigma_i = - \ten{a}_{i}\grad \phi_i,    i=w,n
    
    and
       \psi_{i} = p_{i}/|\vec g|\rho_{w,0}
       b        = \varrho_{n,0}/\varrho_{w,0}
       \hat{mu}_{i} = \mu_i/\mu_{w}

    for the pressure of each phase, p_i, and we have the capillary pressure relation
       \psi_{c} = \psi_{n} - \psi_{w}

    The dependent variables are
      S_w, and \psi_w

    Note S_n = 1-S_w

    needs to be implemented
    r_m = r_m(\vec x,t) 

    TODO:

      Figure out if really need to evaluate potential interpolation points
   """
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_incompressible_split_sd_pressure_het_matType
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nSatModel=1,
                 diagonal_conductivity=True,
                 #for debugging
                 swConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_split_pressure_base.__init__(self,
                                                   nd=nd,
                                                   dimensionless_gravity=dimensionless_gravity,
                                                   density_w=density_w,
                                                   density_n=density_n,
                                                   viscosity_w=viscosity_w,
                                                   viscosity_n=viscosity_n,
                                                   psk_model=psk_model,
                                                   nMaterialTypes=nMaterialTypes,
                                                   nSatModel=nSatModel,
                                                   diagonal_conductivity=diagonal_conductivity,
                                                   swConstant=swConstant,
                                                   capillaryDiffusionScaling=capillaryDiffusionScaling)

        variableNames=['psi_w']
        #these are only nonlinear for compressible flow
        mass      = {}
        advection = {0:{0:'linear'}}
        hamiltonian={}
        diffusion = {0:{0: {0:'linear'}}}
        potential = {0:{0: 'u'}}         
        reaction  = {0:{0:'linear'}}

        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}


        TC_base.__init__(self,
                         self.nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
			 variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
    def evaluate(self,t,c):
        if c[('u',0)].shape == self.q_s_w.shape:
            materialTypes = self.materialTypes_q
            s_w = self.q_s_w
            grad_psic = self.q_grad_psic
            c['psi_n']= numpy.copy(self.q_psic)
            c['psi_n'] += c[('u',0)]
        elif c[('u',0)].shape == self.ebqe_s_w.shape:
            materialTypes = self.materialTypes_ebqe
            s_w = self.ebqe_s_w
            grad_psic = self.ebqe_grad_psic
            c['psi_n']= numpy.copy(self.ebqe_psic)
            c['psi_n'] += c[('u',0)]
        elif c[('u',0)].shape == self.ip_s_w.shape:
            c['psi_n']= numpy.copy(self.ip_psic)
            c['psi_n'] += c[('u',0)]
            #mwf hack 04/03/09 skip
            return
            materialTypes = self.materialTypes_ip
            s_w = self.ip_s_w
            grad_psic = self.ip_grad_psic
        else:
            assert c[('u',0)].shape == self.ebq_s_w.shape
            materialTypes = self.materialTypes_ebq
            s_w = self.ebq_s_w
            grad_psic = self.ebq_grad_psic
            c['psi_n']= numpy.copy(self.ebq_psic)
            c['psi_n'] += c[('u',0)]
        assert self.rwork_psk != None
        self.twophaseDarcy_incompressible_split_sd_pressure_het_matType(self.psk_types[self.psk_model],
                                                                        self.sdInfo[(0,0)][0],
                                                                        self.sdInfo[(0,0)][1],
                                                                        materialTypes,
                                                                        self.muw,
                                                                        self.mun,
                                                                        self.omega_types,
                                                                        self.Ksw_types,
                                                                        self.b,
                                                                        self.capillaryDiffusionScaling,
                                                                        self.rwork_psk,     
                                                                        self.rwork_psk_tolerances,
                                                                        self.rwork_density_w,
                                                                        self.rwork_density_n,
                                                                        self.g[:self.nd],
                                                                        s_w,
                                                                        grad_psic,
                                                                        c[('f',0)],
                                                                        c[('a',0,0)])

        #mwf debug
        if (numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any()):
            import pdb
            pdb.set_trace()
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if c[('u',0)].shape == self.q_s_w.shape:
            c['Ks']=c[('a',0,0)][:,:,0]

#
class TwophaseDarcy_incompressible_split_saturation(TwophaseDarcy_split_saturation_base):
    """
    Aqueous phase mass conservation equation (saturation equation) in
    an incompressible fractional flow formulation

Saturation equation
\begin{eqnarray}
\label{eq:2p-ff-mb-w}
\pd{m_w}{t} + \deld\left(\vec f_w - \ten{a}_{w}\grad \phi_w \right) + r_w &=& 0
\end{eqnarray}
and total flow conservation equation
\begin{eqnarray}
\label{eq:2p-ff-mb-m}
\deld\left(\vec f_m - \ten{a}_m \grad \phi_m \right) + r_m &=& 0 
\end{eqnarray}

\begin{table}[ht]
\caption{Coefficient definitions for Two-phase flow, \eqn{2p-ff-mb-w} and \eqn{2p-ff-mb-m}
\label{tab:2p-ff-coef-1}
}
\begin{tabular}{cc}
\hline
Var. & Def. \\
\hline
$u_w$ & $S_w $ \\
$u_n$ & $\psi_w $ \\
$\phi_w$ & $\psi_c$ \\
$\phi_m$ & $\psi_w$  \\
$m_w$ & $\theta_s \rho_{w}S_w$  \\
$\vec f_w$ & $\gvec{\sigma}_t F_w - \ten{K}_s\lambda_wF_n\left(b\rho_n - \rho_w\right)\vec g_u$  \\
$\vec f_m$ &$-\ten{K}_s\lambda_tF_n\grad \psi_c + \ten{K}_s\vec g\lambda_t\left[\rho_w + F_n\left(b\rho_n - \rho_w\right)\right]$\\
$\ten a_{w,w}$ & $-\lambda_wF_n \ten{K}_{s}$ \\
$\ten a_{m,m}$ & $\lambda_t \ten{K}_{s}$ \\
\hline
$F_i $ & $\lambda_i/\lambda_t $, $i=w,n$  \\
$\grad \psi_c $ & $\od{\psi_c}{S_w}\grad S_w $  \\
$\gvec {\sigma}_t$ & $\gvec \sigma_w + \gvec \sigma_n$\\
\hline 
\end{tabular}
\end{table}

     Here S_i is the saturation for each phase, \varrho_{i} is the density and
     \varrho_{i,0} is a reference value

    (normalized) mass flux for each phase is

       \vec \sigma_i = - \ten{a}_{i}\grad \phi_i,    i=w,n
    
    and
       \psi_{i} = p_{i}/|\vec g|\rho_{w,0}
       b        = \varrho_{n,0}/\varrho_{w,0}
       \hat{mu}_{i} = \mu_i/\mu_{w}

    for the pressure of each phase, p_i, and we have the capillary pressure relation
       \psi_{c} = \psi_{n} - \psi_{w}

    The dependent variables are
      S_w, and \psi_w

    Note S_n = 1-S_w

    needs to be implemented
    r_m = r_m(\vec x,t) 

    TODO:

      Figure out if really need to evaluate potential interpolation points
   """
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_incompressible_split_sd_saturation_het_matType
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nPressModel=0,
                 diagonal_conductivity=True,
                 #for debugging
                 qScalarConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_split_saturation_base.__init__(self,
                                                     nd=nd,
                                                     dimensionless_gravity=dimensionless_gravity,
                                                     density_w=density_w,
                                                     density_n=density_n,
                                                     viscosity_w=viscosity_w,
                                                     viscosity_n=viscosity_n,
                                                     psk_model=psk_model,
                                                     nMaterialTypes=nMaterialTypes,
                                                     nPressModel=nPressModel,
                                                     diagonal_conductivity=diagonal_conductivity,
                                                     qScalarConstant=qScalarConstant,
                                                     capillaryDiffusionScaling=capillaryDiffusionScaling)

        variableNames=['s_w']
        mass      = {0:{0:'linear'}}
        advection = {0:{0:'nonlinear'}}
        hamiltonian={}
        diffusion = {0:{0:{0:'nonlinear'}}}
        potential = {0:{0: 'nonlinear'}} 
        reaction  = {0:{0:'linear'}}

        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}


        TC_base.__init__(self,
                         self.nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
			 variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
    def evaluate(self,t,c):
        if c[('f',0)].shape == self.q_q_t.shape:
            materialTypes = self.materialTypes_q
            q_t = self.q_q_t
            psiw = self.q_psiw
        elif c[('f',0)].shape == self.ebqe_q_t.shape:
            materialTypes = self.materialTypes_ebqe
            q_t = self.ebqe_q_t
            psiw = self.ebqe_psiw
        elif c[('f',0)].shape == self.ip_q_t.shape:
            materialTypes = self.materialTypes_ip
            q_t = self.ip_q_t
            psiw = self.ip_psiw
        else:
            assert c[('f',0)].shape == self.ebq_q_t.shape
            materialTypes = self.materialTypes_ebq
            q_t = self.ebq_q_t
            psiw = self.ebq_psiw
        assert self.rwork_psk != None
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.twophaseDarcy_incompressible_split_sd_saturation_het_matType(self.psk_types[self.psk_model],
                                                                          self.sdInfo[(0,0)][0],
                                                                          self.sdInfo[(0,0)][1],
                                                                          materialTypes,
                                                                          self.muw,
                                                                          self.mun,
                                                                          self.omega_types,
                                                                          self.Ksw_types,
                                                                          self.b,
                                                                          self.capillaryDiffusionScaling,
                                                                          self.rwork_psk,
                                                                          self.rwork_psk_tolerances,
                                                                          self.rwork_density_w,
                                                                          self.rwork_density_n,
                                                                          self.g[:self.nd],
                                                                          q_t,
                                                                          c[('u',0)],
                                                                          c[('m',0)],
                                                                          c[('dm',0,0)],
                                                                          c[('phi',0)],
                                                                          c[('dphi',0,0)],
                                                                          c[('f',0)],
                                                                          c[('df',0,0)],
                                                                          c[('a',0,0)],
                                                                          c[('da',0,0,0)])
        #material types fi
        #mwf debug
        if (numpy.isnan(c[('da',0,0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('phi',0)]).any() or
            numpy.isnan(c[('dphi',0,0)]).any() or
            numpy.isnan(c[('df',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any() or
            numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any()):
            import pdb
            pdb.set_trace()


class TwophaseDarcy_compressible_split_pressure(TwophaseDarcy_split_pressure_base):
    """
    Total flow conservation equation in a split compressible fractional flow formulation
    Right now, the options are 
 
        compressibility for the non-wetting phase (compressibleN) : compressiblityFlag=1
        compressibility for both phases but with the slight       
          compressiblity assumption (spatial density gradients are negligible) compressiblityFlag=2
             
Saturation equation
\begin{eqnarray}
\label{eq:2p-c-ff-mb-w}
\pd{m_w}{t} + \deld\left(\vec f_w - \ten{a}_{w}\grad \phi_w \right) + r_w &=& 0
\end{eqnarray}
and total flow conservation equation
\begin{eqnarray}
\label{eq:2p-c-ff-mb-m}
\pd{m_m}{t} + \deld\left(\vec f_m - \ten{a}_m \grad \phi_m \right) + r_m &=& 0 
\end{eqnarray}

\begin{table}[ht]
\caption{Coefficient definitions for compressible two-phase fractional flow formulation flow, \eqn{2p-c-ff-mb-w} and \eqn{2p-c-ff-mb-m}
\label{tab:2p-c-ff-coef-1}
}
\begin{tabular}{cc}
\hline
Var. & Def. \\
\hline
$u_w$ & $S_w $ \\
$u_n$ & $\psi_w $ \\
$\phi_w$ & $\psi_c$ \\
$\phi_m$ & $\psi_w$  \\
$m_w$ & $\theta_s \rho_{w}S_w$  \\
$m_m$ & $\theta_s \rho_{w}S_w + \theta_s\rho_{n}(1-S_w)$  \\
$\vec f_w^{\dagger}$ & $\gvec{\sigma}_t F_w - \ten{K}_s\lambda_wF_n\left(b\rho_n - \rho_w\right)\vec g_u$  \\
$\vec f_m^{\dagger}$ &$-\ten{K}_s\lambda_tF_n\grad \psi_c + \ten{K}_s\vec g\lambda_t\left[\rho_w + F_n\left(b\rho_n - \rho_w\right)\right]$\\
$\ten a_{w,w}^{\dagger}$ & $-\lambda_wF_n \ten{K}_{s}$ \\
$\ten a_{m,m}^{\dagger}$ & $\lambda_t \ten{K}_{s}$ \\
\hline
$F_i $ & $\lambda_i/\lambda_t $, $i=w,n$  \\
$\grad \psi_c $ & $\od{\psi_c}{S_w}\grad S_w $  \\
$\gvec {\sigma}_t$ & $\gvec \sigma_w + \gvec \sigma_n$\\
\hline 
\multicolumn{2}{l}
{$\dagger$ : $\rho_i \equiv 1$, $\lambda_i= k_{r,i}/\hat{\mu}_{i}$ for slight compressibility assumption}
\end{tabular}
\end{table}
    


     Here S_i is the saturation for each phase, \varrho_{i} is the density and
     \varrho_{i,0} is a reference value

    (normalized) mass flux for each phase is

       \vec \sigma_i = - \ten{a}_{i}\grad \phi_i,    i=w,n
    
    and
       \psi_{i} = p_{i}/|\vec g|\rho_{w,0}
       b        = \varrho_{n,0}/\varrho_{w,0}
       \hat{mu}_{i} = \mu_i/\mu_{w}
       \lambda_i= \rho k_{r,i}/\hat{\mu}_{i}
       \lambda_t= \lambda_n + \lambda_w

    for the pressure of each phase, p_i, and we have the capillary pressure relation
       \psi_{c} = \psi_{n} - \psi_{w}

    The dependent variables are
      S_w, and \psi_w

    Note S_n = 1-S_w

    needs to be implemented
    r_m = r_m(\vec x,t) 

    TODO:

      Figure out if really need to evaluate potential interpolation points
   """
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_slightCompressible_split_sd_pressure_het_matType
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_compressibleN_split_sd_pressure_het_matType
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=TwophaseDarcyFlow_base.default_density_w_parameters,
                 density_n_parameters=TwophaseDarcyFlow_base.default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nSatModel=1,
                 compressibilityFlag=2,
                 diagonal_conductivity=True,
                 #for debugging
                 swConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_split_pressure_base.__init__(self,
                                                   nd=nd,
                                                   dimensionless_gravity=dimensionless_gravity,
                                                   density_w=density_w,
                                                   density_n=density_n,
                                                   viscosity_w=viscosity_w,
                                                   viscosity_n=viscosity_n,
                                                   density_w_parameters=density_w_parameters,
                                                   density_n_parameters=density_n_parameters,
                                                   psk_model=psk_model,
                                                   nMaterialTypes=nMaterialTypes,
                                                   nSatModel=nSatModel,
                                                   diagonal_conductivity=diagonal_conductivity,
                                                   swConstant=swConstant,
                                                   capillaryDiffusionScaling=capillaryDiffusionScaling)

        self.compressibilityFlag=compressibilityFlag
        variableNames=['psi_w']
        #
        mass      = {0:{0:'nonlinear'}}
        advection = {0:{0:'nonlinear'}}
        hamiltonian={}
        diffusion = {0:{0: {0:'nonlinear'}}}
        potential = {0:{0: 'u'}}         
        reaction  = {0:{0:'linear'}}
        if self.compressibilityFlag == 2:
            advection = {0:{0:'linear'}}
            diffusion = {0:{0: {0:'linear'}}}

        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}


        TC_base.__init__(self,
                         self.nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
			 variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
    def evaluate(self,t,c):
        if c[('u',0)].shape == self.q_s_w.shape:
            materialTypes = self.materialTypes_q
            s_w = self.q_s_w
            grad_psic = self.q_grad_psic
            c['psi_n']= numpy.copy(self.q_psic)
            c['psi_n'] += c[('u',0)]
        elif c[('u',0)].shape == self.ebqe_s_w.shape:
            materialTypes = self.materialTypes_ebqe
            s_w = self.ebqe_s_w
            grad_psic = self.ebqe_grad_psic
            c['psi_n']= numpy.copy(self.ebqe_psic)
            c['psi_n'] += c[('u',0)]
        elif c[('u',0)].shape == self.ip_s_w.shape:
            c['psi_n']= numpy.copy(self.ip_psic)
            c['psi_n'] += c[('u',0)]
            #mwf hack 04/03/09 skip
            return
            materialTypes = self.materialTypes_ip
            s_w = self.ip_s_w
            grad_psic = self.ip_grad_psic
        else:
            assert c[('u',0)].shape == self.ebq_s_w.shape
            materialTypes = self.materialTypes_ebq
            s_w = self.ebq_s_w
            grad_psic = self.ebq_grad_psic
            c['psi_n']= numpy.copy(self.ebq_psic)
            c['psi_n'] += c[('u',0)]
        assert self.rwork_psk != None
        if self.compressibilityFlag == 2:
            self.twophaseDarcy_slightCompressible_split_sd_pressure_het_matType(self.psk_types[self.psk_model],
                                                                                  self.density_types[self.density_w_model],
                                                                                  self.density_types[self.density_n_model],
                                                                                  self.sdInfo[(0,0)][0],
                                                                                  self.sdInfo[(0,0)][1],
                                                                                  materialTypes,
                                                                                  self.muw,
                                                                                  self.mun,
                                                                                  self.omega_types,
                                                                                  self.Ksw_types,
                                                                                  self.b,
                                                                                  self.capillaryDiffusionScaling,
                                                                                  self.rwork_psk,     
                                                                                  self.rwork_psk_tolerances,
                                                                                  self.rwork_density_w,
                                                                                  self.rwork_density_n,
                                                                                  self.g[:self.nd],
                                                                                  s_w,
                                                                                  c[('u',0)],
                                                                                  c['psi_n'],
                                                                                  grad_psic,
                                                                                  c[('m',0)],
                                                                                  c[('dm',0,0)],
                                                                                  c[('f',0)],
                                                                                  c[('a',0,0)])

        else:
            self.twophaseDarcy_compressibleN_split_sd_pressure_het_matType(self.psk_types[self.psk_model],
                                                                           self.density_types[self.density_w_model],
                                                                           self.density_types[self.density_n_model],
                                                                           self.sdInfo[(0,0)][0],
                                                                           self.sdInfo[(0,0)][1],
                                                                           materialTypes,
                                                                           self.muw,
                                                                           self.mun,
                                                                           self.omega_types,
                                                                           self.Ksw_types,
                                                                           self.b,
                                                                           self.capillaryDiffusionScaling,
                                                                           self.rwork_psk,     
                                                                           self.rwork_psk_tolerances,
                                                                           self.rwork_density_w,
                                                                           self.rwork_density_n,
                                                                           self.g[:self.nd],
                                                                           s_w,
                                                                           c[('u',0)],
                                                                           c['psi_n'],
                                                                           grad_psic,
                                                                           c[('m',0)],
                                                                           c[('dm',0,0)],
                                                                           c[('f',0)],
                                                                           c[('a',0,0)])

        #mwf debug
        if (numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any()):
            import pdb
            pdb.set_trace()
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if c[('u',0)].shape == self.q_s_w.shape:
            c['Ks']=c[('a',0,0)][:,:,0]


class TwophaseDarcy_compressible_split_saturation(TwophaseDarcy_split_saturation_base):
    """
    Aqueous phase mass conservation equation (saturation equation) in
    a  compressible fractional flow formulation

    Right now, the options are 
 
        compressibility for the non-wetting phase (compressibleN) : compressiblityFlag=1
        compressibility for both phases but with the slight       
          compressiblity assumption (spatial density gradients are negligible) compressiblityFlag=2
Saturation equation
\begin{eqnarray}
\label{eq:2p-c-ff-mb-w}
\pd{m_w}{t} + \deld\left(\vec f_w - \ten{a}_{w}\grad \phi_w \right) + r_w &=& 0
\end{eqnarray}
and total flow conservation equation
\begin{eqnarray}
\label{eq:2p-c-ff-mb-m}
\pd{m_m}{t} + \deld\left(\vec f_m - \ten{a}_m \grad \phi_m \right) + r_m &=& 0 
\end{eqnarray}

\begin{table}[ht]
\caption{Coefficient definitions for compressible two-phase fractional flow formulation flow, \eqn{2p-c-ff-mb-w} and \eqn{2p-c-ff-mb-m}
\label{tab:2p-c-ff-coef-1}
}
\begin{tabular}{cc}
\hline
Var. & Def. \\
\hline
$u_w$ & $S_w $ \\
$u_n$ & $\psi_w $ \\
$\phi_w$ & $\psi_c$ \\
$\phi_m$ & $\psi_w$  \\
$m_w$ & $\theta_s \rho_{w}S_w$  \\
$m_m$ & $\theta_s \rho_{w}S_w + \theta_s\rho_{n}(1-S_w)$  \\
$\vec f_w^{\dagger}$ & $\gvec{\sigma}_t F_w - \ten{K}_s\lambda_wF_n\left(b\rho_n - \rho_w\right)\vec g_u$  \\
$\vec f_m^{\dagger}$ &$-\ten{K}_s\lambda_tF_n\grad \psi_c + \ten{K}_s\vec g\lambda_t\left[\rho_w + F_n\left(b\rho_n - \rho_w\right)\right]$\\
$\ten a_{w,w}^{\dagger}$ & $-\lambda_wF_n \ten{K}_{s}$ \\
$\ten a_{m,m}^{\dagger}$ & $\lambda_t \ten{K}_{s}$ \\
\hline
$F_i $ & $\lambda_i/\lambda_t $, $i=w,n$  \\
$\grad \psi_c $ & $\od{\psi_c}{S_w}\grad S_w $  \\
$\gvec {\sigma}_t$ & $\gvec \sigma_w + \gvec \sigma_n$\\
\hline 
\multicolumn{2}{l}
{$\dagger$ : $\rho_i \equiv 1$, $\lambda_i= k_{r,i}/\hat{\mu}_{i}$ for slight compressibility assumption}
\end{tabular}
\end{table}
    


     Here S_i is the saturation for each phase, \varrho_{i} is the density and
     \varrho_{i,0} is a reference value

    (normalized) mass flux for each phase is

       \vec \sigma_i = - \ten{a}_{i}\grad \phi_i,    i=w,n
    
    and
       \psi_{i} = p_{i}/|\vec g|\rho_{w,0}
       b        = \varrho_{n,0}/\varrho_{w,0}
       \hat{mu}_{i} = \mu_i/\mu_{w}
       \lambda_i= \rho k_{r,i}/\hat{\mu}_{i}
       \lambda_t= \lambda_n + \lambda_w

    for the pressure of each phase, p_i, and we have the capillary pressure relation
       \psi_{c} = \psi_{n} - \psi_{w}

    The dependent variables are
      S_w, and \psi_w

    Note S_n = 1-S_w

    needs to be implemented
    r_m = r_m(\vec x,t) 

    TODO:

      Figure out if really need to evaluate potential interpolation points
    """
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_slightCompressible_split_sd_saturation_het_matType
    from proteus.cTwophaseDarcyCoefficients import twophaseDarcy_compressibleN_split_sd_saturation_het_matType
    def __init__(self,
                 nd=1,
                 dimensionless_gravity=[-1.0],
                 density_w=998.2, #Water kg/m^3
                 density_n=1.205, #Air   kg/m^3
                 viscosity_w=8.9e-4, #Water kg/m s
                 viscosity_n=1.81e-5,#Air kg/ m s
                 density_w_parameters=TwophaseDarcyFlow_base.default_density_w_parameters,
                 density_n_parameters=TwophaseDarcyFlow_base.default_density_n_parameters,
                 psk_model='VGM',
                 nMaterialTypes=1,
                 nPressModel=0,
                 compressibilityFlag=2,
                 diagonal_conductivity=True,
                 #for debugging
                 qScalarConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_split_saturation_base.__init__(self,
                                                     nd=nd,
                                                     dimensionless_gravity=dimensionless_gravity,
                                                     density_w=density_w,
                                                     density_n=density_n,
                                                     viscosity_w=viscosity_w,
                                                     viscosity_n=viscosity_n,
                                                     density_w_parameters=density_w_parameters,
                                                     density_n_parameters=density_n_parameters,
                                                     psk_model=psk_model,
                                                     nMaterialTypes=nMaterialTypes,
                                                     nPressModel=nPressModel,
                                                     diagonal_conductivity=diagonal_conductivity,
                                                     qScalarConstant=qScalarConstant,
                                                     capillaryDiffusionScaling=capillaryDiffusionScaling)

        self.compressibilityFlag=compressibilityFlag
        variableNames=['s_w']
        mass      = {0:{0:'linear'}}
        advection = {0:{0:'nonlinear'}}
        hamiltonian={}
        diffusion = {0:{0:{0:'nonlinear'}}}
        potential = {0:{0: 'nonlinear'}} 
        reaction  = {0:{0:'linear'}}

        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}


        TC_base.__init__(self,
                         self.nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
			 variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
    def evaluate(self,t,c):
        if c[('f',0)].shape == self.q_q_t.shape:
            materialTypes = self.materialTypes_q
            q_t = self.q_q_t
            psiw = self.q_psiw
        elif c[('f',0)].shape == self.ebqe_q_t.shape:
            materialTypes = self.materialTypes_ebqe
            q_t = self.ebqe_q_t
            psiw = self.ebqe_psiw
        elif c[('f',0)].shape == self.ip_q_t.shape:
            materialTypes = self.materialTypes_ip
            q_t = self.ip_q_t
            psiw = self.ip_psiw
        else:
            assert c[('f',0)].shape == self.ebq_q_t.shape
            materialTypes = self.materialTypes_ebq
            q_t = self.ebq_q_t
            psiw = self.ebq_psiw
        assert self.rwork_psk != None
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if self.compressibilityFlag == 2:
            self.twophaseDarcy_slightCompressible_split_sd_saturation_het_matType(self.psk_types[self.psk_model],
                                                                                  self.density_types[self.density_w_model],
                                                                                  self.density_types[self.density_n_model],
                                                                                  self.sdInfo[(0,0)][0],
                                                                                  self.sdInfo[(0,0)][1],
                                                                                  materialTypes,
                                                                                  self.muw,
                                                                                  self.mun,
                                                                                  self.omega_types,
                                                                                  self.Ksw_types,
                                                                                  self.b,
                                                                                  self.capillaryDiffusionScaling,
                                                                                  self.rwork_psk,
                                                                                  self.rwork_psk_tolerances,
                                                                                  self.rwork_density_w,
                                                                                  self.rwork_density_n,
                                                                                  self.g[:self.nd],
                                                                                  q_t,
                                                                                  psiw,
                                                                                  c[('u',0)],
                                                                                  c[('m',0)],
                                                                                  c[('dm',0,0)],
                                                                                  c[('phi',0)],
                                                                                  c[('dphi',0,0)],
                                                                                  c[('f',0)],
                                                                                  c[('df',0,0)],
                                                                                  c[('a',0,0)],
                                                                                  c[('da',0,0,0)])
        else:
            self.twophaseDarcy_compressibleN_split_sd_saturation_het_matType(self.psk_types[self.psk_model],
                                                                             self.density_types[self.density_w_model],
                                                                             self.density_types[self.density_n_model],
                                                                             self.sdInfo[(0,0)][0],
                                                                             self.sdInfo[(0,0)][1],
                                                                             materialTypes,
                                                                             self.muw,
                                                                             self.mun,
                                                                             self.omega_types,
                                                                             self.Ksw_types,
                                                                             self.b,
                                                                             self.capillaryDiffusionScaling,
                                                                             self.rwork_psk,
                                                                             self.rwork_psk_tolerances,
                                                                             self.rwork_density_w,
                                                                             self.rwork_density_n,
                                                                             self.g[:self.nd],
                                                                             q_t,
                                                                             psiw,
                                                                             c[('u',0)],
                                                                             c[('m',0)],
                                                                             c[('dm',0,0)],
                                                                             c[('phi',0)],
                                                                             c[('dphi',0,0)],
                                                                             c[('f',0)],
                                                                             c[('df',0,0)],
                                                                             c[('a',0,0)],
                                                                             c[('da',0,0,0)])
            
        #material types fi
        #mwf debug
        if (numpy.isnan(c[('da',0,0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('phi',0)]).any() or
            numpy.isnan(c[('dphi',0,0)]).any() or
            numpy.isnan(c[('df',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any() or
            numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any()):
            import pdb
            pdb.set_trace()


class IncompressibleFractionalFlowPressureMualemVanGenuchten(TwophaseDarcy_incompressible_split_pressure):
    """
    Total flow equation coefficients for incompressible flow assuming Mualem-Van Genuchten psk's
    """
    def __init__(self,
                 nd,
                 Ksw_types, 
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 dimensionless_gravity,
                 density_w,
                 density_n,
                 viscosity_w,
                 viscosity_n,
                 nSatModel=1,
                 diagonal_conductivity=True,
                 vgm_small_eps=1.0e-16,
                 vgm_ns_del=1.0e-8,
                 #for debugging
                 swConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_incompressible_split_pressure.__init__(self,
                                                             nd=nd,
                                                             dimensionless_gravity=dimensionless_gravity,
                                                             density_w=density_w,
                                                             density_n=density_n,
                                                             viscosity_w=viscosity_w,
                                                             viscosity_n=viscosity_n,
                                                             psk_model='VGM',
                                                             nMaterialTypes=len(thetaR_types),
                                                             nSatModel=nSatModel,
                                                             diagonal_conductivity=diagonal_conductivity,
                                                             swConstant=swConstant,
                                                             capillaryDiffusionScaling=capillaryDiffusionScaling)

        for input in [vgm_n_types,vgm_alpha_types,thetaR_types,thetaSR_types]:
            assert len(input)==self.nMaterialTypes
        
        vgm_m_types    = 1.0-1.0/vgm_n_types
        thetaS_types   = thetaSR_types + thetaR_types
        Sw_max_types   = numpy.ones((self.nMaterialTypes,),'d')
        Sw_min_types   = thetaR_types/thetaS_types
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.psk_tolerances['VGM']['eps_small']=vgm_small_eps
        self.psk_tolerances['VGM']['ns_del']   =vgm_ns_del
        self.setMaterialTypes(Ksw_types=Ksw_types,
                              omega_types=thetaS_types,
                              Sw_max_types=Sw_max_types,
                              Sw_min_types=Sw_min_types,
                              vg_alpha_types=vgm_alpha_types,
                              vg_m_types=vgm_m_types)

#
class IncompressibleFractionalFlowSaturationMualemVanGenuchten(TwophaseDarcy_incompressible_split_saturation):
    """
    Saturation equation coefficients for incompressible flow assuming Mualem-Van Genuchten psk's
    """
    def __init__(self,
                 nd,
                 Ksw_types, 
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 dimensionless_gravity,
                 density_w,
                 density_n,
                 viscosity_w,
                 viscosity_n,
                 nPressModel=1,
                 diagonal_conductivity=True,
                 vgm_small_eps=1.0e-16,
                 vgm_ns_del=1.0e-8,
                 #for debugging
                 qScalarConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_incompressible_split_saturation.__init__(self,
                                                               nd=nd,
                                                               dimensionless_gravity=dimensionless_gravity,
                                                               density_w=density_w,
                                                               density_n=density_n,
                                                               viscosity_w=viscosity_w,
                                                               viscosity_n=viscosity_n,
                                                               psk_model='VGM',
                                                               nMaterialTypes=len(thetaR_types),
                                                               nPressModel=nPressModel,
                                                               diagonal_conductivity=diagonal_conductivity,
                                                               qScalarConstant=qScalarConstant,
                                                               capillaryDiffusionScaling=capillaryDiffusionScaling)

        for input in [vgm_n_types,vgm_alpha_types,thetaR_types,thetaSR_types]:
            assert len(input)==self.nMaterialTypes
        
        vgm_m_types    = 1.0-1.0/vgm_n_types
        thetaS_types   = thetaSR_types + thetaR_types
        Sw_max_types   = numpy.ones((self.nMaterialTypes,),'d')
        Sw_min_types   = thetaR_types/thetaS_types
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.psk_tolerances['VGM']['eps_small']=vgm_small_eps
        self.psk_tolerances['VGM']['ns_del']   =vgm_ns_del
        self.setMaterialTypes(Ksw_types=Ksw_types,
                              omega_types=thetaS_types,
                              Sw_max_types=Sw_max_types,
                              Sw_min_types=Sw_min_types,
                              vg_alpha_types=vgm_alpha_types,
                              vg_m_types=vgm_m_types)

#
class CompressibleFractionalFlowPressureMualemVanGenuchten(TwophaseDarcy_compressible_split_pressure):
    """
    Total flow equation coefficients for slight compressible flow assuming Mualem-Van Genuchten psk's
    """
    def __init__(self,
                 nd,
                 Ksw_types, 
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 dimensionless_gravity,
                 density_w,
                 density_n,
                 viscosity_w,
                 viscosity_n,
                 density_w_parameters,
                 density_n_parameters,
                 nSatModel=1,
                 compressibilityFlag=2,#2 slight-compressiblity for N,W, 1 -- compressiblity in N only
                 diagonal_conductivity=True,
                 vgm_small_eps=1.0e-16,
                 vgm_ns_del=1.0e-8,
                 #for debugging
                 swConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_compressible_split_pressure.__init__(self,
                                                           nd=nd,
                                                           dimensionless_gravity=dimensionless_gravity,
                                                           density_w=density_w,
                                                           density_n=density_n,
                                                           viscosity_w=viscosity_w,
                                                           viscosity_n=viscosity_n,
                                                           density_w_parameters=density_w_parameters,
                                                           density_n_parameters=density_n_parameters,
                                                           psk_model='VGM',
                                                           nMaterialTypes=len(thetaR_types),
                                                           nSatModel=nSatModel,
                                                           compressibilityFlag=compressibilityFlag,
                                                           diagonal_conductivity=diagonal_conductivity,
                                                           swConstant=swConstant,
                                                           capillaryDiffusionScaling=capillaryDiffusionScaling)

        for input in [vgm_n_types,vgm_alpha_types,thetaR_types,thetaSR_types]:
            assert len(input)==self.nMaterialTypes
        
        vgm_m_types    = 1.0-1.0/vgm_n_types
        thetaS_types   = thetaSR_types + thetaR_types
        Sw_max_types   = numpy.ones((self.nMaterialTypes,),'d')
        Sw_min_types   = thetaR_types/thetaS_types
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.psk_tolerances['VGM']['eps_small']=vgm_small_eps
        self.psk_tolerances['VGM']['ns_del']   =vgm_ns_del
        self.setMaterialTypes(Ksw_types=Ksw_types,
                              omega_types=thetaS_types,
                              Sw_max_types=Sw_max_types,
                              Sw_min_types=Sw_min_types,
                              vg_alpha_types=vgm_alpha_types,
                              vg_m_types=vgm_m_types)


#
class CompressibleFractionalFlowSaturationMualemVanGenuchten(TwophaseDarcy_compressible_split_saturation):
    """
    Saturation equation coefficients for slightly compressible flow assuming Mualem-Van Genuchten psk's
    """
    def __init__(self,
                 nd,
                 Ksw_types, 
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 dimensionless_gravity,
                 density_w,
                 density_n,
                 viscosity_w,
                 viscosity_n,
                 density_w_parameters,
                 density_n_parameters,
                 nPressModel=1,
                 compressibilityFlag=2,#2 slight-compressiblity for N,W, 1 -- compressiblity in N only
                 diagonal_conductivity=True,
                 vgm_small_eps=1.0e-16,
                 vgm_ns_del=1.0e-8,
                 #for debugging
                 qScalarConstant=1.0,
                 capillaryDiffusionScaling=1.0):
        TwophaseDarcy_compressible_split_saturation.__init__(self,
                                                             nd=nd,
                                                             dimensionless_gravity=dimensionless_gravity,
                                                             density_w=density_w,
                                                             density_n=density_n,
                                                             viscosity_w=viscosity_w,
                                                             viscosity_n=viscosity_n,
                                                             density_w_parameters=density_w_parameters,
                                                             density_n_parameters=density_n_parameters,
                                                             psk_model='VGM',
                                                             nMaterialTypes=len(thetaR_types),
                                                             nPressModel=nPressModel,
                                                             compressibilityFlag=compressibilityFlag,
                                                             diagonal_conductivity=diagonal_conductivity,
                                                             qScalarConstant=qScalarConstant,
                                                             capillaryDiffusionScaling=capillaryDiffusionScaling)

        for input in [vgm_n_types,vgm_alpha_types,thetaR_types,thetaSR_types]:
            assert len(input)==self.nMaterialTypes
        
        vgm_m_types    = 1.0-1.0/vgm_n_types
        thetaS_types   = thetaSR_types + thetaR_types
        Sw_max_types   = numpy.ones((self.nMaterialTypes,),'d')
        Sw_min_types   = thetaR_types/thetaS_types
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.psk_tolerances['VGM']['eps_small']=vgm_small_eps
        self.psk_tolerances['VGM']['ns_del']   =vgm_ns_del
        self.setMaterialTypes(Ksw_types=Ksw_types,
                              omega_types=thetaS_types,
                              Sw_max_types=Sw_max_types,
                              Sw_min_types=Sw_min_types,
                              vg_alpha_types=vgm_alpha_types,
                              vg_m_types=vgm_m_types)

#
########################################
#Single-phase species transport
########################################

class GroundwaterTransportCoefficientsELLAM(TC_base):
    from proteus.ctransportCoefficients import groundwaterTransportCoefficientsEvaluate
    """
    groundwater advection-dispersion equation with constant coefficients but variable 
    velocity 
    can have a spatially varying 'reaction' term r(x,t) for testing

    differentiates between Darcy velocity and 'adjoint' and 'characteristic' velocities
    anything labeled just velocity corresponds to Darcy velocity

    computes the characteristic and adjoint velocities as 
       q_e = q_e/\bar{\omega}_e   on element \Omega_e where \bar{\omega_e} is the average 
    TODO
        figure out how to get reference point quadrature values here

        allow differerent velocity reps for different components?
 
        move python loops in velocity setup to cython
    """
    def __init__(self,nc=1,
                 omega=0.3,
                 alpha_L=1.0,
                 alpha_T=0.2,
                 d=1.3e-9,
                 nd = 2,
                 velocityFunctions = None,
                 velocitySpaceFlag = 'c0p1',
                 reactionTerms = None,
                 meModelId = 0,
                 flowModelId = None,
                 flowIsTransient=False,
                 flowVelocityComponentMap=None):
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            advection[i] = {i : {i:'linear'}}
            mass[i] = {i : {i:'linear'}}
            reaction[i] = {i : {i:'constant'}}
            potential[i] = {i : 'u'}
        #end i
        sparseDiffusionTensors = {}
        for ci in range(nc):
            sparseDiffusionTensors[(ci,ci)]=(numpy.arange(start=0,stop=nd**2+1,step=nd,dtype='i'),
                                             numpy.array([range(nd) for row in range(nd)],dtype='i'))
        names = ['C_%s' % ci for ci in range(nc)]
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         names,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
        self.omega = omega
        self.d = d
        self.alpha_L = alpha_L
        self.alpha_T = alpha_T
        self.useC = True
        self.velocityFunctions = velocityFunctions
        self.velocitySpaceFlag = velocitySpaceFlag
        self.reactionTerms     = reactionTerms
        self.nd = nd
        self.needToEvaluateVelocity = True
        self.flowIsTransient = flowIsTransient
        self.velocitySpace = {}; self.velocity = {}; self.adjoint_velocity = {};
        self.characteristic_velocitySpace = {}; self.adjoint_velocitySpace = {}
        self.adjoint_velocity_times_last = {};  self.adjoint_velocity_times={}
        self.adjoint_velocity_dofs_last = {};   self.adjoint_velocity_dofs = {}
        self.adjoint_velocity_l2g = {}; self.adjoint_velocity_interpolation_values = {}
        self.velocity_dofs = {}; self.velocity_l2g = {}; self.velocity_interpolation_values = {}
        for ci in range(self.nc):
            self.velocitySpace[ci] = None
            self.characteristic_velocitySpace[ci] = None; self.adjoint_velocitySpace[ci] = None
            self.adjoint_velocity_times_last[ci] = None;  self.adjoint_velocity_times[ci]=None
            self.adjoint_velocity_dofs_last[ci] = None;   self.adjoint_velocity_dofs[ci] = None
            self.adjoint_velocity_l2g[ci] = None ; self.adjoint_velocity_interpolation_values[ci]=None
        self.q={}; self.ebqe={}; self.ebq ={}; self.ebq_global = {}
        self.need_ebq_arrays = False
        self.flowModelId = flowModelId
        self.flowModel   = None
        self.meModelId   = meModelId
        self.flowVelocityComponentMap = flowVelocityComponentMap #map transport components to flow velocity components
        if flowModelId != None and self.flowVelocityComponentMap == None:
            #by default all point to first component of flow
            self.flowVelocityComponentMap = {}
            for ci in range(self.nc):
                self.flowVelocityComponentMap[ci] = 0 
    def initializeMesh(self,mesh):
        self.mesh = mesh
        #adjoint and characteristic velocities are the same for a linear problem
        self.characteristic_velocity_dofs = {}; 
        self.characteristic_velocity_l2g  = {}; 
        #TODO This is not going to be consistent with BDM representation
        if self.velocitySpaceFlag == 'c0p1':
            for ci in range(self.nc):
                self.velocitySpace[ci]  = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh,self.nd)
                self.velocity[ci]       = FemTools.FiniteElementFunction(self.velocitySpace[ci],dim_dof=self.nd,isVector=True,name="velocity_%s" % ci)
                self.velocity_interpolation_values[ci] = numpy.zeros((mesh.nElements_global,
                                                                     self.velocitySpace[ci].referenceFiniteElement.interpolationConditions.nQuadraturePoints,
                                                                     self.nd),'d')
                self.velocity_dofs[ci] = self.velocity[ci].dof
                self.velocity_l2g[ci]  = self.velocitySpace[ci].dofMap.l2g
                #adjoint/characteristic speed might be discontinuous at element boundaries if heterogeneous
                self.adjoint_velocitySpace[ci]= FemTools.DG_AffineLinearOnSimplexWithNodalBasis(mesh,self.nd)
                self.adjoint_velocity[ci]     = FemTools.FiniteElementFunction(self.adjoint_velocitySpace[ci],dim_dof=self.nd,isVector=True,name="adjoint_velocity_%s" % ci)
                self.adjoint_velocity_interpolation_values[ci] = numpy.zeros((mesh.nElements_global,
                                                                              self.adjoint_velocitySpace[ci].referenceFiniteElement.interpolationConditions.nQuadraturePoints,
                                                                              self.nd),'d')
                self.adjoint_velocity_dofs[ci] = self.adjoint_velocity[ci].dof
                self.adjoint_velocity_l2g[ci]  = self.adjoint_velocitySpace[ci].dofMap.l2g
            #
            self.characteristic_velocitySpace = self.adjoint_velocitySpace
            self.characateristic_velocity     = self.adjoint_velocity
            self.characteristic_velocity_interpolation_values = self.adjoint_velocity_interpolation_values
            self.characteristic_velocity_dofs = self.adjoint_velocity_dofs
            self.characteristic_velocity_l2g  = self.adjoint_velocity_l2g

        elif self.velocitySpaceFlag in ['rt0','bdm1']:
            self.nVDOFs_element = mesh.nElementBoundaries_element
            if self.velocitySpaceFlag == 'bdm1':
                self.nVDOFs_element = mesh.nElementBoundaries_element*self.nd
            for ci in range(self.nc):
                self.velocity_dofs[ci] = numpy.zeros((self.mesh.nElements_global*self.nVDOFs_element),'d')
                self.velocity_l2g[ci]  = numpy.arange((self.mesh.nElements_global*self.nVDOFs_element),dtype='i').reshape((mesh.nElements_global,self.nVDOFs_element))
                self.velocity_interpolation_values[ci] = numpy.zeros((mesh.nElements_global,self.nVDOFs_element),'d')
                self.adjoint_velocity_dofs[ci] = numpy.copy(self.velocity_dofs[ci])
            self.need_ebq_arrays = True
            #
            #rt0 can use same representation because just uses local dofs, but doesn't really stay in the same space
            self.adjoint_velocity_interpolation_values = self.velocity_interpolation_values
            self.adjoint_velocity_l2g         = self.velocity_l2g
            self.characteristic_velocity_dofs = self.adjoint_velocity_dofs 
            self.characteristic_velocity_l2g  = self.adjoint_velocity_l2g

    def initializeElementQuadrature(self,t,cq):
        """
        Give the TC object access to the element quadrature storage
        """
        #could reuse trial functions from own transport class in some cases
        for ci in range(self.nc):
            self.q[('velocity',ci)] = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')
            if self.velocitySpace[ci] != None:
                self.q[('v',ci)]= numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],
                                              self.velocitySpace[ci].referenceFiniteElement.localFunctionSpace.dim),
                                             'd')
                #
                self.velocitySpace[ci].getBasisValues(self.elementQuadraturePoints,self.q[('v',ci)])
            #could reuse trial functions from own transport class in some cases
            self.q[('adjoint_velocity',ci)] = numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],self.nd),'d')
            if self.adjoint_velocitySpace[ci] != None:
                self.q[('adjoint_velocity_v',ci)]= numpy.zeros((cq['x'].shape[0],cq['x'].shape[1],
                                                               self.adjoint_velocitySpace[ci].referenceFiniteElement.localFunctionSpace.dim),
                                                              'd')
                #
                self.adjoint_velocitySpace[ci].getBasisValues(self.elementQuadraturePoints,self.q[('adjoint_velocity_v',ci)])
            self.q[('characteristic_velocity',ci)] = self.q[('adjoint_velocity',ci)]
    
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        """
        Give the TC object access to the exterior element boundary quadrature storage
        """
        for ci in range(self.nc):
            self.ebqe[('velocity',ci)] = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')
            if self.velocitySpace[ci] != None:
                self.ebqe[('v',ci)]= numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],
                                                 self.velocitySpace[ci].referenceFiniteElement.localFunctionSpace.dim),
                                                'd')
                self.velocitySpace[ci].getBasisValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,self.ebqe[('v',ci)])
            #
            self.ebqe[('adjoint_velocity',ci)] = numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],self.nd),'d')
            if self.adjoint_velocitySpace[ci] != None:
                self.ebqe[('adjoint_velocity_v',ci)]= numpy.zeros((cebqe['x'].shape[0],cebqe['x'].shape[1],
                                                                  self.adjoint_velocitySpace[ci].referenceFiniteElement.localFunctionSpace.dim),
                                                                 'd')
                self.adjoint_velocitySpace[ci].getBasisValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,self.ebqe[('adjoint_velocity_v',ci)])

    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        """
        Give the TC object access to the element boundary quadrature storage
        doesn't use ebq, ebq_global right now just initialize to zeros so it can run with models that do
        """
        for ci in range(self.nc):
            self.ebq[('velocity',ci)] = numpy.zeros((cebq['x'].shape[0],cebq['x'].shape[1],cebq['x'].shape[2],self.nd),'d')
            self.ebq_global[('velocity',ci)] = numpy.zeros((cebq_global['x'].shape[0],cebq_global['x'].shape[1],self.nd),'d')

    def evaluateVelocity(self,t,q):
        """

        """
        import Quadrature,cfemIntegrals,cpostprocessing
        evaluatedVelocity = False
        tEval = None
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if self.velocityFunctions != None:
            if self.velocitySpaceFlag in ['rt0','bdm1'] and self.need_ebq_arrays:
                self.setup_ebq_arrays()
            tEval = t
            for ci in range(self.nc):
                if self.velocitySpace[ci] != None:
                    for eN in range(self.mesh.nElements_global):
                        for k in range(self.velocitySpace[ci].referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                            self.velocity_interpolation_values[ci][eN,k,:] = self.velocityFunctions[ci](self.velocitySpace[ci].interpolationPoints[eN,k],t)
                    #
                    #evaluate velocity degrees of freedom from its interpolation conditions
                    self.velocity[ci].projectFromInterpolationConditions(self.velocity_interpolation_values[ci])
                    self.velocity[ci].getValues(self.q[('v',ci)],self.q[('velocity',ci)])
                    self.velocity[ci].getValuesGlobalExteriorTrace(self.ebqe[('v',ci)],self.ebqe[('velocity',ci)])
            
                #manual evaluation for now
                elif self.velocitySpaceFlag in ['rt0','bdm1']:
                    if self.velocitySpaceFlag == 'rt0':
                        #rt0 just use element integrals as nodal unknowns
                        #\todo move this loop out of python
                        for eN in range(self.mesh.nElements_global):
                            for ebN in range(self.mesh.nElementBoundaries_element):
                                integral = 0.0
                                for kb in range(self.ebq['x'].shape[-2]):#nElementBoundaryQuadraturePoints_elementBoundary):
                                    v = self.velocityFunctions[ci](self.ebq['x'][eN,ebN,kb],t)
                                    for I in range(self.ebq['n'].shape[-1]):#.nd):
                                        integral += v[I]*self.ebq['n'][eN,ebN,kb,I]*self.ebq['dS'][eN,ebN,kb]
                                    #mwf debug
                                    #print "setup RT0 eN= %s ebN= %s kb=%s v=[%s,%s] n=[%s,%s] integral=%s " % (eN,ebN,kb,v[0],v[1],ebq['n'][eN,ebN,kb,0],ebq['n'][eN,ebN,kb,1],integral)
                                self.velocity_interpolation_values[ci][eN,ebN] = integral
                                self.velocity_dofs[ci][self.velocity_l2g[ci][eN,ebN]] = self.velocity_interpolation_values[ci][eN,ebN]
                        #
                        cpostprocessing.getElementRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                           self.vt.mesh.elementNodesArray,
                                                                           self.vt.q['abs(det(J))'],
                                                                           self.vt.q['x'],
                                                                           self.velocity_dofs[ci],
                                                                           self.q[('velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                                                 self.vt.mesh.elementNodesArray,
                                                                                                 self.vt.mesh.elementBoundaryElementsArray,
                                                                                                 self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                 self.vt.q['abs(det(J))'],
                                                                                                 self.vt.ebqe['x'],
                                                                                                 self.velocity_dofs[ci],
                                                                                                 self.ebqe[('velocity',ci)])

                    elif self.velocitySpaceFlag == 'bdm1':
                       #bdm1 have to use projection
                       #<(\vec v - \vec \hat{v})\cdot \vec n_f,z>_{f} = 0 for faces gamma_f in each element, z in P^1(\gamma_f)
                        #\todo move this loop out of python
                        for eN in range(self.mesh.nElements_global):
                            for ebN in range(self.mesh.nElementBoundaries_element):
                                for kb in range(self.ebq['x'].shape[-2]):#nElementBoundaryQuadraturePoints_elementBoundary):
                                    self.ebq[('velocity',ci)][eN,ebN,kb,:] = self.velocityFunctions[ci](self.ebq['x'][eN,ebN,kb],t)
                        #
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        cpostprocessing.solveLocalBDM1projection(self.BDMprojectionMat_element,
                                                                 self.BDMprojectionMatPivots_element,
                                                                 self.ebq[('w*dS',0)],
                                                                 self.ebq['n'],
                                                                 self.ebq[('velocity',ci)],
                                                                 self.velocity_dofs[ci])
                        cpostprocessing.getElementBDM1velocityValuesLagrangeRep(self.vt.q[('v',ci)],
                                                                                self.velocity_dofs[ci],
                                                                                self.q[('velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(self.vt.mesh.elementBoundaryElementsArray,
                                                                                                      self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                      self.vt.ebqe[('v',ci)],
                                                                                                      self.velocity_dofs[ci],
                                                                                                      self.ebqe[('velocity',ci)])
            evaluatedVelocity = True

        #velocity functions != None
        elif self.flowModel != None and self.velocitySpaceFlag in ['rt0','bdm1']:
            tEval = self.flowModel.timeIntegration.t
            #todo have to make sure conservative flux is using flux representation for tracking or tell tracking
            #to use a different evaluation routine!
            #could/should get these directly from conservative flux?
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #todo add map from ci --> flowModel
            for ci in range(self.nc):
                flow_ci = self.flowVelocityComponentMap[ci]
                assert self.flowModel.conservativeFlux != None and self.flowModel.conservativeFlux.has_key((flow_ci))
                self.velocity_dofs[ci] = self.flowModel.q[('velocity_dofs',flow_ci)]
                self.velocity_l2g[ci] = self.flowModel.q[('velocity_l2g',flow_ci)]
                self.q[('velocity',ci)] = self.flowModel.q[('velocity',flow_ci)]
                self.ebqe[('velocity',ci)] = self.flowModel.ebqe[('velocity',flow_ci)]
                evaluatedVelocity = True
        return evaluatedVelocity
    def evaluateAdjointVelocity(self,t,q):
        import cpostprocessing
        tEval = t
        #mwf debug
        #import pdb
        #pdb.set_trace()
        for ci in range(self.nc):
            if self.adjoint_velocity_times[ci] != None and tEval > self.adjoint_velocity_times[ci]:
                self.adjoint_velocity_times_last[ci] = self.adjoint_velocity_times[ci]
                self.adjoint_velocity_dofs_last[ci].flat[:] = self.adjoint_velocity_dofs[ci]
                self.adjoint_velocity_times[ci]      = tEval
            if self.velocitySpace[ci] != None:
                #now get adjoint/characateristic velocity using average volume fraction over an element
                #assumes spaces are compatible locally, just allowing discontinuity for adjoint velocity
                self.adjoint_velocity_interpolation_values[ci].flat[:] = self.velocity_interpolation_values[ci].flat
                #\todo move this out of python
                for eN in range(self.mesh.nElements_global):
                    omega_e = 0.0
                    vol_e   = 0.0
                    for k in range(q['dV'].shape[1]):
                        omega_e += q['dV'][eN,k]*q[('dm',0,0)][eN,k]
                        vol_e   += q['dV'][eN,k]
                    self.adjoint_velocity_interpolation_values[ci][eN] *= vol_e/(omega_e+1.e-12) 
                #
                #evaluate velocity degrees of freedom from its interpolation conditions
                self.adjoint_velocity[ci].projectFromInterpolationConditions(self.adjoint_velocity_interpolation_values[ci])
                self.adjoint_velocity[ci].getValues(self.q[('adjoint_velocity_v',ci)],self.q[('adjoint_velocity',ci)])
                self.adjoint_velocity[ci].getValuesGlobalExteriorTrace(self.ebqe[('adjoint_velocity_v',ci)],self.ebqe[('adjoint_velocity',ci)])

                self.needToEvaluateVelocity = t <= 0.0 or self.flowIsTransient#True#False #TODO need to allow for transients
            #manual evaluation for now
            elif self.velocitySpaceFlag in ['rt0','bdm1']:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                #
                #\todo move this out of python
                #Warning only works because RT0 velocity dofs stored on each element without any interelement flux continuity assumes
                #can't modify adjoint_velocity_interpolation_values to be consistent with this for now
                self.adjoint_velocity_dofs[ci].flat[:] = self.velocity_dofs[ci].flat
                for eN in range(self.mesh.nElements_global):
                    omega_e = 0.0
                    vol_e   = 0.0
                    for k in range(q['dV'].shape[1]):
                        omega_e += q['dV'][eN,k]*q[('dm',0,0)][eN,k]
                        vol_e   += q['dV'][eN,k]
                    for i in range(self.nVDOFs_element):
                        self.adjoint_velocity_dofs[ci][self.adjoint_velocity_l2g[ci][eN,i]] *= vol_e/(omega_e + 1.0e-12)
                if self.velocityFunctions != None:
                    tEval = t
                    if self.velocitySpaceFlag == 'rt0':
                        cpostprocessing.getElementRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                           self.vt.mesh.elementNodesArray,
                                                                           self.vt.q['abs(det(J))'],
                                                                           self.vt.q['x'],
                                                                           self.adjoint_velocity_dofs[ci],
                                                                           self.q[('adjoint_velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                                                 self.vt.mesh.elementNodesArray,
                                                                                                 self.vt.mesh.elementBoundaryElementsArray,
                                                                                                 self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                 self.vt.q['abs(det(J))'],
                                                                                                 self.vt.ebqe['x'],
                                                                                                 self.adjoint_velocity_dofs[ci],
                                                                                                 self.ebqe[('adjoint_velocity',ci)])
                    elif self.velocitySpaceFlag == 'bdm1':

                        cpostprocessing.getElementBDM1velocityValuesLagrangeRep(self.vt.q[('v',ci)],
                                                                                self.adjoint_velocity_dofs[ci],
                                                                                self.q[('adjoint_velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(self.vt.mesh.elementBoundaryElementsArray,
                                                                                                      self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                      self.vt.ebqe[('v',ci)],
                                                                                                      self.adjoint_velocity_dofs[ci],
                                                                                                      self.ebqe[('adjoint_velocity',ci)])

                elif self.flowModel != None:
                    tEval = self.flowModel.timeIntegration.t
                    flow_ci = self.flowVelocityComponentMap[ci]
                    if self.flowModel.velocityPostProcessor.vpp_algorithms[flow_ci].localVelocityRepresentationFlag == 2:
                        cpostprocessing.getElementRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                           self.vt.mesh.elementNodesArray,
                                                                           self.vt.q['abs(det(J))'],
                                                                           self.vt.q['x'],
                                                                           self.adjoint_velocity_dofs[ci],
                                                                           self.q[('adjoint_velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(self.vt.mesh.nodeArray,
                                                                                                 self.vt.mesh.elementNodesArray,
                                                                                                 self.vt.mesh.elementBoundaryElementsArray,
                                                                                                 self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                 self.vt.q['abs(det(J))'],
                                                                                                 self.vt.ebqe['x'],
                                                                                                 self.adjoint_velocity_dofs[ci],
                                                                                                 self.ebqe[('adjoint_velocity',ci)])
                    elif self.flowModel.velocityPostProcessor.vpp_algorithms[flow_ci].localVelocityRepresentationFlag == 1:
                        cpostprocessing.getElementRT0velocityValues(self.vt.q['x'],
                                                                    self.adjoint_velocity_dofs[ci],
                                                                    self.q[('adjoint_velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryRT0velocityValues(self.vt.mesh.exteriorElementBoundariesArray,
                                                                                          self.vt.mesh.elementBoundaryElementsArray,
                                                                                          self.vt.ebqe['x'],
                                                                                          self.adjoint_velocity_dofs[ci],
                                                                                          self.ebqe[('adjoint_velocity',ci)])

                    elif self.flowModel.velocityPostProcessor.vpp_algorithms[flow_ci].localVelocityRepresentationFlag == 0:
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        cpostprocessing.getElementBDM1velocityValuesLagrangeRep(self.vt.q[('v',ci)],
                                                                                self.adjoint_velocity_dofs[ci],
                                                                                self.q[('adjoint_velocity',ci)])
                        cpostprocessing.getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(self.vt.mesh.elementBoundaryElementsArray,
                                                                                                      self.vt.mesh.exteriorElementBoundariesArray,
                                                                                                      self.vt.ebqe[('v',ci)],
                                                                                                      self.adjoint_velocity_dofs[ci],
                                                                                                      self.ebqe[('adjoint_velocity',ci)])

                else:
                    raise NotImplementedError
                self.needToEvaluateVelocity = t <= 0.0 or self.flowIsTransient#True#False #TODO need to allow for transients
            #mwf debug
            #import pdb
            #pdb.set_trace()
            if self.adjoint_velocity_times_last[ci] == None:
                self.adjoint_velocity_times_last[ci] = tEval
                self.adjoint_velocity_dofs_last[ci]  = numpy.copy(self.adjoint_velocity_dofs[ci])
            if self.adjoint_velocity_times[ci] == None:
                self.adjoint_velocity_times[ci] = tEval
        #nc
    def setup_ebq_arrays(self):
        import Quadrature,cfemIntegrals,cpostprocessing
        if self.need_ebq_arrays:
            boundaryQuadrature = Quadrature.SimplexGaussQuadrature(nd=self.nd-1,order=2)#Quadrature.GaussEdge(2)
            self.elementBoundaryQuadraturePoints = numpy.array([pt for pt in boundaryQuadrature.points],dtype='d')
            self.elementBoundaryQuadratureWeights= numpy.array([wt for wt in boundaryQuadrature.weights],dtype='d')

            nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
            self.ebq['x'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
            self.ebq['inverse(J)'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,self.vt.nSpace_global,self.vt.nSpace_global),'d')
            self.ebq['g'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,max(1,self.vt.nSpace_global-1),max(1,self.vt.nSpace_global-1)),'d')
            self.ebq['sqrt(det(g))'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebq['n'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary,self.vt.nSpace_global),'d')
            self.ebq['dS'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,nElementBoundaryQuadraturePoints_elementBoundary),'d')
            #todo add check that test space is C0P1
            self.vt.testSpace[0].elementMaps.getValuesTrace(self.elementBoundaryQuadraturePoints,
                                                            self.ebq['x'])
            self.vt.testSpace[0].elementMaps.getJacobianValuesTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebq['inverse(J)'],
                                                                    self.ebq['g'],
                                                                    self.ebq['sqrt(det(g))'],
                                                                    self.ebq['n'])
            cfemIntegrals.calculateElementBoundaryIntegrationWeights(self.ebq['sqrt(det(g))'],
                                                                     self.elementBoundaryQuadratureWeights,
                                                                     self.ebq['dS'])
            #need these at least for bdm1
            if self.velocitySpaceFlag == 'bdm1':
                self.ebq[('w*dS',0)] = numpy.zeros((self.mesh.nElements_global,
                                                    self.mesh.nElementBoundaries_element,
                                                    nElementBoundaryQuadraturePoints_elementBoundary,
                                                    self.vt.nDOF_test_element[0]),'d')
                self.ebq[('w',0)]   = numpy.zeros((self.mesh.nElements_global,
                                                   self.mesh.nElementBoundaries_element,
                                                   nElementBoundaryQuadraturePoints_elementBoundary,
                                                   self.vt.nDOF_test_element[0]),'d')
                self.ebq['hat(x)'] = numpy.copy(self.ebq['x'])
                self.vt.testSpace[0].elementMaps.getInverseValuesTrace(self.ebq['inverse(J)'],self.ebq['x'],self.ebq['hat(x)'])
                self.vt.testSpace[0].elementMaps.getPermutations(self.ebq['hat(x)'])
                self.vt.testSpace[0].getBasisValuesTrace(self.vt.testSpace[0].elementMaps.permutations,
                                                         self.ebq['hat(x)'],
                                                         self.ebq[('w',0)])

                cfemIntegrals.calculateWeightedShapeTrace(self.elementBoundaryQuadratureWeights,
                                                          self.ebq['sqrt(det(g))'],
                                                          self.ebq[('w',0)],
                                                          self.ebq[('w*dS',0)])
                self.BDMcomponent=0
                self.BDMprojectionMat_element = numpy.zeros((self.mesh.nElements_global,
                                                             self.nVDOFs_element,
                                                             self.nVDOFs_element),
                                                            'd')
                self.BDMprojectionMatPivots_element = numpy.zeros((self.mesh.nElements_global,
                                                                   self.nVDOFs_element),
                                                                  'i')
                cpostprocessing.buildLocalBDM1projectionMatrices(self.ebq[('w*dS',0)],
                                                                 self.ebq['n'],
                                                                 self.ebq[('w',0)],
                                                                 self.BDMprojectionMat_element)
                cpostprocessing.factorLocalBDM1projectionMatrices(self.BDMprojectionMat_element,
                                                                  self.BDMprojectionMatPivots_element)
                for ci in range(self.nc):
                    #todo could make these ebq_global, call them interpolation values
                    self.ebq[('velocity',ci)] = numpy.zeros((self.mesh.nElements_global,
                                                             self.mesh.nElementBoundaries_element,
                                                             nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.vt.nSpace_global),'d')

            self.need_ebq_arrays = False
        
    def evaluateReactions(self,t,c):
        if self.reactionTerms != None and c.has_key('x') and c.has_key(('r',0)):
            self.reactionTerms.evaluate(t,c)
    def attachModels(self,modelList):
        self.vt = modelList[self.meModelId]
        if self.flowModelId != None:
            self.flowModel = modelList[self.flowModelId]
        if self.velocityFunctions == None:
            assert self.flowModel != None
            #todo for now assume flow model has only 1 velocity, need a mapping to get other components?
            for ci in range(self.nc):
                flow_ci = self.flowVelocityComponentMap[ci]
                self.q[('velocity',ci)] = self.flowModel.q[('velocity',flow_ci)]
                self.ebqe[('velocity',ci)]= self.flowModel.ebqe[('velocity',flow_ci)]
                self.ebq[('velocity',ci)] = self.flowModel.ebq[('velocity',flow_ci)]
                self.ebq_global[('velocity',ci)]= self.flowModel.ebq_global[('velocity',flow_ci)]
    def evaluate(self,t,c):
        """
        TODO
          evaluate velocity is currently setting ebqe when c=q but need to make sure this is done
          before evaluate is called with c=ebqe
        """
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.evaluateReactions(t,c)
        evaluatedVelocity = False
        if self.q[('velocity',0)].shape == c[('df',0,0)].shape:
            if self.needToEvaluateVelocity:
                evaluatedVelocity = self.evaluateVelocity(t,c)
        for ci in range(self.nc):
            if self.q[('velocity',ci)].shape == c[('df',ci,ci)].shape:
                v = self.q[('velocity',ci)]
            elif self.ebqe[('velocity',ci)].shape == c[('df',ci,ci)].shape:
                v = self.ebqe[('velocity',ci)]
            elif self.ebq[('velocity',ci)].shape == c[('df',ci,ci)].shape:
                v = self.ebq[('velocity',ci)]
            elif self.ebq_global[('velocity',ci)].shape == c[('df',ci,ci)].shape:
                v = self.ebq_global[('velocity',ci)]
            else:
                print c[('df',ci,ci)].shape
                print "no v---------------------"
                raise RuntimeError
            if self.useC:
                self.groundwaterTransportCoefficientsEvaluate(self.omega,
                                                              self.d,
                                                              self.alpha_L,
                                                              self.alpha_T,
                                                              v,
                                                              c[('u',ci)],
                                                              c[('m',ci)],
                                                              c[('dm',ci,ci)],
                                                              c[('f',ci)],
                                                              c[('df',ci,ci)],
                                                              c[('a',ci,ci)])

        if evaluatedVelocity:
            self.evaluateAdjointVelocity(t,c)
