#! /usr/bin/env python

"""
common definitions for simple miscible displacement problem

    Flow (pressure equation)
    
    \deld \vec q = s(x,t)

          \vec q =-\frac{K(x,t)}{\mu(c(x,t))} \grad p

    Transport equation      
    \pd{\theta c}{t} + \deld \left[\vec q c - \theta\ten{D}\nabla c\right] = s_c(x,t) 

    where \mu(c) = a*c + b
"""
from proteus import *
from proteus import SubsurfaceTransportCoefficients as STC

import numpy as np

class MiscibleDisplacementCoefficients_Flow(STC.SinglePhaseDarcyCoefficients):
    """
    Simplified transport coefficients for pressure equation in single phase flow with mild viscosity depedence

    Flow (pressure equation)
    
    \deld \vec q = s(x,t)

          \vec q =-\frac{K(x,t)}{\mu(c(x,t))} \grad p

    Base class does most of the work for evaluating coefficients
    """

    def __init__(self,
                 K_types,     #permeability function (x,t) that varies by element type
                 source_types,#source function (x,t) that varies by element type                 
                 nd=2,        #number of space dimensions
                 viscosity_a=0.0, viscosity_b=1.0,#mu = a*c + b
                 concentration_model_id=None,  #which number is the model for the concentration transport equation? 
                 timeVaryingCoefficients=False, #do the coefficients vary in time?
                 materialValuesLocallyConstant=False): #are the material functions constants? e.g., K_j(x,t) = K^0_j ?
        if concentration_model_id != None and viscosity_a <= 1.0e-16:
            print "Warning, specified concentration model with id {0} but no viscosity dependence mu=a*c+b with a={1:10.3e} b={2:10.3e} ".format(concentration_model_id,viscosity_a,viscosity_b)
        if concentration_model_id == None and viscosity_a > 1.0e-16:
            print "Warning, no specified concentration model but have viscosity dependence mu=a*c+b with a={0:10.3e} b={1:10.3e} ".format(viscosity_a,viscosity_b)
        
        self.concentration_model_id = concentration_model_id; self.concentration_model = None
        self.viscosity_a = viscosity_a 
        self.viscosity_b = viscosity_b 
        
        STC.SinglePhaseDarcyCoefficients.__init__(self,K_types,source_types,nc=1,nd=nd,
                                                  timeVaryingCoefficients=timeVaryingCoefficients,
                                                  materialValuesLocallyConstant=materialValuesLocallyConstant)
        
        self.variableNames=['p']
        
        def attachModels(self,modelList):
            if self.concentration_model_id != None: #grab a reference to the model that solves for concentration 
                assert 0 <= self.concentration_model_id and self.concentration_model_id < len(modelList)
                self.concentration_model = modelList[self.concentration_model_id]
                #assumes that the first unknown in the transport equation is the concentration of our species of interest
                #get references to the quadrature dictionaries ...
                #element quadrature points
                self.q_c    = self.concentration_model.q[('u',0)]
                #exterior boundary of domain
                self.ebqe_c = self.concentration_model.ebqe[('u',0)]
                #element boundary points treated as unique per element (e.g., for DG) 
                if self.concentration_model.ebq.has_key(('u',0)):
                    self.ebq_c = self.concentration_model.ebq[('u',0)]
                #element boundary points treated as unique per element boundary
                if self.concentration_model.ebq_global.has_key(('u',0)):
                    self.ebq_global_c = self.concentration_model.ebq_global[('u',0)]

    def initializeElementQuadrature(self,t,cq):
        #call parent classes function, then initialize quadrature point for concentration and viscosity values
        STC.SinglePhaseDarcyCoefficients.initializeElementQuadrature(self,t,cq)
        self.q_c = np.zeros(cq[('u',0)].shape,'d')
        self.q_mu = np.zeros(cq[('u',0)].shape,'d')
        self.q_mu.fill(self.viscosity_b)
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        #call parent classes function, then initialize element boundary quadrature point for concentration and viscosity  values
        STC.SinglePhaseDarcyCoefficients.initializeElementBoundaryQuadrature(self,t,cebq,cebq_global)
        #spatial coordinates have an 'extra' dimension (3). Don't always have a 'u' value
        ebq_shape = cebq['x'].shape[:-1];         ebq_global_shape = cebq_global['x'].shape[:-1]
        self.ebq_c = np.zeros(ebq_shape,'d');  self.ebq_global_c = np.zeros(ebq_global_shape,'d')
        self.ebq_mu= np.zeros(ebq_shape,'d');  self.ebq_global_mu= np.zeros(ebq_global_shape,'d')
        self.ebq_mu.fill(self.viscosity_b);       self.ebq_global_mu.fill(self.viscosity_b)

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        #call parent classes function, then initialize quadrature point for concentration and viscosity values
        STC.SinglePhaseDarcyCoefficients.initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe)
        self.ebqe_c = np.zeros(cebqe[('u',0)].shape,'d')
        self.ebqe_mu = np.zeros(cebqe[('u',0)].shape,'d')
        self.ebqe_mu.fill(self.viscosity_b)

    def evaluate(self,t,c):
        #takes care of K(x,t) and f(x,t)
        STC.SinglePhaseDarcyCoefficients.evaluate(self,t,c)
        #figure out which set of quadrature points this is
        if c['x'].shape == self.q_x_shape:
            mu = self.q_mu; conc = self.q_c;
        elif c['x'].shape == self.ebqe_x_shape:
            mu = self.ebqe_mu; conc = self.ebqe_c;
        elif c['x'].shape == self.ebq_global_x_shape:
            mu = self.ebq_global_mu; conc = self.ebq_global_c;
        elif c['x'].shape == self.ebq_x_shape:
            mu = self.ebq_mu; conc = self.ebq_c;
        else:
            raise NotImplementedError
        #mwf debug
        #import pdb
        #pdb.set_trace()
        #mu = a*c + b
        mu[:] = conc; mu[mu < 0] = 0.0  #trim negative values
        mu *= self.viscosity_a; mu += self.viscosity_b

        #check for divide by zero
        assert np.absolute(mu).min() > 0.0, "problem, viscosity has a zero value about to divide by min |mu| = {}".format(np.absolute(mu).min())
        #
        for i in range(c[('a',0,0)].shape[-1]):
            np.divide(c[('a',0,0)][...,i],mu,c[('a',0,0)][...,i])

    #eval
#Miscible Displacement Class


name = "MiscibleDisplacement"

### Spatial Domain ###
L = (1000,1000,1) #m
nd= 2

#temporal domain
ndays = 10.0
T = ndays*24.0*60.*60. #seconds


#just a simple rectangular domain for now
regular_domain = Domain.RectangularDomain(L[:nd],
                                          units="m")
#generate a .poly file in case we want to use an unstructured grid
regular_domain.writePoly('mesh')
unstructured_mesh = True
#'domain' is the name of the object that Proteus looks for to define geometry 
domain = Domain.PlanarStraightLineGraphDomain(fileprefix='mesh',name='md_mesh')
domain.boundaryTags = regular_domain.boundaryLegend

genMesh = True
refinement_level = 10 #define characteristic length
he = L[0]/float(refinement_level)

triangleOptions = "VApq30Dena{area:8.8f}".format(area=(he**2)/2.0)


### Material Properties ###
#homogeneous for now
nmaterial = 1
k_0 = 1.0e-8 #m^2
f_0 = 0.0    #1/s

Ident = np.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

permeabilities = {1:lambda x,t: k_0*Ident}
sources = {1:lambda x,t: f_0}
#check about this
sources[0] = sources[1]
permeabilities[0] = permeabilities[1]
#constant viscosity
mu_b= 1.0e-3 #Pa*s or kg/(m*s)
mu_a= 0.0


### Boundary Conditions ###
## Flow Equation
#Pressure boundary condition on bottom left and top right
#no flow everywhere else
inflow_length = L[1]*0.1
outflow_length= L[1]*0.1

pressure_inflow = 100.e3 #Pa or kg/(m*s^2)
pressure_outflow= 0.0

## Transport Equation
concentration_inflow = 1.0 # kg/m^3 ?
concentration_background = 0.0 
#
def pressure_bc(x,flag):
    print "flag = %s left= %s " % (flag, domain.boundaryTags['left'])
    if flag == domain.boundaryTags['left'] and x[1] <= inflow_length:
        return lambda x,t: pressure_inflow
    if flag == domain.boundaryTags['right'] and x[1] >= outflow_length:
        return lambda x,t: pressure_outflow

def noflux(x,flag):
    if flag in [domain.boundaryTags['bottom'],domain.boundaryTags['top']]:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['left'] and x[1] > inflow_length:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['right'] and x[1] < outflow_length:
        return lambda x,t: 0.0

#
def concentration_bc(x,flag):
    if flag == domain.boundaryTags['left'] and x[1] <= inflow_length:
        return lambda x,t: concentration_inflow

### Initial conditions ###
class ConstantIC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        return self.val

##Flow Equation
initialConditions_flow = {0:ConstantIC(pressure_outflow)}
##Transport Equation
initialConditions_trans = {0:ConstantIC(concentration_background)}


