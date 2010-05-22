#!/usr/bin/env python

from MeshTools import *
from FemTools import *
from QuadTools import *
import numpy as Numeric

from testFemTools import *

class ProtoQFRKDG1disc:
    """
    An initial effort to implement QFRKDG method for triangular mesh in
    AdhPyUtil. This class uses a discontinuous, locally
    P^1 approximation and RK2 time integration for scalar linear advection.

    This version takes advantage of some of the existing functionality
    in AdhPyUtil, without being fully integrated, since I'm still learing
    a lot about the library and how to implement the QFRKDG algorithm
    efficiently.

    I will begin with native Python loops that are going to be very slow,
    until I make sure I've got things running correctly with the new finite
    element paradigm. Then I will look to moving to loops in c and integrating
    with the ScalarTransportSolver interface.

    """
    def __init__(self,mesh,velTimeDep=False):
        """
        Initialize solver. The first set of information will be 
        quantities holding dimensional information based on a  P^k
        approximation.

        """
        self.nPhysVar = 1
        self.nSpaceDim= 2
        self.polyOrder = 1
        #need this
        self.mesh = mesh
        #global finite element space
        self.dgspace = DG_AffineLinearOnSimplexWithNodalBasis(mesh,
                                                              self.nSpaceDim)
        #number of Runge Kutta Stages: k+1
        self.nRungeKuttaStage= 2
        #dimension of P^k
        self.nElemDof = self.dgspace.referenceFiniteElement.localFunctionSpace.dim
        #dimesion of P^k on an edge
        self.nEdgeDof = self.polyOrder+1
        #number of numerical quadrature points
        self.nElemQuad= 3
        #need to update velocities
        self.velTimeDep = velTimeDep
        self.velUpdated = False
        #build reference element matrices, numerical quadrature info,
        #
        self.setupReferenceInformation()
        #hold information for time integration interface
        self.systemSize = self.dgspace.dim
        #mwf QUESTION: should I make physEdgeNodeArray an entry
        #              in a dict called q although not technically a quadrature
        #              point. Still basic operation is getting values there
        #build trace operators for mesh
        self.traceArray,self.physEdgeNodeArray,self.unitNormalOut0,\
          self.edgeLenFactor,self.iamNeigArray = \
          buildElemEdgeInfoPk(
            self.mesh,self.dgspace,
            self.refNodes,
            self.refEdgeNodes,
            self.polyOrder,
            verboseLevel=0)
        
        #arrays for holding computed quanties
        ngEdges = self.mesh.nElementBoundaries_global
        ngElems = self.mesh.nElements_global
        #mwf QUESTION: should I make physEdgeNodeArray an entry
        #              in a dict called q although not technically a quadrature
        #              point. Still basic operation is getting values there
        #go ahead and get element nodal locations in physical space
        self.physElemNodeArray = Numeric.zeros((ngElems,self.nElemDof,3),
                                               Numeric.Float)
        
        self.dgspace.mapFamily.getValues(self.refNodes,
                                         self.physElemNodeArray)
        
        #get an upwinded flux value for every "edge dof/node" on each edge
        self.edgeFluxArray = Numeric.zeros((ngEdges,self.nEdgeDof),
                                             Numeric.Float)
        #get a physical flux (vector) value for every "element dof/node"
        #on each element 
        self.elemFluxArray = Numeric.zeros((ngElems,self.nElemDof,3),
                                             Numeric.Float)
        #store velocity values at element nodes separately to save
        #evaluations for steady-state computations
        self.physVelAtNodes= Numeric.zeros((ngElems,self.nElemDof,3),
                                             Numeric.Float)
        #default physical coefficients and auxililiary conditions
        self.IC = zeroIC
        self.BC = zeroBC
        self.velocity = constVel
        
        #for cfl calculations
        hxy = [2.0*elem.area/elem.diameter for elem in self.mesh.elementList]
        self.hxyArray = Numeric.array(hxy)
        self.cflFactor = -12345.0

        #mwf DEBUG:
        #print 'QFRDKG1disc init hxy = \n',self.hxyArray
        
        #mwf NOTE: wasting one unknown here. Should clean up RK stage
        #          formulation
        #current solution value always
        #self.sol = ScalarFiniteElementFunction(self.dgspace)
        #solution values for RK stages
        #self.solk = []
        #for k in range(nRungeKuttaStage):
        #    self.solk.append(ScalarFiniteElementFunction(self.dgspace))
        #end i loop
        
    #end init
    def setupReferenceInformation(self):
        """
        Build reference element quantities for the discretization. These are a
        function of the polynomial order and the shape function itself.

        mwf: TODO
          go through this routine and make sure I'm getting all of the
          information possible from the existing AdhPyUtil classes
        """
        #local constants that will be used in initialization
        one12 = 1.0/12.0
        one24 = 0.5*one12
        one3  = 1.0/3.0
        one6  = 0.5*one3
        zero  = 0.0
        eighteen = 18.0
        msix  = -6.0
        
        #reference nodes on an element. Is this in dgspace.referenceSpace?
        #keep as 3d nodes for now even though only 2d discretization
        self.refNodes = Numeric.zeros((self.nElemDof,3),Numeric.Float)
        #number vertex nodes counter clockwise starting at origin
        #for any order P^k. Then go back to origin and repeat in
        #counter clockwise ordering for remaining nodes
        #mwf NOTE: This is transpose of way stored in P2MESH code
        self.refNodes[0,:] = (0.0,0.0,0.0)  #xHat,yHat,zHat
        self.refNodes[1,:] = (1.0,0.0,0.0)
        self.refNodes[2,:] = (0.0,1.0,0.0)

        #nodal locations for edge's reference space
        #keep as 3d nodes even though these are really 1d
        self.refEdgeNodes = Numeric.zeros((self.nEdgeDof,3),Numeric.Float)
        #number vertex nodes from origin for any order P^k and
        #then go back to origin and repeat for higher order nodes
        #mwf NOTE: This is transpose of way stored in P2MESH code
        self.refEdgeNodes[0,:] = (0.0,0.0,0.0)
        self.refEdgeNodes[1,:] = (1.0,0.0,0.0)

        #mass matrix on reference element
        self.M = Numeric.zeros((self.nElemDof,self.nElemDof),Numeric.Float)
        self.M[:,:] = one24
        #diagonal
        for i in range(self.nElemDof):
            self.M[i,i] = one12
        #end diagonal

        #inverse mass matrix on reference element
        self.Minv = Numeric.zeros((self.nElemDof,self.nElemDof),Numeric.Float)
        self.Minv[:,:] = msix
        #diagonal
        for i in range(self.nElemDof):
            self.Minv[i,i] = eighteen
        #end diagonal

        #Dx and Dy (gradient arrays)
        self.Dx = Numeric.zeros((self.nElemDof,self.nElemDof),Numeric.Float)
        self.Dy = Numeric.zeros((self.nElemDof,self.nElemDof),Numeric.Float)

        self.Dx[0,:] = -one6
        self.Dx[1,:] =  one6
        self.Dx[2,:] =  zero

        self.Dy[0,:] = -one6
        self.Dy[1,:] =  zero
        self.Dy[2,:] =  one6

        #products of M^{-1} and D^{x,y}
        self.MinvDx = Numeric.zeros((self.nElemDof,self.nElemDof),
                                    Numeric.Float)
        self.MinvDy = Numeric.zeros((self.nElemDof,self.nElemDof),
                                    Numeric.Float)
        self.MinvDx = Numeric.dot(self.Minv,self.Dx)
        self.MinvDy = Numeric.dot(self.Minv,self.Dy)

        #element-edge mass matrices: nElemDof x nElemDof
        self.Me = []
        nedges = 3
        for i in range(nedges): #number of edges
            self.Me.append(Numeric.zeros((self.nElemDof,self.nElemDof),
                                         Numeric.Float))
        #end initialization for Me
        #mwf NOTE: change convention of edge numbering to match
        #          AdhPyUtil. Edges are numbered according to the
        #          vertex that they are opposite.
        #          Before, they were numbered counter clockwise starting
        #          at the origin
        #edge 0 --> hypo (was edge 1)
        #mwf debug
        #print 'setup ref info Me[0] = ',self.Me[0]
        self.Me[0][0,:] = zero
        self.Me[0][1,:] = (0.0, one3, one6)
        self.Me[0][2,:] = (0.0, one6, one3)

        #edge 1 --> y axis (was edge 2)
        self.Me[1][0,:] = (one3, 0.0, one6)
        self.Me[1][1,:] = zero
        self.Me[1][2,:] = (one6, zero, one3)

        #edge 2 --> x axis (was edge 0)
        self.Me[2][0,:] = (one3, one6, zero)
        self.Me[2][1,:] = (one6, one3, zero)
        self.Me[2][2,:] = zero

        #products M^{-1}*M^{e,i}
        self.MinvMe = []
        for i in range(nedges): #number of edges
            self.MinvMe.append(Numeric.dot(self.Minv,
                                                      self.Me[i]))
        #end end MinvMe loop

        self.elemQuadrature = SimplexGaussQuadrature(2)
        self.elemQuadrature.setOrder(2)
        self.nElemQuad = len(self.elemQuadrature.points)
    #end setup reference information

    def printReferenceInformation(self,filename='refInfoQFRKDG1.txt'):
        """
        print out the reference element information for the QFRKDG1
        discretization.

        
        """
        rout = open(filename,'w')
        rout.write("""
**************************************************
Reference element information
**************************************************
P%d : nElemDof= %d nEdgeDof= %d nRungeKuttaStage= %d 
        """ % (self.polyOrder,
               self.nElemDof,self.nEdgeDof,self.nRungeKuttaStage))

        rout.write("""
Element Nodal Locations:
%s
        """ % self.refNodes)

        rout.write("""
Edge Nodal Locations:
%s
      """ % self.refEdgeNodes)

        rout.write("""
Mass Matrix:
%s
        """ % self.M)
        
        rout.write("""
Inverse Mass Matrix:
%s
        """ % self.Minv)
                   
        rout.write("""
Dx Gradient Matrix:
%s
        """ % self.Dx)
                   
        rout.write("""
Dy Gradient Matrix:
%s
        """ % self.Dy)
                   
        rout.write("""
M^{-1}D^x Matrix:
%s
        """ % self.MinvDx)
                   
        rout.write("""
M^{-1}D^y Matrix:
%s
        """ % self.MinvDy)
                   
        rout.write("""
M^{e,2} Matrix (x axis):
%s
        """ % self.Me[2])
        
        rout.write("""
M^{e,0} Matrix (hypoteneuse):
%s
        """ % self.Me[0])
        
        rout.write("""
M^{e,1} Matrix (y axis):
%s
        """ % self.Me[1])
        
        rout.write("""
M^{-1}M^{e,2} Matrix (x axis):
%s
        """ % self.MinvMe[2])
        
        rout.write("""
M^{-1}M^{e,0} Matrix (hypoteneuse):
%s
        """ % self.MinvMe[0])
        
        rout.write("""
M^{-1}M^{e,1} Matrix (y axis):
%s
        """ % self.MinvMe[1])
        
        rout.write("""
Element Numerical Quadrature order %d:
points:
%s
weights:
%s
        """ % (self.elemQuadrature.order,
               self.elemQuadrature.points,
               self.elemQuadrature.weights))
        
    #end printReferenceInformation

    def setIC(self,icIn):
        """
        set function for initial condition at physical locations
        """
        self.IC = icIn
    #end setIC
    def setBC(self,bcIn):
        """
        set function for boundary conditions at physical locations
        """
        self.BC = bcIn
    #end setBC
    def setVelocity(self,velIn):
        """
        set function for velocity at physical locations
        """
        self.velocity = velIn
    #end setVelocity

    def getMaxTimeStep(self):
        """
        return dt corresponding to cfl=1 for now

        """
        return self.cflFactor
    #
    def updateNumericalFluxes(self,u,t,velRep,boundaryRep,
                              verboseLevel=0):
        """
        compute numerical fluxes at the edge space nodal degrees of
        freedom. u holds the most current solution value

        Physical boundary conditions are given by boundaryRep. 

        Physical velocity is given by velRep

        Values are upwinded using unitNormalOut0 which is
        normal vector from neighbor 0 to neighbor 1
        """
        #for now, loop through all edges and treat differently
        #depending on whether or not they are internal or boundary
        ngEdges = self.mesh.nElementBoundaries_global
        ngElems = self.mesh.nElements_global

        #loop through all edges
        
        #for each edge
        #  get normal vector
        #  for each node on edge
        #    get physical coordinate
        #    get physical velocity
        #    get solution value from each side 
        #    compute upwind numerical flux in (global) normal direction
        #  end node loop
        #end edge loop
        for ie in range(ngEdges):
            nn0= self.unitNormalOut0[ie,:]
            #mwf QUESTION: what impact would using the traceOperator
            #              on each coordinate have?
            ve = velRep(self.physEdgeNodeArray[ie,:,:],t) #velocity at nodes

            #which side to pick 
            veDotn= Numeric.zeros(self.nEdgeDof,Numeric.Float)
            for i in range(self.nEdgeDof):
                veDotn[i] = Numeric.dot(ve[i,:],nn0)
            #end for
            upDir = Numeric.where(veDotn > 0.0, 0,1)

            #global dof indices for 0 and 1 neighbors,
            #element on 0 side
            elem0 = mesh.elementBoundaryElementsArray[ie,0]
            #element on 1 side, may not exist if boundary edge
            elem1 = mesh.elementBoundaryElementsArray[ie,1]


            #mwf NOTE: assumes that they are consecutive
            #          and are numbered 0:nElemDof-1 consistently
            #          Could just as easily use ue[elem0,:] for
            #          ue = reshape(u.dof,(ngElems,nElemDof))
            gdof0 = self.dgspace.dofMap.l2g[elem0,:]
            #solution on 0
            u0 = Numeric.array([u.dof[i] for i in gdof0])
            #now get trace
            utr0= Numeric.dot(self.traceArray[ie,0,:,:],u0)
            #mwf debug
            #print 'update numerical fluxes, tr= \n',self.traceArray[ie,0,:,:]
            #print 'u0= \n',u0
            
            utr1= utr0
            if elem1 > -1:
                gdof1= self.dgspace.dofMap.l2g[elem1,:]
                u1   = Numeric.array([u.dof[i] for i in gdof1])
                utr1 = Numeric.dot(self.traceArray[ie,1,:,:],u1)
            else:
                utr1 = boundaryRep(self.physEdgeNodeArray[ie,:,:],t)
                #this should take care of inflow/outflow boundaries
                #since utr1 will get picked when flow is into element 0
                #and utr0 will get picked when flow is out of element 0
                #which is out of domain
            #end check on boundary
            #get upwind values now
            u4up  = tuple([utr0,utr1])
            uup   = Numeric.zeros(self.nEdgeDof,Numeric.Float)
            for i in range(self.nEdgeDof):
                uup[i] = u4up[upDir[i]][i]
            #end i
            self.edgeFluxArray[ie,:] = veDotn*uup
            if verboseLevel > 0:
                print """
in updateNumericalFluxes ie= %d
ve      =\n%s
nn0     = %s
upDir   = %s
utr0    = %s
utr1    = %s
upw     = %s
numflux = %s
                      """ % (ie,ve,nn0,upDir,utr0,utr1,uup,
                             self.edgeFluxArray[ie,:])
            #end if on print out
        #end ie

    #end updateNumerical fluxes

    def updateElementalFluxes(self,u,t,velRep,verboseLevel=0):
        """
        determine physical flux representation in elemental space by
        nodal interpolation. By this I mean go through and calculate
          \hat{f}^x = \sum (v_x \hat{u})_{i}N_i(x,y,z) and
          \hat{f}^y = \sum (v_y \hat{u})_{i}N_i(x,y,z)

        where
          (v_x u)_i is v_x(\vec p_i)\hat{u}_i for nodes \{\vec p_i\}

        """

        #if I had a finite element representation of velocity I think
        #I could just use something like
        #v.femSpace.mapFamily.getValues(self.refNodes,velValues)

        #try to update velocities here if necessary
        self.updateVelocityInfo(t,velRep,verboseLevel)
        #shorter way that takes advantage of ordering
        ngElems = self.mesh.nElements_global
        ue = Numeric.reshape(u.dof,(ngElems,self.nElemDof))
        
        self.elemFluxArray[:,:,:] = self.physVelAtNodes[:,:,:]
        #mwf DEBUG
        #print 'updateElemFluxes v= \n',self.physVelAtNodes[:,:,:]
        #print 'ue= \n',ue
        for e in range(ngElems):
            for i in range(self.nElemDof):
                self.elemFluxArray[e,i,:] *= ue[e,i]
            #end i
            #print 'updateElemFluxes elem= ',e
            #print ' v= \n',self.physVelAtNodes[e,:,:]
            #print ' u= \n',ue[e,:]
            #print 'uv= \n',self.elemFluxArray[e,:,:]
        #end e
        #recall physElemNodeArray is dimensioned like
        #  physElemNodeArray(elem #, node #, 0:2)
        #try to evaluate velocity at all physical locations at once?
        #go ahead and try orig evaluation to double check
        #for elem in range(mesh.nElements_global):
        #    ve = velRep(self.physElemNodeArray[elem,:,:],t) #localDim x 3
        #    #dof for u on element
        #    gdof = self.dgspace.dofMap.l2g[elem,:]
        #    ue2  = Numeric.array([u.dof[i] for i in gdof])
        #    self.elemFluxArray[elem,:,:]  = ue2*ve
        #end original calculation      
        if verboseLevel > 2:
            print 'updateElemFluxes physElemNodeArray= ',
            print self.physElemNodeArray
            for elem in range(mesh.nElements_global):
                ve = velRep(self.physElemNodeArray[elem,:,:],t) #localDim x 3
                #dof for u on element
                gdof = self.dgspace.dofMap.l2g[elem,:]
                ue2  = Numeric.array([u.dof[i] for i in gdof])
                elemFluxArrayOrig  = ue2*ve
                print 'orig elemFlux= \n',elemFluxArrayOrig
                print '\n new elemFluxArray= \n',self.elemFluxArray[elem,:,:]
            #end for
        #end if verbose
    #end updateElementalFluxes

    def updateVelocityInfo(self,t,velRep,verboseLevel=0):
        """
        determine physical velocity representation in elemental space by
        nodal interpolation. By this I mean go through and calculate
          \hat{f}^x = \sum (v_x )_{i}N_i(x,y,z) and
          \hat{f}^y = \sum (v_y )_{i}N_i(x,y,z)

        where
          (v_x )_i is v_x(\vec p_i) for nodes \{\vec p_i\}

        """
        #skip if steady state velocity
        if self.velUpdated and not self.velTimeDep:
            return
        #
        #if I had a finite element representation of velocity I think
        #I could just use something like
        #v.femSpace.mapFamily.getValues(self.refNodes,velValues)

        #recall physElemNodeArray is dimensioned like
        #  physElemNodeArray(elem #, node #, 0:2)
        #try to evaluate velocity at all physical locations at once?
        if verboseLevel > 1:
            print 'updateVelocityInfo physElemNodeArray= ',
            print self.physElemNodeArray
        #end verboseLevel
        ngElems = self.mesh.nElements_global
        ve = [velRep(self.physElemNodeArray[i,:,:],t) for i in range(ngElems)]
        self.physVelAtNodes = Numeric.array(ve)
        if verboseLevel > 1:
            print 'updateVelocityInfo physVel= \n',self.physVelAtNodes
        #end if
        self.cflFactor = 12345.0
        for e in range(ngElems):
            vdot = [Numeric.dot(self.physVelAtNodes[e,i,:],
                                self.physVelAtNodes[e,i,:])
                    for i in range(self.nElemDof)]
            vnrm = max(Numeric.sqrt(vdot))
            self.cflFactor = min(self.cflFactor,self.hxyArray[e]/(vnrm+1.0e-12))
        #end e
        #end elem loop
        self.velUpdated = True
        if verboseLevel > 1:
            print 'updateVelocityInfo self.cflFactor= ',self.cflFactor
        #end if verbose
    #end updateElementalFluxes

    def calculateRHS(self,u,t,verboseLevel=0):
        """
        calculate QFDG discretization for conservative linear
        transport system as the right hand side of the ode

          \od{\vec u}{t} = \mat{M}^{-1}\vec f

        where
          \vec f is a combination of elemental flux and boundary flux
           terms.

        Assumes
            u        --- the most recent solution value
            t        --- physical time for evaluation
            
        """
        self.updateElementalFluxes(u,t,self.velocity,
                                   verboseLevel)
        self.updateNumericalFluxes(u,t,self.velocity,self.BC,
                                   verboseLevel)

        rhs = Numeric.zeros(self.dgspace.dim,Numeric.Float)

        ngElems = self.mesh.nElements_global
        ngEdges = self.mesh.nElementBoundaries_global
        nEdgeLoc= 3 #nd+1

        #do one big loop for now ...

        #loop through elements
        #  loop through edges
        #     get trace operator from edge to element
        #     compute edge flux in elemental space using trace transpose
        #     accumulate product of MinvMe*edgeFlux for edge into rhs
        #  end edge loop
        #  adjust elemental flux for reference space
        #  accumulate Minv*Dx and Minv*Dy for elemental flux into rhs
        #end loop
        for elem in range(ngElems):
            #hope this is right element
            invJac  = self.dgspace.mapFamily.getMap(elem).getInverseJacobian()
            rhsLoc  = Numeric.zeros(self.dgspace.localDim,Numeric.Float)
            edgeOut = Numeric.zeros(self.dgspace.localDim,
                                    Numeric.Float)

            for edge in range(nEdgeLoc):
                globedge = self.mesh.elementBoundariesArray[elem,edge]
                #brute force. Figure out how to get safely from chris
                #p0 = self.mesh.nodeArray[mesh.edgeNodesArray[globedge,0],:]
                #p1 = self.mesh.nodeArray[mesh.edgeNodesArray[globedge,1],:]
                #dx = p1-p0
                #elen = Numeric.sqrt(Numeric.dot(dx,dx))
                lenfact = self.edgeLenFactor[elem,edge]
                #is current element neighbor 0 or 1?
                iamneig = self.iamNeigArray[elem,edge] 
                #end if
                traceLoc = self.traceArray[globedge,iamneig,:,:]
                numLoc   = self.edgeFluxArray[globedge,:]
                edge2elem= Numeric.dot(Numeric.transpose(traceLoc),
                                                  numLoc)
                #mwf DEBUG:
                if verboseLevel > 2:
                    e2elmRaw = Numeric.zeros(self.nElemDof,Numeric.Float)
                    trt = Numeric.transpose(traceLoc)
                    for i in range(self.nElemDof):
                        for j in range(self.nEdgeDof):
                            e2elmRaw[i] += trt[i,j]*numLoc[j]
                        #end j
                    #end i
                    print 'compRHS e= ',elem,' edge= ',edge
                    print '\t edge2elem= \n',edge2elem
                    print '\t e2elmRaw=  \n',e2elmRaw
                    print 'lenfact= ',lenfact
                #end if
                edge2elem *= lenfact   #scale by length/|det J|
                edgeOut   += Numeric.dot(self.MinvMe[edge],
                                                    edge2elem)
            #end ege loop
            #
            #compute element flux vector g on reference element
            #mwf NOTE: QUESTION: what effect does orientation leading
            #          to negative determinants have here???
            gv = Numeric.zeros((self.nElemDof,3),Numeric.Float)
            fv = self.elemFluxArray[elem,:,:]
            for i in range(self.nElemDof):
                gv[i,0] = invJac[0,0]*fv[i,0]+invJac[0,1]*fv[i,1]
                gv[i,1] = invJac[1,0]*fv[i,0]+invJac[1,1]*fv[i,1]
            #end i
            if verboseLevel > 2:
                print 'rhs: elem= ',elem,'Jinv= \n',invJac
                for id in range(self.nElemDof):
                    print '\t xi= ',self.physElemNodeArray[elem,id,:]
                    print '\t vi= ',self.physVelAtNodes[elem,id,:]
                    print '\t ui= ',u.dof[elem*self.nElemDof+id]
                    print '\t fvi= ',fv[id,:],'\n\t gvi= ',gv[id,:]
                #end id
            #end verbose
            #now get derivative matrix contributions
            rhsLoc += Numeric.dot(self.MinvDx,gv[:,0])
            rhsLoc += Numeric.dot(self.MinvDy,gv[:,1])
            rhsLoc -= edgeOut #ege numerical flux contributiosn
            #mwf DEBUG:
            if verboseLevel > 2:
                print 'in computeRHS e= ',elem,'\n\tMinvDx*gvx= '
                print Numeric.dot(self.MinvDx,gv[:,0])
                print '\n\tMinvDy*gvy= '
                print Numeric.dot(self.MinvDy,gv[:,1])
            #end another local edge loop
            #embed in global unknown
            for i in range(self.nElemDof):
                ig = self.dgspace.dofMap.l2g[elem,i]
                rhs[ig] += rhsLoc[i]
            #end embedding
        #end global element loop
        return rhs
    #end calculateRHS

    def computeError(self,u,t,exactRep):
        """
        compute discrete L1,L2, and L_inf norm assuming exact solution
        is given pointwise by exactRep

        Uses numerical quadrature for L1

        Uses mass matrix as inner product to compute L2 error assuming
          representation of error in trial space

        Computes L_inf by taking max at nodal locations
        """

        #do I assume that exactRep is a ScalarFiniteElementFunction
        #already so that I can use the getValues functionality?

        L1err = -12345.0
        L2err = -12345.0
        LIerr = -12345.0
        Umass = 0.0
        ngElems = self.mesh.nElements_global
        #hold shape function values at quadrature points over
        #whole domain
        shapeVals = Numeric.zeros((ngElems,self.nElemQuad,self.nElemDof),
                                  Numeric.Float)
        self.dgspace.getBasisValues(self.elemQuadrature.points,shapeVals)
        
        #practice taking advantage of u unknown ordering
        uglob = Numeric.reshape(u.dof,(ngElems,self.nElemDof))

        L1sum = 0.0 ; L2sum = 0.0 ; 
        for elem in range(ngElems):
            det = Numeric.absolute(
                self.dgspace.mapFamily.getMap(elem).getJacobianDeterminant())
            uVals  = uglob[elem,:]
            exVals = exactRep(self.physElemNodeArray[elem,:,:],t)
            errVals= uVals-exVals
            l2tmp  = Numeric.dot(self.M,errVals)
            L2sum += Numeric.dot(errVals,l2tmp)*det
            #get error values at numerical quadrature points for L1 error
            errQuad= Numeric.dot(shapeVals[elem,:,:],errVals)
            l1tmp  = Numeric.absolute(errQuad)
            L1sum += Numeric.dot(l1tmp,self.elemQuadrature.weights)*det
            LIerr  = max(max(l1tmp),LIerr)
            #get mass too
            uQuad  = Numeric.dot(shapeVals[elem,:,:],uVals)
            Umass += Numeric.dot(uQuad,self.elemQuadrature.weights)*det
        #end elem loop
        L2err = math.sqrt(L2sum)
        L1err = L1sum
        return (L1err,L2err,LIerr,Umass)
    #end compute error
#end protoQFRKDG1

class ProtoSSPRKintegrator:
    """
    Crude implementation of Explicit RK schemes for ODE that are SSP for
    linear operators (aka) spatial discretizations in 

    \od{\vec y}{t} = \mat{L}\vec y

    See Gottlieb, Shu, Tadmor siam review article and notes

    This is implemented following Gottlieb and Shu etal formulation rather
    than more traditional one used in p2mesh code
    """

    def __init__(self,order,cfl):
        """
        basically just setup stage coefficient arrays based on order

        could use recursive formula if want to fill in entire
        alpha matrix (lower tridiagonal)
        """
        self.order = order
        self.cfl   = cfl
        alpha = Numeric.zeros(self.order,Numeric.Float)
        if order == 1:
            alpha[0] = 1.0
        elif order == 2:
            alpha[0] = 0.5;     alpha[1] = 0.5
        elif order == 3:
            alpha[0] = 1./3.;   alpha[1] = 1./2.; alpha[2] = 1./6.
        elif order == 4:
            alpha[0] = 3./8.;   alpha[1] = 1./3.; alpha[2] = 1./4.;
            alpha[3] = 1./24.;
        elif order == 5:
            alpha[0] = 11./30.; alpha[1] = 3./8.; alpha[2] = 1./6.;
            alpha[3] = 1./12.;  alpha[4] = 1./120.
        else:
            print 'order= ',order,' not implemented. Quitting !!!'
            sys.exit(1)
        #end if on order
        self.alpha= alpha
        self.rhs = None
    #end init
    def setRHS(self,rhsIn):
        """
        set a right hand side operator for integration
        """
        self.rhs = rhsIn
    def step(self,u,t,dt,verboseLevel=0):
        """
        take a step from t to t+dt assuming starting with u = u(t)
        returns tt holds value of time reached
        here u is a scalar finite element function
        """
        if verboseLevel > -1:
            print 'RK',self.order,' taking step from (',t,' --> ',t+dt,')'
        #end if
        if verboseLevel > 1:
            print 'RK',self.order,' u(',t,')= ',u.dof
        #end if
        sol = []
        for k in range(self.order+1):
            sol.append(ScalarFiniteElementFunction(u.femSpace))
        #end
        sol[0].dof[:] = u.dof[:]
        #take order-1 forward euler steps
        for k in range(1,self.order):
            teval    = t  #weight for time level for evaluating each stage,
                          #all zeros for the linear SSP RK discretizations 
            f        = self.rhs.calculateRHS(sol[k-1],teval,verboseLevel)
            sol[k].dof[:]   = sol[k-1].dof[:]
            sol[k].dof[:]  += dt*f[:]
        #end loop from 0 to s-2
        for k in range(self.order):
            sol[self.order].dof[:] += self.alpha[k]*sol[k].dof[:]
        #end k
        f  = self.rhs.calculateRHS(sol[self.order-1],teval)
        sol[self.order].dof[:] += self.alpha[self.order-1]*dt*f[:]

        if verboseLevel > 1:
            print 'RK',self.order,' u(',t+dt,')= ',sol[self.order].dof
        #end if
        tt = t+dt
        return sol[self.order],tt  
    #end step
    def calculateSolution(self,u,tin,tout,maxSteps=10000,
                          verboseLevel=0):
        """
        integrate from tin to tout 
        """
        #mwf FIX: hardwire max time step for now
        failed= 0
        dtMax = self.rhs.getMaxTimeStep()
        dtCFL = self.cfl*dtMax
        t = tin
        its = 0
        unew = u
        while t < tout and its < maxSteps:
            dt = dtCFL
            if t+dt > tout:
                dt = tout-t
            #end if
            tt=-1  #time reached
            unew,tt  = self.step(unew,t,dt,verboseLevel)
            t = tt
            its += 1
        #end while
        return failed,unew,t
# #
def buildElemEdgeInfoPk(mesh,dgspace,elemRefNodes,edgeRefNodes,
                        polyOrder,verboseLevel=0):
    """
    build a trace operator array for going from element P^k space
    to an edge P^k space, where both are defined in terms of the usual
    Lagrangian shape functions

    Also,...

    Returns physical locations of edge nodal degrees of freedom

    Builds an array that tells me what the normal
     vector for each edge is assuming that it is oriented from
     neighbor 0 to neighbor 1 according to the same convention
     that is in mesh.elementBoundaryElementsArray and traceArray.

    Builds an array telling me what neighbor (0 or 1) an element is for
    its local edges.
    
    Builds dumb array for looking up len/det factor for each edge
     that depends on the element and edge and whether or not the
     element is neighbor 0 of the edge or not
   
    dgspace is the elemental space, I m not sure yet how the edge/face
    space will be defined

    format for trace array is

       traceArray(edge,0:1,0:edgeDim-1,0:polyDim-1)

    where
       edge    -- global edge number
       0,1     -- denotes the element neigbhor
       edgeDim -- the dimension of the edge space (k+1) for P^k in 2d
       polyDim -- the dimension of the element space P^k in 2d

    format for physEdgeNodes is
       physEdgeNodes(edge,i,:) = p_i

    where
       p_i = (x_i,y_i,z_i) -- physical location
       i = 0:edgeDim-1     -- local id for edge space 

    format for unitNormalOut0 is
       unitNormalOut0(edge,:) -- normal vector at edge from neighbor 0

    format for iamNeigArray is
       iamNeigArray(e,i) = 1 if element e is neighbor 1 for edge i
                             wher i=0,2 is a  local edge numbering
                           0 otherwise
    format for edgeLenFactor is
       edgeLenFactor(e,i) -- +/- length_i/|det J_e|

    for global element e with local edge i. value is positive if
       element e is neighbor 0 for edge i
    """
    #if not (polyOrder == dgspace.polyOrder):
    #    print "buildTraceOpPk order mismatch ",polyOrder," ",dgspace.polyOrder
    #    return
    #
    edgeDim    = polyOrder+1
    elemDim    = dgspace.referenceFiniteElement.localFunctionSpace.dim
    #global number of edges
    ngEdges = mesh.nElementBoundaries_global
    ngElems = mesh.nElements_global
    nEdgeLoc= mesh.nElementBoundaries_element #nd+1
    traceArray = Numeric.zeros((ngEdges,2,edgeDim,elemDim),Numeric.Float)

    #edge space nodal location is (ie,k,:) for ie global edge id, k polydim id
    physEdgeNodes = Numeric.zeros((ngEdges,edgeDim,3),Numeric.Float)

    #normal vector from neighbor 0 to neighbor 1 (or outside domain)
    unitNormalOut0 = Numeric.zeros((ngEdges,3),Numeric.Float)

    #signed length / |det Jac| for element edge pairs
    edgeLenFactor  = Numeric.zeros((ngElems,nEdgeLoc),Numeric.Float)

    #array for keeping track of element-->edge relationship
    #essentially the inverse of elementBoundaryElementsArray
    iamNeigArray   = Numeric.zeros((ngElems,nEdgeLoc),Numeric.Int)
    
    #loop through mesh edges and neighboring elements that are defined
    for ie in range(ngEdges):
        #map edge reference nodes to physical space
        #assume a unique parameterization for edge from node 0 to
        # node 1 in mesh.edgeNodesArray
        #edge space nodal locations in physical space,
        #shape is (edgeDim,3)
        #mwf NOTE: need to synchronize this with refEdgeNodes in discretization
        p0 = mesh.nodeArray[mesh.edgeNodesArray[ie,0],:]
        p1 = mesh.nodeArray[mesh.edgeNodesArray[ie,1],:]
        physEdgeNodes[ie,:,:]= getEdgeNodesInPhysicalSpace(
            p0,
            p1,
            polyOrder) 

        #normal relative to orientation p0-->p1
        nx = p1[1]-p0[1]        #y^1-y^0
        ny = p0[0]-p1[0]        #x^0-x^1
        nn = sqrt(nx*nx+ny*ny)  #normalize by length
        nnx= nx/nn
        nny= ny/nn
        nnz= 0.0
        #midpoint of edge
        pmx = 0.5*(p0[0]+p1[0])
        pmy = 0.5*(p0[1]+p1[1])
        pmz = 0.5*(p0[2]+p1[2])
        
        #very clumsy as usual ....
        #always want to compute upwind value relative to normal
        #going from neighbor 0 to neighbor 1
        #element on 0 side
        elem0 = mesh.elementBoundaryElementsArray[ie,0]
        #element on 1 side, may not exist if boundary edge
        elem1 = mesh.elementBoundaryElementsArray[ie,1]
        #can I check here to see if normal is consistent with
        #outer normal going from 0 to 1
        #mwf CHECK: how to make sure get barycenter of correct element?
        #           can I index elementList safely?
        e0 = mesh.nodeArray[mesh.elementNodesArray[elem0,0]] 
        e1 = mesh.nodeArray[mesh.elementNodesArray[elem0,1]] 
        e2 = mesh.nodeArray[mesh.elementNodesArray[elem0,2]]
        #mwf NOTE: not always correct assumption on barycenter?
        xc0= (e0[X]+e1[X]+e2[X])/3.0
        yc0= (e0[Y]+e1[Y]+e2[Y])/3.0
        zc0= (e0[Z]+e1[Z]+e2[Z])/3.0
        bc2pm = (pmx-xc0,pmy-yc0,pmz-zc0)
        midDotNorm = bc2pm[0]*nnx + bc2pm[1]*nny + bc2pm[2]*nnz
        if midDotNorm < 0.0:
            #switch normal orientation
            nnx = -nnx
            nny = -nny
            nnz = -nnz
        #end
        unitNormalOut0[ie,:] = (nnx,nny,nnz)
    #end ie loop
    if verboseLevel > 0:
        print 'in buildTraceOpPk '
        for ie in range(ngEdges):
            for k in range(edgeDim):
                print 'edge ',ie,' edge dof ',k,' = ',physEdgeNodes[ie,k,:]
        #end printout
    #end if verbose

    #for evaluating determinants
    baryRef = Numeric.Array([1./3.,1./3.])
    jacArray    = Numeric.zeros((ngElems,1,2,2),Numeric.Float)
    jacInvArray = Numeric.zeros((ngElems,1,2,2),Numeric.Float)
    jacDetArray = Numeric.zeros((ngElems,1),Numeric.Float)
    dgspace.elementMaps.getJacobianValues(xiArray,jacArray,
                                          jacInvArray,jacDetArray)
    #brute force compute inverse values of edge nodes on each element
    #holds refCoords for edges
    refCoords  = Numeric.zeros((ngElems,edgeDim,3),Numeric.Float) 
    dgspace.elementMaps.getInverseValues(,refCoords)
    #holds basis values for edges
    basisValues= Numeric.zeros((ngElems,edgeDim,elemDim),Numeric.Float)
   
    #easiest if loop through all of the elements
    for elid in range(ngElems):
        det = Numeric.absolute(jacDetArray[elid])
        #loop through edges,
        #  get physical coords,
        #  map them back to ref space
        #  get basis values at them
        #  store basis values in trace op
        #end loop
        for i in range(nEdgeLoc):
            globedge = mesh.elementBoundariesArray[elid,i]
            #still have edgeList's around?
            elen     = mesh.edgeList[globedge].length
            #physLocal is local nodes in physical space, size is
            # polyOrder+1 x 3
            physLocal = physEdgeNodes[globedge,:,:]
            iamneig = 0 #which trace op setting (0 or 1)
            lenfact = elen/det
            if elid == mesh.elementBoundaryElementsArray[globedge,1]:
                iamneig = 1
                lenfact = -elen/det  #neighbor 1 has flux stored as inner flux
                iamNeigArray[elid,i] = 1
            #end if
            edgeLenFactor[elid,i]=lenfact
            for j in range(edgeDim): 
                refCoords[j,:] = dgspace.mapFamily.getMap(elid).getInverseValue(physLocal[j,:])
                for k,w in enumerate(dgspace.referenceSpace.basis):
                    basisValues[j,k] = w(tuple(refCoords[j,:]))
                    #if globedge == 47:
                    #    print 'edge ',globedge,' j= ',j,' k= ',k
                    #    print 'w(1.,0.,0.)= ',w((1.0,0.0,0.0))
                    #    print 'physLocal= ',physLocal[j,:]
                    #    print 'refCoords= ',refCoords[j,:]
                    #    print 'basisValues= ',basisValues[j,k]
                    #end if debug
                #end k
            #end j loop 
            traceArray[globedge,iamneig,:,:] = basisValues
            #if globedge == 47: #mwf debug
            #    print 'edge ',globedge,' physLocal= ',physLocal
            #    print '\t iamneig= ',iamneig,' elid = ',elid
            #    print 'refCoords=\n',refCoords
            #    print 'basisValues=\n',basisValues
            #end if 
        #end loop through local edges
    #end loop through elements
    #mwf debug
    #print out traceArray
    if (verboseLevel > 0):
        print 'in buildTraceArray, result is '
        for ie in range(ngEdges):
            print 'glob edge ',ie,'\n 0 neig= \n',traceArray[ie,0,:,:]
            print 'glob edge ',ie,'\n 1 neig= \n',traceArray[ie,1,:,:]
            elem0 = mesh.elementBoundaryElementsArray[ie,0]
            print '\t unit normal 0 (',elem0,') = ',unitNormalOut0[ie,:]
        #end ie
    #end if verbose
        
    return (traceArray,physEdgeNodes,unitNormalOut0,edgeLenFactor,iamNeigArray)
#end buildTraceOpArrayPk

def getScalarFemFunctionFromFunction(dgspace,refNodes,icRep,t0=0.0):
    """
    build an initial condition by evaluating function icRep at
    nodal values for dgspace
    """
    u = ScalarFiniteElementFunction(dgspace)
    physNodes = Numeric.zeros((dgspace.mesh.nElements_global,
                               dgspace.localDim,3),Numeric.Float)
    u.femSpace.mapFamily.getValues(refNodes,physNodes)
    #set value of u manually 
    for ie in range(mesh.nElements_global):
        px    = physNodes[ie,:,:]
        icVal = icRep(px,t0)
        ig    = u.femSpace.dofMap.l2g[ie,:] #assumes consecutive for slice?
        for i in range(len(ig)):
            u.dof[ig[i]] = icVal[i]
        #end i
    #end ie loop
    return u
def constVel(xArray,t):
    """
    return a velocity value at array of physical coordinates
    input is an array dimensioned
       nNodes x 3
    """
    
    velArray = Numeric.zeros(xArray.shape,Numeric.Float)
    velArray[:,:] = (1.0,0.0,0.0)
    #velArray[:,:] = (0.0,-1.0,0.0)
    return velArray
#end constantVelRep

def zeroBC(xArray,t):
    """
    return a boundary value at array of physical coordinates
    input is an array dimensioned
       nNodes x 1
    """
    bndArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    bndArray[:] = 0.0
    return bndArray

#end boundary rep

def constBC(xArray,t):
    """
    return a boundary value at array of physical coordinates
    input is an array dimensioned
       nNodes x 1
    """
    bndArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    bndArray[:] = 1.0
    return bndArray

#end boundary rep

def zeroIC(xArray,t0):
    """
    scalar analytical function 
    input is an array dimensioned
       nNodes x 3

    """
    icArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    icArray[:] = 0.0
    return icArray
#end 

def constIC(xArray,t0):
    """
    scalar analytical function 
    input is an array dimensioned
       nNodes x 3

    """
    icArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    icArray[:] = 1.0
    return icArray
#end 

def chrisCone(xArray,t0):
    """
    cone on [0,1]x[0,1]
    """
    radius = 1.0/8.0
    centerX= 0.25*math.sin(2.0*math.pi*t0)+0.5
    centerY= 0.25*math.cos(2.0*math.pi*t0)+0.5
    uArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    coneX  = xArray[:,0]-centerX
    coneY  = xArray[:,1]-centerY
    for i in range(len(coneX)):
        if coneX[i]*coneX[i] + coneY[i]*coneY[i] < radius*radius:
            uArray[i] = 0.25*(1.0+Numeric.cos(math.pi*coneX[i]/radius))*\
             (1.0+Numeric.cos(math.pi*coneY[i]/radius))
        #end if
    #end for

    return uArray
#end chrisCone

def chrisRotVel(xArray,t):
    """
    rotating velocity on unit square
    """
    velArray = Numeric.zeros(xArray.shape,Numeric.Float)
    velArray[:,0] = 2.0*math.pi*(xArray[:,1]-0.5)
    velArray[:,1] = 2.0*math.pi*(0.5 - xArray[:,0])
    return velArray
#end
def sinXsinY(xArray,t0):
    """
    scalar analytical function 
    input is an array dimensioned
       nNodes x 3

    """
    icArray = Numeric.zeros((xArray.shape[0]),Numeric.Float)
    icArray[:] = Numeric.sin(xArray[:,0]*math.pi)*Numeric.sin(xArray[:,1]*math.pi)
    return icArray
#end 

def readme():
    """
    helpful documentation string?
    """
    usage = """
    The protoQFRKDG1 module contains an initial attempt to implement a
    quadrature free Runge-Kutta discontinuous Galerkin discretization
    for the conservative form of the level set equation using P^1 on
    triangles in 2d. In essence, this is just conservative linear advection. 

    So far, it has been tested for the rotating cone test problem defined
    in the chrisCone and chrisRotVel functions. To run this problem on
    a mesh defined by 51 nodes in the x and y directions at a target cfl
    number of 0.3 (the max should be 1/3 = 1/(2p+1)) try

    python protoQFRKDG1.py --nx 51 --ny 51 --cfl 0.3 --maxIts 10001 > & outputPy51 &

    and wait. The output I got from this was

    At  t = 0.5
    L1err = 0.000686934
    L2err = 0.00275261
    LIerr = 0.0319275
    Umass = 0.0155379

    The solution is output using the saveScalarFEMfunctionMatlab
    function defined in testFemTools.py. To visualize the results in
    matlab, one should be able to type something like

      mesh ; sol0 ; u=u0;  exact ; uex=u; sol

    This loads in a matlab mesh representation in arrays p,e,t and a
    P^1, C0 version of the solution in the vector u. The DG solution
    is made continuous by just averaging the nodal solution values
    from neighboring elements.
    
    The reference element information is printed out in the file
    refInfoQFRKDG1.txt by default
       
    The ProtoQFRKDG1disc class defines the vast majority of the work
    and content of the algorithm. The actual DG femSpace is
      DG_AffineOnSimplexLagrangeBasis
    defined in testFemTools.py

    The function buildElemEdgeInfoPk does most of the work setting up the
    array data structures for negotiating the edge-element relationships
    in the assembly process.
    """
#end readme
if __name__ == '__main__':
    #test basic QFRKDG1 functionality
    from optparse import OptionParser
    import sys,shutil,os,re

    from testMesh import *
    #from testFemTools import *
    
    parser = OptionParser()
    #options controlling input details
    parser.add_option('-C','--cfl',
                      default=0.1,
                      help="""target Courant number
                           """)
    parser.add_option('-x','--nx',
                      default=11,
                      help="""number of nodes in x direction
                           """)
    parser.add_option('-y','--ny',
                      default=11,
                      help="""number of nodes in y direction
                           """)
    parser.add_option('-m','--maxIts',
                      default=10001,
                      help="""max number of time steps allowed 
                           """)
    parser.add_option('-P','--problem',
                      default='cone',
                      help="""type of problem to run
                           """)
    parser.add_option('-T','--tend',
                      default=0.5,
                      help="""final time
                           """)
    parser.add_option('-v','--verbose',default=0,
                      help="""verbosity level for simulation""")

    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args
    verb = int(options.verbose)
    #number of nodes for rectangular grid upon which triangular mesh
    #will be built should get 2 triangles for each rectangle
    #(nx-1)(ny-1) in the original grid
    cfl  = float(options.cfl)
    nx   = int(options.nx) 
    ny   = int(options.ny)
    prob = options.problem
    tend = float(options.tend)
    maxits= int(options.maxIts)
    #initial time
    t0 = 0.0
    #first just create a simple triangular mesh
    Lx = 1.0   #domain length in x and y
    Ly = 1.0

    #flag for viewing mesh in construction
    #0 -- do nothing (default)
    #1 -- gnuplot
    #2 -- matlab
    viewMesh = 2
    mesh = constructTriangularMeshOnRectangle(Lx,Ly,nx,ny,viewMesh)

    print 'mesh Info says \n',mesh.meshInfo()
    

    #create a discretization

    dgdisc = ProtoQFRKDG1disc(mesh)
    
    dgdisc.printReferenceInformation()

    exactRep = zeroIC
    if prob == 'cone':
        dgdisc.setIC(chrisCone)
        dgdisc.setBC(zeroBC)
        dgdisc.setVelocity(chrisRotVel)
        exactRep = chrisCone
    elif prob == 'cone-constVel':
        dgdisc.setIC(chrisCone)
        dgdisc.setBC(zeroBC)
        dgdisc.setVelocity(constVel)
        exactRep = chrisCone
    else:
        dgdisc.setIC(constIC)
        dgdisc.setBC(constBC)
        dgdisc.setVelocity(constVel)
    #end prob tag
    #dgdisc.setIC(sinXsinY)

    
    u0 = getScalarFemFunctionFromFunction(dgdisc.dgspace,
                                          dgdisc.refNodes,
                                          dgdisc.IC,t0)

    saveScalarFEMfunctionMatlab(u0,'sol0','mesh0')

    dgdisc.updateNumericalFluxes(u0,t0,dgdisc.velocity,dgdisc.BC,
                                 verboseLevel=verb)

    dgdisc.updateElementalFluxes(u0,t0,dgdisc.velocity,
                                 verboseLevel=verb)

    dgdisc.calculateRHS(u0,t0,verboseLevel=verb)


    rkint =  ProtoSSPRKintegrator(2,cfl)
    rkint.setRHS(dgdisc)

    failed,u,tf  = rkint.calculateSolution(u0,t0,tend,
                                        maxSteps=maxits,verboseLevel=verb)

    saveScalarFEMfunctionMatlab(u,'sol','mesh')

    errVals = dgdisc.computeError(u,tf,exactRep)
    errout = """
    At  t = %g
    L1err = %g
    L2err = %g
    LIerr = %g
    Umass = %g
    """ % (tf,errVals[0],errVals[1],errVals[2],errVals[3])
    print errout

    #look ate exact solution too
    uex = getScalarFemFunctionFromFunction(dgdisc.dgspace,
                                           dgdisc.refNodes,
                                           exactRep,tf)

    saveScalarFEMfunctionMatlab(uex,'exact')


#end main
