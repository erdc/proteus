#!/usr/bin/env python
from MeshTools import *
from FemTools import *
import numpy


def getEdgeNodesInPhysicalSpace(p0,p1,order):
    """
    get Lagrange nodes for P^{order} assuming parameterization between
    2d nodes p0 and p1 ,
      p = p0*t + (1-t)*p1

    return is array that's (k+1) x 3
    """
    if order == 2:
        edgeRefNodes = (0.0,1.0,0.5)
    elif order == 3:
        edgeRefNodes = (0.0,1.0,1.0/3.0,2.0/3.0)
    elif order == 4:
        edgeRefNodes = (0.0,1.0,0.25,0.5,0.75)
    else: #default is k=1
        edgeRefNodes = (0.0,1.0)
    #end elif
    physNodes = Numeric.zeros((order+1,3),Numeric.Float)
    for i in range(order+1):
        physNodes[i,:] = p0*(1.0-edgeRefNodes[i]) + p1*edgeRefNodes[i]
    #end i loop

    return physNodes
def buildTraceOpArray(mesh,dgspace,verboseLevel=0):
    """
    build a proto trace operator array for going from element P^k space
    to an edge P^k space, where both are defined in terms of the usual
    Lagrangian shape functions

    dgspace is the elemental space, I'm not sure yet how the edge/face
    space will be defined

    format for output array is

       traceArray(edge,0:1,0:edgeDim-1,0:polyDim-1)

    where
       0,1     -- denotes the element neigbhor
       edgeDim -- the dimension of the edge space (k+1) for P^k in 2d
       polyDim -- the dimension of the element space P^k in 2d

    """
    polyOrder  = dgspace.polyOrder
    edgeDim    = polyOrder+1
    elemDim    = dgspace.localDim
    #global number of edges
    ngEdges = mesh.nElementBoundaries_global

    traceArray = Numeric.zeros((ngEdges,2,edgeDim,elemDim),Numeric.Float)

    #edge space nodal location is (ie,k,:) for ie global edge id, k polydim id
    physEdgeNodes = Numeric.zeros((ngEdges,edgeDim,3),Numeric.Float)

    #loop through mesh edges and neighboring elements that are defined
    for ie in range(ngEdges):
        #map edge reference nodes to physical space
        #assume a unique parameterization for edge from node 0 to
        # node 1 in mesh.edgeNodesArray
        #edge space nodal locations in physical space,
        #shape is (edgeDim,3)
        p0 = mesh.nodeArray[mesh.edgeNodesArray[ie,0],:]
        p1 = mesh.nodeArray[mesh.edgeNodesArray[ie,1],:]
        physEdgeNodes[ie,:,:]= getEdgeNodesInPhysicalSpace(
            p0,
            p1,
            polyOrder)

    #end ie loop
    if verboseLevel > 0:
        print('in buildTraceOp ')
        for ie in range(ngEdges):
            for k in range(edgeDim):
                print('edge ',ie,' edge dof ',k,' = ',physEdgeNodes[ie,k,:])
        #end printout
    #end if verbose

    if verboseLevel > 0:
        refCoordsDebug  = Numeric.zeros((elemDim,3),Numeric.Float)
        refCoordsDebug[0,:] = (0.0,0.0,0.0)
        refCoordsDebug[1,:] = (1.0,0.0,0.0)
        refCoordsDebug[2,:] = (0.0,1.0,0.0)

        print('in buildTraceOp ')
        for j in range(3):
            print('node ',j,' =\n',refCoordsDebug[j,:])
            for k,w in enumerate(dgspace.referenceSpace.basis):
                print('dg basis ', k,' value: ')
                print(w(refCoordsDebug[j,:]))
            #end k loop
        #end j loop
    #end if verbose

    #brute force compute inverse values of edge nodes on each element
    #holds refCoords for edges
    refCoords  = Numeric.zeros((edgeDim,3),Numeric.Float)
    #holds basis values for edges
    basisValues= Numeric.zeros((edgeDim,elemDim),Numeric.Float)
    #easiest if loop through all of the elements
    for elid in range(mesh.nElements_global):
        #loop through edges,
        #  get physical coords,
        #  map them back to ref space
        #  get basis values at them
        #  store basis values in trace op
        #end loop
        for i in range(mesh.nElementBoundaries_element):
            globedge = mesh.elementBoundariesArray[elid,i]
            #physLocal is local nodes in physical space, size is
            # polyOrder+1 x 3
            physLocal = physEdgeNodes[globedge,:,:]
            iamneig = 0 #which trace op setting (0 or 1)
            if elid == mesh.elementBoundaryElementsArray[globedge,1]:
                iamneig = 1
            #end if
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
        print('in buildTraceArray, result is ')
        for ie in range(ngEdges):
            print('glob edge ',ie,'\n 0 neig= \n',traceArray[ie,0,:,:])
            print('glob edge ',ie,'\n 1 neig= \n',traceArray[ie,1,:,:])
        #end ie
    #end if verbose

    return traceArray

#end buildTraceOp

def checkFunctionTrace(u,mesh,dgspace,traceArray,jumpTol=1.0e-8):
    """
    take a dg finite element function and see if values
    differ from left and right neighbors on an edge
    Basically, it's computing a jump term.
    For a continuous solution, should get no difference.
    """

    #global number of edges and elements
    ngEdges = mesh.nElementBoundaries_global
    ngElems = mesh.nElements_global
    #
    polyOrder  = dgspace.polyOrder
    edgeDim    = polyOrder+1
    elemDim    = dgspace.localDim
    #hold local element dofs
    locDofValues = Numeric.zeros((elemDim),Numeric.Float)
    for globedge in range(ngEdges):
        jumpValues   = Numeric.zeros((edgeDim),Numeric.Float)
        for neig in range(2):
            elid = mesh.elementBoundaryElementsArray[globedge,neig]
            if elid > -1: #element neighbor exists
                for k in range(elemDim):
                    ig = u.femSpace.dofMap.l2g[elid,k]
                    locDofValues[k] = u.dof[ig]
                #end k
                for j in range(edgeDim):
                    for k in range(elemDim):
                        jumpValues[j] += (1.0-2.*neig)*traceArray[globedge,neig,j,k]*locDofValues[k]
                    #end k
                #end j
            else: #what to do about edges on bounary
                jumpValues *= 0.0
            #end if
        #end neig loop
        if Numeric.sum(abs(jumpValues)) > jumpTol:
            print('Warning edge ',globedge,' jumpValues big= \n',jumpValues)
        #end if big jump
    #end for loop on edges
#end checkFunctionTrace
if __name__ == '__main__':
    #test out some of the FemTools functionality here
    from testMesh import *
    import math

    verboseLevel = 1
    #first just create a simple triangular mesh and look at it in a
    #couple of different ways
    Lx = 1.0   #domain length in x and y
    Ly = 1.0

    #number of nodes for rectangular grid upon which triangular mesh
    #will be built should get 2 triangles for each rectangle
    #(nx-1)(ny-1) in the original grid
    nx = 21
    ny = 21

    #flag for viewing mesh in construction
    #0 -- do nothing (default)
    #1 -- gnuplot
    #2 -- matlab
    viewMesh = 2
    mesh = constructTriangularMeshOnRectangle(Lx,Ly,nx,ny,viewMesh)

    print('mesh Info says \n',mesh.meshInfo())

    dgDofMap = DG1LagrangeDofMap(mesh)

    if verboseLevel > 0:
        print("DG1 dofMap l2g is ")
        for i in range(dgDofMap.l2g.shape[0]):
            sn = '\t %d %d %d' % tuple(dgDofMap.l2g[i,:])
            print(sn)
        #end for i
    #end if verbose

    #now create a DG1 finite element space
    globalDGspace = DG1_AffineOnSimplexLagrangeBasis(mesh)

    ngdim = globalDGspace.dim
    nelem = mesh.nElements_global
    polydim = 3
    #what about a member of the DG1 space
    u = ScalarFiniteElementFunction(globalDGspace)

    #try to just calculate u as an analytical function of x
    nodalRefPoints = Numeric.zeros((polydim,3),Numeric.Float)
    nodalRefPoints[0,:] = (0.0,0.0,0.0)
    nodalRefPoints[1,:] = (1.0,0.0,0.0)
    nodalRefPoints[2,:] = (0.0,1.0,0.0)

    if verboseLevel > 0:
        print('nodalRefPoints = \n',nodalRefPoints)
    #end verbose check

    print('finite element space says num in map family is ')
    print('\t', u.femSpace.mapFamily.nMaps)
    physNodes = Numeric.zeros((nelem,polydim,3),Numeric.Float)

    u.femSpace.mapFamily.getValues(nodalRefPoints,physNodes)

    if verboseLevel > 0:
        print('physical Nodes = \n',physNodes)
    #end verbose check

    #set value of u manually
    for ie in range(mesh.nElements_global):
        for iv in range(polydim):
            px = physNodes[ie,iv,0]
            py = physNodes[ie,iv,1]
            ig = u.femSpace.dofMap.l2g[ie,iv]
            u.dof[ig] = math.sin(math.pi*px)*math.cos(math.pi*py)
        #end local shape loop
    #end element loop

    if verboseLevel > 0:
        print('u dof values = \n',u.dof)
    #end verbose check

    saveScalarFEMfunctionMatlab(u,'sol','mesh')

    if verboseLevel > 0:
        #now look at edge --> element loops for dg algorithms
        testEdgeToElementMapping(mesh)
    #end verbose
    #try to build trace operators in crude way
    traceOp = buildTraceOpArray(mesh,globalDGspace,verboseLevel)

    #check if contiuous solution has no jumps?
    checkFunctionTrace(u,mesh,globalDGspace,traceOp)
