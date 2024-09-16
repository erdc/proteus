"""
Module for diagnostic utilities

.. inheritance-diagram:: proteus.DiagUtils
   :parts: 1
"""
from .EGeometry import *
from .MeshTools import *
from .FemTools import *
from .LinearAlgebraTools import *
from .LinearSolvers import *
from .Transport import *
from .Norms import *
from .Profiling import logEvent

def L2errorFEMvsAF(analyticalFunction,quadraturePointArray,quadratureWeightArray,
                    functionValueArray,T=None):
    """
    supposed to be L2 norm of error in vector quantity
    I think just using dot would cover both scalar and vector case
    """
    error=0.0
    range_nQuadraturePoints_element = list(range(quadraturePointArray.shape[1]))
    for eN in range(quadraturePointArray.shape[0]):
        for k in range_nQuadraturePoints_element:
            AF = analyticalFunction.uOfXT(quadraturePointArray[eN,k],T)
            eloc = functionValueArray[eN,k]-AF
            error += numpy.dot(eloc,eloc)*quadratureWeightArray[eN,k]
    error = sqrt(abs(error))
    return error

def getQuadraturePhysPointsAndWeights(mesh,femSpace,quadrature,verbose=0):
    """
    for convenience, hide steps for generating quadrature points and
    weights, with Jacobians, on physical mesh based on points and
    weights on reference element

    returns points array that's nelem x nquadloc x 3
            weight array that's nelem x nquadloc
    """
    nd = femSpace.referenceFiniteElement.referenceElement.dim
    nquad    = len(quadrature.points)
    qpoints  = numpy.zeros((nquad,3),'d')
    qweights = numpy.zeros(nquad,'d')
    for k,p in enumerate(quadrature.points):
        qpoints[k][:] = p
    for k,w in enumerate(quadrature.weights):
        qweights[k] = w

    quadX = numpy.zeros((mesh.nElements_global,nquad,3),'d')
    quadW = numpy.zeros((mesh.nElements_global,nquad),'d')
    jacTmp   = numpy.zeros((mesh.nElements_global,nquad,nd,nd),'d')
    jInvTmp  = numpy.zeros((mesh.nElements_global,nquad,nd,nd),'d')
    detJTmp  = numpy.zeros((mesh.nElements_global,nquad),'d')

    femSpace.elementMaps.getValues(qpoints,quadX)
    femSpace.elementMaps.getJacobianValues(qpoints,jacTmp,
                                           jInvTmp,detJTmp)

    for eN in range(mesh.nElements_global):
        for k in range(nquad):
            quadW[eN,k] = abs(detJTmp[eN,k])*qweights[k]
        #end k
    #end eN

    return quadX,quadW,qpoints,qweights

def getFEMvals(u,xiArray,verbose=0):
    """
    for convenience, hide steps for generating finite element solution
    at physical points corresponding to reference points held in
    xiArray

    returns array that's nelem x npointloc
    """
    nelems = u.femSpace.elementMaps.mesh.nElements_global
    ndofs  = u.femSpace.referenceFiniteElement.localFunctionSpace.dim
    nploc  = xiArray.shape[0]
    bvals = numpy.zeros((nelems,nploc,ndofs),'d')
    uvals = numpy.zeros((nelems,nploc),'d')
    u.femSpace.getBasisValues(xiArray,bvals)

    for eN in range(nelems):
        for k in range(nploc):
            for j in range(ndofs):
                J = u.femSpace.dofMap.l2g[eN,j]
                uvals[eN,k] += bvals[eN,k,j]*u.dof[J]
                logEvent("""getFemValues eN=%d xiArray[%d]= %s
                    jloc=%d J=%d u.dof[%d]= %g
                    uvals[%d,%d]= %g
                    """ % (eN,k,xiArray[k],j,J,J,u.dof[J],eN,k,uvals[eN,k]),level=3)
                #end verbose
            #end j
        #end k
    #end eN
    return uvals

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#stuff for running test problems
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
identity tensor in 1d,2d,3d
"""
Ident1 = numpy.ones((1,1),'d')

Ident2 = numpy.zeros((2,2),'d')
for k in range(2):
    Ident2[k,k] = 1.0
#end k
Ident3 = numpy.zeros((3,3),'d')
for k in range(3):
    Ident3[k,k] = 1.0
#end k

# # # # # # # # # # # # # # # # # # # # # # # # #
#some useful test routines
# # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # #
#examples for testing a new finite element space
# # # # # # # # # # # # # # # # # # # # # # # # #
def testCrRavNodalBasis(nd,verbose=0):
    """
    test local Crouzeix-Raviart element space

    """
    if verbose > -1:
        print('creating CrouzeixRaviartWithNodalBasis space dim= ',nd)
    #end
    #look at values at element barycenter, and face barycenters
    npoints = 2 + nd
    xiArray = numpy.zeros((npoints,3),'d')
    if nd == 1:
        xiArray[:,0] = [0.5, 1., 0.]
    elif nd == 2:
        xiArray[:,0:2] = [[1./3., 1./3.],
                          [0.5, 0.5],
                          [0., .5],
                          [0.5, 0.]]
    elif nd == 3:
        xiArray[:,:] = [[1./4., 1./4., 1./4.],
                        [1./3., 1./3., 1./3.],
                        [0., 1./3., 1./3.],
                        [1./3., 0., 1./3.],
                        [1./3., 1./3., 0.]]
    #end if
    if verbose > 1:
        print('trying to get values at points \n',xiArray)
    space = CrouzeixRaviartWithNodalBasis(nd)
    #space = LinearOnSimplexWithNodalBasis(nd)
    if verbose > -1:
        print('number of dofs= \n',space.dim)
        bvals = numpy.zeros((npoints,space.dim),'d')
        bgrads= numpy.zeros((npoints,space.dim,nd),'d')
        for j in range(npoints):
            for k in space.range_dim:
                bvals[j,k] = space.basis[k](xiArray[j])
                bgrads[j,k,:]= space.basisGradients[k](xiArray[j])
            #end k
        #end j
        print('basis values at \n',xiArray)
        print('are \n',bvals)
        print('basis gradients are \n',bgrads)
    #end if

    #look at values on faces as mapping from lower dimensional space
    exiArray = numpy.zeros((1,max(nd-1,1)),'d')
    if nd == 1:
        exiArray[0,0] = 0.
    elif nd == 2:
        exiArray[0,0] = 0.5
    else:
        exiArray[0,0:2] = [1./3., 1./3.]
    #end else
    if verbose > -1:
        nElementBoundaries = nd+1
        bvals = numpy.zeros((nElementBoundaries,space.dim),'d')
        bgrads= numpy.zeros((nElementBoundaries,space.dim,nd),'d')
        for k in range(nElementBoundaries):
            for j in space.range_dim:
                bvals[k,j] = space.basisTrace[k][j](exiArray[0])
                bgrads[k,j,:] = space.basisGradientsTrace[k][j](exiArray[0])
            #end j
        #end k
        print('trace basis values at ',exiArray,' on edges 0:nd+1 are ')
        print(bvals)
        print('trace basis gradients are ')
        print(bgrads)
    #end if
#end testCr
def testQuadNodalBasis(nd,verbose=0):
    """
    test local P2 nodal finite element space

    """
    if verbose > -1:
        print('creating QuadraticOnSimplexWithNodal space dim= ',nd)
    #end
    #look at values at element barycenter, and face barycenters
    tdim = '1d'
    if nd == 2:
        tdim= '2d'
    #end if
    if nd == 3:
        tdim= '3d'
    #end if
    #npoints = nd+2
    #xiArray = numpy.zeros((npoints,nd),'d')
    xiArray = p2refNodes[nd-1]
    npoints = xiArray.shape[0]

    #end if
    if verbose > 1:
        print('trying to get values at points ',xiArray)
    space = QuadraticOnSimplexWithNodalBasis(nd)
    if verbose > -1:
        print('number of dofs= ',space.dim)
        bvals = numpy.zeros((npoints,space.dim),'d')
        bgrads= numpy.zeros((npoints,space.dim,nd),'d')
        if verbose > 6:
            for k in range(nd+1):
                print('baryCoord ',k,'(',xiArray[0],')=',baryCoords[tdim][k](xiArray[0]))
            #end k
            for k in space.range_dim:
                print('basis func ',k,'(',xiArray[0],')=',space.basis[k](xiArray[0]))
            #end k
        #end verbose
        for j in range(npoints):
            for k in space.range_dim:
                bvals[j,k] = space.basis[k](xiArray[j])
                bgrads[j,k,:]= space.basisGradients[k](xiArray[j])
            #end k
        #end j
        print('basis values at \n',xiArray)
        print('are \n',bvals)
        print('basis gradients are \n',bgrads)
    #end if

    #look at values on faces as mapping from lower dimensional space
    exiArray = numpy.zeros((1,max(nd-1,1)),'d')
    if nd == 1:
        exiArray[0,0] = 0.
    elif nd == 2:
        exiArray[0,0] = 0.5
    else:
        exiArray[0,0:2] = [1./3., 1./3.]
    #end else
    if verbose > -1:
        nElementBoundaries = nd+1
        bvals = numpy.zeros((nElementBoundaries,space.dim),'d')
        bgrads= numpy.zeros((nElementBoundaries,space.dim,nd),'d')
        for k in range(nElementBoundaries):
            for j in space.range_dim:
                bvals[k,j] = space.basisTrace[k][j](exiArray[0])
                bgrads[k,j,:] = space.basisGradientsTrace[k][j](exiArray[0])
            #end j
        #end k
        print('trace basis values at ',exiArray,' on edges 0:nd+1 are')
        print('are \n',bvals)
        print('trace basis gradients are \n',bgrads)
    #end if
#end testQuad


def testEdgeDOFMap(mesh,nd):
    """
    test edge dof map to see what its doing
    """

    #dofMap = EdgeDOFMap(mesh)
    dofMap = NodalDOFMap(mesh)
    if nd == 1:
        ndofLoc= 1
        #try to do a proto loop over elements and assemble local stiffness matrix

        stiffMat = numpy.array([[1.0,-1.0],
                                  [-1.0,1.0]])
    #end 1d
    elif nd == 2:
        ndofLoc= 3
        #try to do a proto loop over elements and assemble local stiffness matrix
        #what I'm getting out of diffusion jacobian for p1c0
        stiffMat = numpy.array([[0.5, 0., -0.5],
                                  [0., 0., 0.],
                                  [-0.5, 0., 0.5]])
        #what I'm getting out of diffusion jacobian for p1nc
        #stiffMat = numpy.array([[2.0, 0., -2.0],
        #                          [0., 0., 0.],
        #                          [-2.0, 0., 2.0]])

        stiffMat = numpy.array([[4.0, -2., -2.],
                                  [-2., 2., 0.],
                                  [-2., 0., 2.]])
    #end 2d

    A = Mat(dofMap.nDOF,dofMap.nDOF)
    for eN in range(mesh.nElements_global):
        for i in range(ndofLoc):
            ig = dofMap.l2g[eN,i]
            for j in range(ndofLoc):
                jg = dofMap.l2g[eN,j]
                print('loc(',i,',',j,') = ',stiffMat[i,j],' --> A(',ig,',',jg,')= ',A[ig,jg])
                A[ig,jg] += stiffMat[i,j]
            #end j
        #end i
    #end eN

    print('leaving testEdgeDofMap A= \n',A)

def testQuadRefMats(nd,verbose=0):
    """
    test quad reference matrices to see what its doing
    """

    lspace = QuadraticOnSimplexWithNodalBasis(nd)
    ndofLoc= lspace.dim
    volWeights = [1.0,0.5,1.0/6.0]

    #compute mass matrix numerically
    quadRule = SimplexGaussQuadrature(nd)
    quadRule.setOrder(4)

    stiffMat = numpy.zeros((lspace.dim,lspace.dim),'d')
    massMat = numpy.zeros((lspace.dim,lspace.dim),'d')
    for p,w in zip(quadRule.points,quadRule.weights):
        for i in lspace.range_dim:
            for j in lspace.range_dim:
                stiffMat[i,j] += numpy.dot(lspace.basisGradients[i](p),
                                             lspace.basisGradients[j](p))*w*volWeights[nd-1]
                massMat[i,j] += lspace.basis[i](p)*lspace.basis[j](p)*w*volWeights[nd-1]
            #end j
        #end i
    #end p,w
    print('P2 localStiffMat = \n',stiffMat)
    print('P2 localMassMat = \n',massMat)


#end testQuadDofMap

def testQuadDOFMap(mesh,nd,verbose=0):
    """
    test quad dof map to see what its doing
    """

    lspace = QuadraticOnSimplexWithNodalBasis(nd)
    #dofMap = NodalDOFMap(mesh)
    dofMap = QuadraticLagrangeDOFMap(mesh,lspace,nd)
    ndofLoc= lspace.dim
    volWeights = [1.0,0.5,1.0/6.0]

    #compute mass matrix numerically
    quadRule = SimplexGaussQuadrature(nd)
    quadRule.setOrder(4)

    stiffMat = numpy.zeros((lspace.dim,lspace.dim),'d')
    massMat = numpy.zeros((lspace.dim,lspace.dim),'d')
    for p,w in zip(quadRule.points,quadRule.weights):
        for i in lspace.range_dim:
            for j in lspace.range_dim:
                stiffMat[i,j] += numpy.dot(lspace.basisGradients[i](p),
                                             lspace.basisGradients[j](p))*w*volWeights[nd-1]
                massMat[i,j] += lspace.basis[i](p)*lspace.basis[j](p)*w*volWeights[nd-1]
            #end j
        #end i
    #end p,w
    if verbose > -1:
        print('P2 localStiffMat = \n',stiffMat)
        print('P2 localMassMat = \n',massMat)
    #end verbose

    if verbose > 2:
        print('testQuadNodalDOF locDof= ',ndofLoc,' global nDof=',dofMap.nDOF)
    A = Mat(dofMap.nDOF,dofMap.nDOF)
    for eN in range(mesh.nElements_global):
        for i in range(ndofLoc):
            ig = dofMap.l2g[eN,i]
            for j in range(ndofLoc):
                jg = dofMap.l2g[eN,j]
                print('loc(',i,',',j,') = ',stiffMat[i,j],' --> A(',ig,',',jg,')= ',A[ig,jg])
                A[ig,jg] += stiffMat[i,j]
            #end j
        #end i
    #end eN

    print('leaving testQuadDofMap A= \n',A)

#end testQuadDofMap

## @}