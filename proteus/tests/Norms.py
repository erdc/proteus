import numpy as np
from past.utils import old_div
from proteus import MeshTools
from proteus import Quadrature
from proteus import FemTools

#Assumption is that simplical meshes are used

def getDummyMesh(IEN,nodeIDs):
    dummyMesh = MeshTools.Mesh()
    dummyMesh.nodeOffsets_subdomain_owned = [0,0]
    dummyMesh.globalMesh = dummyMesh
    dummyMesh.max_nNodeNeighbors_node = 0
    dummyMesh.nElements_global = IEN.shape[0]
    if(IEN.shape[1] == 3): #2D triangle has 3 vertices
        dummyMesh.dim = 2
    else: #3D
        dummyMesh.dim = 3

    dummyMesh.nodeArray=nodeIDs
    dummyMesh.elementNodesArray=IEN

    return dummyMesh

def get_L2_norm(h5file,field):

    IEN = np.array(h5file.root.elementsSpatial_Domain1)
    nodeIDs = np.array(h5file.root.Mesh_Spatial_Domain_1)
    nodeCoords = np.array(h5file.root.nodesSpatial_Domain1)

    dummyMesh = getDummyMesh(IEN,nodeIDs)

    elementQuadrature = Quadrature.SimplexGaussQuadrature(dummyMesh.dim, 3)
    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh=dummyMesh,nd=dummyMesh.dim)

    L2_norm_total = 0

    rowOfOnes = []
    referenceVolumeFactor = 0
    if(dummyMesh.dim == 2):
        rowOfOnes = [1.0,1.0,1.0]
        referenceVolumeFactor = 2.0
    else:
        rowOfOnes = [1.0,1.0,1.0,1.0]
        referenceVolumeFactor = 6.0

    for eID,ele in enumerate(IEN):
        vertices = nodeCoords[ele]
        
        if(dummyMesh.dim==2): #remove column of zeros
            vertices = np.delete(vertices,2,1) 
        volumeMatrix = np.vstack([vertices.transpose(),rowOfOnes])
        volume = abs(np.linalg.det(volumeMatrix)/referenceVolumeFactor)
        localBasis= femSpace.getBasisValuesRef(np.array(elementQuadrature.points))
        #jacobian
        J = volume/(old_div(1.0,referenceVolumeFactor))
        scalar =  np.matmul(field[dummyMesh.elementNodesArray[eID]].transpose(),localBasis.transpose())
        L2_norm = np.dot(np.square(scalar),np.array(elementQuadrature.weights))
        L2_norm = L2_norm*J
        L2_norm_total+=L2_norm

    L2_norm_total = np.sqrt(L2_norm_total)

    return L2_norm_total

def get_L2_vectorNorm(h5file,field):

    IEN = np.array(h5file.root.elementsSpatial_Domain1)
    nodeIDs = np.array(h5file.root.Mesh_Spatial_Domain_1)
    nodeCoords = np.array(h5file.root.nodesSpatial_Domain1)
    dummyMesh = getDummyMesh(IEN,nodeIDs)

    elementQuadrature = Quadrature.SimplexGaussQuadrature(dummyMesh.dim, 2)
    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh=dummyMesh,nd=dummyMesh.dim)

    L2_norm_total = np.array([0.0,0.0,0.0])

    rowOfOnes = []
    referenceVolumeFactor = 0
    if(dummyMesh.dim == 2):
        rowOfOnes = [1.0,1.0,1.0]
        referenceVolumeFactor = 2.0
    else:
        rowOfOnes = [1.0,1.0,1.0,1.0]
        referenceVolumeFactor = 6.0

    for eID,ele in enumerate(IEN):
        vertices = nodeCoords[ele]
        if(dummyMesh.dim==2): #remove column of zeros
            vertices = np.delete(vertices,2,1) 
        volumeMatrix = np.vstack([vertices.transpose(),rowOfOnes])
        volume = abs(np.linalg.det(volumeMatrix)/referenceVolumeFactor)
        localBasis=femSpace.getBasisValuesRef(np.array(elementQuadrature.points))
        J = volume/(old_div(1.0,referenceVolumeFactor))
        for i in range(dummyMesh.dim):
            vectorComp = np.dot(field[dummyMesh.elementNodesArray[eID],i],localBasis)
            L2_norm = np.dot(np.square(vectorComp),np.array(elementQuadrature.weights))
            L2_norm = L2_norm*J
            L2_norm_total[i]+=L2_norm

    L2_norm_total = np.sqrt(L2_norm_total)
    return L2_norm_total



