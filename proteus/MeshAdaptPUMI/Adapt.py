from __future__ import division
from builtins import str
from builtins import range
import proteus
import sys
import numpy
from proteus import Profiling


def reconstructMesh(domain,mesh):

   if hasattr(domain,"PUMIMesh") and not isinstance(domain,proteus.Domain.PUMIDomain) :

     Profiling.logEvent("Reconstruct based on Proteus, convert PUMI mesh to Proteus")

     nd = domain.nd
     from scipy import spatial
     meshVertexTree = spatial.cKDTree(mesh.nodeArray)
     meshVertex2Model= [0]*mesh.nNodes_owned

     assert domain.vertices, "model vertices (domain.vertices) were not specified"
     assert domain.vertexFlags, "model classification (domain.vertexFlags) needs to be specified"

     for idx,vertex in enumerate(domain.vertices):
       if(nd==2 and len(vertex) == 2): #there might be a smarter way to do this
         vertex.append(0.0) #need to make a 3D coordinate
       closestVertex = meshVertexTree.query(vertex)
       meshVertex2Model[closestVertex[1]] = 1

     isModelVert = numpy.asarray(meshVertex2Model).astype("i")

     meshBoundaryConnectivity = numpy.zeros((mesh.nExteriorElementBoundaries_global,2+nd),dtype=numpy.int32)
     for elementBdyIdx in range(len(mesh.exteriorElementBoundariesArray)):
       exteriorIdx = mesh.exteriorElementBoundariesArray[elementBdyIdx]
       meshBoundaryConnectivity[elementBdyIdx][0] = mesh.elementBoundaryMaterialTypes[exteriorIdx]
       meshBoundaryConnectivity[elementBdyIdx][1] = mesh.elementBoundaryElementsArray[exteriorIdx][0]
       meshBoundaryConnectivity[elementBdyIdx][2] = mesh.elementBoundaryNodesArray[exteriorIdx][0]
       meshBoundaryConnectivity[elementBdyIdx][3] = mesh.elementBoundaryNodesArray[exteriorIdx][1]
       if(nd==3):
         meshBoundaryConnectivity[elementBdyIdx][4] = mesh.elementBoundaryNodesArray[exteriorIdx][2]

     domain.PUMIMesh.reconstructFromProteus2(mesh.cmesh,isModelVert,meshBoundaryConnectivity)


