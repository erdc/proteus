#!/usr/bin/env python

#from pyadh import *
#from proteus.iproteus import *
from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus import MeshAdaptPUMI
from proteus import Archiver
from tables import *


from petsc4py import PETSc
#from mpi4py import MPI
#mpi = CDLL('libmpi.so.0', RTLD_GLOBAL)
#print dir(MeshAdaptPUMI)

MeshAdaptPUMI = MeshAdaptPUMI.MeshAdaptPUMI()
MeshAdaptPUMI.helloworld('hello: Entering MeshAdaptPUMI')
MeshAdaptPUMI.readGeomModel('yy')
MeshAdaptPUMI.readPUMIMesh('geom.sms')

#mesh = MeshTools.Mesh()
mesh = MeshTools.TetrahedralMesh()
mesh.cmesh = cmeshTools.CMesh()

#mesh2.buildFromC(mesh2.cmesh)

#mesh = MeshAdaptPUMIDrvr.mesh_proteus
MeshAdaptPUMI.initProteusMesh(mesh)
MeshAdaptPUMI.ConstructFromSerialPUMIMesh(mesh) #not working
#cmeshTools.CMesh_FromMesh(mesh)

print "Done reading in mesh" 

#cmeshTools.write3dmFiles(mesh,"mesh",1)
#mesh.printMesh(mesh)
mesh.writeMeshEnsight("mesh","n")
#mesh.meshInfo()
#mesh.writeMeshXdmf(cmesh,ar,"mesh",0.0,False,False,0,False)

#cmeshTools.writeTetgenFiles(mesh,"tetgen",1)
#cmeshTools.writeTetgenFiles(mesh,"tetgen",1)
#mesh.writeTetgenFiles(mesh,"tetgen",1)
MeshAdaptPUMI.helloworld('Done MeshAdaptPUMI')

#cmeshTools.cmeshToolsDeleteMeshDataStructures(mesh2.cmesh)
#from pyadh import vtkViewers as vtkViewers
#vtkViewers.viewMesh(mesh2)

