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
MeshAdaptPUMI.loadModelAndMesh('???','geom.sms')

mesh = MeshTools.TetrahedralMesh()
mesh.cmesh = cmeshTools.CMesh()

#MeshAdaptPUMI.initProteusMesh(mesh)
MeshAdaptPUMI.ConstructFromSerialPUMIMesh(mesh.cmesh)

#mesh.cmesh2 = cmesh
#cmeshTools.CMesh_FromMesh(mesh)

mesh.buildFromC(mesh.cmesh)
print "Done reading in mesh"

#cmeshTools.write3dmFiles(cmesh,"mesh",1)
#mesh.printMesh()
print "meshInfo says : \n", mesh.meshInfo()
mesh.writeMeshEnsight("mesh","n")
#mesh.writeMeshXdmf(cmesh,ar,"mesh",0.0,False,False,0,False)

cmeshTools.writeTetgenFiles(mesh.cmesh,"tetgen",1)
#cmeshTools.writeTetgenFiles(mesh,"tetgen",1)
#mesh.writeTetgenFiles(mesh,"tetgen",1)
MeshAdaptPUMI.helloworld('Done MeshAdaptPUMI')

#cmeshTools.cmeshToolsDeleteMeshDataStructures(mesh2.cmesh)
#from pyadh import vtkViewers as vtkViewers
#vtkViewers.viewMesh(mesh2)
