#!/usr/bin/env python

#from proteus.iproteus import *
from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus import MeshAdaptPUMI

#import mpi4py
#mpi = CDLL('libmpi.so.0', RTLD_GLOBAL)
#print dir(MeshAdaptPUMI)

MeshAdaptPUMI = MeshAdaptPUMI.MeshAdaptPUMI()
MeshAdaptPUMI.helloworld('hello: Entering MeshAdaptPUMI')
MeshAdaptPUMI.readGeomModel('yy')
MeshAdaptPUMI.readPUMIMesh('geom.sms')

mesh2 = MeshTools.Mesh()
mesh2.cmesh = cmeshTools.CMesh()

#mesh2.buildFromC(mesh2.cmesh)

#mesh = MeshAdaptPUMIDrvr.mesh_proteus
MeshAdaptPUMI.initProteusMesh(mesh2)
#MeshAdaptPUMIDrvr.ConstructFromPUMIMesh(mesh2.cmesh) #not working
#cmeshTools.CMesh_FromMesh(mesh)

#cmeshTools.writeTetgenFiles(cmesh,"tetgen",1)
#MeshAdaptPUMIDrvr.deleteProteusMesh(mesh2)


