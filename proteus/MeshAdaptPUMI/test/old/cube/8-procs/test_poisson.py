#!/usr/bin/env python

from proteus.iproteus import *
from proteus.Profiling import logEvent
from ctypes import *
from proteus import Domain
from proteus import MeshTools
from proteus import cmeshTools
from proteus import MeshAdaptPUMI
from proteus import Archiver
from tables import *

#from proteus import iproteus
import poisson_3d_p
import poisson_3d_c0p1_n

Profiling.logLevel=7
Profiling.verbose=True

from petsc4py import PETSc
from proteus import Comm
Comm.init()
comm = Comm.get()
MeshAdaptPUMI = MeshAdaptPUMI.MeshAdaptPUMI()
MeshAdaptPUMI.helloworld('hello: Entering MeshAdaptPUMI')
#MeshAdaptPUMI.readGeomModel('cube.smd')
#MeshAdaptPUMI.readPUMIMesh('cube.sms')

#domain = Domain.PUMIDomain(fileprefix="cube.sms",'cube.smd')
mesh = MeshTools.TetrahedralMesh()
mesh.cmesh = cmeshTools.CMesh()
#cmesh = cmeshTools.CMesh()

#MeshAdaptPUMI.ConstructFromSerialPUMIMesh(mesh.cmesh) #not working
#mesh.generateFromPUMI()

pList = [poisson_3d_p]
nList = [poisson_3d_c0p1_n]
so = default_so
so.name = pList[0].name = "poisson_3d_c0p1"+"pe"+`comm.size()`
so.sList=[default_s]
opts.logLevel=7
opts.verbose=True
opts.profile=True
opts.gatherArchive=True
nList[0].linearSolver=default_n.KSP_petsc4py
nList[0].multilevelLinearSolver=default_n.KSP_petsc4py
ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
ns.calculateSolution('poisson_3d_c0p1')

#mesh.cmesh2 = cmesh
#cmeshTools.CMesh_FromMesh(mesh)

mesh.buildFromC(mesh.cmesh)
print "Done reading in mesh" 

#print "meshInfo says : \n", mesh.meshInfo()
#mesh.writeMeshEnsight("mesh","n")

MeshAdaptPUMI.helloworld('Done MeshAdaptPUMI')

#from pyadh import vtkViewers as vtkViewers
#vtkViewers.viewMesh(mesh2)

