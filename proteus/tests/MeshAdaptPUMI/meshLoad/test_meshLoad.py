#!/usr/bin/env python

#from pyadh import *
#from proteus.iproteus import *
from ctypes import *
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdaptPUMI
from proteus import Archiver
from tables import *

from petsc4py import PETSc

import os
print os.getcwd()

testDir='./proteus/MeshAdaptPUMI/test/meshLoad/'
cubeMdl=testDir + 'cube.dmg'
cube670p1=testDir + 'pumi670/cube.smb'
MeshAdaptPUMI = MeshAdaptPUMI.MeshAdaptPUMI()
MeshAdaptPUMI.loadModelAndMesh(cubeMdl, cube670p1)

mesh = MeshTools.TetrahedralMesh()
mesh.cmesh = cmeshTools.CMesh()
MeshAdaptPUMI.constructFromSerialPUMIMesh(mesh.cmesh)
mesh.buildFromC(mesh.cmesh)

print "Done reading in mesh"
print "meshInfo says : \n", mesh.meshInfo()

mesh.writeMeshEnsight("mesh","n")
cmeshTools.writeTetgenFiles(mesh.cmesh,"tetgen",1)
