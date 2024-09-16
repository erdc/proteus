#!/usr/bin/env python
from pyadh import *
comm = Comm.init()
cmeshTools.commInit()
import vtk
vtkMesh = vtk.vtkUnstructuredGrid()
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName('mymesh.vtu')
reader.SetOutput(vtkMesh)
reader.Update()
mesh = MeshTools.Mesh()
mesh.cmesh = cmeshTools.CMesh()
from pyadh import cvtkviewers as cvtkviewers
cvtkviewers.getMeshFromVTKUnstructuredGrid(vtkMesh,mesh.cmesh)

if mesh.nNodes_element == 2:
    cmeshTools.computeGeometricInfo_edge(mesh.cmesh)
elif  mesh.nNodes_element == 3:
    cmeshTools.computeGeometricInfo_triangle(mesh.cmesh)
else:
    cmeshTools.computeGeometricInfo_tetrahedron(mesh.cmesh)
mesh.buildFromC(mesh.cmesh)
from pyadh import vtkViewers as vtkViewers
vtkViewers.viewMesh(mesh)

cmeshTools.commDestroy()
