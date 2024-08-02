import numpy
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdapt
from proteus import Domain
from proteus import Comm
from petsc4py import PETSc
import os
import pytest

def test_gmshLoadAndAdapt(verbose=0):
    """Test for loading gmsh mesh through PUMI, estimating error and adapting for 
    a Couette flow case"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    Model=testDir + '/Couette.null'
    Mesh=testDir + '/Couette.msh'

    domain = Domain.PUMIDomain(manager=MeshAdapt.AdaptManager()) #initialize the domain

    domain.AdaptManager.modelDict = {'flow':0}
    domain.AdaptManager.sizeInputs = [b'error_erm']
    domain.AdaptManager.adapt = 1
    domain.AdaptManager.hmax = 0.01
    domain.AdaptManager.hmin= 0.008
    domain.AdaptManager.hphi= 0.008
    domain.AdaptManager.numIterations= 1
    domain.AdaptManager.targetError= 1


    domain.AdaptManager.PUMIAdapter.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    domain.AdaptManager.PUMIAdapter.setAdaptProperties(domain.AdaptManager)

    domain.faceList=[[80],[76],[42],[24],[82],[78]]
    domain.boundaryLabels=[1,2,3,4,5,6]

    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    comm = Comm.init()

    nElements_initial = mesh.nElements_global
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter,
                         domain.faceList,domain.regList, parallel = comm.size() > 1, dim = domain.nd)

    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"coordinates",mesh.nodeArray)

    rho = numpy.array([998.2,998.2])
    nu = numpy.array([1.004e-6, 1.004e-6])
    g = numpy.asarray([0.0,0.0,0.0])
    deltaT = 1.0 #dummy number
    epsFact = 1.0 #dummy number
    domain.AdaptManager.PUMIAdapter.transferPropertiesToPUMI(rho,nu,g,deltaT,deltaT,deltaT,epsFact)


    #Couette Flow
    Lz = 0.05
    Uinf = 2e-3
    #hard code solution
    vector=numpy.zeros((mesh.nNodes_global,3),'d')
    dummy = numpy.zeros((mesh.nNodes_global,),'d') 
    vector[:,0] = dummy
    vector[:,1] = Uinf*mesh.nodeArray[:,2]/Lz #v-velocity
    vector[:,2] = dummy
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"velocity", vector)
    del vector
    del dummy

    scalar=numpy.zeros((mesh.nNodes_global,1),'d')
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"p", scalar)
    
    scalar[:,0] = mesh.nodeArray[:,2]
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"phi", scalar)
    del scalar

    scalar = numpy.ones((mesh.nNodes_global,1),'d')
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"vof", scalar)
    del scalar
    
    errorTotal=domain.AdaptManager.PUMIAdapter.get_local_error()
    assert errorTotal<1e-14

    #assert domain.AdaptManager.PUMIAdapter.willAdapt()
    
    domain.AdaptManager.PUMIAdapter.adaptPUMIMesh(b"")
    
    mesh = MeshTools.TetrahedralMesh()
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter,
                     domain.faceList,
                     domain.regList,
                     parallel = comm.size() > 1,
                     dim = domain.nd)
    nElements_final = mesh.nElements_global
    assert nElements_final>nElements_initial

def test2DgmshLoadAndAdapt(verbose=0):
    """Test for loading gmsh mesh through PUMI, estimating error and adapting for 
    a 2D Couette flow case"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    Model=testDir + '/Couette2D.null'
    Mesh=testDir + '/Couette2D.msh'
    domain = Domain.PUMIDomain(dim=2,manager=MeshAdapt.AdaptManager()) #initialize the domain

    modelDict = {'flow':0}
    domain.AdaptManager.modelDict = modelDict
    domain.AdaptManager.sizeInputs = [b'error_erm']
    domain.AdaptManager.adapt = 1
    domain.AdaptManager.hmax = 0.01
    domain.AdaptManager.hmin= 0.008
    domain.AdaptManager.hphi= 0.008
    domain.AdaptManager.numIterations= 1
    domain.AdaptManager.targetError= 1

    domain.AdaptManager.PUMIAdapter.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    domain.AdaptManager.PUMIAdapter.setAdaptProperties(domain.AdaptManager)

    domain.faceList=[[14],[12],[11],[13]]
    domain.boundaryLabels=[1,2,3,4]

    mesh = MeshTools.TriangularMesh()
    mesh.cmesh = cmeshTools.CMesh()
    comm = Comm.init()

    nElements_initial = mesh.nElements_global
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter, domain.faceList,domain.regList, parallel = comm.size() > 1, dim = domain.nd)

    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"coordinates",mesh.nodeArray)

    rho = numpy.array([998.2,998.2])
    nu = numpy.array([1.004e-6, 1.004e-6])
    g = numpy.asarray([0.0,0.0])
    deltaT = 1.0 #dummy number
    epsFact = 1.0 #dummy number
    domain.AdaptManager.PUMIAdapter.transferPropertiesToPUMI(rho,nu,g,deltaT,deltaT,deltaT,epsFact)

    #Couette Flow
    Lz = 0.05
    Uinf = 2e-3
    #hard code solution
    vector=numpy.zeros((mesh.nNodes_global,3),'d')
    dummy = numpy.zeros(mesh.nNodes_global); 
    vector[:,0] = Uinf*mesh.nodeArray[:,1]/Lz #v-velocity
    vector[:,1] = dummy
    vector[:,2] = dummy
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"velocity", vector)
    del vector
    del dummy

    scalar=numpy.zeros((mesh.nNodes_global,1),'d')
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"p", scalar)

    scalar[:,0] = mesh.nodeArray[:,1]
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"phi", scalar)
    del scalar

    scalar = numpy.zeros((mesh.nNodes_global,1),'d')+1.0
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"vof", scalar)

    errorTotal=domain.AdaptManager.PUMIAdapter.get_local_error()
    assert errorTotal<1e-14

    #assert domain.AdaptManager.PUMIAdapter.willAdapt()
    
    domain.AdaptManager.PUMIAdapter.adaptPUMIMesh(b"")
    
    mesh = MeshTools.TriangularMesh()
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter,
                     domain.faceList,
                     domain.regList,
                     parallel = comm.size() > 1,
                     dim = domain.nd)
    nElements_final = mesh.nElements_global
    assert (nElements_final>nElements_initial)

def t2DmultiRegion(verbose=0):
    """Test for loading gmsh mesh through PUMI with multiple-regions"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    Model=testDir + '/TwoQuads.dmg'
    Mesh=testDir + '/TwoQuads.smb'
    domain = Domain.PUMIDomain(dim=2,manager=MeshAdapt.AdaptManager()) #initialize the domain

    domain.AdaptManager.reconstructedFlag=0 #this is used to indicate that no mesh reconstruction is being done.
    domain.AdaptManager.PUMIAdapter.loadModelAndMesh(bytes(Model,'utf-8'), bytes(Mesh,'utf-8'))
    domain.faceList=[[14],[12],[11],[13],[15],[16]]
    domain.boundaryLabels=[1,2,3,4,5,6]
    domain.regList=[[41],[42]]

    mesh = MeshTools.TriangularMesh()
    mesh.cmesh = cmeshTools.CMesh()
    comm = Comm.init()
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter, domain.faceList,domain.regList, parallel = comm.size() > 1, dim = domain.nd)
    assert mesh.elementMaterialTypes[0]==1
    assert mesh.elementMaterialTypes[-1]==2
