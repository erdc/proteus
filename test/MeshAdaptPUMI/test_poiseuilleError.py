import numpy
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdapt
from proteus import Domain
from proteus import Comm
from petsc4py import PETSc
from proteus import Quadrature
from proteus.default_n import *
import os

def test_poiseuilleError(verbose=0):
    """Test for loading gmsh mesh through PUMI, estimating error for 
    a Poiseuille flow case. The estimated error should be larger than the
    exact error in the seminorm"""
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
    domain.boundaryLabels = [1,2,3,4,5,6]

    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    comm = Comm.init()

    nElements_initial = mesh.nElements_global
    mesh.convertFromPUMI(domain,domain.AdaptManager.PUMIAdapter, domain.faceList, domain.regList,parallel = comm.size() > 1, dim = domain.nd)

    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"coordinates",mesh.nodeArray)


    rho = numpy.array([998.2,998.2])
    nu = numpy.array([1.004e-6, 1.004e-6])
    g = numpy.asarray([0.0,0.0,0.0])
    deltaT = 1.0 #dummy number 
    epsFact = 1.0 #dummy number 
    domain.AdaptManager.PUMIAdapter.transferPropertiesToPUMI(rho,nu,g,deltaT,deltaT,deltaT,epsFact)

    #Poiseuille Flow
    Ly=0.2
    Lz = 0.05
    Re = 100
    Umax = Re*nu[0]/Lz

    def vOfX(x):
        return 4*Umax/(Lz**2)*(x[2])*(Lz-x[2])
    def dvOfXdz(x):
        return 4*Umax/(Lz**2)*(Lz-2*x[2])

    #hard code solution
    vector=numpy.zeros((mesh.nNodes_global,3),'d')
    dummy = numpy.zeros(mesh.nNodes_global); 

    vector[:,0] = dummy
    vector[:,1] = 4*Umax/(Lz**2)*(mesh.nodeArray[:,2])*(Lz-mesh.nodeArray[:,2]) #v-velocity
    vector[:,2] = dummy
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"velocity", vector)

    scalar=numpy.zeros((mesh.nNodes_global,1),'d')
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"p", scalar)

    scalar[:,0] = mesh.nodeArray[:,2]
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"phi", scalar)
    del scalar

    scalar = numpy.zeros((mesh.nNodes_global,1),'d')+1.0
    domain.AdaptManager.PUMIAdapter.transferFieldToPUMI(b"vof", scalar)

    errorTotal=domain.AdaptManager.PUMIAdapter.get_local_error()


    # load the femspace with linear basis and get the quadrature points on a reference element
    elementQuadrature = Quadrature.SimplexGaussQuadrature(domain.nd,3)

    assert mesh.nNodes_element == 4 #confirm all of the elements have 4 nodes

    #hard code computation for H1 seminorm; ideally will be reformatted using the classes within proteus
    derivativeArrayRef = [[1,0,0],[0,1,0],[0,0,1],[-1,-1,-1]]
    error = 0
    for eID in range(mesh.nElements_global):
        nodes=mesh.elementNodesArray[eID]
        coords=[];
        for i in range(mesh.nNodes_element):
            coords.append(mesh.nodeArray[nodes[i]])
        J = numpy.matrix([[coords[0][0]-coords[3][0],coords[1][0]-coords[3][0],coords[2][0]-coords[3][0]],
                          [coords[0][1]-coords[3][1],coords[1][1]-coords[3][1],coords[2][1]-coords[3][1]],
                          [coords[0][2]-coords[3][2],coords[1][2]-coords[3][2],coords[2][2]-coords[3][2]]])
        invJ = J.I
        detJ = numpy.linalg.det(J)
        gradPhi_h = 0
        for k in range(len(elementQuadrature.points)):
            tempQpt = 0 
            zCoord = elementQuadrature.points[k][0]*coords[0][2] \
                +elementQuadrature.points[k][1]*coords[1][2] \
                +elementQuadrature.points[k][2]*coords[2][2] \
                +(1-elementQuadrature.points[k][0]-elementQuadrature.points[k][1]-elementQuadrature.points[k][2])*coords[3][2]
            for i in range(mesh.nNodes_element):
                temp = 0
                for j in range(domain.nd):
                    temp = temp + derivativeArrayRef[i][j]*invJ[j,2]
                tempQpt = tempQpt + vector[nodes[i]][1]*temp
            exactgradPhi = dvOfXdz([0,0,zCoord])
            gradPhi_h = gradPhi_h + tempQpt
            error = error + (exactgradPhi-gradPhi_h)**2*elementQuadrature.weights[k]*abs(detJ)

    error = sqrt(error)
    assert error < errorTotal
