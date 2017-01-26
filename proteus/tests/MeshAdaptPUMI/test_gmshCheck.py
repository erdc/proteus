import numpy
from proteus import MeshTools
from proteus import cmeshTools
from proteus.MeshAdaptPUMI import MeshAdaptPUMI
from proteus import Domain
from proteus import Comm
from petsc4py import PETSc
from nose.tools import eq_ as eq
from nose.tools import ok_ as ok
import os

def test_gmshLoadAndAdapt(verbose=0):
    """Test for loading gmsh mesh through PUMI, estimating error and adapting for 
    a Couette flow case"""
    testDir=os.path.dirname(os.path.abspath(__file__))
    Model=testDir + '/Couette.null'
    Mesh=testDir + '/Couette.msh'

    domain = Domain.PUMIDomain() #initialize the domain
    domain.PUMIMesh=MeshAdaptPUMI.MeshAdaptPUMI(hmax=0.01, hmin=0.008, numIter=1,sfConfig='ERM',maType='isotropic')
    domain.PUMIMesh.loadModelAndMesh(Model, Mesh)
    domain.faceList=[[80],[76],[42],[24],[82],[78]]

    mesh = MeshTools.TetrahedralMesh()
    mesh.cmesh = cmeshTools.CMesh()
    comm = Comm.init()

    nElements_initial = mesh.nElements_global
    mesh.convertFromPUMI(domain.PUMIMesh, domain.faceList, parallel = comm.size() > 1, dim = domain.nd)

    domain.PUMIMesh.transferFieldToPUMI("coordinates",mesh.nodeArray)

    rho = numpy.array([998.2,998.2])
    nu = numpy.array([1.004e-6, 1.004e-6])
    g = numpy.asarray([0.0,0.0,0.0])
    domain.PUMIMesh.transferPropertiesToPUMI(rho,nu,g)

    #Couette Flow
    Lz = 0.05
    Uinf = 2e-3
    #hard code solution
    vector=numpy.zeros((mesh.nNodes_global,3),'d')
    dummy = numpy.zeros(mesh.nNodes_global); 
    vector[:,0] = dummy
    vector[:,1] = Uinf*mesh.nodeArray[:,2]/Lz #v-velocity
    vector[:,2] = dummy
    domain.PUMIMesh.transferFieldToPUMI("velocity", vector)
    del vector
    del dummy

    scalar=numpy.zeros((mesh.nNodes_global,1),'d')
    domain.PUMIMesh.transferFieldToPUMI("p", scalar)

    scalar[:,0] = mesh.nodeArray[:,2]
    domain.PUMIMesh.transferFieldToPUMI("phi", scalar)
    del scalar

    scalar = numpy.zeros((mesh.nNodes_global,1),'d')+1.0
    domain.PUMIMesh.transferFieldToPUMI("vof", scalar)

    domain.PUMIMesh.get_local_error()

    ok(domain.PUMIMesh.willAdapt(),1)

    domain.PUMIMesh.adaptPUMIMesh()
    
    mesh = MeshTools.TetrahedralMesh()
    mesh.convertFromPUMI(domain.PUMIMesh,
                     domain.faceList,
                     parallel = comm.size() > 1,
                     dim = domain.nd)
    nElements_final = mesh.nElements_global
    ok(nElements_final>nElements_initial)

if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_gmshCheck:test_gmshLoadAndAdapt')

