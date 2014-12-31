import numpy as np
#from .mock_model import MockModel

#variableNames = ['u0']

#nodeArray = np.array([[0, 0, 0],
#                      [1, 0, 0],
#                      [0, 1, 0],
#                      [1, 1, 0]])


# The majority of PointGauges.buildQuantityRow should be moved to either
# femTools or meshTools.  This will simplify test writing
# as that function can then be safely mocked.


# u = TODO
#t = 0
# mesh = TODO
# coefficients = TODO

# model = MockModel(variableNames, nodeArray, u, t, mesh, coefficients)

#cek mock setup of multilevel model
from proteus import Profiling, Comm, Domain, Archiver, TransportCoefficients, MeshTools, FemTools, Transport, Quadrature, NumericalFlux
from proteus import default_p as p
from proteus import default_n as n
from math import ceil,sqrt,pow
import xml.etree.ElementTree as ElementTree
Profiling.openLog("test_gauges.log",5)
Profiling.verbose=False
comm = Comm.init()
#simple physics
p.name="test_gauges"
p.nd = 3
class LinearSolution:
    def uOfXT(self,x,t):
        return (x[0]+x[1]+x[2])*(t+1.0);
p.initialConditions = {0:LinearSolution()}
p.dirichletConditions = {0: lambda x,flag: None}
p.domain = Domain.RectangularDomain(name="test_gauges_domain")
p.coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX = [lambda x: np.eye(p.nd,p.nd)], fOfX = [lambda x: 0],nc=1,nd=p.nd)
#simple/incomplete numerics
n.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis}
n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,3)
n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,3)
n.numericalFluxType = NumericalFlux.NoFlux
#build multilevel model
archive = Archiver.XdmfArchive(dataDir=".",filename="test_gauges_archive")
if p.nd == 1:
    nn = comm.size()+1
    mlMesh = MeshTools.MultilevelEdgeMesh(nn,1,1,
                                          p.domain.L[0],1,1,
                                          refinementLevels=1,
                                          nLayersOfOverlap=0,#1
                                          parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)#element
elif p.nd == 2:
    nnx = nny = int(ceil(sqrt(comm.size())))+1
    mlMesh = MeshTools.MultilevelTriangularMesh(nnx,nny,1,
                                                 p.domain.L[0],p.domain.L[1],1,
                                                 refinementLevels=1,
                                                 nLayersOfOverlap=0,
                                                 parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)
elif p.nd == 3:
    nnx = nny = nnz = int(ceil(pow(comm.size(),1.0/3.0)))+1
    mlMesh = MeshTools.MultilevelTetrahedralMesh(nnx,nny,nnz,
                                                  p.domain.L[0],p.domain.L[1],p.domain.L[2],
                                                  refinementLevels=1,
                                                  nLayersOfOverlap=0,
                                                  parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)
model = Transport.MultilevelTransport(p,n,mlMesh)
tnList=[0.0,1.0,2.0]
#
#gauges.attachModel(model,archive)
#

#step through tnList
m = model.levelModelList[-1]
m.setInitialConditions(p.initialConditions,tnList[0])
#gauges.calculate()
archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
tCount = 0
m.archiveFiniteElementSolutions(archive,tnList[0],tCount,initialPhase=True,writeVectors=True,meshChanged=True)
for t in tnList[1:]:
    tCount +=1
    m.setInitialConditions(p.initialConditions,t)
    m.archiveFiniteElementSolutions(archive,t,tCount,initialPhase=False,writeVectors=True,meshChanged=True)
    #gauges.calculate()
archive.close()
