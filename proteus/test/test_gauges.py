from math import ceil,sqrt,pow
import xml.etree.ElementTree as ElementTree

import numpy as np
from proteus import (Comm,
                     Domain,
                     Archiver,
                     TransportCoefficients,
                     MeshTools,
                     FemTools,
                     Transport,
                     Quadrature,
                     NumericalFlux)
from proteus import default_p as p
from proteus import default_n as n

from proteus.Gauges import PointGauges

def build1DMesh(p, nnx):
    return MeshTools.MultilevelEdgeMesh(nnx, 1, 1,
                                        p.domain.L[0], 1, 1,
                                        refinementLevels=1,
                                        nLayersOfOverlap=0,
                                        parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)
def build2DMesh(p, nnx, nny):
    return MeshTools.MultilevelTriangularMesh(nnx,nny,1,
                                              p.domain.L[0], p.domain.L[1],1,
                                              refinementLevels=1,
                                              nLayersOfOverlap=0,
                                              parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)

def build3DMesh(p, nnx, nny, nnz):
    return MeshTools.MultilevelTetrahedralMesh(nnx,nny,nnz,
                                               p.domain.L[0], p.domain.L[1], p.domain.L[2],
                                               refinementLevels=1,
                                               nLayersOfOverlap=0,
                                               parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)

def gauge_setup():
    comm = Comm.init()

    #Simplified Physics
    p.name="test_gauges"
    p.nd = 3

    class LinearSolution:
        def uOfXT(self,x,t):
            return (x[0]+10*x[1]+100*x[2])*(t+1.0);

    p.initialConditions = {0:LinearSolution()}
    p.dirichletConditions = {0: lambda x,flag: None}
    p.domain = Domain.RectangularDomain(name="test_gauges_domain")
    p.coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX = [lambda x: np.eye(p.nd, p.nd)],
                                                                       fOfX = [lambda x: 0], nc=1, nd=p.nd)
    #Simplified and Incomplete Numerics
    n.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis}
    n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,3)
    n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,3)
    n.numericalFluxType = NumericalFlux.NoFlux

    #build multilevel model
    archive = Archiver.XdmfArchive(dataDir=".",filename="test_gauges_archive")

    total_nodes = comm.size()

    if p.nd == 1:
        mlMesh = build1DMesh(p, total_nodes+1)
    elif p.nd == 2:
        nnx = nny = int(ceil(sqrt(total_nodes)))+1
        mlMesh = build2DMesh(p, nnx, nny)
    elif p.nd == 3:
        nnx = nny = nnz = int(ceil(pow(total_nodes,1.0/3.0)))+1
        mlMesh = build3DMesh(p, nnx, nny, nnz)

    model = Transport.MultilevelTransport(p,n,mlMesh)

    return model, archive

def test_gauges():
    tnList=[0.0,1.0,2.0]
    model, archive = gauge_setup()

    p = PointGauges(gauges=((('u0',), ((0.5, 0.5, 0), (1, 0.5, 0)))),
                    activeTime=(0, 2.5),
                    sampleRate=0.2,
                    fileName='combined_gauge_0_0.5_sample_all.csv')

    p.attachModel(model,archive)

    m = model.levelModelList[-1]
    m.setInitialConditions(p.initialConditions,tnList[0])
    p.calculate()

    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    tCount = 0
    m.archiveFiniteElementSolutions(archive,tnList[0],tCount,initialPhase=True,writeVectors=True,meshChanged=True)
    for t in tnList[1:]:
        tCount +=1
        m.setInitialConditions(p.initialConditions,t)
        m.archiveFiniteElementSolutions(archive,t,tCount,initialPhase=False,writeVectors=True,meshChanged=True)
        p.calculate()
    archive.close()


if __name__ == '__main__':
    import nose
    nose.main()