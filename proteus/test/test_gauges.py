from math import ceil, sqrt, pow

import numpy as np
import numpy.testing as npt

from proteus import (Comm,
                     Domain,
                     TransportCoefficients,
                     MeshTools,
                     FemTools,
                     Transport,
                     Quadrature,
                     NumericalFlux)
from proteus import default_p as p
from proteus import default_n as n

from proteus.Gauges import PointGauges

from nose.tools import eq_

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
            return (x[0]+10*x[1]+100*x[2])*(t+1.0)

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

    total_nodes = comm.size()

    if p.nd == 1:
        mlMesh = build1DMesh(p, total_nodes+1)
    elif p.nd == 2:
        nnx = nny = int(ceil(sqrt(total_nodes)))+1
        mlMesh = build2DMesh(p, nnx, nny)
    elif p.nd == 3:
        nnx = nny = nnz = int(ceil(pow(total_nodes, 1.0/3.0)))+1
        mlMesh = build3DMesh(p, nnx, nny, nnz)

    model = Transport.MultilevelTransport(p, n, mlMesh)

    return model, p.initialConditions

def run_gauge(p, time_list):
    model, initialConditions = gauge_setup()

    p.attachModel(model, None)

    m = model.levelModelList[-1]
    m.setInitialConditions(initialConditions, time_list[0])
    tCount = 0
    p.calculate()

    for t in time_list[1:]:
        tCount +=1
        m.setInitialConditions(initialConditions, t)
        p.calculate()


def parse_gauge_output(filename):
    with open(filename) as f:
        header = f.readline().split(',')
        eq_(header[0].strip(), 'time')
        gauge_names = [gauge_name.strip() for gauge_name in header[1:]]
        f.seek(0)
        data = np.genfromtxt(f, delimiter=",", skip_header=1)
    return gauge_names, data

def test_gauge_output():
    filename = 'test_gauge_output.csv'

    p = PointGauges(gauges=((('u0',), ((0, 0, 0), (1, 1, 1))),),
                    fileName=filename)
    time_list=[0.0, 1.0, 2.0]
    run_gauge(p, time_list)

    correct_gauge_names = ['u0 [        0         0         0]', 'u0 [        1         1         1]']
    correct_data = np.asarray([[   0.,    0.,  111.],
                               [   1.,    0.,  222.],
                               [   2.,    0.,  333.]])

    # synchronize processes before attempting to read file

    from proteus import Comm
    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)

    eq_(correct_gauge_names, gauge_names)
    npt.assert_equal(correct_data, data)

if __name__ == '__main__':
    test_gauge_output()