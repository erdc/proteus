from math import ceil, sqrt, pow
import os
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
from importlib import reload
reload(p)
reload(n)

from proteus.Gauges import PointGauges, LineGauges, LineIntegralGauges

from proteus.test_utils.util import setup_profiling, silent_rm

def build1DMesh(p, nnx):
    return MeshTools.MultilevelEdgeMesh(nnx, 1, 1,
                                        p.domain.x[0], 0.0, 0.0,
                                        p.domain.L[0], 1.0, 1.0,
                                        refinementLevels=1,
                                        nLayersOfOverlap=0,
                                        parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)
def build2DMesh(p, nnx, nny):
    return MeshTools.MultilevelTriangularMesh(nnx,nny,1,
                                              p.domain.x[0], p.domain.x[1], 1.0,
                                              p.domain.L[0], p.domain.L[1], 1.0,
                                              refinementLevels=1,
                                              nLayersOfOverlap=0,
                                              parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)

def build3DMesh(p, nnx, nny, nnz):
    return MeshTools.MultilevelTetrahedralMesh(nnx,nny,nnz,
                                               p.domain.x[0], p.domain.x[1], p.domain.x[2],
                                               p.domain.L[0], p.domain.L[1], p.domain.L[2],
                                               refinementLevels=1,
                                               nLayersOfOverlap=0,
                                               parallelPartitioningType=MeshTools.MeshParallelPartitioningTypes.node)

def gauge_setup(nd, total_nodes=None):
    comm = Comm.get()

    #Simplified Physics
    p.name="test_gauges"

    p.nd = nd

    class LinearSolution(object):
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
    n.cfluxtag = None
    n.conservativeFlux = None

    if total_nodes is None:
        total_nodes = 2*comm.size()

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

def run_gauge(p, time_list, nd=3, total_nodes=None):
    model, initialConditions = gauge_setup(nd, total_nodes)

    p.attachModel(model, None)

    m = model.levelModelList[-1]
    m.setInitialConditions(initialConditions, time_list[0])
    # behave like a time-stepping model, set tLast to this time
    m.timeIntegration.tLast = time_list[0]
    tCount = 0
    p.calculate()

    for t in time_list[1:]:
        tCount +=1
        m.setInitialConditions(initialConditions, t)
        m.timeIntegration.tLast = t
        p.calculate()

def parse_gauge_output(filename):
    with open(filename) as f:
        header = f.readline().split(',')
        assert header[0].strip() == 'time'
        gauge_names = [gauge_name.strip() for gauge_name in header[1:]]
        f.seek(0)
        data = np.genfromtxt(f, delimiter=",", skip_header=1)
    return gauge_names, data

def test_2D_point_gauge_output():
    filename = 'test_2D_gauge_output.csv'
    silent_rm(filename)

    p = PointGauges(gauges=((('u0',), ((0, 0, 0), (1, 1, 0))),),
                    fileName=filename)
    time_list=[0.0, 1.0, 2.0]
    run_gauge(p, time_list)

    correct_gauge_names = ['u0 [        0         0         0]', 'u0 [        1         1         0]']
    correct_data = np.asarray([[   0.,    0.,  11.],
                               [   1.,    0.,  22.],
                               [   2.,    0.,  33.]])

    # synchronize processes before attempting to read file

    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)

    assert correct_gauge_names == gauge_names
    npt.assert_allclose(correct_data, data)
    delete_file(filename)

def test_point_gauge_output():
    filename = 'test_gauge_output.csv'
    silent_rm(filename)

    p = PointGauges(gauges=((('u0',), ((0, 0, 0), (1, 1, 1))),),
                    fileName=filename)
    time_list=[0.0, 1.0, 2.0]
    run_gauge(p, time_list)

    correct_gauge_names = ['u0 [        0         0         0]', 'u0 [        1         1         1]']
    correct_data = np.asarray([[   0.,    0.,  111.],
                               [   1.,    0.,  222.],
                               [   2.,    0.,  333.]])

    # synchronize processes before attempting to read file

    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)

    assert correct_gauge_names == gauge_names
    npt.assert_allclose(correct_data, data)
    delete_file(filename)


def test_point_gauge_output_2():
    filename = 'test_gauge_output_2.csv'
    silent_rm(filename)

    p = PointGauges(gauges=((('u0',), ((0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1))),),
                    fileName=filename)
    time_list=[0.0, 1.0, 2.0]
    run_gauge(p, time_list)

    correct_gauge_names = ['u0 [        0         0         0]',
                           'u0 [      0.5       0.5       0.5]',
                           'u0 [        1         1         1]']
    correct_data = np.asarray([[   0.,    0.,  55.5, 111.],
                               [   1.,    0.,  111., 222.],
                               [   2.,    0.,  166.5, 333.]])

    # synchronize processes before attempting to read file
    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)

    assert correct_gauge_names == gauge_names
    npt.assert_allclose(correct_data, data)
    delete_file(filename)

def test_line_integral_gauge_output():
    filename = 'test_line_integral_gauge_output.csv'
    silent_rm(filename)

    lines = (((0, 0, 0), (1, 1, 1)),)
    fields = ('u0', )

    l = LineIntegralGauges(gauges=((fields, lines),),
                   fileName=filename)
    time_list=[0.0, 1.0, 2.0]
    run_gauge(l, time_list)

    correct_gauge_names = ['u0 [        0         0         0] - [        1         1         1]']
    correct_data = np.asarray([[   0., 96.128819820072678],
                               [   1., 192.25763964014536],
                               [   2., 288.38645946021808]])

    # synchronize processes before attempting to read file
    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)
    assert correct_gauge_names == gauge_names
    npt.assert_allclose(correct_data, data)
    delete_file(filename)

def test_2D_line_integral_gauge_output():
    filename = 'test_2D_line_integral_gauge_output.csv'
    silent_rm(filename)

    lines = (((0, 0, 0), (1, 1, 0)),
             ((0, 0, 0), (1, 0, 0)),
             ((0, 0, 0), (0, 1, 0)),
             ((0, 0.5, 0), (1, 0.5, 0)),
             ((0, 1, 0), (1, 1, 0)),
             ((0.5, 0, 0), (0.5, 1, 0)))
    fields = ('u0', )

    l = LineIntegralGauges(gauges=((fields, lines),),
                   fileName=filename)

    time_list=[0.0, 1.0, 2.0, 2.5]
    run_gauge(l, time_list)

    correct_gauge_names = ['u0 [        0         0         0] - [        1         1         0]',
                           'u0 [        0         0         0] - [        1         0         0]',
                           'u0 [        0         0         0] - [        0         1         0]',
                           'u0 [        0       0.5         0] - [        1       0.5         0]',
                           'u0 [        0         1         0] - [        1         1         0]',
                           'u0 [      0.5         0         0] - [      0.5         1         0]']

    correct_data = np.asarray([[0.,   7.77817459,  0.5,   5.,   5.5,  10.5,   5.5],
                               [1.,  15.55634919,  1.,   10.,  11.,   21.,   11.],
                               [2.,  23.33452378,  1.5,  15.,  16.5,  31.5,  16.5],
                               [2.5, 27.22361108,  1.75, 17.5, 19.25, 36.75, 19.25]])

    # synchronize processes before attempting to read file

    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)
    assert correct_gauge_names == gauge_names

    npt.assert_allclose(correct_data, data)
    delete_file(filename)

def test_line_gauge_output():
    filename = 'test_line_output.csv'
    silent_rm(filename)

    lines = (((0, 0, 0), (1, 1, 1)),)
    fields = ('u0', )

    l = LineGauges(gauges=((fields, lines),),
                   fileName=filename)

    time_list=[0.0, 1.0, 2.0]
    # good hard-code for up to 1,024 processes or so, but slow (about 8 seconds on single process).
    run_gauge(l, time_list, total_nodes=2048)
    correct_gauge_names = ['u0 [        0         0         0]', 'u0 [ 0.076923  0.076923  0.076923]',
                           'u0 [  0.15385   0.15385   0.15385]', 'u0 [  0.23077   0.23077   0.23077]',
                           'u0 [  0.30769   0.30769   0.30769]', 'u0 [  0.38462   0.38462   0.38462]',
                           'u0 [  0.46154   0.46154   0.46154]', 'u0 [  0.53846   0.53846   0.53846]',
                           'u0 [  0.61538   0.61538   0.61538]', 'u0 [  0.69231   0.69231   0.69231]',
                           'u0 [  0.76923   0.76923   0.76923]', 'u0 [  0.84615   0.84615   0.84615]',
                           'u0 [  0.92308   0.92308   0.92308]', 'u0 [        1         1         1]']
    correct_data = np.asarray([[0., 0., 8.53846154, 17.07692308, 25.61538462,
                                34.15384615, 42.69230769, 51.23076923, 59.76923077, 68.30769231,
                                76.84615385, 85.38461538, 93.92307692, 102.46153846, 111.],
                               [1., 0., 17.07692308, 34.15384615, 51.23076923,
                                68.30769231, 85.38461538, 102.46153846, 119.53846154, 136.61538462,
                                153.69230769, 170.76923077, 187.84615385, 204.92307692, 222.],
                               [2., 0., 25.61538462, 51.23076923, 76.84615385,
                                102.46153846, 128.07692308, 153.69230769, 179.30769231, 204.92307692,
                                230.53846154, 256.15384615, 281.76923077, 307.38461538, 333.]])

    # synchronize processes before attempting to read file

    Comm.get().barrier()

    gauge_names, data = parse_gauge_output(filename)

    assert correct_gauge_names == gauge_names
    npt.assert_allclose(correct_data, data)
    delete_file(filename)

def delete_file(filename):
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except OSError as e:
            print ("Error: %s - %s" %(e.filename,e.strerror))
        else:
            pass
    

if __name__ == '__main__':
    setup_profiling()
    test_line_gauge_output()
    test_point_gauge_output()
    test_point_gauge_output_2()
    test_line_integral_gauge_output()
