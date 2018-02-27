from proteus import Comm, Profiling
from collections import namedtuple
import numpy.testing as npt
import numpy as np
from nose.tools import eq_, ok_
import os
import math
import pytest
import xml.etree.ElementTree as ElementTree
from proteus.EGeometry import (EVec,
                               ETen)


from proteus.MeshTools import (Node,
                               Element,
                               Edge,
                               Polygon,
                               getNodesFromEdges,
                               getEdgesFromPolygons,
                               Triangle,
                               Quadrilateral,
                               Tetrahedron,
                               Hexahedron,
                               PointMesh,
                               EdgeGrid,
                               RectangularGrid,
                               QuadrilateralGrid,
                               EdgeMesh,
                               TriangularMesh,
                               TetrahedralMesh,
                               MeshParallelPartitioningTypes,
                               triangleVerticesToNormals,
                               tetrahedronVerticesToNormals,
                               intersectPoints,
                               intersectEdges,
                               intersectPolyhedron,
                               getMeshIntersections,
                               MultilevelEdgeMesh,
                               MultilevelTriangularMesh,
                               MultilevelTetrahedralMesh,
                               MultilevelHexahedralMesh,
                               InterpolatedBathymetryMesh,)

comm = Comm.init()
Profiling.procID = comm.rank()

GNUPLOT=False

@pytest.mark.MeshTools
class TestMeshTools():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        pass

    def teardown_method(self,method):
        filenames = []
        filenames.extend(['tetgen'+'.'+post for post in ['face','ele','node']])
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError, e:
                    print("Error: %s - %s." %(e.filename,e.strerror))
            else:
                    pass

    def test_Node(self):
        origin_default = Node()
        origin = Node(nodeNumber=0,x=0.0,y=0.0,z=0.0)
        ppp = Node(nodeNumber=1, x= 0.5, y= 0.5, z= 0.5)
        ppm = Node(nodeNumber=2, x= 0.5, y= 0.5, z=-0.5)
        pmp = Node(nodeNumber=3, x= 0.5, y=-0.5, z= 0.5)
        pmm = Node(nodeNumber=4, x= 0.5, y=-0.5, z=-0.5)
        mpp = Node(nodeNumber=5, x=-0.5, y= 0.5, z= 0.5)
        mpm = Node(nodeNumber=6, x=-0.5, y= 0.5, z=-0.5)
        mmp = Node(nodeNumber=7, x=-0.5, y=-0.5, z= 0.5)
        mmm = Node(nodeNumber=8, x=-0.5, y=-0.5, z=-0.5)
        lexicographic_ordering = [ppp, ppm, pmp, pmm, mpp, mpm, mmp, mmm]
        ok_(origin.N == 0)
        ok_(origin.p[0] == 0.0)
        ok_(origin.p[1] == 0.0)
        ok_(origin.p[2] == 0.0)
        ok_(origin.length == 1.0)
        ok_(origin.diameter == 1.0)
        ok_((origin.unitNormal == Node.xUnitVector).all())
        ok_(origin == origin_default)
        ok_(str(origin) == str(origin_default))
        ok_(hash(origin) == hash(origin_default))
        ok_(ppp > origin)
        for  i in range(8):
            p = lexicographic_ordering[i]
            for pg in lexicographic_ordering[:i]:
                ok_(pg > p)
                ok_(pg >= p)
                ok_(pg != p)
            for pl in lexicographic_ordering[i+1:]:
                ok_(pl < p)
                ok_(pl <= p)
                ok_(pl != p)
        v = EVec(0.25,0.35,0.45)
        ntest = Node(29,0.0,0.0,0.0)
        ntest.p +=v
        ok_((ntest.p == (0.25, 0.35, 0.45)).all())

    def test_Element(self):
        nodesCube= [Node(0,0.0,0.0,0.0),
                    Node(1,0.0,1.0,0.0),
                    Node(2,1.0,1.0,0.0),
                    Node(3,1.0,0.0,0.0),
                    Node(4,0.0,0.0,1.0),
                    Node(5,0.0,1.0,1.0),
                    Node(6,1.0,1.0,1.0),
                    Node(7,1.0,0.0,1.0)]
        e0 = Element(elementNumber=1, nodes = nodesCube)
        nodesCube.sort()
        ok_(e0.N == 1)
        for N, n   in enumerate(nodesCube):
            ok_(e0.nodes[N] == n)

    def test_Edge(self):
        edge0 = Edge(0,nodes=[Node(0,0.0,0.0,0.0),
                              Node(1,1.0,1.0,1.0)])
        edge1 = Edge(1,nodes=[Node(1,1.0,1.0,1.0),
                              Node(2,2.0,1.0,1.0)])
        edgesOrdered = [edge0, edge1]
        edgesDisordered = {edge1.nodes: edge1, edge0.nodes: edge0}
        edgeKeys = edgesDisordered.keys()
        edgeKeys.sort()
        for e1,e2_key in zip(edgesOrdered, edgeKeys):
            ok_(e1.nodes == edgesDisordered[e2_key].nodes)
        edge0.computeGeometricInfo()
        ok_((edge0.barycenter == EVec(0.5,0.5,0.5)).all())
        ok_(edge0.length == math.sqrt(3.0))
        ok_(edge0.diameter == math.sqrt(3.0))
        ok_(edge0.innerDiameter == math.sqrt(3.0))
        ok_(edge0.hasGeometricInfo)
        nodeSet = set(edge0.nodes + edge1.nodes)
        ok_(nodeSet == set(getNodesFromEdges([edge0,edge1])))

    def test_Polygon(self):
        nodeList = [Node(nodeNumber=1, x= 0.5, y= 0.5, z= 0.5),
                    Node(nodeNumber=2, x= 0.5, y= 0.5, z=-0.5),
                    Node(nodeNumber=3, x= 0.5, y=-0.5, z= 0.5)]
        polygon = Polygon(12,nodes=nodeList)
        ok_(polygon.N == 12)
        ok_(len(polygon.nodes) == 3)
        nodeList.sort()
        ok_(tuple(nodeList) == polygon.nodes)
        polygon.edges = [Edge(0,[nodeList[0], nodeList[1]]),
                         Edge(1,[nodeList[1], nodeList[2]]),
                         Edge(2,[nodeList[2], nodeList[0]])]
        edges = getEdgesFromPolygons([polygon])
        ok_(len(edges) == 3)
        polygon.edges = [Edge(0,[nodeList[0], nodeList[1]]),
                         Edge(1,[nodeList[1], nodeList[2]]),
                         Edge(3,[nodeList[1], nodeList[2]]),
                         Edge(2,[nodeList[2], nodeList[0]])]
        edges = getEdgesFromPolygons([polygon])
        ok_(len(edges) == 3)

    def test_Triangle(self):
        t0 = Triangle(0,nodes=[Node(0, 0.0, 0.0, 0.0),
                               Node(1, 0.0, 1.0, 0.0),
                               Node(2, 1.0, 0.0, 0.0)])
        t1 = Triangle(0,nodes=[Node(3, 1.0, 1.0, 0.0),
                               Node(1, 0.0, 1.0, 0.0),
                               Node(2, 1.0, 0.0, 0.0)])
        ok_(not t0.hasGeometricInfo)
        ok_(t0.nodes < t1.nodes)
        t0.computeGeometricInfo()
        ok_((t0.barycenter == EVec(1.0/3.0,1.0/3.0,0.0)).all())
        #needs more

    def test_Quadrilateral(self):
        nodes=[Node(0, 0.0, 0.0, 0.0),
               Node(1, 0.0, 1.0, 0.0),
               Node(2, 1.0, 1.0, 0.0),
               Node(3, 1.0, 0.0, 0.0),
               Node(4, 2.0, 1.0, 0.0),
               Node(5, 2.0, 0.0, 0.0)]
        edges0 = [Edge(0,[nodes[0],nodes[1]]),
                  Edge(1,[nodes[1],nodes[2]]),
                  Edge(2,[nodes[2],nodes[3]]),
                  Edge(3,[nodes[3],nodes[0]])]
        edges1 = [Edge(2,[nodes[3],nodes[2]]),
                  Edge(4,[nodes[2],nodes[4]]),
                  Edge(5,[nodes[4],nodes[5]]),
                  Edge(6,[nodes[5],nodes[3]])]
        q0 = Quadrilateral(0,edges0)
        q1 = Quadrilateral(1,edges1)
        ok_(q0.nodes < q1.nodes)
        ok_(q0.nodes != q1.nodes)
        ok_( not q0.nodes > q1.nodes)
        ok_(not q0.hasGeometricInfo)
        q0.computeGeometricInfo()
        ok_(q0.area == 1.0)

    def test_Tetrahedron(self):
        nodes0=[Node(0, 0.0, 0.0, 0.0),
               Node(1, 0.0, 1.0, 0.0),
               Node(2, 1.0, 1.0, 0.0),
               Node(3, 0.0, 0.0, 1.0)]
        nodes1=nodes0[1:] + [Node(4, 1.0, 0.0, 1.0)]
        T0 = Tetrahedron(0,nodes0)
        T1 = Tetrahedron(1,nodes1)
        ok_(T0.nodes < T1.nodes)
        ok_(T0.nodes != T1.nodes)
        ok_(not T0.nodes > T1.nodes)
        ok_(not T0.hasGeometricInfo)
        T0.computeGeometricInfo()
        ok_(T0.volume == 1.0/6.0)
        triangleDict={}
        for t in T0.triangles:
            triangleDict[t.nodes] = t
        edgeDict={}
        for e in T0.edges:
            edgeDict[e.nodes] = e
        T0_1 = Tetrahedron(0,nodes0, edgeDict=edgeDict)
        T0_2 = Tetrahedron(0,nodes0, triangleDict=triangleDict)
        T0_3 = Tetrahedron(0,nodes0, edgeDict=edgeDict, triangleDict=triangleDict)
        ok_(T0.nodes == T0_1.nodes == T0_2.nodes)

    def test_Hexahedron(self):
        hexGrid = RectangularGrid(3,2,2,
                                  2.0,1.0,1.0)
        H0 = hexGrid.hexahedronList[0]
        H1 = hexGrid.hexahedronList[1]
        ok_(H0.nodes < H1.nodes)
        ok_(H0.nodes != H1.nodes)
        ok_(not H0.nodes > H1.nodes)
        ok_(not H0.hasGeometricInfo)

    def test_MeshParallelPartitioningTypes(self):
        ok_(MeshParallelPartitioningTypes.element == 0)
        ok_(MeshParallelPartitioningTypes.node == 1)
        ok_(MeshParallelPartitioningTypes.node !=
            MeshParallelPartitioningTypes.element)

    def test_intersect_points(self):
        npt.assert_equal(intersectPoints(([0, 1], [0, 2]), [[0, 1], [0., 1.5]]),
                                              [([0, 1]), ([0., 1.5])])

        npt.assert_equal(intersectPoints(([0, 1], [0, 1]), [[0, 1], [1., 2.]]),
                                              [([0, 1]), None])

    def test_intersect_edges(self):
        # check that overlaps grab furthest point on edge
        npt.assert_equal(intersectEdges([(0, 1), (0, 2)], [[(0, 1), (0, 2)],]),
                         [[0, 2]])

        # check 3D case
        npt.assert_equal(intersectEdges([(3, 1, 1), (3, 1, 2)], [[(3, 1, 1), (3, 1, 4)],],),
                         [[3, 1, 2]])

        # check a proper 3D intersection
        npt.assert_equal(intersectEdges([(1, 1, 0), (1, 1, 2)], [[(5, 5, 5), (-1, -1, -1)],]),
                         [[1, 1, 1]])

    def test_intersect_polyhedron(self):
        # check 2-D triangle intersection
        l = ((-1, 0.5), (1, 0.5))
        p = (((-1, 0), (0, 0)), ((0, -1), (0, 0)), ((1, 1), (0.5, 0.5)))
        npt.assert_equal(intersectPolyhedron(l, p),
                         ((0., 0.5), ((0.5, 0.5))))

        # check 2-D square intersection in 3 space
        l = ((0, 0, 0), (2, 2, 0))
        p = (((-1, 0, 0), (0, 0, 0)),
             ((0, -1, 0), (0, 0, 0)),
             ((1, 0, 0), (1, 1, 0)),
             ((0, 1, 0), (1, 1, 0)))
        npt.assert_equal(intersectPolyhedron(l, p),
                         ((0., 0., 0.), ((1., 1., 0.))))


    def test_triangle_vertices_to_normals(self):
        npt.assert_equal(triangleVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0])),
                         [(([-1.,  0.,  0.]), ([0, 0, 0])),
                          (([-0., -1., -0.]), ([0, 0, 0])),
                          (([ 1.,  1.,  0.]), ([0, 1, 0]))])

    def test_tetrahedron_vertices_to_normals(self):
        npt.assert_equal(tetrahedronVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1])),
                                                           [(([ 0,  0, -1]), ([0, 0, 0])),
                                                            (([-1,  0,  0]), ([0, 0, 0])),
                                                            (([ 0, -1,  0]), ([0, 0, 0])),
                                                            (([1, 1, 1]), ([0, 1, 0]))])

    def test_mesh_intersections(self):
        """
        Test correctness of intersecting a line and a tetrahedral mesh.
        """

        endpoints = np.asarray(((0, 0, 0), (1, 1, 1)))
        mesh_type = namedtuple('mesh', 'nodeArray elementNodesArray')
        nodeArray = np.asarray(([[ 0. ,  0. ,  0. ],
           [ 0.5,  0. ,  0. ],
           [ 1. ,  0. ,  0. ],
           [ 0. ,  0.5,  0. ],
           [ 0.5,  0.5,  0. ],
           [ 1. ,  0.5,  0. ],
           [ 0. ,  1. ,  0. ],
           [ 0.5,  1. ,  0. ],
           [ 1. ,  1. ,  0. ],
           [ 0. ,  0. ,  0.5],
           [ 0.5,  0. ,  0.5],
           [ 1. ,  0. ,  0.5],
           [ 0. ,  0.5,  0.5],
           [ 0.5,  0.5,  0.5],
           [ 1. ,  0.5,  0.5],
           [ 0. ,  1. ,  0.5],
           [ 0.5,  1. ,  0.5],
           [ 1. ,  1. ,  0.5],
           [ 0. ,  0. ,  1. ],
           [ 0.5,  0. ,  1. ],
           [ 1. ,  0. ,  1. ],
           [ 0. ,  0.5,  1. ],
           [ 0.5,  0.5,  1. ],
           [ 1. ,  0.5,  1. ],
           [ 0. ,  1. ,  1. ],
           [ 0.5,  1. ,  1. ],
           [ 1. ,  1. ,  1. ]]))

        elementNodesArray = np.asarray([[ 0,  1,  4, 13],
           [ 0,  1, 10, 13],
           [ 0,  3,  4, 13],
           [ 0,  3, 12, 13],
           [ 0,  9, 10, 13],
           [ 0,  9, 12, 13],
           [ 1,  2,  5, 14],
           [ 1,  2, 11, 14],
           [ 1,  4,  5, 14],
           [ 1,  4, 13, 14],
           [ 1, 10, 11, 14],
           [ 1, 10, 13, 14],
           [ 3,  4,  7, 16],
           [ 3,  4, 13, 16],
           [ 3,  6,  7, 16],
           [ 3,  6, 15, 16],
           [ 3, 12, 13, 16],
           [ 3, 12, 15, 16],
           [ 4,  5,  8, 17],
           [ 4,  5, 14, 17],
           [ 4,  7,  8, 17],
           [ 4,  7, 16, 17],
           [ 4, 13, 14, 17],
           [ 4, 13, 16, 17],
           [ 9, 10, 13, 22],
           [ 9, 10, 19, 22],
           [ 9, 12, 13, 22],
           [ 9, 12, 21, 22],
           [ 9, 18, 19, 22],
           [ 9, 18, 21, 22],
           [10, 11, 14, 23],
           [10, 11, 20, 23],
           [10, 13, 14, 23],
           [10, 13, 22, 23],
           [10, 19, 20, 23],
           [10, 19, 22, 23],
           [12, 13, 16, 25],
           [12, 13, 22, 25],
           [12, 15, 16, 25],
           [12, 15, 24, 25],
           [12, 21, 22, 25],
           [12, 21, 24, 25],
           [13, 14, 17, 26],
           [13, 14, 23, 26],
           [13, 16, 17, 26],
           [13, 16, 25, 26],
           [13, 22, 23, 26],
           [13, 22, 25, 26]])

        mesh = mesh_type(nodeArray=nodeArray, elementNodesArray=elementNodesArray)
        toPolyhedron = tetrahedronVerticesToNormals
        eq_(getMeshIntersections(mesh, toPolyhedron, endpoints),
            set([((0.5, 0.5, 0.5), (1, 1, 1)), ((0, 0, 0), (0.5, 0.5, 0.5))]))


    def test_PointMesh(self):
        points = np.array([[0.0,0.0,0.0],[1.0,0.0,0.0]])
        pg=PointMesh(points)
        npt.assert_almost_equal(pg.nodeArray, np.array([[ 0., 0., 0.],
                                                        [ 1., 0., 0.]]))
        ok_(pg.nNodes_global == 2)
        ok_((pg.elementNodesArray == np.array([0,1])).all())
        ok_(pg.nElements_global == 2)

    def test_EdgeGrid(self):
         eg=EdgeGrid(nx=3,Lx=1.0)
         if GNUPLOT:
             eg.writeEdgesGnuplot2('edgegrid')
             eg.viewMeshGnuplotPipe('edgegrid')
         npt.assert_almost_equal(eg.nodeArray,np.array([[ 0.,   0.,   0. ],
                                                       [ 0.5,  0.,   0. ],
                                                       [ 1.,   0.,   0. ]]))
         ok_(eg.nNodes_global == 3)
         ok_((eg.elementNodesArray == np.array([[0, 1],
                                                [1, 2]])).all())
         ok_(eg.nElements_global == 2)
         ok_((eg.elementBoundariesArray == [[ 0.,   0.,   0. ],
                                            [ 0.5,  0.,   0. ],
                                            [ 1.,   0.,   0. ]]).all())

    def test_QuadrilateralGrid(self):
         qg=QuadrilateralGrid(nx=3,ny=3,Lx=1.0,Ly=1.0)
         if GNUPLOT:
             qg.writeEdgesGnuplot2('quadrilateralgrid')
             qg.viewMeshGnuplotPipe('quadrilateralgrid')
         npt.assert_almost_equal(qg.nodeArray, np.array([[ 0.,   0.,   0., ],
                                                         [ 0.,   0.5,  0., ],
                                                         [ 0.,   1.,   0., ],
                                                         [ 0.5,  0.,   0., ],
                                                         [ 0.5,  0.5,  0., ],
                                                         [ 0.5,  1.,   0., ],
                                                         [ 1.,   0.,   0., ],
                                                         [ 1.,   0.5,  0., ],
                                                         [ 1.,   1.,   0., ]]))
         ok_(qg.nNodes_global == 9)
         ok_((qg.elementNodesArray == np.array([[0, 1, 3, 4],
                                                [1, 2, 4, 5],
                                                [3, 4, 6, 7],
                                                [4, 5, 7, 8]])).all())
         ok_(qg.nElements_global == 4)
         edges = np.array([[0, 1],
                           [1, 2],
                           [0, 3],
                           [1, 4],
                           [2, 5],
                           [3, 4],
                           [4, 5],
                           [3, 6],
                           [4, 7],
                           [5, 8],
                           [6, 7],
                           [7, 8]])
         ok_((qg.elementBoundariesArray == edges).all())
         ok_(qg.nElementBoundaries_global == 12)
         ok_((qg.edgeNodesArray == edges).all())
         ok_(qg.nEdges_global == 12)

    def test_RectangularGrid_1D(self):
         grid1d = RectangularGrid(3,1,1,1.0,1.0,1.0)
         if GNUPLOT:
             grid1d.writeEdgesGnuplot('grid1d')
             grid1d.viewMeshGnuplotPipe('grid1d')

    def test_RectangularGrid_2D(self):
        grid2d = RectangularGrid(3,3,1,1.0,1.0,1.0)
        if GNUPLOT:
            grid2d.writeEdgesGnuplot('grid2d')
            grid2d.viewMeshGnuplotPipe('grid2d')

    def test_RectangularGrid_3D(self):
        grid3d = RectangularGrid(3,3,3,1.0,1.0,1.0)
        if GNUPLOT:
            grid3d.writeEdgesGnuplot('grid3d')
            grid3d.viewMeshGnuplotPipe('grid3d')

    def test_EdgeMesh_1D(self):
        mesh1d = EdgeMesh()
        mesh1d.generateEdgeMeshFromRectangularGrid(3,1.0)
        if GNUPLOT:
           mesh1d.buildLists()
           mesh1d.writeEdgesGnuplot('mesh1d')
           mesh1d.viewMeshGnuplotPipe('mesh1d')

    def test_TriangularMesh(self):
        mesh2d = TriangularMesh()
        mesh2d.generateTriangularMeshFromRectangularGrid(3,3,1.0,1.0)
        if GNUPLOT:
           mesh2d.buildLists()
           mesh2d.writeEdgesGnuplot('mesh2d')
           mesh2d.viewMeshGnuplotPipe('mesh2d')

    def test_TetrahedralMesh(self):
        mesh3d = TetrahedralMesh()
        mesh3d.generateTetrahedralMeshFromRectangularGrid(3,3,3,1.0,1.0,1.0)
        if GNUPLOT:
           mesh3d.buildLists()
           mesh3d.writeEdgesGnuplot('mesh3d')
           mesh3d.viewMeshGnuplotPipe('mesh3d')

    def test_Refine_1D(self):
        grid1d = RectangularGrid(3,1,1,1.0,1.0,1.0)
        grid1dFine = RectangularGrid()
        children = grid1dFine.refine(grid1d,2)
        parent=0
        child=0
        for pN,cL in children.iteritems():
            ok_(parent == pN)
            for c in cL:
                ok_(child == c.N)
                child +=1
            parent +=1
        if GNUPLOT:
            grid1dFine.writeEdgesGnuplot('grid1dFine')
            grid1dFine.viewMeshGnuplotPipe('grid1dFine')

    def test_Refine_2D(self):
        grid2d = RectangularGrid(3,3,1,1.0,1.0,1.0)
        grid2dFine = RectangularGrid()
        children = grid2dFine.refine(grid2d,2,2)
        childElements = [[0, 4, 1, 5],
                         [2, 6, 3, 7],
                         [8, 12, 9, 13],
                         [10, 14, 11, 15]]
        parent = 0
        for pN,cL in children.iteritems():
            ok_(parent == pN)
            for ci,c in enumerate(cL):
                ok_(childElements[pN][ci] == c.N)
            parent += 1

    def test_Refine_3D(self):
        grid3d = RectangularGrid(3,3,3,1.0,1.0,1.0)
        grid3dFine = RectangularGrid()
        children = grid3dFine.refine(grid3d,2,2,2)
        childElements = [[0 ,16,4 ,20,1 ,17,5 ,21],
                         [2 ,18,6 ,22,3 ,19,7 ,23],
                         [8 ,24,12,28,9 ,25,13,29],
                         [10,26,14,30,11,27,15,31],
                         [32,48,36,52,33,49,37,53],
                         [34,50,38,54,35,51,39,55],
                         [40,56,44,60,41,57,45,61],
                         [42,58,46,62,43,59,47,63]]
        parent = 0
        for pN,cL in children.iteritems():
            ok_(parent == pN)
            for ci,c in enumerate(cL):
                ok_(childElements[pN][ci] == c.N)
            parent += 1

    def test_MultilevelEdgeMesh(self):
        n = 3
        for ptype in [MeshParallelPartitioningTypes.node,
                       MeshParallelPartitioningTypes.element]:
            mlMesh = MultilevelEdgeMesh(3,1,1,refinementLevels=n,parallelPartitioningType=ptype)
            elementChildren= [np.array([0, 1, 2, 3]),
                              np.array([0, 1, 2, 3, 4, 5, 6, 7])]
            elementChildrenOffsets= [np.array([0, 2, 4]),
                                     np.array([0, 2, 4, 6, 8])]
            elementParents= [np.array([]),
                             np.array([0, 0, 1, 1]),
                             np.array([0, 0, 1, 1, 2, 2, 3, 3])]
            for l in range(n):
                if l < n-1:
                    ok_((elementChildren[l] == mlMesh.elementChildrenArrayList[l]).all())
                    ok_((elementChildrenOffsets[l] == mlMesh.elementChildrenOffsetsList[l]).all())
                ok_((elementParents[l] == mlMesh.elementParentsArrayList[l]).all())

    def test_MultilevelTriangularMesh(self):
        n = 3
        for ptype in [MeshParallelPartitioningTypes.element,
                      MeshParallelPartitioningTypes.node]:
            mlMesh = MultilevelTriangularMesh(3,3,1,refinementLevels=n,parallelPartitioningType=ptype)
            elementChildren= [np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
                                         17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], dtype=np.int32), np.array([  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
                                                                                                                                   13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
                                                                                                                                   26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
                                                                                                                                   39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
                                                                                                                                   52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
                                                                                                                                   65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
                                                                                                                                   78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
                                                                                                                                   91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
                                                                                                                                   104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                                                                                                                                   117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127], dtype=np.int32)]
            elementChildrenOffsets= [np.array([ 0,  4,  8, 12, 16, 20, 24, 28, 32], dtype=np.int32), np.array([  0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,
                                                                                                                 52,  56,  60,  64,  68,  72,  76,  80,  84,  88,  92,  96, 100,
                                                                                                                 104, 108, 112, 116, 120, 124, 128], dtype=np.int32)]
            elementParents= [np.array([], dtype=np.int32), np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5,
                                                                     5, 6, 6, 6, 6, 7, 7, 7, 7], dtype=np.int32), np.array([ 0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,
                                                                                                                             4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,
                                                                                                                             8,  8,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12,
                                                                                                                             12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16,
                                                                                                                             17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21,
                                                                                                                             21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 25, 25,
                                                                                                                             25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29,
                                                                                                                             29, 30, 30, 30, 30, 31, 31, 31, 31], dtype=np.int32)]
            for l in range(n):
                if l < n-1:
                    #pass
                    ok_((elementChildren[l] == mlMesh.elementChildrenArrayList[l]).all())
                    ok_((elementChildrenOffsets[l] == mlMesh.elementChildrenOffsetsList[l]).all())
                ok_((elementParents[l] == mlMesh.elementParentsArrayList[l]).all())
            mlMesh = MultilevelTriangularMesh(3,3,1,refinementLevels=n,triangleFlag=1)
            mlMesh = MultilevelTriangularMesh(3,3,1,refinementLevels=n,triangleFlag=2)
            mlMesh2 = MultilevelTriangularMesh(0,0,0,skipInit=True)
            mlMesh2.generateFromExistingCoarseMesh(mlMesh.meshList[0],
                                                   refinementLevels=n)


    def test_MultilevelTetrahedralMesh(self):
        n = 2
        for ptype in [MeshParallelPartitioningTypes.element,
                      MeshParallelPartitioningTypes.node]:
            mlMesh = MultilevelTetrahedralMesh(3,3,3,refinementLevels=n,parallelPartitioningType=ptype)
            elementChildren= [np.array([  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
                                          13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
                                          26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
                                          39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
                                          52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
                                          65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
                                          78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
                                          91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
                                          104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                                          117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
                                          130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
                                          143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                                          156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
                                          169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
                                          182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
                                          195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
                                          208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
                                          221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
                                          234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
                                          247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,
                                          260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
                                          273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285,
                                          286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298,
                                          299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311,
                                          312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
                                          325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,
                                          338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350,
                                          351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363,
                                          364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376,
                                          377, 378, 379, 380, 381, 382, 383], dtype=np.int32)]
            elementChildrenOffsets= [np.array([  0,   8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88,  96,
                                                 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 200,
                                                 208, 216, 224, 232, 240, 248, 256, 264, 272, 280, 288, 296, 304,
                                                 312, 320, 328, 336, 344, 352, 360, 368, 376, 384], dtype=np.int32)]
            elementParents= [np.array([], dtype=np.int32), np.array([ 0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2,
                                                                      2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4,  4,
                                                                      4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,
                                                                      6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,
                                                                      8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,
                                                                      10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12,
                                                                      12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14,
                                                                      14, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16,
                                                                      17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19,
                                                                      19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21,
                                                                      21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23,
                                                                      23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25,
                                                                      25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27,
                                                                      27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29,
                                                                      29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31,
                                                                      31, 32, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33,
                                                                      34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 36,
                                                                      36, 36, 36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38,
                                                                      38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 40, 40, 40,
                                                                      40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42,
                                                                      42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44,
                                                                      44, 44, 44, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46,
                                                                      46, 46, 47, 47, 47, 47, 47, 47, 47, 47], dtype=np.int32)]
            for l in range(n):
                if l < n-1:
                    ok_((elementChildren[l] == mlMesh.elementChildrenArrayList[l]).all())
                    ok_((elementChildrenOffsets[l] == mlMesh.elementChildrenOffsetsList[l]).all())
                ok_((elementParents[l] == mlMesh.elementParentsArrayList[l]).all())
            mlMesh2 = MultilevelTetrahedralMesh(0,0,0,skipInit=True)
            mlMesh2.generateFromExistingCoarseMesh(mlMesh.meshList[0],
                                                   refinementLevels=n)

    def test_MultilevelHexahedralMesh(self):
        n = 1
        mlMesh = MultilevelHexahedralMesh(4,4,4,
                                          refinementLevels=n)


if __name__ == '__main__':
    pass
