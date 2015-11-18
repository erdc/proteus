from collections import namedtuple
import numpy.testing as npt
import numpy as np
from nose.tools import eq_, ok_
import math
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
                               getMeshIntersections)

GNUPLOT=False

def test_Node():
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

def test_Element():
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

def test_Edge():
    edge0 = Edge(0,nodes=[Node(0,0.0,0.0,0.0),Node(1,1.0,1.0,1.0)])
    edge1 = Edge(0,nodes=[Node(0,1.0,1.0,1.0),Node(1,2.0,1.0,1.0)])
    edgesOrdered = [edge0, edge1]
    sortTest = [edge1, edge0]
    sortTest.sort()
    for eO,eS in zip(edgesOrdered, sortTest):
        ok_(eO == eS)
    edge0.computeGeometricInfo()
    ok_((edge0.barycenter == EVec(0.5,0.5,0.5)).all())
    ok_(edge0.length == math.sqrt(3.0))
    ok_(edge0.diameter == math.sqrt(3.0))
    ok_(edge0.innerDiameter == math.sqrt(3.0))
    ok_(edge0.hasGeometricInfo)
    nodeSet = set(edge0.nodes + edge1.nodes)
    ok_(nodeSet == set(getNodesFromEdges([edge0,edge1])))

def test_Polygon():
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

def test_Triangle():
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

def test_Quadrilateral():
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
    #todo  more tests...

#def test_Tetrahedron():
#def test_Hexahedron():

def test_MeshParallelPartitioningTypes():
    ok_(MeshParallelPartitioningTypes.element == 0)
    ok_(MeshParallelPartitioningTypes.node == 1)
    ok_(MeshParallelPartitioningTypes.node !=
        MeshParallelPartitioningTypes.element)

def test_intersect_points():
    npt.assert_equal(intersectPoints(([0, 1], [0, 2]), [[0, 1], [0., 1.5]]),
                                          [([0, 1]), ([0., 1.5])])

    npt.assert_equal(intersectPoints(([0, 1], [0, 1]), [[0, 1], [1., 2.]]),
                                          [([0, 1]), None])

def test_intersect_edges():
    # check that overlaps grab furthest point on edge
    npt.assert_equal(intersectEdges([(0, 1), (0, 2)], [[(0, 1), (0, 2)],]),
                     [[0, 2]])

    # check 3D case
    npt.assert_equal(intersectEdges([(3, 1, 1), (3, 1, 2)], [[(3, 1, 1), (3, 1, 4)],],),
                     [[3, 1, 2]])

    # check a proper 3D intersection
    npt.assert_equal(intersectEdges([(1, 1, 0), (1, 1, 2)], [[(5, 5, 5), (-1, -1, -1)],]),
                     [[1, 1, 1]])

def test_intersect_polyhedron():
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


def test_triangle_vertices_to_normals():
    npt.assert_equal(triangleVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0])),
                     [(([-1.,  0.,  0.]), ([0, 0, 0])),
                      (([-0., -1., -0.]), ([0, 0, 0])),
                      (([ 1.,  1.,  0.]), ([0, 1, 0]))])

def test_tetrahedron_vertices_to_normals():
    npt.assert_equal(tetrahedronVerticesToNormals(([0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1])),
                                                       [(([ 0,  0, -1]), ([0, 0, 0])),
                                                        (([-1,  0,  0]), ([0, 0, 0])),
                                                        (([ 0, -1,  0]), ([0, 0, 0])),
                                                        (([1, 1, 1]), ([0, 1, 0]))])

def test_mesh_intersections():
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


def test_PointMesh():
    points = np.array([[0.0,0.0,0.0],[1.0,0.0,0.0]])
    pg=PointMesh(points)
    npt.assert_almost_equal(pg.nodeArray, np.array([[ 0., 0., 0.],
                                                    [ 1., 0., 0.]]))
    ok_(pg.nNodes_global == 2)
    ok_((pg.elementNodesArray == np.array([0,1])).all())
    ok_(pg.nElements_global == 2)

def test_EdgeGrid():
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

def test_QuadrilateralGrid():
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

def test_RectangularGrid_1D():
     print "Testing 1D Rectangular Grid"
     grid1d = RectangularGrid(3,1,1,1.0,1.0,1.0)
     if GNUPLOT:
         grid1d.writeEdgesGnuplot('grid1d')
         grid1d.viewMeshGnuplotPipe('grid1d')

def test_RectangularGrid_2D():
    grid2d = RectangularGrid(3,3,1,1.0,1.0,1.0)
    if GNUPLOT:
        grid2d.writeEdgesGnuplot('grid2d')
        grid2d.viewMeshGnuplotPipe('grid2d')

def test_RectangularGrid_3D():
    grid3d = RectangularGrid(3,3,3,1.0,1.0,1.0)
    if GNUPLOT:
        grid3d.writeEdgesGnuplot('grid3d')
        grid3d.viewMeshGnuplotPipe('grid3d')

def test_EdgeMesh_1D():
    mesh1d = EdgeMesh()
    mesh1d.generateEdgeMeshFromRectangularGrid(3,1.0)
    #if GNUPLOT:
    #    mesh1d.buildLists()
    #    mesh1d.writeEdgesGnuplot('mesh1d')
    #    mesh1d.viewMeshGnuplotPipe('mesh1d')

def test_TriangularMesh():
    mesh2d = TriangularMesh()
    mesh2d.generateTriangularMeshFromRectangularGrid(3,3,1.0,1.0)
    #if GNUPLOT:
    #    mesh2d.buildLists()
    #    mesh2d.writeEdgesGnuplot('mesh2d')
    #    mesh2d.viewMeshGnuplotPipe('mesh2d')

def test_TetrahedralMesh():
    mesh3d = TetrahedralMesh()
    mesh3d.generateTetrahedralMeshFromRectangularGrid(3,3,3,1.0,1.0,1.0)
    #if GNUPLOT:
    #    mesh3d.buildLists()
    #    mesh3d.writeEdgesGnuplot('mesh3d')
    #    mesh3d.viewMeshGnuplotPipe('mesh3d')

def test_Refine_1D():
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

#def test_Refine_2D():
    #ok_(True)
    #grid2d = RectangularGrid(3,3,1,1.0,1.0,1.0)
    #grid2dFine = RectangularGrid()
    #children = grid2dFine.refine(grid2d,2,2)
    # childElements = [[0, 4, 1, 5],
    #                  [2, 6, 3, 7],
    #                  [8, 12, 9, 13],
    #                  [10, 14, 11, 15]]
    # parent = 0
    # for pN,cL in children.iteritems():
    #     ok_(parent == pN)
    #     for ci,c in enumerate(cL):
    #         ok_(childElements[pN][ci] == c.N)
    #     parent += 1
    #if GNUPLOT:
    #    grid2dFine.writeEdgesGnuplot('grid2dFine')
    #    grid2dFine.viewMeshGnuplotPipe('grid2dFine')
#      print "Testing 3D Rectangular Grid"
#      grid3dFine = RectangularGrid()
#      children = grid3dFine.refine(grid3d,2,2,2)
#      for pN,cL in children.iteritems():
#          print "Parent Element "+str(pN)
#          print "Child Elements "
#          for c in cL:
#              print str(c.N)
#      grid3dFine.writeEdgesGnuplot('grid3dFine')
#      grid3dFine.viewMeshGnuplotPipe('grid3dFine')
#      print "Testing 1D Edge Mesh Refinement"
#      mesh1dFine = EdgeMesh()
#      children = mesh1dFine.refine(mesh1d)
#      for pN,cL in children.iteritems():
#          print "Parent Element "+str(pN)
#          print "Child Elements "
#          for c in cL:
#              print str(c.N)
#      mesh1dFine.writeEdgesGnuplot('mesh1dFine')
#      mesh1dFine.viewMeshGnuplotPipe('mesh1dFine')
#      print "Testing 2D Triangular Mesh Refinement"
#      mesh2dFine = TriangularMesh()
#      children = mesh2dFine.refine(mesh2d)
#      for pN,cL in children.iteritems():
#          print "Parent Element "+str(pN)
#          print "Child Elements "
#          for c in cL:
#              print str(c.N)
#      mesh2dFine.writeEdgesGnuplot('mesh2dFine')
#      mesh2dFine.viewMeshGnuplotPipe('mesh2dFine')
#      print "Testing 3D Tetrahedral Mesh Refinement"
#      mesh3dFine = TetrahedralMesh()
#      children = mesh3dFine.refine(mesh3d)
#      for pN,cL in children.iteritems():
#          print "Parent Element "+str(pN)
#          print "Child Elements "
#          for c in cL:
#              print str(c.N)
#      mesh3dFine.writeEdgesGnuplot('mesh3dFine')
#      mesh3dFine.viewMeshGnuplotPipe('mesh3dFine')
#      print "Testing writeMeshADH"
#      mesh3d.writeMeshADH('mesh')
#      print "Testing MultilevelMeshes"
#      print "Testing MultilevelEdgeMesh"
#      mlMesh = MultilevelEdgeMesh(3,1,1,refinementLevels=3)
#      for l in range(len(mlMesh.meshList)):
#          meshFile="mesh"+str(l)
#          mlMesh.meshList[l].writeEdgesGnuplot(meshFile)
#          mlMesh.meshList[l].viewMeshGnuplotPipe(meshFile)
#          print "++++++++++++++++Level "+str(l)+"++++++++++++++++++"
#          for e in mlMesh.meshList[l].elementList:
#              print "Parent Element is "+str(e.N)
#              print e.nodes
#              if l < len(mlMesh.meshList)-1:
#                  for ec in mlMesh.elementChildren[l][e.N]:
#                      print "Child Element is "+str(ec.N)
#                      print ec.nodes
#      print "Testing MultilevelTriangularMesh"
#      mlMesh = MultilevelTriangularMesh(3,3,1,refinementLevels=2)
#      level=0
#      for l in range(len(mlMesh.meshList)):
#          meshFile="mesh"+str(l)
#          mlMesh.meshList[l].writeEdgesGnuplot(meshFile)
#          mlMesh.meshList[l].viewMeshGnuplotPipe(meshFile)
#          print "++++++++++++++++Level "+str(l)+"++++++++++++++++++"
#          for e in mlMesh.meshList[l].elementList:
#              print "Parent Element is "+str(e.N)
#              print e.nodes
#              if l < len(mlMesh.meshList)-1:
#                  for ec in mlMesh.elementChildren[l][e.N]:
#                      print "Child Element is "+str(ec.N)
#                      print ec.nodes
#      print "Testing MultiLevlTetrahedralMesh"
#      mlMesh = MultilevelTetrahedralMesh(3,3,3,refinementLevels=3)
#      level=0
#      for l in range(len(mlMesh.meshList)):
#          meshFile="mesh"+str(l)
#          mlMesh.meshList[l].writeEdgesGnuplot(meshFile)
#          mlMesh.meshList[l].viewMeshGnuplotPipe(meshFile)
#          print "++++++++++++++++Level "+str(l)+"++++++++++++++++++"
#          for e in mlMesh.meshList[l].elementList:
#              print "Parent Element is "+str(e.N)
#              print e.nodes
#              if l < len(mlMesh.meshList)-1:
#                  for ec in mlMesh.elementChildren[l][e.N]:
#                      print "Child Element is "+str(ec.N)
#                      print ec.nodes
#     #
#     # debuggin code from mwf
#     #
#     #how much junk to print out
#     verboseLevel = 2
#     #first just create a simple triangular mesh and look at it in a
#     #couple of different ways
#     Lx = 1.0   #domain length in x and y
#     Ly = 1.0

#     #number of nodes for rectangular grid upon which triangular mesh
#     #will be built should get 2 triangles for each rectangle
#     #(nx-1)(ny-1) in the original grid
#     nx = 3
#     ny = 3

#     #flag for viewing mesh in construction
#     #0 -- do nothing (default)
#     #1 -- gnuplot
#     #2 -- matlab
#     viewMesh = 2
#     meshFileBase='mesh2d'
#     nz = 1
#     Lz = 1.0
#     grid = RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
#     #grid2d.writeEdgesGnuplot('grid2d')
#     #grid2d.viewMeshGnuplotPipe('grid2d')

#     mesh = TriangularMesh()

#     mesh.rectangularToTriangular(grid)

#     if viewMesh == 1:
#         #print mesh in gnuplot format
#         mesh.writeEdgesGnuplot(meshFileBase)
#         #can view with
#         #mesh.viewMeshGnuplotPipe(meshFileBase)
#     elif viewMesh == 2:
#         mesh.writeEdgesMatlab(meshFileBase)
#         #view in matlab with meshFileBase.m
#     #end else

#     print 'mesh Info says \n',mesh.meshInfo()
#     fileName2 = 'meshV2'
#     mp,me,mt = mesh.buildMatlabMeshDataStructures(fileName2)

#     if verboseLevel > 1:
#         #do brute force loop through array to look at it
#         print 'matlab node array is '
#         for j in xrange(mp.shape[1]): #number of columns is number of nodes
#             print '\t',mp[0,j],' ',mp[1,j]
#         #end for

#         #do brute force loop through edge array too
#         print 'matlab edge array holds (matlab edge id, node 0, node 1)'
#         print 'note base 0'
#         for j in xrange(me.shape[1]): #number of columns is number of edges
#             print '\t',me[4,j]-1,' ',me[0,j]-1,' ',me[1,j]-1
#         #end for

#         #do brute force loop through element array too
#         print 'matlab elem array (matlab elem id, node 0, node 1, node 3)'
#         print 'note base 0'
#         for j in xrange(mt.shape[1]): #number of columns is number of edges
#             print '\t',j,' ',mt[0,j]-1,' ',mt[1,j]-1,' ',mt[2,j]-1
#         #end for
#     #end verbose print out for mesh
# #def testEdgeToElementMapping(mesh):
#     """
#     test mesh interface for going from a global edge identifier to its 2
#     neighboring elements.

#       globElem = mesh.elementBoundaryElementsArray[globEdge,neigId]

#     where
#       globEdge is a global edge identifier, neigId is 0,1 for interior edges
#       and 0 for boundary edges (I think). globElem is the global element id
#       for the element on local side neigId.

#     I'm not sure about what I can deduce from the value of neigId in
#     terms of the orientation of the edge and neighboring elements.


#       mesh.exteriorBoundaryElementsArray holds the list of edges on
#       the physical boundary and similarly,
#       mesh.interiorBoundaryElementsArray holds the interior edges.


#     """
#     print "printing mesh edges and neighboring elements"
#     print "format is globEdgeId : locId ---> element Id "
#     for ie in range(mesh.elementBoundaryElementsArray.shape[0]):
#         for neig in range(len(mesh.elementBoundaryElementsArray[ie,:])):
#             elid = mesh.elementBoundaryElementsArray[ie,neig]
#             print "\t ",ie," : ",neig," ---> ",elid
#         #end loop through local element neigs
#     #end loop through global edges
#     print "printing mesh edges and neighboring elements that are defined"
#     for ie in range(mesh.elementBoundaryElementsArray.shape[0]):
#         for neig in range(len(mesh.elementBoundaryElementsArray[ie,:])):
#             elid = mesh.elementBoundaryElementsArray[ie,neig]
#             if (elid > -1):
#                 print "\t ",ie," : ",neig," ---> ",elid
#             #end check if valid index
#         #end loop through local element neigs
#     #end loop through global edges


#     print "print element neighbors for interior edges"
#     for ieI in range(mesh.nInteriorElementBoundaries_global):
#         ie = mesh.interiorElementBoundariesArray[ieI]
#         for neig in range(len(mesh.elementBoundaryElementsArray[ie,:])):
#             elid = mesh.elementBoundaryElementsArray[ie,neig]
#             print "\t ",ie," : ",neig," ---> ",elid
#         #end loop through local element neigs
#     #end loop through global edges

#     print "print element neighbors for exterior edges"
#     for ieE in range(mesh.nExteriorElementBoundaries_global):
#         ie = mesh.exteriorElementBoundariesArray[ieE]
#         for neig in range(len(mesh.elementBoundaryElementsArray[ie,:])):
#             elid = mesh.elementBoundaryElementsArray[ie,neig]
#             print "\t ",ie," : ",neig," ---> ",elid
#         #end loop through local element neigs
#     #end loop through global edges
#end testEdgeToElementMapping
