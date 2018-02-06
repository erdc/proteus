from numpy import dot, sum
import numpy.testing as npt
from proteus.EGeometry import (X,Y,Z)

from proteus.Quadrature import (SimplexGaussQuadrature,
                                SimplexLobattoQuadrature,
                                CubeGaussQuadrature,
                                LobattoEdgeAlt,
                                CompositeTrapezoidalEdge,
                                FaceBarycenterEdge,
                                CompositeTrapezoidalTriangle,
                                CompositeTriangle,
                                FaceBarycenterTriangle)

gaussPoint=SimplexGaussQuadrature(nd=0, order=1)
gaussEdge=SimplexGaussQuadrature(nd=1, order=1)
compositeTrapezoidalEdge = CompositeTrapezoidalEdge(order=1)
faceBarycenterEdge = FaceBarycenterEdge(order=1)
gaussTriangle=SimplexGaussQuadrature(nd=2, order=1)
compositeTrapezoidalTriangle = CompositeTrapezoidalTriangle(order=1)
faceBarycenterTriangle = FaceBarycenterTriangle(order=1)
gaussTetrahedron=SimplexGaussQuadrature(nd=3, order=1)
gaussSquare=CubeGaussQuadrature(nd=2, order=1)
gaussCube=CubeGaussQuadrature(nd=3, order=1)

#define some simple functions to integrate
a=1.1
b=0.92
c=1.34

def f0(x):
    return map(lambda y: 1.0, x)

def f1(x):
    return map(lambda y: 1.0 + a*y[X] + b*y[Y] + c*y[Z], x)

def f2(x):
    return map(lambda y: 1.0 + a*y[X]**2 + c*y[Y]**2 + b*y[Z]**2, x)

def f3(x):
    return map(lambda y: 1.0 + b*y[X]**3 + a*y[Y]**3 + c*y[Z]**3, x)

def f4(x):
    return map(lambda y: 1.0 + c*y[X]**4 + b*y[Y]**4 + a*y[Z]**4, x)

def test_gauss_point():
    gaussPoint.setOrder(1)
    int0_f4 = dot(f4(gaussPoint.points),gaussPoint.weights)
    print int0_f4
    gaussPoint.setOrder(2)
    int1_f4 = dot(f4(gaussPoint.points),gaussPoint.weights)
    print int1_f4
    assert(int0_f4 == int1_f4)
    
def test_gauss_tri4():
    print "4th Order Polynomial"
    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int0_f4
    gaussTriangle.setOrder(2)
    int1_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int1_f4
    gaussTriangle.setOrder(3)
    int2_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int2_f4
    gaussTriangle.setOrder(4)
    int3_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int3_f4
    gaussTriangle.setOrder(5)
    int4_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int4_f4
    gaussTriangle.setOrder(6)
    int5_f4 = dot(f4(gaussTriangle.points),gaussTriangle.weights)
    print int5_f4
    npt.assert_almost_equal(int3_f4,int4_f4)
    npt.assert_almost_equal(int4_f4,int5_f4)

def test_gauss_tet4():
    print "4th Order Polynomial"
    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f4
    gaussTetrahedron.setOrder(2)
    int1_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f4
    gaussTetrahedron.setOrder(3)
    int2_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f4
    gaussTetrahedron.setOrder(4)
    int3_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int3_f4
    gaussTetrahedron.setOrder(5)
    int4_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int4_f4
    gaussTetrahedron.setOrder(6)
    int5_f4 = dot(f4(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int5_f4
    npt.assert_almost_equal(int3_f4,int4_f4)
    npt.assert_almost_equal(int4_f4,int5_f4)

def test_gauss_edge3():
    print "3rd Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int0_f3
    gaussEdge.setOrder(2)
    int1_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int1_f3
    gaussEdge.setOrder(3)
    int2_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int2_f3
    gaussEdge.setOrder(4)
    int3_f3 = dot(f3(gaussEdge.points),gaussEdge.weights)
    print int3_f3
    npt.assert_almost_equal(int2_f3,int3_f3)

def test_gauss_tri3():
    print "3rd Order Polynomial"
    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int0_f3
    gaussTriangle.setOrder(2)
    int1_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int1_f3
    gaussTriangle.setOrder(3)
    int2_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int2_f3
    gaussTriangle.setOrder(4)
    int3_f3 = dot(f3(gaussTriangle.points),gaussTriangle.weights)
    print int3_f3
    npt.assert_almost_equal(int2_f3,int3_f3)

def test_gauss_tet3():
    print "3rd Order Polynomial"
    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f3
    gaussTetrahedron.setOrder(2)
    int1_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f3
    gaussTetrahedron.setOrder(3)
    int2_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f3
    gaussTetrahedron.setOrder(4)
    int3_f3 = dot(f3(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int3_f3
    npt.assert_almost_equal(int2_f3,int3_f3)

def test_gauss_edge2():
    print "2nd Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int0_f2
    gaussEdge.setOrder(2)
    int1_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int1_f2
    gaussEdge.setOrder(3)
    int2_f2 = dot(f2(gaussEdge.points),gaussEdge.weights)
    print int2_f2
    npt.assert_almost_equal(int1_f2,int2_f2)

def test_gauss_tri2():
    print "2nd Order Polynomial"
    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int0_f2
    gaussTriangle.setOrder(2)
    int1_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int1_f2
    gaussTriangle.setOrder(3)
    int2_f2 = dot(f2(gaussTriangle.points),gaussTriangle.weights)
    print int2_f2
    npt.assert_almost_equal(int1_f2,int2_f2)
    
def test_gauss_tet2():
    print "2nd Order Polynomial"
    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f2
    gaussTetrahedron.setOrder(2)
    int1_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f2
    gaussTetrahedron.setOrder(3)
    int2_f2 = dot(f2(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f2
    npt.assert_almost_equal(int1_f2,int2_f2)
    
def test_gauss_edge1():
    print "1st Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int0_f1
    gaussEdge.setOrder(2)
    int1_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int1_f1
    gaussEdge.setOrder(3)
    int2_f1 = dot(f1(gaussEdge.points),gaussEdge.weights)
    print int1_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)

def test_gauss_tri1():
    print "1st Order Polynomial"
    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int0_f1
    gaussTriangle.setOrder(2)
    int1_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int1_f1
    gaussTriangle.setOrder(3)
    int2_f1 = dot(f1(gaussTriangle.points),gaussTriangle.weights)
    print int1_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)

def test_gauss_tet1():
    print "1st Order Polynomial"
    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f1
    gaussTetrahedron.setOrder(2)
    int1_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f1
    gaussTetrahedron.setOrder(3)
    int2_f1 = dot(f1(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)
    
def test_gauss_edge0():
    print "0th Order Polynomial"
    print "Edge"
    gaussEdge.setOrder(1)
    int0_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int0_f0
    gaussEdge.setOrder(2)
    int1_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int1_f0
    gaussEdge.setOrder(3)
    int2_f0 = dot(f0(gaussEdge.points),gaussEdge.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_gauss_tri0():
    print "0th Order Polynomial"
    print "Triangle"
    gaussTriangle.setOrder(1)
    int0_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int0_f0
    gaussTriangle.setOrder(2)
    int1_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int1_f0
    gaussTriangle.setOrder(3)
    int2_f0 = dot(f0(gaussTriangle.points),gaussTriangle.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_gauss_tet0():
    print "0th Order Polynomial"
    print "Tetrahedron"
    gaussTetrahedron.setOrder(1)
    int0_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int0_f0
    gaussTetrahedron.setOrder(2)
    int1_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int1_f0
    gaussTetrahedron.setOrder(3)
    int2_f0 = dot(f0(gaussTetrahedron.points),gaussTetrahedron.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_gauss_square4():
    print "4th Order Polynomial"
    print "Square"
    gaussSquare.setOrder(1)
    int0_f4 = dot(f4(gaussSquare.points),gaussSquare.weights)
    print int0_f4
    gaussSquare.setOrder(2)
    int1_f4 = dot(f4(gaussSquare.points),gaussSquare.weights)
    print int1_f4
    gaussSquare.setOrder(3)
    int2_f4 = dot(f4(gaussSquare.points),gaussSquare.weights)
    print int2_f4
    gaussSquare.setOrder(4)
    int3_f4 = dot(f4(gaussSquare.points),gaussSquare.weights)
    print int3_f4
    gaussSquare.setOrder(5)
    int4_f4 = dot(f4(gaussSquare.points),gaussSquare.weights)
    print int4_f4
    npt.assert_almost_equal(int3_f4,int4_f4)

def test_gauss_cube4():
    print "4th Order Polynomial"
    print "Cube"
    gaussCube.setOrder(1)
    int0_f4 = dot(f4(gaussCube.points),gaussCube.weights)
    print int0_f4
    gaussCube.setOrder(2)
    int1_f4 = dot(f4(gaussCube.points),gaussCube.weights)
    print int1_f4
    gaussCube.setOrder(3)
    int2_f4 = dot(f4(gaussCube.points),gaussCube.weights)
    print int2_f4
    gaussCube.setOrder(4)
    int3_f4 = dot(f4(gaussCube.points),gaussCube.weights)
    print int3_f4
    gaussCube.setOrder(5)
    int4_f4 = dot(f4(gaussCube.points),gaussCube.weights)
    print int4_f4
    npt.assert_almost_equal(int3_f4,int4_f4)

def test_gauss_square3():
    print "3rd Order Polynomial"
    print "Square"
    gaussSquare.setOrder(1)
    int0_f3 = dot(f3(gaussSquare.points),gaussSquare.weights)
    print int0_f3
    gaussSquare.setOrder(2)
    int1_f3 = dot(f3(gaussSquare.points),gaussSquare.weights)
    print int1_f3
    gaussSquare.setOrder(3)
    int2_f3 = dot(f3(gaussSquare.points),gaussSquare.weights)
    print int2_f3
    gaussSquare.setOrder(4)
    int3_f3 = dot(f3(gaussSquare.points),gaussSquare.weights)
    print int3_f3
    npt.assert_almost_equal(int2_f3,int3_f3)

def test_gauss_cube3():
    print "3rd Order Polynomial"
    print "Cube"
    gaussCube.setOrder(1)
    int0_f3 = dot(f3(gaussCube.points),gaussCube.weights)
    print int0_f3
    gaussCube.setOrder(2)
    int1_f3 = dot(f3(gaussCube.points),gaussCube.weights)
    print int1_f3
    gaussCube.setOrder(3)
    int2_f3 = dot(f3(gaussCube.points),gaussCube.weights)
    print int2_f3
    gaussCube.setOrder(4)
    int3_f3 = dot(f3(gaussCube.points),gaussCube.weights)
    print int3_f3
    npt.assert_almost_equal(int2_f3,int3_f3)

def test_gauss_square2():
    print "2nd Order Polynomial"
    print "Square"
    gaussSquare.setOrder(1)
    int0_f2 = dot(f2(gaussSquare.points),gaussSquare.weights)
    print int0_f2
    gaussSquare.setOrder(2)
    int1_f2 = dot(f2(gaussSquare.points),gaussSquare.weights)
    print int1_f2
    gaussSquare.setOrder(3)
    int2_f2 = dot(f2(gaussSquare.points),gaussSquare.weights)
    print int2_f2
    npt.assert_almost_equal(int1_f2,int2_f2)
    
def test_gauss_cube2():
    print "2nd Order Polynomial"
    print "Cube"
    gaussCube.setOrder(1)
    int0_f2 = dot(f2(gaussCube.points),gaussCube.weights)
    print int0_f2
    gaussCube.setOrder(2)
    int1_f2 = dot(f2(gaussCube.points),gaussCube.weights)
    print int1_f2
    gaussCube.setOrder(3)
    int2_f2 = dot(f2(gaussCube.points),gaussCube.weights)
    print int2_f2
    npt.assert_almost_equal(int1_f2,int2_f2)
    
def test_gauss_square1():
    print "1st Order Polynomial"
    print "Square"
    gaussSquare.setOrder(1)
    int0_f1 = dot(f1(gaussSquare.points),gaussSquare.weights)
    print int0_f1
    gaussSquare.setOrder(2)
    int1_f1 = dot(f1(gaussSquare.points),gaussSquare.weights)
    print int1_f1
    gaussSquare.setOrder(3)
    int2_f1 = dot(f1(gaussSquare.points),gaussSquare.weights)
    print int1_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)

def test_gauss_cube1():
    print "1st Order Polynomial"
    print "Cube"
    gaussCube.setOrder(1)
    int0_f1 = dot(f1(gaussCube.points),gaussCube.weights)
    print int0_f1
    gaussCube.setOrder(2)
    int1_f1 = dot(f1(gaussCube.points),gaussCube.weights)
    print int1_f1
    gaussCube.setOrder(3)
    int2_f1 = dot(f1(gaussCube.points),gaussCube.weights)
    print int2_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)
    
def test_gauss_square0():
    print "0th Order Polynomial"
    print "Square"
    gaussSquare.setOrder(1)
    int0_f0 = dot(f0(gaussSquare.points),gaussSquare.weights)
    print int0_f0
    gaussSquare.setOrder(2)
    int1_f0 = dot(f0(gaussSquare.points),gaussSquare.weights)
    print int1_f0
    gaussSquare.setOrder(3)
    int2_f0 = dot(f0(gaussSquare.points),gaussSquare.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_gauss_cube0():
    print "0th Order Polynomial"
    print "Cube"
    gaussCube.setOrder(1)
    int0_f0 = dot(f0(gaussCube.points),gaussCube.weights)
    print int0_f0
    gaussCube.setOrder(2)
    int1_f0 = dot(f0(gaussCube.points),gaussCube.weights)
    print int1_f0
    gaussCube.setOrder(3)
    int2_f0 = dot(f0(gaussCube.points),gaussCube.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_compositeTrapezoidal_edge1():
    print "1st Order Polynomial"
    print "Edge"
    compositeTrapezoidalEdge.setOrder(1)
    int0_f1 = dot(f1(compositeTrapezoidalEdge.points),compositeTrapezoidalEdge.weights)
    print int0_f1
    compositeTrapezoidalEdge.setOrder(2)
    int1_f1 = dot(f1(compositeTrapezoidalEdge.points),compositeTrapezoidalEdge.weights)
    print int1_f1
    compositeTrapezoidalEdge.setOrder(3)
    int2_f1 = dot(f1(compositeTrapezoidalEdge.points),compositeTrapezoidalEdge.weights)
    print int2_f1
    compositeTrapezoidalEdge.setOrder(4)
    int3_f1 = dot(f1(compositeTrapezoidalEdge.points),compositeTrapezoidalEdge.weights)
    print int3_f1
    compositeTrapezoidalEdge.setOrder(5)
    int4_f1 = dot(f1(compositeTrapezoidalEdge.points),compositeTrapezoidalEdge.weights)
    print int4_f1
    npt.assert_almost_equal(int3_f1,int4_f1)
    
def test_faceBarycenter_edge1():
    print "1st Order Polynomial"
    print "Edge"
    faceBarycenterEdge.setOrder(1)
    int0_f1 = dot(f1(faceBarycenterEdge.points),faceBarycenterEdge.weights)
    print int0_f1
    faceBarycenterEdge.setOrder(2)
    int1_f1 = dot(f1(faceBarycenterEdge.points),faceBarycenterEdge.weights)
    print int1_f1
    faceBarycenterEdge.setOrder(3)
    int2_f1 = dot(f1(faceBarycenterEdge.points),faceBarycenterEdge.weights)
    print int2_f1
    faceBarycenterEdge.setOrder(4)
    int3_f1 = dot(f1(faceBarycenterEdge.points),faceBarycenterEdge.weights)
    print int3_f1
    faceBarycenterEdge.setOrder(5)
    int4_f1 = dot(f1(faceBarycenterEdge.points),faceBarycenterEdge.weights)
    print int4_f1
    npt.assert_almost_equal(int3_f1,int4_f1)
    
def test_compositeTrapezoidal_triangle1():
    print "1st Order Polynomial"
    print "Triangle"
    compositeTrapezoidalTriangle.setOrder(1)
    int0_f1 = dot(f1(compositeTrapezoidalTriangle.points),compositeTrapezoidalTriangle.weights)
    print sum(compositeTrapezoidalTriangle.weights)
    print int0_f1
    compositeTrapezoidalTriangle.setOrder(2)
    int1_f1 = dot(f1(compositeTrapezoidalTriangle.points),compositeTrapezoidalTriangle.weights)
    print sum(compositeTrapezoidalTriangle.weights)
    print int1_f1
    compositeTrapezoidalTriangle.setOrder(3)
    int2_f1 = dot(f1(compositeTrapezoidalTriangle.points),compositeTrapezoidalTriangle.weights)
    print sum(compositeTrapezoidalTriangle.weights)
    print int2_f1
    compositeTrapezoidalTriangle.setOrder(4)
    int3_f1 = dot(f1(compositeTrapezoidalTriangle.points),compositeTrapezoidalTriangle.weights)
    print sum(compositeTrapezoidalTriangle.weights)
    print int3_f1
    compositeTrapezoidalTriangle.setOrder(5)
    int4_f1 = dot(f1(compositeTrapezoidalTriangle.points),compositeTrapezoidalTriangle.weights)
    print sum(compositeTrapezoidalTriangle.weights)
    print int4_f1
    npt.assert_almost_equal(int3_f1,int4_f1)
    
def test_faceBarycenter_triangle1():
    print "1st Order Polynomial"
    print "Triangle"
    faceBarycenterTriangle.setOrder(1)
    int0_f1 = dot(f1(faceBarycenterTriangle.points),faceBarycenterTriangle.weights)
    print sum(faceBarycenterTriangle.weights)
    print int0_f1
    faceBarycenterTriangle.setOrder(2)
    int1_f1 = dot(f1(faceBarycenterTriangle.points),faceBarycenterTriangle.weights)
    print sum(faceBarycenterTriangle.weights)
    print int1_f1
    faceBarycenterTriangle.setOrder(3)
    int2_f1 = dot(f1(faceBarycenterTriangle.points),faceBarycenterTriangle.weights)
    print sum(faceBarycenterTriangle.weights)
    print int2_f1
    faceBarycenterTriangle.setOrder(4)
    int3_f1 = dot(f1(faceBarycenterTriangle.points),faceBarycenterTriangle.weights)
    print sum(faceBarycenterTriangle.weights)
    print int3_f1
    faceBarycenterTriangle.setOrder(5)
    int4_f1 = dot(f1(faceBarycenterTriangle.points),faceBarycenterTriangle.weights)
    print sum(faceBarycenterTriangle.weights)
    print int4_f1
    npt.assert_almost_equal(int3_f1,int4_f1)
    

lobattoPoint=SimplexLobattoQuadrature(nd=0, order=1)
lobattoEdge=SimplexLobattoQuadrature(nd=1, order=1)
lobattoEdgeAlt=LobattoEdgeAlt(order=1)
lobattoTriangle=SimplexLobattoQuadrature(nd=2, order=1)
lobattoTetrahedron=SimplexLobattoQuadrature(nd=3, order=1)

def test_lobatto_point():
    lobattoPoint.setOrder(1)
    int0_f4 = dot(f4(lobattoPoint.points),lobattoPoint.weights)
    print int0_f4
    lobattoPoint.setOrder(2)
    int1_f4 = dot(f4(lobattoPoint.points),lobattoPoint.weights)
    print int1_f4
    assert(int0_f4 == int1_f4)

def test_lobatto_edge1():
    print "1st Order Polynomial"
    print "Edge"
    lobattoEdge.setOrder(1)
    int0_f1 = dot(f1(lobattoEdge.points),lobattoEdge.weights)
    print int0_f1
    lobattoEdge.setOrder(2)
    int1_f1 = dot(f1(lobattoEdge.points),lobattoEdge.weights)
    print int1_f1
    lobattoEdge.setOrder(3)
    int2_f1 = dot(f1(lobattoEdge.points),lobattoEdge.weights)
    print int1_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)

def test_lobatto_edgeAlt1():
    print "1st Order Polynomial"
    print "Edge"
    lobattoEdgeAlt.setOrder(1)
    int0_f1 = dot(f1(lobattoEdgeAlt.points),lobattoEdgeAlt.weights)
    print int0_f1
    lobattoEdgeAlt.setOrder(2)
    int1_f1 = dot(f1(lobattoEdgeAlt.points),lobattoEdgeAlt.weights)
    print int1_f1
    lobattoEdgeAlt.setOrder(3)
    int2_f1 = dot(f1(lobattoEdgeAlt.points),lobattoEdgeAlt.weights)
    print int1_f1
    lobattoEdgeAlt.setOrder(4)
    int3_f1 = dot(f1(lobattoEdgeAlt.points),lobattoEdgeAlt.weights)
    print int3_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)
    npt.assert_almost_equal(int3_f1,int3_f1)

def test_lobatto_tri1():
    print "1st Order Polynomial"
    print "Triangle"
    lobattoTriangle.setOrder(1)
    int0_f1 = dot(f1(lobattoTriangle.points),lobattoTriangle.weights)
    print int0_f1
    lobattoTriangle.setOrder(2)
    int1_f1 = dot(f1(lobattoTriangle.points),lobattoTriangle.weights)
    print int1_f1
    lobattoTriangle.setOrder(3)
    int2_f1 = dot(f1(lobattoTriangle.points),lobattoTriangle.weights)
    print int1_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)

def test_lobatto_tet1():
    print "1st Order Polynomial"
    print "Tetrahedron"
    lobattoTetrahedron.setOrder(1)
    int0_f1 = dot(f1(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int0_f1
    lobattoTetrahedron.setOrder(2)
    int1_f1 = dot(f1(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int1_f1
    lobattoTetrahedron.setOrder(3)
    int2_f1 = dot(f1(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int2_f1
    npt.assert_almost_equal(int0_f1,int1_f1)
    npt.assert_almost_equal(int1_f1,int2_f1)
    
def test_lobatto_edge0():
    print "0th Order Polynomial"
    print "Edge"
    lobattoEdge.setOrder(1)
    int0_f0 = dot(f0(lobattoEdge.points),lobattoEdge.weights)
    print int0_f0
    lobattoEdge.setOrder(2)
    int1_f0 = dot(f0(lobattoEdge.points),lobattoEdge.weights)
    print int1_f0
    lobattoEdge.setOrder(3)
    int2_f0 = dot(f0(lobattoEdge.points),lobattoEdge.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_lobatto_tri0():
    print "0th Order Polynomial"
    print "Triangle"
    lobattoTriangle.setOrder(1)
    int0_f0 = dot(f0(lobattoTriangle.points),lobattoTriangle.weights)
    print int0_f0
    lobattoTriangle.setOrder(2)
    int1_f0 = dot(f0(lobattoTriangle.points),lobattoTriangle.weights)
    print int1_f0
    lobattoTriangle.setOrder(3)
    int2_f0 = dot(f0(lobattoTriangle.points),lobattoTriangle.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

def test_lobatto_tet0():
    print "0th Order Polynomial"
    print "Tetrahedron"
    lobattoTetrahedron.setOrder(1)
    int0_f0 = dot(f0(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int0_f0
    lobattoTetrahedron.setOrder(2)
    int1_f0 = dot(f0(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int1_f0
    lobattoTetrahedron.setOrder(3)
    int2_f0 = dot(f0(lobattoTetrahedron.points),lobattoTetrahedron.weights)
    print int2_f0
    npt.assert_almost_equal(int0_f0,int1_f0)
    npt.assert_almost_equal(int1_f0,int2_f0)

if __name__ == '__main__':
    unittest.main(verbosity=2)
