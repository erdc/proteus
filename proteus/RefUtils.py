#!/usr/bin/env python
"""
Classes, functions, and some global data that are useful
for doing FEM calculations on reference elements, etc

.. inheritance-diagram:: proteus.RefUtils
   :parts: 1
"""
from .EGeometry import *
from .Quadrature import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#global data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

baryCoords = {}
baryCoords['1d'] = []
baryCoords['1d'].append(lambda xi: 1.0-xi[0])
baryCoords['1d'].append(lambda xi: xi[0])
#
baryCoords['2d'] = []
baryCoords['2d'].append(lambda xi: 1.0-xi[0]-xi[1])
baryCoords['2d'].append(lambda xi: xi[0])
baryCoords['2d'].append(lambda xi: xi[1])
#
baryCoords['3d'] = []
baryCoords['3d'].append(lambda xi: 1.0-xi[0]-xi[1]-xi[2])
baryCoords['3d'].append(lambda xi: xi[0])
baryCoords['3d'].append(lambda xi: xi[1])
baryCoords['3d'].append(lambda xi: xi[2])

#
#
baryGrads = {}
baryGrads['1d'] = []
baryGrads['1d'].append(numpy.array([-1.0]))
baryGrads['1d'].append(numpy.array([ 1.0]))
#
baryGrads['2d'] = []
baryGrads['2d'].append(numpy.array([-1.0,-1.0]))
baryGrads['2d'].append(numpy.array([ 1.0, 0.0]))
baryGrads['2d'].append(numpy.array([ 0.0, 1.0]))
#
baryGrads['3d'] = []
baryGrads['3d'].append(numpy.array([-1.0,-1.0,-1.0]))
baryGrads['3d'].append(numpy.array([ 1.0, 0.0, 0.0]))
baryGrads['3d'].append(numpy.array([ 0.0, 1.0, 0.0]))
baryGrads['3d'].append(numpy.array([ 0.0, 0.0, 1.0]))
#


#keep around nodal points for P2 lagrange shape functions on
#unit simplex
p2refNodes = []
p2refNodes.append(numpy.array([[0.0],[1.0],[0.5]]))
p2refNodes.append(numpy.array([[0.0, 0.0],
                                 [1.0, 0.0],
                                 [0.0, 1.0],
                                 [0.5, 0.0],
                                 [0.5, 0.5],
                                 [0.0, 0.5]]))
p2refNodes.append(numpy.array([[0.0,   0.0,   0.0],
                                 [1.0,   0.0,   0.0],
                                 [0.0,   1.0,   0.0],
                                 [0.0,   0.0,   1.0],
                                 [0.5,   0.0,   0.0],  #(0,1)
                                 [0.5,   0.5,   0.0],  #(1,2)
                                 [0.0,   0.5,   0.5],  #(2,3)
                                 [0.0,   0.5,   0.0],  #(0,2)
                                 [0.5,   0.0,   0.5],  #(1,3)
                                 [0.0,   0.0,   0.5]]))#(0,3)
#todo generate these for use in interpolation conditions
q2refNodes = []
q2refNodes.append(numpy.array([[-1.0],
                               [ 1.0],
                               [ 0.5]]))

q2refNodes.append(numpy.array([[-1.0, -1.0],
                               [-1.0,  1.0],
                               [ 1.0,  1.0],
                               [ 1.0, -1.0],
                               [-1.0,  0.0],
                               [ 0.0,  1.0],
                               [ 1.0,  0.0],
                               [ 0.0, -1.0],
                               [ 0.0,  0.0]]))

q2refNodes.append(numpy.array([[-1.0,   -1.0,   -1.0],#nodes on bottom
                               [ 1.0,   -1.0,   -1.0],
                               [ 1.0,    1.0,   -1.0],
                               [-1.0,    1.0,   -1.0],

                               [-1.0,   -1.0,    1.0],#nodes on top
                               [ 1.0,   -1.0,    1.0],
                               [ 1.0,    1.0,    1.0],
                               [-1.0,    1.0,    1.0],

                               [ 0.0,   -1.0,   -1.0],#nodes on bottom edges
                               [ 1.0,    0.0,   -1.0],
                               [ 0.0,    1.0,   -1.0],
                               [-1.0,    0.0,   -1.0],

                               [-1.0,   -1.0,    0.0],#nodes on side edges
                               [ 1.0,   -1.0,    0.0],
                               [ 1.0,    1.0,    0.0],
                               [-1.0,    1.0,    0.0],

                               [ 0.0,   -1.0,    1.0],#nodes on top edges
                               [ 1.0,    0.0,    1.0],
                               [ 0.0,    1.0,    1.0],
                               [-1.0,    0.0,    1.0],

                               [ 0.0,    0.0,   -1.0],#nodes on face centers
                               [ 0.0,   -1.0,    0.0],
                               [ 1.0,    0.0,    0.0],
                               [ 0.0,    1.0,    0.0],
                               [-1.0,    0.0,    0.0],
                               [ 0.0,    0.0,    1.0],

                               [ 0.0,    0.0,    0.0]])) #node on element center


#which local boundaries are ref nodes "on"
p2tetrahedronLocalBoundaryLookup = {0:[1,2,3],
                                    1:[0,2,3],
                                    2:[0,1,3],
                                    3:[0,1,2],
                                    4:[2,3],
                                    5:[0,3],
                                    6:[0,1],
                                    7:[1,3],
                                    8:[0,2],
                                    9:[1,2]}

quadrilateralLocalBoundaryLookup = {0:[3,0],
                                    1:[0,1],
                                    2:[1,2],
                                    3:[2,3]}

hexahedronLocalBoundaryLookup = {0:[0,1,4],
                                 1:[0,1,2],
                                 2:[0,2,3],
                                 3:[0,3,4],
                                 4:[1,4,5],
                                 5:[1,2,5],
                                 6:[2,3,5],
                                 7:[3,4,5]}

q2quadrilateralLocalBoundaryLookup = {0:[3,0],
                                      1:[0,1],
                                      2:[1,2],
                                      3:[2,3],
                                      4:[0],
                                      5:[1],
                                      6:[2],
                                      7:[3],
                                      8:[]}

q2hexahedronLocalBoundaryLookup = {0:[0,1,4],#corner nodes
                                   1:[0,1,2],
                                   2:[0,2,3],
                                   3:[0,3,4],
                                   4:[1,4,5],
                                   5:[1,2,5],
                                   6:[2,3,5],
                                   7:[3,4,5],
                                   8:[3,4,5],#edge nodes
                                   9:[3,4,5],
                                   10:[3,4,5],
                                   11:[3,4,5],
                                   12:[3,4,5],
                                   13:[3,4,5],
                                   14:[3,4,5],
                                   15:[3,4,5],
                                   16:[3,4,5],
                                   17:[3,4,5],
                                   18:[3,4,5],
                                   19:[3,4,5],
                                   20:[3,4,5],# face nodes
                                   21:[3,4,5],
                                   22:[3,4,5],
                                   23:[3,4,5],
                                   24:[3,4,5],
                                   25:[3,4,5],
                                   26:[]}#center node
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#local convenience functions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def fact(n):
    if n == 0: return 1
    return n*fact(n-1)
#end fact

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#test functions for internal data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def testBarycentricCoords(verbose=0):
    """

    make sure barycentric coordinates working as expected


    """
    #1d
    nd=1

    #xiArray = numpy.array([[0.0],[0.25],[0.5],[0.75],[1.0]])
    xiArray = p2refNodes[nd-1]
    nxi = xiArray.shape[0]
    if verbose > 2:
        print('1d creating xiArray= ',xiArray)
    #end if
    lamVals = numpy.zeros((nd+1,nxi),'d')
    dlamVals= numpy.zeros((nd+1,),'d')
    for i,lam in enumerate(baryCoords['1d']):
        for j in range(nxi):
            lamVals[i,j]=lam(xiArray[j])
        #end j
    #end i
    for i,dlam in enumerate(baryGrads['1d']):
        dlamVals[i]=dlam
        #end j
    #end i
    out = """
dim = %d
xi  = %s
la0 = %s
la1 = %s
dla = %s
""" % (nd,xiArray,lamVals[0,:],lamVals[1,:],dlamVals)

    print(out)

    #2d
    nd=2
    #nxi     = 6
    #xiArray = numpy.array([[0.0, 0.0],
    #                         [0.5, 0.0],
    #                         [1.0, 0.0],
    #                         [0.5, 0.5],
    #                         [0.0, 1.0],
    #                         [0.0, 0.5]])
    xiArray = p2refNodes[nd-1]
    nxi = xiArray.shape[0]
    if verbose > 2:
        print('2d creating xiArray= ',xiArray)
    #end if
    lamVals = numpy.zeros((nd+1,nxi),'d')
    dlamVals= numpy.zeros((nd+1,nd),'d')
    for i,lam in enumerate(baryCoords['2d']):
        for j in range(nxi):
            lamVals[i,j]=lam(xiArray[j])
        #end j
    #end i
    for i,dlam in enumerate(baryGrads['2d']):
        dlamVals[i]=dlam
        #end j
    #end i
    out = """
dim = %d
xi  =\n%s
la0 = %s
la1 = %s
la2 = %s
dla =\n%s
""" % (nd,xiArray,lamVals[0,:],lamVals[1,:],lamVals[2,:],dlamVals)

    print(out)

    #3d
    nd=3
    #nxi     = 8
    #xiArray = numpy.array([[0.0,   0.0,   0.0],
    #                         [1.0,   0.0,   0.0],
    #                         [0.0,   1.0,   0.0],
    #                         [0.0,   0.0,   1.0],
    #                         [1./3., 1./3., 0.0],
    #                         [1./3., 0.0,   1./3.],
    #                         [0.0,   1./3., 1./3.],
    #                         [1./3., 1./3., 1./3.]])
    xiArray  = p2refNodes[nd-1]
    nxi = xiArray.shape[0]

    if verbose > 2:
        print('3d creating xiArray= \n',xiArray)
    #end if
    lamVals = numpy.zeros((nd+1,nxi),'d')
    dlamVals= numpy.zeros((nd+1,nd),'d')
    for i,lam in enumerate(baryCoords['3d']):
        for j in range(nxi):
            lamVals[i,j]=lam(xiArray[j])
        #end j
    #end i
    for i,dlam in enumerate(baryGrads['3d']):
        dlamVals[i]=dlam
        #end j
    #end i
    out = """
dim = %d
xi  =\n%s
la0 = %s
la1 = %s
la2 = %s
la3 = %s
dla =\n%s
""" % (nd,xiArray,lamVals[0,:],lamVals[1,:],lamVals[2,:],lamVals[3,:],dlamVals)

    print(out)


#end testBarycentricCoords

## @}
