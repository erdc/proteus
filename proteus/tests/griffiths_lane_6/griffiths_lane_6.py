#! /usr/bin/env python
from math import *
from proteus import Domain

def gl_6_3d(width):
    """
    A 3d levee profile for the sacramento levee modeling
    """
    boundaryLegend = {'left':1,
                      'right':2,
                      'front':3,
                      'back':4,
                      'bottom':5,
                      'top':6,
                      'leftTop':7,
                      'rightTop':8}
    vertices_front = [[0.0,0.0,0.0],#0
                      [0.0,0.0,7.3],#1
                      [33.5,0.0,7.3],#2
                      [33.5 + 21.3/tan(2.0*pi*18.0/360.0),0.0,7.3+21.3],#3
                      [33.5 + 21.3/tan(2.0*pi*18.0/360.0)+7.3,0.0,7.3+21.3],#4
                      [33.5 + 124.4,0.0,7.3],#5
                      [2*33.5+124.4,0.0,7.3],#6
                      [2*33.5+124.4,0.0,0.0]]#7
    vertexFlags_front = [boundaryLegend['left'],#0
                         boundaryLegend['left'],#1
                         boundaryLegend['leftTop'],#2
                         boundaryLegend['leftTop'],#3
                         boundaryLegend['rightTop'],#4
                         boundaryLegend['rightTop'],#5
                         boundaryLegend['right'],#6
                         boundaryLegend['right']]#7
    vertices_back = [[v[0],width,v[2]] for v in vertices_front]
    vertexFlags_back = vertexFlags_front
    vertices = vertices_front + vertices_back
    vertexFlags = vertexFlags_front + vertexFlags_back
    facets = [[[0,1,2,3,4,5,6,7]],#front
              [[8,9,10,11,12,13,14,15]],#back
              [[0,1,9,8]],#left
              [[1,9,10,2]],#leftTop1
              [[2,10,11,3]],#leftTop2
              [[3,11,12,4]],#top
              [[4,12,13,5]],#rightTop2
              [[5,13,14,6]],#rightTop1
              [[6,14,15,7]],#right
              [[0,8,15,7]]]#bottom
    facetFlags = [boundaryLegend['front'],
                  boundaryLegend['back'],
                  boundaryLegend['left'],
                  boundaryLegend['leftTop'],
                  boundaryLegend['leftTop'],
                  boundaryLegend['top'],
                  boundaryLegend['rightTop'],
                  boundaryLegend['rightTop'],
                  boundaryLegend['right'],
                  boundaryLegend['bottom']]
    regions = [[0.001,0.001,0.001]]
    regionFlags = [1]
    print vertices,vertexFlags,facets,facetFlags,regions,regionFlags
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 regions=regions,
                                                 regionFlags=regionFlags)
    domain.boundaryFlags=boundaryLegend
    domain.regionLegend = {'levee':1,'default':0}
    return domain

if __name__=='__main__':
    import os
    domain =  gl_6_3d(5.0)
    domain.writeAsymptote("gl_6_3d")
    domain.writePoly("gl_6_3d")
    domain.writePLY("gl_6_3d")
    print domain.boundaryFlags
    #os.system("asy -V gl_6_3d")
    os.system("tetgen -KVApfen gl_6_3d.poly")
