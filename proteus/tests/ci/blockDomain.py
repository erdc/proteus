#!/usr/bin/env python
import numpy

def genPoly(polyfileBase="blockDomain",
            nx=4,ny=4,
            Lx=1.0,Ly=1.0):
    """
    create a simple block domain in 2d
    """
    dx=Lx/nx; dy = Ly/ny

    vertices = []
    for j in range(ny+1):
        for i in range(nx+1):
            vertices.append((0.0+dx*i,0.0+dy*j))
        #
    #
    nVertices = len(vertices)
    assert nVertices == (nx+1)*(ny+1)

    #
    boundaryTags = {'bottom':1,
                    'right':2,
                    'top':3,
                    'left':5,
                    'interior':0}
    
    base = 1
    #write a segment for each edge
    nSegmentsTotal = ny*(nx+1) + nx*(ny+1)
    segments = {}
    nSegments = 0
    #xfaces
    #i=0, left, i = Nx right
    for i in range(nx+1):
        for j in range(ny):
            #vertical edges
            vb = i + j*(nx+1); vt = i + (j+1)*(nx+1)
            tag = boundaryTags['interior']
            if i == 0: tag = boundaryTags['left']
            if i == nx: tag = boundaryTags['right']
            #segment number, start vertex, end vertex, id
            segments[nSegments] = (nSegments,vb,vt,tag)
            nSegments += 1
    #
    #yfaces
    #j=0, bottom, j=Ny, top
    for j in range(ny+1):
        for i in range(nx):
            vl = i + j*(nx+1); vr = i+1 + j*(nx+1)
            tag = boundaryTags['interior']
            if j == 0: tag = boundaryTags['bottom']
            if j == ny: tag = boundaryTags['top']
            segments[nSegments] = (nSegments,vl,vr,tag)
            nSegments += 1
    #
    assert nSegments == nSegmentsTotal

    #return a table to identify regions by a unique flag too
    regions = {}
    curRegion = 0
    for j in range(ny):
        for i in range(nx):
            #region number, x,y, region id
            regions[(i,j)] =(curRegion + base,0.+(i+0.5)*dx,0.+(j+0.5)*dy,curRegion)
            curRegion += 1
        #
    #
    polyfile = open(polyfileBase+'.poly','w')
    polyfile.write('#poly file for [%s,%s] domain with %s x %s blocks \n' % (Lx,Ly,nx,ny))
    polyfile.write('%d %d %d %d \n' % (nVertices,2,1,0))
    polyfile.write('#vertices \n')
    for iv in range(len(vertices)):
        polyfile.write('%d %12.5e %12.5e %d \n' % (iv+base,vertices[iv][0],vertices[iv][1],1))
    #

    #write a segment for each edge
    polyfile.write('%d %d \n' % (nSegments,1))
    polyfile.write('#segments \n')
    for seg in range(nSegments):
        polyfile.write('%d %d %d %d \n ' % (segments[seg][0]+base,
                                            segments[seg][1]+base,
                                            segments[seg][2]+base,
                                            segments[seg][3]))
    polyfile.write('#holes\n 0\n')
    polyfile.write('#regions\n')
    nRegions = nx*ny
    polyfile.write('%d \n' % nRegions)
    curRegion = 0
    for j in range(ny):
        for i in range(nx):
            polyfile.write('%d %12.5e %12.5e %d \n' % (regions[(i,j)][0],
                                                       regions[(i,j)][1],
                                                       regions[(i,j)][2],
                                                       regions[(i,j)][3]))
    #
    polyfile
    return (Lx,Ly,1.0),boundaryTags,regions
if __name__ == '__main__':
    Lx = 1.0; Ly = 1.0; nx = 4; ny = 2;

    L,bt,reg = genPoly("blockDomain",
                       nx=nx,ny=ny,
                       Lx=Lx,Ly=Ly)

    print("domain= %s nx=%s ny=%s boundaryTags=%s regions=%s " % (L,
                                                                  nx,ny,
                                                                  bt,
                                                                  reg))
